// Copyright (c) 2026 jk9576. All rights reserved.
// Halite: Deterministic Lattice Boltzmann Fluid / Gas Automaton (D2Q9)
// Licensed under the MIT License.
//
// Physical and integer-exact validation scenarios, shared by the test suite
// and the `halite validate` command.

const std = @import("std");
const lbm = @import("lbm.zig");

pub const PoiseuilleResult = struct {
    collision: lbm.Collision,
    tau: f64,
    height: usize,
    steps: u64,
    /// sqrt(sum (u - u_a)^2 / sum u_a^2) over the channel cross-section.
    l2_error: f64,
    /// Numerical and analytical centreline velocity.
    u_max: f64,
    u_max_analytic: f64,
};

/// Body-force-driven channel flow between halfway bounce-back walls, periodic
/// in x. The analytic profile is u(y) = F / (2 nu) (y + 1/2)(H - 1/2 - y) with
/// nu = (tau - 1/2) / 3, using the quantised force actually applied.
pub fn poiseuille(gpa: std.mem.Allocator, collision: lbm.Collision, tau_q16: i32, height: usize, steps: u64) !PoiseuilleResult {
    var l = try lbm.Lattice.create(gpa, 4, height);
    defer l.destroy(gpa);
    try l.setTau(tau_q16);
    const tau = lbm.fromQ(tau_q16);
    const nu = (tau - 0.5) / 3.0;
    const h: f64 = @floatFromInt(height);
    // Choose the force for a centreline speed near 0.05.
    const force_q = @max(1, lbm.toQ(0.05 * 8.0 * nu / (h * h)));
    try l.setConfig(.{
        .collision = collision,
        .x_boundary = .periodic,
        .y_boundary = .wall,
        .force_x_q16 = force_q,
    });
    l.resetFlow();
    for (0..steps) |_| l.step();

    const force = lbm.fromQ(force_q);
    var num: f64 = 0;
    var den: f64 = 0;
    var u_max: f64 = 0;
    for (0..height) |y| {
        const yy: f64 = @floatFromInt(y);
        const ua = force / (2.0 * nu) * (yy + 0.5) * (h - 0.5 - yy);
        // Average the four columns; they are identical for a converged run.
        var u: f64 = 0;
        for (0..4) |x| u += lbm.fromQ(l.velocity(y * 4 + x)[0]);
        u /= 4.0;
        num += (u - ua) * (u - ua);
        den += ua * ua;
        u_max = @max(u_max, u);
    }
    return .{
        .collision = collision,
        .tau = tau,
        .height = height,
        .steps = steps,
        .l2_error = @sqrt(num / den),
        .u_max = u_max,
        .u_max_analytic = force / (8.0 * nu) * h * h,
    };
}

/// Deterministic pseudo-random obstacle field (about 1 in 12 cells).
pub fn scatterObstacles(l: *lbm.Lattice, seed: u64) void {
    var rng = std.Random.DefaultPrng.init(seed);
    const r = rng.random();
    for (0..l.height) |y| {
        for (0..l.width) |x| {
            if (r.uintLessThan(u32, 12) == 0) l.setObstacleCell(x, y, true);
        }
    }
}

pub const ConservationResult = struct {
    steps: u64,
    mass_drift: i64,
    momentum_x_drift: i64,
    momentum_y_drift: i64,
};

/// Fully periodic box with scattered obstacles and a disturbed initial state.
/// Mass is exactly conserved; momentum is not (obstacles absorb it), so the
/// result reports the mass drift and the momentum drift of an obstacle-free run.
pub fn conservation(gpa: std.mem.Allocator, collision: lbm.Collision, steps: u64) !ConservationResult {
    var l = try lbm.Lattice.create(gpa, 48, 32);
    defer l.destroy(gpa);
    try l.setConfig(.{ .collision = collision, .x_boundary = .periodic, .y_boundary = .periodic });
    disturb(&l);
    scatterObstacles(&l, 0x5eed);
    const m0 = l.totals();
    for (0..steps) |_| l.step();
    const m1 = l.totals();

    var free = try lbm.Lattice.create(gpa, 48, 32);
    defer free.destroy(gpa);
    try free.setConfig(.{ .collision = collision, .x_boundary = .periodic, .y_boundary = .periodic });
    disturb(&free);
    const p0 = free.totals();
    for (0..steps) |_| free.step();
    const p1 = free.totals();
    return .{
        .steps = steps,
        .mass_drift = m1.rho - m0.rho,
        .momentum_x_drift = p1.jx - p0.jx,
        .momentum_y_drift = p1.jy - p0.jy,
    };
}

/// A shear wave plus a density bump, built from exact equilibria.
fn disturb(l: *lbm.Lattice) void {
    for (0..l.height) |y| {
        for (0..l.width) |x| {
            const ux: i64 = if ((y / 4) % 2 == 0) 3000 else -2000;
            const uy: i64 = if ((x / 6) % 2 == 0) 1500 else -700;
            const rho: i64 = if (x < 8 and y < 8) lbm.Q + 2000 else lbm.Q;
            var eq: [9]i64 = undefined;
            lbm.equilibrium(rho, ux, uy, &eq);
            for (0..9) |i| l.f[y * l.width + x][i] = lbm.sat32(eq[i]);
        }
    }
}

pub const BalanceResult = struct {
    steps: u64,
    /// Momentum change minus (body force input - obstacle force - wall force), summed over all steps.
    residual_x: i64,
    residual_y: i64,
    /// Obstacle drag in the last step and the body force on the fluid.
    final_drag: i64,
    body_force_total: i64,
};

/// Periodic-x channel with walls and a cylinder, driven by a body force. In
/// integer arithmetic the momentum-exchange forces account for every unit of
/// momentum: P(t+1) - P(t) = N F - F_obstacle - F_wall exactly.
pub fn momentumBalance(gpa: std.mem.Allocator, collision: lbm.Collision, steps: u64) !BalanceResult {
    var l = try lbm.Lattice.create(gpa, 64, 33);
    defer l.destroy(gpa);
    try l.setTau(lbm.Q);
    try l.setConfig(.{
        .collision = collision,
        .x_boundary = .periodic,
        .y_boundary = .wall,
        .force_x_q16 = 40,
        .force_y_q16 = 3,
    });
    l.setCircleObstacle(20, 16, 5, true);
    l.resetFlow();
    var prev = l.totals();
    var res_x: i64 = 0;
    var res_y: i64 = 0;
    const n: i64 = @intCast(l.cellCount() - l.solidCount());
    for (0..steps) |_| {
        l.step();
        const now = l.totals();
        res_x += (now.jx - prev.jx) - (n * l.config.force_x_q16 - l.stats.obstacle_force_x - l.stats.wall_force_x);
        res_y += (now.jy - prev.jy) - (n * l.config.force_y_q16 - l.stats.obstacle_force_y - l.stats.wall_force_y);
        prev = now;
    }
    return .{
        .steps = steps,
        .residual_x = res_x,
        .residual_y = res_y,
        .final_drag = l.stats.obstacle_force_x,
        .body_force_total = n * l.config.force_x_q16,
    };
}

pub const SymmetryResult = struct {
    steps: u64,
    /// Largest |lift| over all steps (Q16.16); exactly 0 for a symmetric run.
    max_abs_lift: i64,
    /// Cells whose distributions differ from their mirror image.
    asymmetric_cells: usize,
    drag: f64,
    drag_coefficient: f64,
};

/// Channel with a cylinder on the centreline of an odd-height domain and a
/// parabolic Zou-He inlet. Sign-symmetric rounding keeps the flow exactly
/// mirror symmetric, so the lift is exactly zero.
pub fn mirrorSymmetry(gpa: std.mem.Allocator, collision: lbm.Collision, steps: u64) !SymmetryResult {
    const height = 41;
    var l = try lbm.Lattice.create(gpa, 120, height);
    defer l.destroy(gpa);
    try l.setTau(lbm.toQ(0.8));
    try l.setConfig(.{
        .collision = collision,
        .x_boundary = .channel,
        .y_boundary = .wall,
        .inflow_q16 = lbm.toQ(0.06),
        .inlet_profile = .parabolic,
    });
    l.resetFlow();
    l.setCircleObstacle(30, 20, 5, true);
    var max_lift: i64 = 0;
    for (0..steps) |_| {
        l.step();
        max_lift = @max(max_lift, @as(i64, @intCast(@abs(l.stats.obstacle_force_y))));
    }
    const mirror = [9]u4{ 0, 1, 4, 3, 2, 8, 7, 6, 5 };
    var asym: usize = 0;
    for (0..height) |y| {
        for (0..l.width) |x| {
            const a = l.f[y * l.width + x];
            const b = l.f[(height - 1 - y) * l.width + x];
            for (0..9) |i| {
                if (a[i] != b[mirror[i]]) {
                    asym += 1;
                    break;
                }
            }
        }
    }
    const drag = lbm.fromQ(l.stats.obstacle_force_x);
    const u_mean = 2.0 / 3.0 * lbm.fromQ(l.config.inflow_q16);
    const d: f64 = @floatFromInt(l.obstacleExtent().rows);
    return .{
        .steps = steps,
        .max_abs_lift = max_lift,
        .asymmetric_cells = asym,
        .drag = drag,
        .drag_coefficient = 2.0 * drag / (u_mean * u_mean * d),
    };
}

pub const ChannelResult = struct {
    steps: u64,
    inlet_flux: f64,
    /// Largest relative deviation of the column mass flux from the inlet flux.
    max_flux_deviation: f64,
    outlet_density: f64,
};

/// Open channel with a uniform Zou-He velocity inlet and a pressure outlet:
/// at steady state every column carries the same mass flux.
pub fn channelFlux(gpa: std.mem.Allocator, collision: lbm.Collision, steps: u64) !ChannelResult {
    var l = try lbm.Lattice.create(gpa, 96, 24);
    defer l.destroy(gpa);
    try l.setConfig(.{ .collision = collision, .x_boundary = .channel, .y_boundary = .wall, .inflow_q16 = lbm.toQ(0.05) });
    l.resetFlow();
    for (0..steps) |_| l.step();
    var inlet: f64 = 0;
    var worst: f64 = 0;
    var out_rho: f64 = 0;
    for (0..l.width) |x| {
        var flux: f64 = 0;
        for (0..l.height) |y| flux += lbm.fromQ(l.velocity(y * l.width + x)[0]);
        if (x == 0) inlet = flux;
        worst = @max(worst, @abs(flux - inlet) / inlet);
    }
    for (0..l.height) |y| out_rho += lbm.fromQ(lbm.momentsOf(l.f[y * l.width + l.width - 1]).rho);
    return .{
        .steps = steps,
        .inlet_flux = inlet,
        .max_flux_deviation = worst,
        .outlet_density = out_rho / @as(f64, @floatFromInt(l.height)),
    };
}
