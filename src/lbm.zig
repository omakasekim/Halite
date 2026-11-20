// Copyright (c) 2026 jk9576. All rights reserved.
// Halite: Deterministic Lattice Boltzmann Fluid / Gas Automaton (D2Q9)
// Licensed under the MIT License.
//
// Integer-only D2Q9 solver. Distributions are stored in Q16.16 fixed point and
// every arithmetic step uses 64-bit intermediates with a single, sign-symmetric
// rounding, so a run is bit-for-bit reproducible on every target and mirror
// symmetric geometries stay exactly mirror symmetric.
//
// Model: incompressible He-Luo equilibrium (rho0 = 1, u = j), BGK or TRT
// collision, Guo body forcing, halfway bounce-back on obstacles, Zou-He velocity
// inlet and pressure outlet, periodic or no-slip domain edges, and
// momentum-exchange force measurement on obstacles and walls. The collision
// restores the exact cell mass and momentum after rounding, so mass is conserved
// to the last bit in closed and periodic domains.

const std = @import("std");

pub const WIDTH: usize = 160;
pub const HEIGHT: usize = 80;
pub const NUM_CELLS: usize = WIDTH * HEIGHT;

// Q16.16 fixed-point math helpers
pub const Q: i32 = 65536; // 1.0 in q16.16

/// Truncating Q16.16 multiply (floor of the exact product). Kept for API
/// compatibility; the solver uses the rounding helpers below.
pub inline fn qMul(a: i32, b: i32) i32 {
    const p = @as(i64, a) * @as(i64, b);
    return @truncate(p >> 16);
}

/// Truncating Q16.16 divide; returns 0 when dividing by zero.
pub inline fn qDiv(a: i32, b: i32) i32 {
    if (b == 0) return 0;
    const p = @as(i64, a) << 16;
    return @truncate(@divTrunc(p, b));
}

/// num / den rounded to nearest, ties away from zero. `den` must be positive.
/// The rounding is odd-symmetric: divRound(-n, d) == -divRound(n, d).
pub inline fn divRound(num: i64, den: i64) i64 {
    const half = @divTrunc(den, 2);
    // Truncation is odd-symmetric, so shifting by +-half before it rounds both signs alike.
    return @divTrunc(num + (if (num >= 0) half else -half), den);
}

/// Rounded Q16.16 multiply.
pub inline fn qMulRound(a: i32, b: i32) i32 {
    return sat32(divRound(@as(i64, a) * @as(i64, b), Q));
}

/// Rounded Q16.16 divide; returns 0 when dividing by zero.
pub inline fn qDivRound(a: i32, b: i32) i32 {
    if (b == 0) return 0;
    const n = @as(i64, a) << 16;
    return sat32(if (b > 0) divRound(n, b) else divRound(-n, -@as(i64, b)));
}

pub inline fn sat32(v: i64) i32 {
    return @intCast(std.math.clamp(v, std.math.minInt(i32), std.math.maxInt(i32)));
}

/// Converts a real number to Q16.16 with rounding, saturating at the i32 range.
pub fn toQ(x: f64) i32 {
    if (std.math.isNan(x)) return 0;
    const scaled = @round(x * 65536.0);
    if (scaled >= 2147483647.0) return std.math.maxInt(i32);
    if (scaled <= -2147483648.0) return std.math.minInt(i32);
    return @intFromFloat(scaled);
}

pub fn fromQ(v: i64) f64 {
    return @as(f64, @floatFromInt(v)) / 65536.0;
}

// D2Q9 lattice velocity vectors: (x, y)
pub const CX = [9]i8{ 0, 1, 0, -1, 0, 1, -1, -1, 1 };
pub const CY = [9]i8{ 0, 0, 1, 0, -1, 1, 1, -1, -1 };

// Opposite direction indices for bounce-back
pub const OPPOSITE = [9]u4{ 0, 3, 4, 1, 2, 7, 8, 5, 6 };

// D2Q9 weights in Q16.16 (rounded), kept for API compatibility:
// w0 = 4/9 (~0.44444) -> 29127
// w1..4 = 1/9 (~0.11111) -> 7282
// w5..8 = 1/36 (~0.02778) -> 1820
pub const W = [9]i32{ 29127, 7282, 7282, 7282, 7282, 1820, 1820, 1820, 1820 };

/// Exact D2Q9 weights as multiples of 1/36, used by the solver.
pub const W36 = [9]i64{ 16, 4, 4, 4, 4, 1, 1, 1, 1 };

pub const TAU_DEFAULT: i32 = Q + Q / 5; // 1.2
/// Smallest accepted relaxation time (about 0.505).
pub const TAU_MIN: i32 = 33096;
pub const TAU_MAX: i32 = 16 * Q;
/// TRT magic parameter 3/16: places halfway bounce-back walls exactly for Poiseuille flow.
pub const MAGIC_POISEUILLE: i32 = 12288;
pub const MAGIC_MAX: i32 = 4 * Q;
pub const INFLOW_DEFAULT: i32 = 5243; // 0.08
pub const INFLOW_MAX: i32 = 19661; // 0.3
pub const FORCE_MAX: i32 = 1024; // 1/64 per step
pub const DENSITY_MIN: i32 = Q / 2;
pub const DENSITY_MAX: i32 = 2 * Q;
pub const MIN_DIM: usize = 3;
pub const MAX_DIM: usize = 8192;
pub const MAX_CELLS: usize = 1 << 24;
/// Velocity clamp used inside the collision so diverging runs cannot overflow.
const U_CLAMP: i64 = 2 * Q;
/// A cell faster than 0.5 lattice units per step is flagged as unstable.
const UNSTABLE_SPEED_SQ: i64 = (Q / 2) * (Q / 2);

pub const Collision = enum(u8) {
    /// Single relaxation time (Bhatnagar-Gross-Krook).
    bgk = 0,
    /// Two relaxation times with a tunable magic parameter.
    trt = 1,
};

pub const XBoundary = enum(u8) {
    /// Zou-He velocity inlet at x = 0 and Zou-He pressure outlet at x = width - 1.
    channel = 0,
    /// Periodic in x.
    periodic = 1,
    /// No-slip bounce-back walls at both x ends (sealed box).
    closed = 2,
};

pub const YBoundary = enum(u8) {
    /// No-slip halfway bounce-back walls at y = -0.5 and y = height - 0.5.
    wall = 0,
    periodic = 1,
};

pub const InletProfile = enum(u8) {
    uniform = 0,
    /// Fully developed parabola with peak `inflow` between the y walls.
    parabolic = 1,
};

pub const RenderMode = enum(u8) {
    /// Legacy speed color map.
    speed = 0,
    /// Signed vorticity: red counter-clockwise, blue clockwise.
    vorticity = 1,
    /// Density deviation from 1 as a gray ramp.
    density = 2,
};

pub const ConfigError = error{
    TauOutOfRange,
    MagicOutOfRange,
    InflowOutOfRange,
    ForceOutOfRange,
    DensityOutOfRange,
};

pub const Config = struct {
    collision: Collision = .bgk,
    magic_q16: i32 = MAGIC_POISEUILLE,
    x_boundary: XBoundary = .channel,
    y_boundary: YBoundary = .wall,
    inflow_q16: i32 = INFLOW_DEFAULT,
    inlet_profile: InletProfile = .uniform,
    outlet_density_q16: i32 = Q,
    force_x_q16: i32 = 0,
    force_y_q16: i32 = 0,

    pub fn validate(self: Config) ConfigError!void {
        if (self.magic_q16 <= 0 or self.magic_q16 > MAGIC_MAX) return error.MagicOutOfRange;
        if (self.inflow_q16 < -INFLOW_MAX or self.inflow_q16 > INFLOW_MAX) return error.InflowOutOfRange;
        if (@abs(self.force_x_q16) > FORCE_MAX or @abs(self.force_y_q16) > FORCE_MAX) return error.ForceOutOfRange;
        if (self.outlet_density_q16 < DENSITY_MIN or self.outlet_density_q16 > DENSITY_MAX) return error.DensityOutOfRange;
    }
};

pub fn validateTau(tau_q16: i32) ConfigError!void {
    if (tau_q16 < TAU_MIN or tau_q16 > TAU_MAX) return error.TauOutOfRange;
}

/// Per-step telemetry. Forces are momentum exchanged during the last streaming
/// step, in Q16.16 lattice units, acting on the solid (positive x = drag).
pub const Stats = struct {
    obstacle_force_x: i64 = 0,
    obstacle_force_y: i64 = 0,
    wall_force_x: i64 = 0,
    wall_force_y: i64 = 0,
    /// Total post-collision mass and momentum over fluid cells (Q16.16).
    mass: i64 = 0,
    momentum_x: i64 = 0,
    momentum_y: i64 = 0,
    fluid_cells: u64 = 0,
    /// Largest |u|^2 over fluid cells, Q32.32.
    max_speed_sq: i64 = 0,
    /// Distributions clamped to the i32 range in the last step.
    saturated: u64 = 0,
    /// Fluid cells whose density was zero or negative in the last step.
    nonpositive_density: u64 = 0,

    pub fn stable(self: Stats) bool {
        return self.saturated == 0 and self.nonpositive_density == 0 and self.max_speed_sq < UNSTABLE_SPEED_SQ;
    }

    pub fn maxSpeed(self: Stats) f64 {
        return @sqrt(@as(f64, @floatFromInt(self.max_speed_sq))) / 65536.0;
    }
};

/// Equilibrium distribution with exact moments: sum = rho, sum c = (ux, uy).
/// All inputs are Q16.16.
pub fn equilibrium(rho: i64, ux: i64, uy: i64, out: *[9]i64) void {
    const usq = ux * ux + uy * uy; // Q32
    var sum: i64 = 0;
    var mx: i64 = 0;
    var my: i64 = 0;
    inline for (1..9) |i| {
        const cu = @as(i64, CX[i]) * ux + @as(i64, CY[i]) * uy;
        // Twice w-free equilibrium in Q32: 2 (rho + 3 cu) + 9 cu^2 - 3 u^2
        const t2 = ((rho + 3 * cu) << 17) + 9 * cu * cu - 3 * usq;
        const v = divRound(W36[i] * t2, 36 << 17);
        out[i] = v;
        sum += v;
        mx += CX[i] * v;
        my += CY[i] * v;
    }
    fixMoments(out, rho, ux - mx, uy - my, sum);
}

/// Adds the momentum residual to the direction that moves momentum the right
/// way (sign-symmetric), then gives the rest particle whatever mass is left.
inline fn fixMoments(out: *[9]i64, rho: i64, rx: i64, ry: i64, sum_1_8: i64) void {
    if (rx > 0) out[1] += rx else out[3] -= rx;
    if (ry > 0) out[2] += ry else out[4] -= ry;
    out[0] = rho - (sum_1_8 + @as(i64, @intCast(@abs(rx))) + @as(i64, @intCast(@abs(ry))));
}

pub const Moments = struct { rho: i64, jx: i64, jy: i64 };

pub inline fn momentsOf(f: [9]i32) Moments {
    var m = Moments{ .rho = 0, .jx = 0, .jy = 0 };
    inline for (0..9) |i| {
        m.rho += f[i];
        m.jx += @as(i64, CX[i]) * f[i];
        m.jy += @as(i64, CY[i]) * f[i];
    }
    return m;
}

/// A D2Q9 lattice over caller-owned buffers. Cell (x, y) is at index y * width + x.
pub const Lattice = struct {
    width: usize,
    height: usize,
    f: [][9]i32,
    f_next: [][9]i32,
    solid: []bool,
    tau_q16: i32 = TAU_DEFAULT,
    config: Config = .{},
    stats: Stats = .{},
    steps: u64 = 0,

    pub const DimensionError = error{InvalidDimensions};

    pub fn checkDimensions(width: usize, height: usize) DimensionError!void {
        if (width < MIN_DIM or height < MIN_DIM or width > MAX_DIM or height > MAX_DIM) return error.InvalidDimensions;
        if (width * height > MAX_CELLS) return error.InvalidDimensions;
    }

    /// Allocates a lattice filled with fluid at rest (rho = 1) and no obstacles.
    pub fn create(gpa: std.mem.Allocator, width: usize, height: usize) !Lattice {
        try checkDimensions(width, height);
        const n = width * height;
        const f = try gpa.alloc([9]i32, n);
        errdefer gpa.free(f);
        const f_next = try gpa.alloc([9]i32, n);
        errdefer gpa.free(f_next);
        const solid = try gpa.alloc(bool, n);
        var l = Lattice{ .width = width, .height = height, .f = f, .f_next = f_next, .solid = solid };
        @memset(l.solid, false);
        l.fillEquilibrium(Q, 0, 0);
        return l;
    }

    pub fn destroy(self: *Lattice, gpa: std.mem.Allocator) void {
        gpa.free(self.f);
        gpa.free(self.f_next);
        gpa.free(self.solid);
        self.* = undefined;
    }

    pub inline fn cellCount(self: *const Lattice) usize {
        return self.width * self.height;
    }

    pub fn setTau(self: *Lattice, tau_q16: i32) ConfigError!void {
        try validateTau(tau_q16);
        self.tau_q16 = tau_q16;
    }

    pub fn setConfig(self: *Lattice, config: Config) ConfigError!void {
        try config.validate();
        self.config = config;
    }

    /// Sets every cell, solid or not, to the equilibrium for (rho, ux, uy).
    pub fn fillEquilibrium(self: *Lattice, rho: i32, ux: i32, uy: i32) void {
        var eq: [9]i64 = undefined;
        equilibrium(rho, ux, uy, &eq);
        var cell: [9]i32 = undefined;
        for (0..9) |i| cell[i] = sat32(eq[i]);
        @memset(self.f, cell);
        @memset(self.f_next, cell);
        self.steps = 0;
        self.stats = .{};
    }

    /// Initial flow for the configured boundaries: the inlet velocity for a
    /// channel, fluid at rest otherwise.
    pub fn resetFlow(self: *Lattice) void {
        const ux: i32 = if (self.config.x_boundary == .channel) self.config.inflow_q16 else 0;
        self.fillEquilibrium(Q, ux, 0);
    }

    pub fn clearObstacles(self: *Lattice) void {
        for (0..self.cellCount()) |idx| {
            if (self.solid[idx]) self.setCellState(idx, false);
        }
    }

    pub fn setObstacleCell(self: *Lattice, x: usize, y: usize, is_solid: bool) void {
        if (x < self.width and y < self.height) self.setCellState(y * self.width + x, is_solid);
    }

    /// A cell released from an obstacle starts as fluid at rest with rho = 1.
    fn setCellState(self: *Lattice, idx: usize, is_solid: bool) void {
        if (self.solid[idx] and !is_solid) {
            var eq: [9]i64 = undefined;
            equilibrium(Q, 0, 0, &eq);
            for (0..9) |i| {
                self.f[idx][i] = sat32(eq[i]);
                self.f_next[idx][i] = sat32(eq[i]);
            }
        }
        self.solid[idx] = is_solid;
    }

    pub fn setCircleObstacle(self: *Lattice, cx: i32, cy: i32, radius: i32, is_solid: bool) void {
        if (radius < 0) return;
        const r: i64 = radius;
        const x0: i64 = @max(0, @as(i64, cx) - r);
        const x1: i64 = @min(@as(i64, @intCast(self.width)) - 1, @as(i64, cx) + r);
        const y0: i64 = @max(0, @as(i64, cy) - r);
        const y1: i64 = @min(@as(i64, @intCast(self.height)) - 1, @as(i64, cy) + r);
        var y = y0;
        while (y <= y1) : (y += 1) {
            var x = x0;
            while (x <= x1) : (x += 1) {
                const dx = x - cx;
                const dy = y - cy;
                if (dx * dx + dy * dy <= r * r) self.setObstacleCell(@intCast(x), @intCast(y), is_solid);
            }
        }
    }

    pub fn solidCount(self: *const Lattice) usize {
        var n: usize = 0;
        for (self.solid) |s| n += @intFromBool(s);
        return n;
    }

    /// Number of rows and columns that contain at least one solid cell
    /// (the frontal height and streamwise length of the obstacles).
    pub fn obstacleExtent(self: *const Lattice) struct { rows: usize, cols: usize } {
        var rows: usize = 0;
        for (0..self.height) |y| {
            for (0..self.width) |x| {
                if (self.solid[y * self.width + x]) {
                    rows += 1;
                    break;
                }
            }
        }
        var cols: usize = 0;
        for (0..self.width) |x| {
            for (0..self.height) |y| {
                if (self.solid[y * self.width + x]) {
                    cols += 1;
                    break;
                }
            }
        }
        return .{ .rows = rows, .cols = cols };
    }

    /// Inlet velocity for row y.
    pub fn inletVelocity(self: *const Lattice, y: usize) i64 {
        const u = @as(i64, self.config.inflow_q16);
        if (self.config.inlet_profile == .uniform or self.config.y_boundary == .periodic) return u;
        const h: i64 = @intCast(self.height);
        const yy: i64 = @intCast(y);
        // 4 u (y + 1/2)(H - 1/2 - y) / H^2, walls at y = -1/2 and H - 1/2.
        return divRound(u * (2 * yy + 1) * (2 * h - 1 - 2 * yy), h * h);
    }

    /// Fluid velocity of a cell after the last step (Q16.16): the Guo velocity
    /// j + F/2 that the collision relaxed towards.
    pub fn velocity(self: *const Lattice, idx: usize) [2]i64 {
        const m = momentsOf(self.f[idx]);
        return .{
            m.jx - divRound(self.config.force_x_q16, 2),
            m.jy - divRound(self.config.force_y_q16, 2),
        };
    }

    /// Advances one stream, boundary and collide cycle.
    pub fn step(self: *Lattice) void {
        self.stream();
        if (self.config.x_boundary == .channel) self.applyChannelBoundaries();
        self.collide();
        self.steps +%= 1;
    }

    fn stream(self: *Lattice) void {
        const w = self.width;
        const h = self.height;
        const x_periodic = self.config.x_boundary == .periodic;
        const y_periodic = self.config.y_boundary == .periodic;
        const count_x_walls = self.config.x_boundary == .closed;
        var ofx: i64 = 0;
        var ofy: i64 = 0;
        var wfx: i64 = 0;
        var wfy: i64 = 0;
        const wi: isize = @intCast(w);
        const hi: isize = @intCast(h);

        for (0..h) |y| {
            for (0..w) |x| {
                const idx = y * w + x;
                if (self.solid[idx]) continue;
                const cell = self.f[idx];
                self.f_next[idx][0] = cell[0];
                for (1..9) |i| {
                    var nx: isize = @as(isize, @intCast(x)) + CX[i];
                    var ny: isize = @as(isize, @intCast(y)) + CY[i];
                    var edge = false;
                    var edge_counts = false;
                    if (nx < 0 or nx >= wi) {
                        if (x_periodic) {
                            nx = if (nx < 0) wi - 1 else 0;
                        } else {
                            edge = true;
                            edge_counts = count_x_walls;
                        }
                    }
                    if (ny < 0 or ny >= hi) {
                        if (y_periodic) {
                            ny = if (ny < 0) hi - 1 else 0;
                        } else {
                            edge = true;
                            edge_counts = true;
                        }
                    }
                    const fi = cell[i];
                    if (edge) {
                        self.f_next[idx][OPPOSITE[i]] = fi;
                        if (edge_counts) {
                            wfx += 2 * @as(i64, CX[i]) * fi;
                            wfy += 2 * @as(i64, CY[i]) * fi;
                        }
                        continue;
                    }
                    const n_idx: usize = @intCast(ny * wi + nx);
                    if (self.solid[n_idx]) {
                        self.f_next[idx][OPPOSITE[i]] = fi;
                        ofx += 2 * @as(i64, CX[i]) * fi;
                        ofy += 2 * @as(i64, CY[i]) * fi;
                    } else {
                        self.f_next[n_idx][i] = fi;
                    }
                }
            }
        }
        self.stats.obstacle_force_x = ofx;
        self.stats.obstacle_force_y = ofy;
        self.stats.wall_force_x = wfx;
        self.stats.wall_force_y = wfy;
    }

    /// Zou-He velocity inlet (x = 0) and pressure outlet (x = width - 1),
    /// written for the incompressible model (rho0 = 1).
    fn applyChannelBoundaries(self: *Lattice) void {
        const w = self.width;
        const rho_out: i64 = self.config.outlet_density_q16;
        for (0..self.height) |y| {
            const in_idx = y * w;
            if (!self.solid[in_idx]) {
                const g = &self.f_next[in_idx];
                const jx = self.inletVelocity(y);
                const half_d = @divTrunc(@as(i64, g[2]) - g[4], 2);
                const f1 = @as(i64, g[3]) + divRound(2 * jx, 3);
                const f5 = @as(i64, g[7]) - half_d + divRound(jx, 6);
                const f8 = @as(i64, g[6]) + half_d + divRound(jx, 6);
                g[1] = sat32(f1);
                g[5] = sat32(f5);
                g[8] = sat32(f8);
            }
            const out_idx = y * w + (w - 1);
            if (!self.solid[out_idx]) {
                const g = &self.f_next[out_idx];
                const known = @as(i64, g[0]) + g[2] + g[4] + 2 * (@as(i64, g[1]) + g[5] + g[8]);
                const jx = known - rho_out;
                const half_d = @divTrunc(@as(i64, g[2]) - g[4], 2);
                g[3] = sat32(@as(i64, g[1]) - divRound(2 * jx, 3));
                g[7] = sat32(@as(i64, g[5]) + half_d - divRound(jx, 6));
                g[6] = sat32(@as(i64, g[8]) - half_d - divRound(jx, 6));
            }
        }
    }

    /// Relaxation rates (Q16.16): omega+ = 1/tau, and for TRT omega- from
    /// Lambda = (tau+ - 1/2)(tau- - 1/2).
    pub fn relaxationRates(self: *const Lattice) struct { plus: i64, minus: i64 } {
        const tau = std.math.clamp(self.tau_q16, TAU_MIN, TAU_MAX);
        const plus: i64 = qDivRound(Q, tau);
        if (self.config.collision == .bgk) return .{ .plus = plus, .minus = plus };
        const magic = std.math.clamp(self.config.magic_q16, 1, MAGIC_MAX);
        const tau_minus = @as(i64, qDivRound(magic, tau - Q / 2)) + Q / 2;
        return .{ .plus = plus, .minus = qDivRound(Q, sat32(tau_minus)) };
    }

    fn collide(self: *Lattice) void {
        const rates = self.relaxationRates();
        const sum_w = rates.plus + rates.minus;
        const diff_w = rates.plus - rates.minus;
        const fx: i64 = self.config.force_x_q16;
        const fy: i64 = self.config.force_y_q16;
        const forced = fx != 0 or fy != 0;
        const hx = divRound(fx, 2);
        const hy = divRound(fy, 2);
        const k_plus = 2 * Q - rates.plus; // (1 - omega+/2) * 2Q
        const k_minus = 2 * Q - rates.minus;

        var st = self.stats;
        st.mass = 0;
        st.momentum_x = 0;
        st.momentum_y = 0;
        st.fluid_cells = 0;
        st.max_speed_sq = 0;
        st.saturated = 0;
        st.nonpositive_density = 0;

        for (0..self.cellCount()) |idx| {
            if (self.solid[idx]) continue;
            const g = self.f_next[idx];
            const m = momentsOf(g);
            const ux = std.math.clamp(m.jx + hx, -U_CLAMP, U_CLAMP);
            const uy = std.math.clamp(m.jy + hy, -U_CLAMP, U_CLAMP);

            var eq: [9]i64 = undefined;
            equilibrium(m.rho, ux, uy, &eq);
            var neq: [9]i64 = undefined;
            inline for (0..9) |i| neq[i] = @as(i64, g[i]) - eq[i];

            var post: [9]i64 = undefined;
            var sum: i64 = 0;
            var mx: i64 = 0;
            var my: i64 = 0;
            const u_dot_f = ux * fx + uy * fy; // Q32
            inline for (1..9) |i| {
                const relax = divRound(sum_w * neq[i] + diff_w * neq[OPPOSITE[i]], 2 * Q);
                var v = @as(i64, g[i]) - relax;
                if (forced) {
                    const cu = @as(i64, CX[i]) * ux + @as(i64, CY[i]) * uy;
                    const cf = @as(i64, CX[i]) * fx + @as(i64, CY[i]) * fy;
                    const even = k_plus * (9 * cu * cf - 3 * u_dot_f); // Q48
                    const odd = k_minus * ((3 * cf) << 16); // Q48
                    v += divRound(W36[i] * (even + odd), 72 << 32);
                }
                post[i] = v;
                sum += v;
                mx += CX[i] * v;
                my += CY[i] * v;
            }
            fixMoments(&post, m.rho, m.jx + fx - mx, m.jy + fy - my, sum);

            const out = &self.f[idx];
            inline for (0..9) |i| {
                const s = sat32(post[i]);
                st.saturated += @intFromBool(s != post[i]);
                out[i] = s;
            }

            st.mass += m.rho;
            st.momentum_x += m.jx + fx;
            st.momentum_y += m.jy + fy;
            st.fluid_cells += 1;
            st.nonpositive_density += @intFromBool(m.rho <= 0);
            const sp = ux * ux + uy * uy;
            if (sp > st.max_speed_sq) st.max_speed_sq = sp;
        }
        self.stats = st;
    }

    /// Recomputes mass and momentum totals from the current distributions.
    pub fn totals(self: *const Lattice) Moments {
        var t = Moments{ .rho = 0, .jx = 0, .jy = 0 };
        for (0..self.cellCount()) |idx| {
            if (self.solid[idx]) continue;
            const m = momentsOf(self.f[idx]);
            t.rho += m.rho;
            t.jx += m.jx;
            t.jy += m.jy;
        }
        return t;
    }

    /// FNV-1a 64 over the distributions (little-endian) and the solid mask.
    pub fn digest(self: *const Lattice) u64 {
        var h = std.hash.Fnv1a_64.init();
        var dims: [16]u8 = undefined;
        std.mem.writeInt(u64, dims[0..8], self.width, .little);
        std.mem.writeInt(u64, dims[8..16], self.height, .little);
        h.update(&dims);
        for (self.f) |cell| {
            var buf: [36]u8 = undefined;
            for (0..9) |i| std.mem.writeInt(i32, buf[i * 4 ..][0..4], cell[i], .little);
            h.update(&buf);
        }
        for (self.solid) |s| h.update(&[_]u8{@intFromBool(s)});
        return h.final();
    }

    pub fn digest32(self: *const Lattice) u32 {
        const d = self.digest();
        return @truncate(d ^ (d >> 32));
    }

    /// Vorticity du_y/dx - du_x/dy at a fluid cell from central differences of
    /// the momentum (Q16.16, times two). Solid or out-of-domain neighbours use
    /// the centre value.
    pub fn vorticity2(self: *const Lattice, x: usize, y: usize) i64 {
        const c = momentsOf(self.f[y * self.width + x]);
        const xl = self.neighbourMoments(x, y, -1, 0) orelse c;
        const xr = self.neighbourMoments(x, y, 1, 0) orelse c;
        const yd = self.neighbourMoments(x, y, 0, -1) orelse c;
        const yu = self.neighbourMoments(x, y, 0, 1) orelse c;
        return (xr.jy - xl.jy) - (yu.jx - yd.jx);
    }

    fn neighbourMoments(self: *const Lattice, x: usize, y: usize, dx: isize, dy: isize) ?Moments {
        const nx = @as(isize, @intCast(x)) + dx;
        const ny = @as(isize, @intCast(y)) + dy;
        if (nx < 0 or ny < 0 or nx >= @as(isize, @intCast(self.width)) or ny >= @as(isize, @intCast(self.height))) return null;
        const idx: usize = @intCast(ny * @as(isize, @intCast(self.width)) + nx);
        if (self.solid[idx]) return null;
        return momentsOf(self.f[idx]);
    }

    /// Renders the lattice into an RGBA buffer of width * height * 4 bytes.
    pub fn render(self: *const Lattice, rgba: []u8, mode: RenderMode) void {
        std.debug.assert(rgba.len >= self.cellCount() * 4);
        for (0..self.height) |y| {
            for (0..self.width) |x| {
                const idx = y * self.width + x;
                const px = rgba[idx * 4 ..][0..4];
                if (self.solid[idx]) {
                    // Obstacle: cybernetic dark gray/slate
                    px.* = .{ 30, 35, 45, 255 };
                    continue;
                }
                switch (mode) {
                    .speed => {
                        const m = momentsOf(self.f[idx]);
                        const speed: u64 = @abs(m.jx) + @abs(m.jy);
                        const val: u64 = @min(255, speed / 32);
                        // Cybernetic neon fluid color map (cyan/blue/magenta)
                        px.* = .{ @intCast(@min(255, val * 2)), @intCast(@min(255, 100 + val)), 220, 255 };
                    },
                    .vorticity => {
                        const v = std.math.clamp(@divTrunc(self.vorticity2(x, y), 4), -215, 215);
                        const a: u8 = @intCast(@abs(v));
                        px.* = if (v >= 0) .{ 20 + a, 24, 40, 255 } else .{ 20, 24, 40 + a, 255 };
                    },
                    .density => {
                        const m = momentsOf(self.f[idx]);
                        const g: u8 = @intCast(std.math.clamp(128 + @divTrunc(m.rho - Q, 8), 0, 255));
                        px.* = .{ g, g, g, 255 };
                    },
                }
            }
        }
    }
};

/// The fixed 160 x 80 lattice with static storage used by the WebAssembly build.
pub const Grid = struct {
    f: [NUM_CELLS][9]i32,
    f_next: [NUM_CELLS][9]i32,
    solid: [NUM_CELLS]bool,
    rgba: [NUM_CELLS * 4]u8,
    tau_q16: i32, // Relaxation time, default 1.2
    config: Config,
    stats: Stats,
    steps: u64,
    render_mode: RenderMode,

    /// Default solver settings, fluid flowing at the inlet velocity, and the
    /// default cylinder obstacle.
    pub fn init() Grid {
        var grid: Grid = undefined;
        grid.initInPlace();
        return grid;
    }

    /// Same as `init` without returning the large struct by value.
    pub fn initInPlace(self: *Grid) void {
        self.tau_q16 = TAU_DEFAULT;
        self.config = .{};
        self.render_mode = .speed;
        self.reset();
    }

    /// Restarts the flow and restores the default cylinder at (40, 40),
    /// radius 10. Solver settings (tau, collision, boundaries, forcing) are kept.
    pub fn reset(self: *Grid) void {
        @memset(&self.solid, false);
        @memset(&self.rgba, 0);
        var l = self.view();
        l.resetFlow();
        self.sync(l);
        self.setCircleObstacle(40, 40, 10, true);
    }

    pub fn view(self: *Grid) Lattice {
        return .{
            .width = WIDTH,
            .height = HEIGHT,
            .f = &self.f,
            .f_next = &self.f_next,
            .solid = &self.solid,
            .tau_q16 = self.tau_q16,
            .config = self.config,
            .stats = self.stats,
            .steps = self.steps,
        };
    }

    fn sync(self: *Grid, l: Lattice) void {
        self.tau_q16 = l.tau_q16;
        self.config = l.config;
        self.stats = l.stats;
        self.steps = l.steps;
    }

    pub fn setCircleObstacle(self: *Grid, cx: i32, cy: i32, radius: i32, is_solid: bool) void {
        var l = self.view();
        l.setCircleObstacle(cx, cy, radius, is_solid);
    }

    pub fn setObstacleCell(self: *Grid, x: usize, y: usize, is_solid: bool) void {
        var l = self.view();
        l.setObstacleCell(x, y, is_solid);
    }

    pub fn clearObstacles(self: *Grid) void {
        var l = self.view();
        l.clearObstacles();
    }

    pub fn setTau(self: *Grid, tau_q16: i32) ConfigError!void {
        try validateTau(tau_q16);
        self.tau_q16 = tau_q16;
    }

    pub fn setConfig(self: *Grid, config: Config) ConfigError!void {
        try config.validate();
        self.config = config;
    }

    pub fn step(self: *Grid) void {
        var l = self.view();
        l.step();
        self.sync(l);
    }

    /// Renders current fluid state into RGBA linear buffer for Canvas ImageData.
    pub fn renderToRgba(self: *Grid) void {
        var l = self.view();
        l.render(&self.rgba, self.render_mode);
    }

    pub fn digest32(self: *Grid) u32 {
        var l = self.view();
        return l.digest32();
    }
};
