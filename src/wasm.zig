// Copyright (c) 2026 jk9576. All rights reserved.
// Halite: Deterministic Lattice Boltzmann Fluid / Gas Automaton (D2Q9)
// Licensed under the MIT License.
//
// Freestanding WebAssembly exports. The module imports nothing. Every value
// crossing the boundary is an integer or a Q16.16 integer, except the
// telemetry getters, which convert to f64 only on the way out.

pub const lbm = @import("lbm.zig");

/// Encoded as major << 16 | minor << 8 | patch.
pub const VERSION: u32 = 0x000200;

var global_grid: lbm.Grid = undefined;
var initialized: bool = false;
/// Scratch obstacle mask written by the host: one byte per cell, nonzero = solid.
var mask_buffer: [lbm.NUM_CELLS]u8 = @splat(0);

fn grid() *lbm.Grid {
    if (!initialized) halite_init();
    return &global_grid;
}

pub export fn halite_init() void {
    global_grid.initInPlace();
    initialized = true;
}

/// Restarts the flow and restores the default cylinder; keeps solver settings.
pub export fn halite_reset() void {
    _ = grid();
    global_grid.reset();
}

pub export fn halite_step() void {
    grid().step();
}

pub export fn halite_render() void {
    grid().renderToRgba();
}

pub export fn halite_set_obstacle(x: u32, y: u32, is_solid: u32) void {
    grid().setObstacleCell(x, y, is_solid != 0);
}

pub export fn halite_get_rgba_ptr() [*]u8 {
    return &grid().rgba;
}

pub export fn halite_get_width() u32 {
    return lbm.WIDTH;
}

pub export fn halite_get_height() u32 {
    return lbm.HEIGHT;
}

// ---- 0.2.0 additions ----------------------------------------------------

pub export fn halite_get_version() u32 {
    return VERSION;
}

/// Runs `count` steps in one call (capped at 100000). Returns 1 while stable.
pub export fn halite_step_n(count: u32) u32 {
    const g = grid();
    var i: u32 = 0;
    while (i < @min(count, 100_000)) : (i += 1) g.step();
    return @intFromBool(g.stats.stable());
}

/// 0 = BGK, 1 = TRT. Returns 1 if accepted.
pub export fn halite_set_collision(op: u32) u32 {
    const g = grid();
    var cfg = g.config;
    cfg.collision = switch (op) {
        0 => .bgk,
        1 => .trt,
        else => return 0,
    };
    g.setConfig(cfg) catch return 0;
    return 1;
}

/// TRT magic parameter in Q16.16, in (0, 4]. 12288 (3/16) is the default.
pub export fn halite_set_trt_magic(magic_q16: i32) u32 {
    const g = grid();
    var cfg = g.config;
    cfg.magic_q16 = magic_q16;
    g.setConfig(cfg) catch return 0;
    return 1;
}

/// Relaxation time in Q16.16, accepted range [33096, 1048576] (about 0.505 to 16).
pub export fn halite_set_tau(tau_q16: i32) u32 {
    grid().setTau(tau_q16) catch return 0;
    return 1;
}

pub export fn halite_get_tau() i32 {
    return grid().tau_q16;
}

/// x_mode: 0 channel (Zou-He inlet/outlet), 1 periodic, 2 closed walls.
/// y_mode: 0 no-slip walls, 1 periodic. Returns 1 if accepted.
pub export fn halite_set_boundary(x_mode: u32, y_mode: u32) u32 {
    const g = grid();
    var cfg = g.config;
    cfg.x_boundary = switch (x_mode) {
        0 => .channel,
        1 => .periodic,
        2 => .closed,
        else => return 0,
    };
    cfg.y_boundary = switch (y_mode) {
        0 => .wall,
        1 => .periodic,
        else => return 0,
    };
    g.setConfig(cfg) catch return 0;
    return 1;
}

/// Inlet velocity in Q16.16 (|u| <= 0.3) and profile (0 uniform, 1 parabolic).
pub export fn halite_set_inflow(ux_q16: i32, profile: u32) u32 {
    const g = grid();
    var cfg = g.config;
    cfg.inflow_q16 = ux_q16;
    cfg.inlet_profile = switch (profile) {
        0 => .uniform,
        1 => .parabolic,
        else => return 0,
    };
    g.setConfig(cfg) catch return 0;
    return 1;
}

/// Outlet density in Q16.16, accepted range [0.5, 2].
pub export fn halite_set_outlet_density(rho_q16: i32) u32 {
    const g = grid();
    var cfg = g.config;
    cfg.outlet_density_q16 = rho_q16;
    g.setConfig(cfg) catch return 0;
    return 1;
}

/// Uniform body force per cell and step in Q16.16 (|F| <= 1024 each).
pub export fn halite_set_body_force(fx_q16: i32, fy_q16: i32) u32 {
    const g = grid();
    var cfg = g.config;
    cfg.force_x_q16 = fx_q16;
    cfg.force_y_q16 = fy_q16;
    g.setConfig(cfg) catch return 0;
    return 1;
}

pub export fn halite_set_circle(cx: i32, cy: i32, radius: i32, is_solid: u32) void {
    grid().setCircleObstacle(cx, cy, radius, is_solid != 0);
}

pub export fn halite_clear_obstacles() void {
    grid().clearObstacles();
}

/// Host-writable mask of width * height bytes; see halite_apply_mask.
pub export fn halite_get_mask_ptr() [*]u8 {
    return &mask_buffer;
}

/// Replaces all obstacles with the mask (nonzero byte = solid). Returns the
/// number of solid cells.
pub export fn halite_apply_mask() u32 {
    const g = grid();
    var n: u32 = 0;
    for (0..lbm.HEIGHT) |y| {
        for (0..lbm.WIDTH) |x| {
            const s = mask_buffer[y * lbm.WIDTH + x] != 0;
            g.setObstacleCell(x, y, s);
            n += @intFromBool(s);
        }
    }
    return n;
}

/// 0 speed (default), 1 vorticity, 2 density. Returns 1 if accepted.
pub export fn halite_set_render_mode(mode: u32) u32 {
    grid().render_mode = switch (mode) {
        0 => .speed,
        1 => .vorticity,
        2 => .density,
        else => return 0,
    };
    return 1;
}

/// Momentum-exchange force on the obstacles during the last step (lattice units).
pub export fn halite_get_force_x() f64 {
    return lbm.fromQ(grid().stats.obstacle_force_x);
}

pub export fn halite_get_force_y() f64 {
    return lbm.fromQ(grid().stats.obstacle_force_y);
}

/// Mean fluid density after the last step.
pub export fn halite_get_mean_density() f64 {
    const s = grid().stats;
    if (s.fluid_cells == 0) return 0;
    return lbm.fromQ(s.mass) / @as(f64, @floatFromInt(s.fluid_cells));
}

pub export fn halite_get_max_speed() f64 {
    return grid().stats.maxSpeed();
}

/// 1 if the last step had no clamped values, no non-positive density and |u| < 0.5.
pub export fn halite_is_stable() u32 {
    return @intFromBool(grid().stats.stable());
}

/// Steps since the last init or reset (wraps at 2^32).
pub export fn halite_get_step_count() u32 {
    return @truncate(grid().steps);
}

/// 32-bit digest of the full lattice state, identical to `halite run --json`.
pub export fn halite_get_digest() u32 {
    return grid().digest32();
}
