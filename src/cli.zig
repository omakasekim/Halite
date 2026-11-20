// Copyright (c) 2026 jk9576. All rights reserved.
// Halite: Deterministic Lattice Boltzmann Fluid / Gas Automaton (D2Q9)
// Licensed under the MIT License.
//
// Command-line front end: run (with JSON telemetry and field export),
// bench (MLUPS) and validate (analytic and integer-exact checks).

const std = @import("std");
const lbm = @import("lbm.zig");
const fieldio = @import("fieldio.zig");
const validate = @import("validate.zig");

pub const VERSION = "0.2.0";

pub const usage =
    \\usage: halite [command] [options]
    \\
    \\commands:
    \\  (none)      run the 50-step demo on the default 160x80 lattice
    \\  run         simulate and report telemetry
    \\  bench       measure throughput in MLUPS for BGK and TRT
    \\  validate    run the Poiseuille, conservation, momentum-balance,
    \\              mirror-symmetry and channel-flux checks (exit 1 on failure)
    \\  version     print the version
    \\  help        print this text
    \\
    \\simulation options (run, bench):
    \\  --width N --height N        lattice size (default 160 x 80)
    \\  --steps N                   steps to run (run: 1000, bench: 500)
    \\  --tau X                     relaxation time, 0.505..16 (default 1.2)
    \\  --collision bgk|trt         collision operator (default bgk)
    \\  --magic X                   TRT magic parameter, 0..4 (default 0.1875)
    \\  --x-boundary channel|periodic|closed   (default channel)
    \\  --y-boundary wall|periodic  (default wall)
    \\  --inflow X                  inlet velocity, |X| <= 0.3 (default 0.08)
    \\  --profile uniform|parabolic inlet profile (default uniform)
    \\  --outlet-density X          outlet density, 0.5..2 (default 1)
    \\  --force-x X --force-y X     body force per cell and step, |F| <= 1/64
    \\  --cylinder CX,CY,R          add a cylinder (repeatable)
    \\  --obstacle FILE.pgm         obstacle mask (dark = solid), resampled
    \\  --invert                    treat light pixels of the mask as solid
    \\  --no-obstacle               no default cylinder at (W/4, H/2, H/8)
    \\
    \\run output options:
    \\  --json                      print a JSON summary
    \\  --every N                   with --json, print a JSON line every N steps
    \\  --pgm FILE --field F        write an 8-bit PGM of speed|vorticity|density|solid
    \\  --vtk FILE                  write a legacy VTK file (density, velocity, vorticity, solid)
    \\
    \\validate options:
    \\  --json                      print results as JSON
    \\
;

pub const CliError = error{
    UnknownCommand,
    UnknownOption,
    MissingValue,
    BadValue,
    TooManyCylinders,
};

const Cylinder = struct { cx: i32, cy: i32, r: i32 };

pub const RunOptions = struct {
    width: usize = lbm.WIDTH,
    height: usize = lbm.HEIGHT,
    steps: ?u64 = null,
    tau_q16: i32 = lbm.TAU_DEFAULT,
    config: lbm.Config = .{},
    cylinders: [16]Cylinder = undefined,
    cylinder_count: usize = 0,
    obstacle_path: ?[]const u8 = null,
    invert: bool = false,
    no_obstacle: bool = false,
    json: bool = false,
    every: u64 = 0,
    pgm_path: ?[]const u8 = null,
    field: fieldio.Field = .speed,
    vtk_path: ?[]const u8 = null,
};

fn parseUint(s: []const u8, max: u64) CliError!u64 {
    const v = std.fmt.parseInt(u64, s, 10) catch return error.BadValue;
    if (v > max) return error.BadValue;
    return v;
}

fn parseReal(s: []const u8) CliError!f64 {
    const v = std.fmt.parseFloat(f64, s) catch return error.BadValue;
    if (!std.math.isFinite(v)) return error.BadValue;
    return v;
}

fn parseQ(s: []const u8) CliError!i32 {
    const v = try parseReal(s);
    if (@abs(v) > 30000.0) return error.BadValue;
    return lbm.toQ(v);
}

fn parseCylinder(s: []const u8) CliError!Cylinder {
    var it = std.mem.splitScalar(u8, s, ',');
    var vals: [3]i32 = undefined;
    for (&vals) |*v| {
        const part = it.next() orelse return error.BadValue;
        v.* = std.fmt.parseInt(i32, std.mem.trim(u8, part, " "), 10) catch return error.BadValue;
    }
    if (it.next() != null or vals[2] < 0 or vals[2] > lbm.MAX_DIM) return error.BadValue;
    if (@abs(vals[0]) > 1 << 20 or @abs(vals[1]) > 1 << 20) return error.BadValue;
    return .{ .cx = vals[0], .cy = vals[1], .r = vals[2] };
}

/// Parses options after the command name. `err_arg` receives the offending argument.
pub fn parseRunOptions(args: []const []const u8, err_arg: *[]const u8) CliError!RunOptions {
    var o = RunOptions{};
    var i: usize = 0;
    while (i < args.len) : (i += 1) {
        const a = args[i];
        err_arg.* = a;
        if (std.mem.eql(u8, a, "--json")) {
            o.json = true;
            continue;
        }
        if (std.mem.eql(u8, a, "--invert")) {
            o.invert = true;
            continue;
        }
        if (std.mem.eql(u8, a, "--no-obstacle")) {
            o.no_obstacle = true;
            continue;
        }
        const known = [_][]const u8{
            "--width",      "--height",     "--steps",    "--tau",     "--collision",      "--magic",
            "--x-boundary", "--y-boundary", "--inflow",   "--profile", "--outlet-density", "--force-x",
            "--force-y",    "--cylinder",   "--obstacle", "--every",   "--pgm",            "--field",
            "--vtk",
        };
        var is_known = false;
        for (known) |k| is_known = is_known or std.mem.eql(u8, a, k);
        if (!is_known) return error.UnknownOption;
        if (i + 1 >= args.len) return error.MissingValue;
        i += 1;
        const v = args[i];
        err_arg.* = v;
        if (std.mem.eql(u8, a, "--width")) {
            o.width = @intCast(try parseUint(v, lbm.MAX_DIM));
        } else if (std.mem.eql(u8, a, "--height")) {
            o.height = @intCast(try parseUint(v, lbm.MAX_DIM));
        } else if (std.mem.eql(u8, a, "--steps")) {
            o.steps = try parseUint(v, 1 << 40);
        } else if (std.mem.eql(u8, a, "--tau")) {
            o.tau_q16 = try parseQ(v);
        } else if (std.mem.eql(u8, a, "--collision")) {
            o.config.collision = std.meta.stringToEnum(lbm.Collision, v) orelse return error.BadValue;
        } else if (std.mem.eql(u8, a, "--magic")) {
            o.config.magic_q16 = try parseQ(v);
        } else if (std.mem.eql(u8, a, "--x-boundary")) {
            o.config.x_boundary = std.meta.stringToEnum(lbm.XBoundary, v) orelse return error.BadValue;
        } else if (std.mem.eql(u8, a, "--y-boundary")) {
            o.config.y_boundary = std.meta.stringToEnum(lbm.YBoundary, v) orelse return error.BadValue;
        } else if (std.mem.eql(u8, a, "--inflow")) {
            o.config.inflow_q16 = try parseQ(v);
        } else if (std.mem.eql(u8, a, "--profile")) {
            o.config.inlet_profile = std.meta.stringToEnum(lbm.InletProfile, v) orelse return error.BadValue;
        } else if (std.mem.eql(u8, a, "--outlet-density")) {
            o.config.outlet_density_q16 = try parseQ(v);
        } else if (std.mem.eql(u8, a, "--force-x")) {
            o.config.force_x_q16 = try parseQ(v);
        } else if (std.mem.eql(u8, a, "--force-y")) {
            o.config.force_y_q16 = try parseQ(v);
        } else if (std.mem.eql(u8, a, "--cylinder")) {
            if (o.cylinder_count == o.cylinders.len) return error.TooManyCylinders;
            o.cylinders[o.cylinder_count] = try parseCylinder(v);
            o.cylinder_count += 1;
        } else if (std.mem.eql(u8, a, "--obstacle")) {
            o.obstacle_path = v;
        } else if (std.mem.eql(u8, a, "--every")) {
            o.every = try parseUint(v, 1 << 40);
        } else if (std.mem.eql(u8, a, "--pgm")) {
            o.pgm_path = v;
        } else if (std.mem.eql(u8, a, "--field")) {
            o.field = fieldio.Field.parse(v) orelse return error.BadValue;
        } else if (std.mem.eql(u8, a, "--vtk")) {
            o.vtk_path = v;
        }
    }
    return o;
}

/// Builds the lattice described by the options. The default obstacle is the
/// cylinder at (W/4, H/2, H/8), which for 160 x 80 is the WebAssembly default.
pub fn setupLattice(gpa: std.mem.Allocator, o: RunOptions, mask: ?fieldio.Image) !lbm.Lattice {
    var l = try lbm.Lattice.create(gpa, o.width, o.height);
    errdefer l.destroy(gpa);
    try l.setTau(o.tau_q16);
    try l.setConfig(o.config);
    l.resetFlow();
    if (mask) |img| _ = fieldio.applyMask(&l, img, o.invert);
    for (o.cylinders[0..o.cylinder_count]) |c| l.setCircleObstacle(c.cx, c.cy, c.r, true);
    if (mask == null and o.cylinder_count == 0 and !o.no_obstacle) {
        l.setCircleObstacle(@intCast(o.width / 4), @intCast(o.height / 2), @intCast(o.height / 8), true);
    }
    return l;
}

/// Reference speed for force coefficients: the inlet speed, or the mean
/// speed 2/3 of the peak for a parabolic inlet between walls.
pub fn referenceSpeed(l: *const lbm.Lattice) f64 {
    const u = @abs(lbm.fromQ(l.config.inflow_q16));
    const parabolic = l.config.inlet_profile == .parabolic and l.config.y_boundary == .wall;
    return if (parabolic) u * 2.0 / 3.0 else u;
}

fn writeOptionalFloat(w: *std.Io.Writer, v: ?f64) !void {
    if (v) |x| {
        if (std.math.isFinite(x)) return w.print("{d}", .{x});
    }
    try w.writeAll("null");
}

pub fn writeTelemetry(w: *std.Io.Writer, l: *const lbm.Lattice) !void {
    const s = l.stats;
    const n: f64 = @floatFromInt(@max(s.fluid_cells, 1));
    const u_ref = referenceSpeed(l);
    const d: f64 = @floatFromInt(l.obstacleExtent().rows);
    const fx = lbm.fromQ(s.obstacle_force_x);
    const fy = lbm.fromQ(s.obstacle_force_y);
    const coef: ?f64 = if (u_ref > 0 and d > 0) 2.0 / (u_ref * u_ref * d) else null;
    try w.print("{{\"step\":{d},\"fluid_cells\":{d},\"mass\":{d},\"mean_density\":{d},\"momentum_x\":{d},\"momentum_y\":{d},\"max_speed\":{d},", .{
        l.steps,
        s.fluid_cells,
        s.mass,
        lbm.fromQ(s.mass) / n,
        s.momentum_x,
        s.momentum_y,
        s.maxSpeed(),
    });
    try w.print("\"drag\":{d},\"lift\":{d},\"wall_force_x\":{d},\"wall_force_y\":{d},\"drag_coefficient\":", .{
        fx,
        fy,
        lbm.fromQ(s.wall_force_x),
        lbm.fromQ(s.wall_force_y),
    });
    try writeOptionalFloat(w, if (coef) |c| c * fx else null);
    try w.writeAll(",\"lift_coefficient\":");
    try writeOptionalFloat(w, if (coef) |c| c * fy else null);
    try w.print(",\"stable\":{},\"digest\":\"{x:0>16}\",\"digest32\":{d}}}", .{ s.stable(), l.digest(), l.digest32() });
}

fn now(io: std.Io) i96 {
    return std.Io.Timestamp.now(io, .awake).toNanoseconds();
}

fn writeFile(io: std.Io, path: []const u8, l: *const lbm.Lattice, kind: enum { pgm, vtk }, field: fieldio.Field) !void {
    const file = try std.Io.Dir.cwd().createFile(io, path, .{});
    defer file.close(io);
    var buf: [64 * 1024]u8 = undefined;
    var fw = file.writer(io, &buf);
    switch (kind) {
        .pgm => try fieldio.writePgm(&fw.interface, l, field),
        .vtk => try fieldio.writeVtk(&fw.interface, l),
    }
    try fw.interface.flush();
}

fn loadMask(io: std.Io, gpa: std.mem.Allocator, path: []const u8) !fieldio.Image {
    const bytes = try std.Io.Dir.cwd().readFileAlloc(io, path, gpa, .limited(fieldio.MAX_IMAGE_PIXELS * 2 + 4096));
    defer gpa.free(bytes);
    return fieldio.parsePgm(gpa, bytes);
}

pub fn cmdRun(io: std.Io, gpa: std.mem.Allocator, w: *std.Io.Writer, o: RunOptions) !void {
    var mask: ?fieldio.Image = if (o.obstacle_path) |p| try loadMask(io, gpa, p) else null;
    defer if (mask) |*m| m.deinit(gpa);
    var l = try setupLattice(gpa, o, mask);
    defer l.destroy(gpa);

    const steps = o.steps orelse 1000;
    const t0 = now(io);
    var s: u64 = 0;
    while (s < steps) : (s += 1) {
        l.step();
        if (o.json and o.every > 0 and l.steps % o.every == 0) {
            try writeTelemetry(w, &l);
            try w.writeByte('\n');
            try w.flush();
        }
    }
    const elapsed_ns = now(io) - t0;
    const secs = @as(f64, @floatFromInt(elapsed_ns)) / 1e9;
    const mlups: ?f64 = if (secs > 0) @as(f64, @floatFromInt(l.stats.fluid_cells * steps)) / secs / 1e6 else null;

    if (o.pgm_path) |p| try writeFile(io, p, &l, .pgm, o.field);
    if (o.vtk_path) |p| try writeFile(io, p, &l, .vtk, o.field);

    if (o.json) {
        try w.print("{{\"version\":\"{s}\",\"width\":{d},\"height\":{d},\"steps\":{d},\"collision\":\"{s}\",\"tau\":{d},\"magic\":{d},", .{
            VERSION,                       l.width,                      l.height,
            steps,                         @tagName(l.config.collision), lbm.fromQ(l.tau_q16),
            lbm.fromQ(l.config.magic_q16),
        });
        try w.print("\"x_boundary\":\"{s}\",\"y_boundary\":\"{s}\",\"inflow\":{d},\"profile\":\"{s}\",\"force_x\":{d},\"force_y\":{d},\"solid_cells\":{d},\"elapsed_ms\":{d:.3},\"mlups\":", .{
            @tagName(l.config.x_boundary),    @tagName(l.config.y_boundary),   lbm.fromQ(l.config.inflow_q16),
            @tagName(l.config.inlet_profile), lbm.fromQ(l.config.force_x_q16), lbm.fromQ(l.config.force_y_q16),
            l.solidCount(),                   secs * 1000.0,
        });
        try writeOptionalFloat(w, mlups);
        try w.writeAll(",\"telemetry\":");
        try writeTelemetry(w, &l);
        try w.writeAll("}\n");
    } else {
        const st = l.stats;
        try w.print("halite {s}: {d} steps on {d}x{d} ({s}, tau {d:.4}, x {s}, y {s})\n", .{
            VERSION,                      steps,                l.width,                       l.height,
            @tagName(l.config.collision), lbm.fromQ(l.tau_q16), @tagName(l.config.x_boundary), @tagName(l.config.y_boundary),
        });
        try w.print("  mean density {d:.6}  max speed {d:.6}  stable {}\n", .{
            lbm.fromQ(st.mass) / @as(f64, @floatFromInt(@max(st.fluid_cells, 1))),
            st.maxSpeed(),
            st.stable(),
        });
        try w.print("  obstacle force ({d:.6}, {d:.6})  solid cells {d}\n", .{ lbm.fromQ(st.obstacle_force_x), lbm.fromQ(st.obstacle_force_y), l.solidCount() });
        try w.print("  digest {x:0>16}  ({d:.3} ms", .{ l.digest(), secs * 1000.0 });
        if (mlups) |m| try w.print(", {d:.2} MLUPS", .{m});
        try w.writeAll(")\n");
    }
}

pub fn cmdBench(io: std.Io, gpa: std.mem.Allocator, w: *std.Io.Writer, base: RunOptions) !void {
    const steps = base.steps orelse 500;
    var mask: ?fieldio.Image = if (base.obstacle_path) |p| try loadMask(io, gpa, p) else null;
    defer if (mask) |*m| m.deinit(gpa);
    const ops = [_]lbm.Collision{ .bgk, .trt };
    if (base.json) try w.print("{{\"version\":\"{s}\",\"width\":{d},\"height\":{d},\"steps\":{d},\"results\":[", .{ VERSION, base.width, base.height, steps });
    for (ops, 0..) |op, k| {
        var o = base;
        o.config.collision = op;
        var l = try setupLattice(gpa, o, mask);
        defer l.destroy(gpa);
        for (0..@min(steps / 10 + 1, 50)) |_| l.step(); // warm-up
        const t0 = now(io);
        for (0..steps) |_| l.step();
        const ns = now(io) - t0;
        const secs = @as(f64, @floatFromInt(@max(ns, 1))) / 1e9;
        const lups = @as(f64, @floatFromInt(l.stats.fluid_cells * steps)) / secs;
        if (base.json) {
            if (k > 0) try w.writeByte(',');
            try w.print("{{\"collision\":\"{s}\",\"fluid_cells\":{d},\"seconds\":{d:.6},\"mlups\":{d:.3},\"stable\":{}}}", .{ @tagName(op), l.stats.fluid_cells, secs, lups / 1e6, l.stats.stable() });
        } else {
            try w.print("{s}: {d} steps x {d} fluid cells in {d:.3} s = {d:.2} MLUPS\n", .{ @tagName(op), steps, l.stats.fluid_cells, secs, lups / 1e6 });
        }
    }
    if (base.json) try w.writeAll("]}\n");
}

pub const Check = struct {
    name: []const u8,
    passed: bool,
    detail: []const u8,
};

/// Runs every validation scenario. Returns true if all passed.
pub fn cmdValidate(gpa: std.mem.Allocator, w: *std.Io.Writer, json: bool) !bool {
    var arena_state = std.heap.ArenaAllocator.init(gpa);
    defer arena_state.deinit();
    const arena = arena_state.allocator();
    var checks: std.ArrayList(Check) = .empty;

    const tau2 = 2 * lbm.Q;
    const h = 24;
    const pb = try validate.poiseuille(gpa, .bgk, tau2, h, 6000);
    const pt = try validate.poiseuille(gpa, .trt, tau2, h, 6000);
    const pt1 = try validate.poiseuille(gpa, .trt, lbm.toQ(0.8), h, 20000);
    try checks.append(arena, .{
        .name = "poiseuille_trt_tau2",
        .passed = pt.l2_error < 0.01,
        .detail = try std.fmt.allocPrint(arena, "L2 error {e:.3} (limit 1e-2), u_max {d:.6} vs {d:.6}", .{ pt.l2_error, pt.u_max, pt.u_max_analytic }),
    });
    try checks.append(arena, .{
        .name = "poiseuille_trt_tau0.8",
        .passed = pt1.l2_error < 0.01,
        .detail = try std.fmt.allocPrint(arena, "L2 error {e:.3} (limit 1e-2), u_max {d:.6} vs {d:.6}", .{ pt1.l2_error, pt1.u_max, pt1.u_max_analytic }),
    });
    try checks.append(arena, .{
        .name = "poiseuille_trt_beats_bgk_tau2",
        .passed = pt.l2_error * 4 < pb.l2_error,
        .detail = try std.fmt.allocPrint(arena, "TRT {e:.3} vs BGK {e:.3}", .{ pt.l2_error, pb.l2_error }),
    });
    for ([_]lbm.Collision{ .bgk, .trt }) |op| {
        const c = try validate.conservation(gpa, op, 2000);
        try checks.append(arena, .{
            .name = try std.fmt.allocPrint(arena, "exact_conservation_{s}", .{@tagName(op)}),
            .passed = c.mass_drift == 0 and c.momentum_x_drift == 0 and c.momentum_y_drift == 0,
            .detail = try std.fmt.allocPrint(arena, "mass drift {d}, momentum drift ({d}, {d}) over {d} steps", .{ c.mass_drift, c.momentum_x_drift, c.momentum_y_drift, c.steps }),
        });
        const b = try validate.momentumBalance(gpa, op, 1500);
        try checks.append(arena, .{
            .name = try std.fmt.allocPrint(arena, "momentum_exchange_balance_{s}", .{@tagName(op)}),
            .passed = b.residual_x == 0 and b.residual_y == 0,
            .detail = try std.fmt.allocPrint(arena, "residual ({d}, {d}); drag {d} of body force {d}", .{ b.residual_x, b.residual_y, b.final_drag, b.body_force_total }),
        });
        const m = try validate.mirrorSymmetry(gpa, op, 1500);
        try checks.append(arena, .{
            .name = try std.fmt.allocPrint(arena, "mirror_symmetry_{s}", .{@tagName(op)}),
            .passed = m.max_abs_lift == 0 and m.asymmetric_cells == 0 and m.drag > 0,
            .detail = try std.fmt.allocPrint(arena, "max |lift| {d}, asymmetric cells {d}, Cd {d:.4}", .{ m.max_abs_lift, m.asymmetric_cells, m.drag_coefficient }),
        });
        const ch = try validate.channelFlux(gpa, op, 6000);
        try checks.append(arena, .{
            .name = try std.fmt.allocPrint(arena, "zou_he_channel_flux_{s}", .{@tagName(op)}),
            .passed = ch.max_flux_deviation < 0.01 and @abs(ch.outlet_density - 1.0) < 1e-3,
            .detail = try std.fmt.allocPrint(arena, "max flux deviation {e:.3} (limit 1e-2), outlet density {d:.6}", .{ ch.max_flux_deviation, ch.outlet_density }),
        });
    }

    var all = true;
    for (checks.items) |c| all = all and c.passed;
    if (json) {
        try w.print("{{\"version\":\"{s}\",\"passed\":{},\"checks\":[", .{ VERSION, all });
        for (checks.items, 0..) |c, i| {
            if (i > 0) try w.writeByte(',');
            try w.print("{{\"name\":\"{s}\",\"passed\":{},\"detail\":\"{s}\"}}", .{ c.name, c.passed, c.detail });
        }
        try w.writeAll("]}\n");
    } else {
        for (checks.items) |c| try w.print("{s} {s}: {s}\n", .{ if (c.passed) "PASS" else "FAIL", c.name, c.detail });
        try w.print("{s}\n", .{if (all) "all checks passed" else "some checks failed"});
    }
    return all;
}

/// The original 50-step demo, now printing the digest it computes.
pub fn cmdDemo(gpa: std.mem.Allocator, w: *std.Io.Writer) !void {
    var l = try setupLattice(gpa, .{}, null);
    defer l.destroy(gpa);
    try w.print("Halite: Running 50 deterministic LBM simulation steps ({d}x{d} grid)...\n", .{ l.width, l.height });
    for (0..50) |_| l.step();
    try w.print("Halite: 50 steps completed. Bit-for-bit state digest {x:0>16}.\n", .{l.digest()});
}

/// Dispatches a command line (without the program name). Returns the exit code.
pub fn dispatch(io: std.Io, gpa: std.mem.Allocator, out: *std.Io.Writer, err: *std.Io.Writer, args: []const []const u8) !u8 {
    if (args.len == 0) {
        try cmdDemo(gpa, out);
        return 0;
    }
    const cmd = args[0];
    const rest = args[1..];
    if (std.mem.eql(u8, cmd, "help") or std.mem.eql(u8, cmd, "--help") or std.mem.eql(u8, cmd, "-h")) {
        try out.writeAll(usage);
        return 0;
    }
    if (std.mem.eql(u8, cmd, "version") or std.mem.eql(u8, cmd, "--version")) {
        try out.print("halite {s}\n", .{VERSION});
        return 0;
    }
    const is_run = std.mem.eql(u8, cmd, "run");
    const is_bench = std.mem.eql(u8, cmd, "bench");
    const is_validate = std.mem.eql(u8, cmd, "validate");
    if (!is_run and !is_bench and !is_validate) {
        try err.print("halite: unknown command '{s}'\n\n{s}", .{ cmd, usage });
        return 2;
    }
    var bad: []const u8 = "";
    const o = parseRunOptions(rest, &bad) catch |e| {
        try err.print("halite: {s}: '{s}'\n", .{ switch (e) {
            error.UnknownOption => "unknown option",
            error.MissingValue => "missing value for",
            error.BadValue => "invalid value",
            error.TooManyCylinders => "at most 16 cylinders",
            error.UnknownCommand => "unknown command",
        }, bad });
        return 2;
    };
    if (is_validate) {
        return if (try cmdValidate(gpa, out, o.json)) 0 else 1;
    }
    lbm.Lattice.checkDimensions(o.width, o.height) catch {
        try err.print("halite: lattice must be between {d} and {d} cells per side and at most {d} cells\n", .{ lbm.MIN_DIM, lbm.MAX_DIM, lbm.MAX_CELLS });
        return 2;
    };
    lbm.validateTau(o.tau_q16) catch {
        try err.writeAll("halite: --tau must be between 0.505 and 16\n");
        return 2;
    };
    o.config.validate() catch |e| {
        try err.print("halite: invalid configuration: {s}\n", .{@errorName(e)});
        return 2;
    };
    if (is_run) {
        cmdRun(io, gpa, out, o) catch |e| {
            try err.print("halite: run failed: {s}\n", .{@errorName(e)});
            return 1;
        };
    } else {
        cmdBench(io, gpa, out, o) catch |e| {
            try err.print("halite: bench failed: {s}\n", .{@errorName(e)});
            return 1;
        };
    }
    return 0;
}
