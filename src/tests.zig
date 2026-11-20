// Copyright (c) 2026 jk9576. All rights reserved.
// Halite: Deterministic Lattice Boltzmann Fluid / Gas Automaton (D2Q9)
// Licensed under the MIT License.

const std = @import("std");
const lbm = @import("lbm.zig");
const fieldio = @import("fieldio.zig");
const validate = @import("validate.zig");
const cli = @import("cli.zig");

const testing = std.testing;
const gpa = testing.allocator;

// One static grid shared by the 160 x 80 tests (it is too large for the test stack twice over).
var shared_grid: lbm.Grid = undefined;

test "Deterministic LBM simulation produces consistent state" {
    var a = try lbm.Lattice.create(gpa, lbm.WIDTH, lbm.HEIGHT);
    defer a.destroy(gpa);
    var b = try lbm.Lattice.create(gpa, lbm.WIDTH, lbm.HEIGHT);
    defer b.destroy(gpa);
    for ([_]*lbm.Lattice{ &a, &b }) |l| {
        l.resetFlow();
        l.setCircleObstacle(40, 40, 10, true);
    }
    for (0..10) |_| {
        a.step();
        b.step();
    }
    const rgba_a = try gpa.alloc(u8, lbm.NUM_CELLS * 4);
    defer gpa.free(rgba_a);
    const rgba_b = try gpa.alloc(u8, lbm.NUM_CELLS * 4);
    defer gpa.free(rgba_b);
    a.render(rgba_a, .speed);
    b.render(rgba_b, .speed);
    try testing.expectEqualSlices(u8, rgba_a, rgba_b);
    try testing.expectEqual(a.digest(), b.digest());
}

test "Obstacle cell setting works" {
    shared_grid.initInPlace();
    shared_grid.setObstacleCell(10, 10, true);
    try testing.expect(shared_grid.solid[10 * lbm.WIDTH + 10]);
    shared_grid.setObstacleCell(10, 10, false);
    try testing.expect(!shared_grid.solid[10 * lbm.WIDTH + 10]);
    // Out-of-range coordinates are ignored.
    shared_grid.setObstacleCell(lbm.WIDTH, 0, true);
    shared_grid.setObstacleCell(0, lbm.HEIGHT, true);
}

test "default workbench grid stays stable where the 0.1.0 solver overflowed" {
    // 0.1.0 overflowed i32 within about 50 steps of this setup (a closed right
    // wall and a non-conservative collision); Debug builds trap on overflow.
    shared_grid.initInPlace();
    shared_grid.clearObstacles();
    shared_grid.setCircleObstacle(35, 40, 8, true);
    for (0..600) |_| shared_grid.step();
    try testing.expect(shared_grid.stats.stable());
    const mean = lbm.fromQ(shared_grid.stats.mass) / @as(f64, @floatFromInt(shared_grid.stats.fluid_cells));
    try testing.expect(@abs(mean - 1.0) < 0.05);
    try testing.expect(shared_grid.stats.obstacle_force_x > 0);
    shared_grid.renderToRgba();
    try testing.expectEqual(@as(u8, 255), shared_grid.rgba[3]);
}

test "static grid and heap lattice produce identical states" {
    shared_grid.initInPlace();
    var l = try lbm.Lattice.create(gpa, lbm.WIDTH, lbm.HEIGHT);
    defer l.destroy(gpa);
    l.resetFlow();
    l.setCircleObstacle(40, 40, 10, true);
    for (0..40) |_| {
        shared_grid.step();
        l.step();
    }
    try testing.expectEqual(l.digest32(), shared_grid.digest32());
    try testing.expectEqual(@as(u64, 40), shared_grid.steps);
}

test "reset keeps solver settings and restores the default cylinder" {
    shared_grid.initInPlace();
    try shared_grid.setTau(lbm.toQ(0.9));
    var cfg = shared_grid.config;
    cfg.collision = .trt;
    try shared_grid.setConfig(cfg);
    shared_grid.clearObstacles();
    for (0..5) |_| shared_grid.step();
    shared_grid.reset();
    try testing.expectEqual(lbm.toQ(0.9), shared_grid.tau_q16);
    try testing.expectEqual(lbm.Collision.trt, shared_grid.config.collision);
    try testing.expect(shared_grid.solid[40 * lbm.WIDTH + 40]);
    try testing.expectEqual(@as(u64, 0), shared_grid.steps);
}

test "divRound rounds to nearest and is odd-symmetric" {
    try testing.expectEqual(@as(i64, 3), lbm.divRound(5, 2));
    try testing.expectEqual(@as(i64, -3), lbm.divRound(-5, 2));
    try testing.expectEqual(@as(i64, 2), lbm.divRound(7, 4));
    try testing.expectEqual(@as(i64, 1), lbm.divRound(4, 3));
    try testing.expectEqual(@as(i64, 2), lbm.divRound(5, 3));
    var rng = std.Random.DefaultPrng.init(1);
    const r = rng.random();
    for (0..20000) |_| {
        const num = r.intRangeAtMost(i64, -(1 << 50), 1 << 50);
        const den = r.intRangeAtMost(i64, 1, 1 << 30);
        const q = lbm.divRound(num, den);
        try testing.expectEqual(-q, lbm.divRound(-num, den));
        // |num - q den| <= den / 2
        try testing.expect(2 * @abs(num - q * den) <= den);
    }
}

test "equilibrium has exact density and momentum for any input" {
    var rng = std.Random.DefaultPrng.init(2);
    const r = rng.random();
    for (0..20000) |_| {
        const rho = r.intRangeAtMost(i64, lbm.Q / 4, 4 * lbm.Q);
        const ux = r.intRangeAtMost(i64, -lbm.Q / 3, lbm.Q / 3);
        const uy = r.intRangeAtMost(i64, -lbm.Q / 3, lbm.Q / 3);
        var eq: [9]i64 = undefined;
        lbm.equilibrium(rho, ux, uy, &eq);
        var s: i64 = 0;
        var mx: i64 = 0;
        var my: i64 = 0;
        for (0..9) |i| {
            s += eq[i];
            mx += lbm.CX[i] * eq[i];
            my += lbm.CY[i] * eq[i];
        }
        try testing.expectEqual(rho, s);
        try testing.expectEqual(ux, mx);
        try testing.expectEqual(uy, my);
        // Rounding moves each population by at most a few units from the real value.
        const fr = lbm.fromQ(rho);
        const fx = lbm.fromQ(ux);
        const fy = lbm.fromQ(uy);
        for (1..9) |i| {
            const cu = @as(f64, @floatFromInt(lbm.CX[i])) * fx + @as(f64, @floatFromInt(lbm.CY[i])) * fy;
            const w = @as(f64, @floatFromInt(lbm.W36[i])) / 36.0;
            const exact = w * (fr + 3 * cu + 4.5 * cu * cu - 1.5 * (fx * fx + fy * fy)) * 65536.0;
            try testing.expect(@abs(@as(f64, @floatFromInt(eq[i])) - exact) <= 2.5);
        }
    }
}

test "mass and momentum are conserved exactly in a periodic box (BGK and TRT)" {
    for ([_]lbm.Collision{ .bgk, .trt }) |op| {
        const c = try validate.conservation(gpa, op, 400);
        try testing.expectEqual(@as(i64, 0), c.mass_drift);
        try testing.expectEqual(@as(i64, 0), c.momentum_x_drift);
        try testing.expectEqual(@as(i64, 0), c.momentum_y_drift);
    }
}

test "closed box with walls and obstacles conserves mass exactly" {
    var l = try lbm.Lattice.create(gpa, 30, 20);
    defer l.destroy(gpa);
    try l.setConfig(.{ .collision = .trt, .x_boundary = .closed, .y_boundary = .wall, .force_x_q16 = 300, .force_y_q16 = -200 });
    validate.scatterObstacles(&l, 99);
    const m0 = l.totals().rho;
    for (0..500) |_| l.step();
    try testing.expectEqual(m0, l.totals().rho);
    try testing.expect(l.stats.stable());
}

test "momentum exchange accounts for every unit of momentum" {
    for ([_]lbm.Collision{ .bgk, .trt }) |op| {
        const b = try validate.momentumBalance(gpa, op, 300);
        try testing.expectEqual(@as(i64, 0), b.residual_x);
        try testing.expectEqual(@as(i64, 0), b.residual_y);
        try testing.expect(b.final_drag > 0);
    }
}

test "symmetric channel flow has exactly zero lift" {
    for ([_]lbm.Collision{ .bgk, .trt }) |op| {
        const m = try validate.mirrorSymmetry(gpa, op, 300);
        try testing.expectEqual(@as(i64, 0), m.max_abs_lift);
        try testing.expectEqual(@as(usize, 0), m.asymmetric_cells);
        try testing.expect(m.drag > 0);
    }
}

test "Poiseuille profile matches the analytic parabola; TRT walls are exact" {
    const trt2 = try validate.poiseuille(gpa, .trt, 2 * lbm.Q, 24, 6000);
    const bgk2 = try validate.poiseuille(gpa, .bgk, 2 * lbm.Q, 24, 6000);
    const trt1 = try validate.poiseuille(gpa, .trt, lbm.Q, 16, 4000);
    try testing.expect(trt2.l2_error < 0.005);
    try testing.expect(trt1.l2_error < 0.005);
    // BGK bounce-back walls drift with tau; at tau = 2 the slip error is large.
    try testing.expect(bgk2.l2_error > 0.02 and bgk2.l2_error < 0.06);
    try testing.expect(trt2.l2_error * 4 < bgk2.l2_error);
    try testing.expect(@abs(trt2.u_max - trt2.u_max_analytic) / trt2.u_max_analytic < 0.005);
}

test "Zou-He channel carries a uniform mass flux at steady state" {
    const ch = try validate.channelFlux(gpa, .trt, 3000);
    try testing.expect(ch.max_flux_deviation < 0.01);
    try testing.expect(@abs(ch.outlet_density - 1.0) < 1e-3);
    try testing.expect(@abs(ch.inlet_flux - 24 * 0.05) < 0.01);
}

test "parabolic inlet profile is symmetric and peaks at the inflow speed" {
    var l = try lbm.Lattice.create(gpa, 8, 21);
    defer l.destroy(gpa);
    try l.setConfig(.{ .inflow_q16 = 6000, .inlet_profile = .parabolic });
    for (0..21) |y| try testing.expectEqual(l.inletVelocity(y), l.inletVelocity(20 - y));
    try testing.expect(l.inletVelocity(10) <= 6000 and l.inletVelocity(10) > 5980);
    try testing.expect(l.inletVelocity(0) > 0);
}

test "configuration validation rejects out-of-range values" {
    var l = try lbm.Lattice.create(gpa, 3, 3);
    defer l.destroy(gpa);
    try testing.expectError(error.TauOutOfRange, l.setTau(lbm.Q / 2));
    try testing.expectError(error.TauOutOfRange, l.setTau(17 * lbm.Q));
    try testing.expectError(error.MagicOutOfRange, l.setConfig(.{ .magic_q16 = 0 }));
    try testing.expectError(error.InflowOutOfRange, l.setConfig(.{ .inflow_q16 = lbm.Q }));
    try testing.expectError(error.ForceOutOfRange, l.setConfig(.{ .force_y_q16 = -2000 }));
    try testing.expectError(error.DensityOutOfRange, l.setConfig(.{ .outlet_density_q16 = 0 }));
    try testing.expectError(error.InvalidDimensions, lbm.Lattice.create(gpa, 2, 10));
    try testing.expectError(error.InvalidDimensions, lbm.Lattice.create(gpa, 8192, 8192));
    try testing.expectEqual(lbm.TAU_DEFAULT, l.tau_q16);
}

test "released obstacle cells restart as fluid at rest" {
    var l = try lbm.Lattice.create(gpa, 10, 10);
    defer l.destroy(gpa);
    l.resetFlow();
    l.setObstacleCell(4, 4, true);
    for (0..20) |_| l.step();
    l.setObstacleCell(4, 4, false);
    const m = lbm.momentsOf(l.f[44]);
    try testing.expectEqual(@as(i64, lbm.Q), m.rho);
    try testing.expectEqual(@as(i64, 0), m.jx);
    try testing.expectEqual(@as(i64, 0), m.jy);
}

test "hostile configurations never overflow and replay identically" {
    var rng = std.Random.DefaultPrng.init(0xbad5eed);
    const r = rng.random();
    for (0..24) |_| {
        const w = r.intRangeAtMost(usize, 3, 24);
        const h = r.intRangeAtMost(usize, 3, 24);
        const cfg = lbm.Config{
            .collision = if (r.boolean()) .bgk else .trt,
            .magic_q16 = r.intRangeAtMost(i32, 1, lbm.MAGIC_MAX),
            .x_boundary = r.enumValue(lbm.XBoundary),
            .y_boundary = r.enumValue(lbm.YBoundary),
            .inflow_q16 = r.intRangeAtMost(i32, -lbm.INFLOW_MAX, lbm.INFLOW_MAX),
            .inlet_profile = r.enumValue(lbm.InletProfile),
            .outlet_density_q16 = r.intRangeAtMost(i32, lbm.DENSITY_MIN, lbm.DENSITY_MAX),
            .force_x_q16 = r.intRangeAtMost(i32, -lbm.FORCE_MAX, lbm.FORCE_MAX),
            .force_y_q16 = r.intRangeAtMost(i32, -lbm.FORCE_MAX, lbm.FORCE_MAX),
        };
        const tau = if (r.boolean()) lbm.TAU_MIN else r.intRangeAtMost(i32, lbm.TAU_MIN, lbm.TAU_MAX);
        const seed = r.int(u64);
        var digests: [2]u64 = undefined;
        for (&digests) |*d| {
            var l = try lbm.Lattice.create(gpa, w, h);
            defer l.destroy(gpa);
            try l.setTau(tau);
            try l.setConfig(cfg);
            l.resetFlow();
            validate.scatterObstacles(&l, seed);
            for (0..150) |_| l.step();
            d.* = l.digest();
        }
        try testing.expectEqual(digests[0], digests[1]);
    }
}

test "corrupted distributions are clamped and reported as unstable" {
    var l = try lbm.Lattice.create(gpa, 6, 6);
    defer l.destroy(gpa);
    try l.setConfig(.{ .x_boundary = .periodic, .y_boundary = .periodic });
    for (l.f, 0..) |*cell, i| {
        for (cell) |*v| v.* = if (i % 2 == 0) std.math.maxInt(i32) else std.math.minInt(i32);
    }
    for (0..5) |_| l.step();
    try testing.expect(!l.stats.stable());
    // Exporting and rendering a corrupted state must not trap either.
    var aw: std.Io.Writer.Allocating = .init(gpa);
    defer aw.deinit();
    for ([_]fieldio.Field{ .speed, .vorticity, .density, .solid }) |field| try fieldio.writePgm(&aw.writer, &l, field);
    try fieldio.writeVtk(&aw.writer, &l);
    var rgba: [6 * 6 * 4]u8 = undefined;
    for ([_]lbm.RenderMode{ .speed, .vorticity, .density }) |mode| l.render(&rgba, mode);
}

test "TRT with Lambda = 3/16 has omega- from the magic parameter" {
    var l = try lbm.Lattice.create(gpa, 3, 3);
    defer l.destroy(gpa);
    try l.setTau(lbm.Q); // tau+ = 1 -> tau- = 0.1875 / 0.5 + 0.5 = 0.875
    try l.setConfig(.{ .collision = .trt });
    const rates = l.relaxationRates();
    try testing.expectEqual(@as(i64, lbm.Q), rates.plus);
    try testing.expectEqual(lbm.divRound(@as(i64, lbm.Q) << 16, lbm.toQ(0.875)), rates.minus);
    try l.setConfig(.{ .collision = .bgk });
    const b = l.relaxationRates();
    try testing.expectEqual(b.plus, b.minus);
}

test "render modes fill every pixel with opaque colour" {
    var l = try lbm.Lattice.create(gpa, 32, 16);
    defer l.destroy(gpa);
    l.resetFlow();
    l.setCircleObstacle(8, 8, 3, true);
    for (0..50) |_| l.step();
    const rgba = try gpa.alloc(u8, 32 * 16 * 4);
    defer gpa.free(rgba);
    for ([_]lbm.RenderMode{ .speed, .vorticity, .density }) |mode| {
        @memset(rgba, 0);
        l.render(rgba, mode);
        var i: usize = 3;
        while (i < rgba.len) : (i += 4) try testing.expectEqual(@as(u8, 255), rgba[i]);
    }
}

// ---- field export and obstacle import ----------------------------------------------

test "PGM export writes a header and one byte per cell" {
    var l = try lbm.Lattice.create(gpa, 20, 10);
    defer l.destroy(gpa);
    l.resetFlow();
    l.setCircleObstacle(5, 5, 2, true);
    for (0..30) |_| l.step();
    for ([_]fieldio.Field{ .speed, .vorticity, .density, .solid }) |field| {
        var aw: std.Io.Writer.Allocating = .init(gpa);
        defer aw.deinit();
        try fieldio.writePgm(&aw.writer, &l, field);
        var img = try fieldio.parsePgm(gpa, aw.written());
        defer img.deinit(gpa);
        try testing.expectEqual(@as(usize, 20), img.width);
        try testing.expectEqual(@as(usize, 10), img.height);
        try testing.expectEqual(@as(u16, 255), img.maxval);
        if (field == .solid) {
            try testing.expectEqual(@as(u16, 0), img.pixels[5 * 20 + 5]);
            try testing.expectEqual(@as(u16, 255), img.pixels[0]);
        }
    }
}

test "PGM solid export round-trips through mask import" {
    var a = try lbm.Lattice.create(gpa, 24, 12);
    defer a.destroy(gpa);
    validate.scatterObstacles(&a, 7);
    var aw: std.Io.Writer.Allocating = .init(gpa);
    defer aw.deinit();
    try fieldio.writePgm(&aw.writer, &a, .solid);
    var img = try fieldio.parsePgm(gpa, aw.written());
    defer img.deinit(gpa);
    var b = try lbm.Lattice.create(gpa, 24, 12);
    defer b.destroy(gpa);
    const n = fieldio.applyMask(&b, img, false);
    try testing.expectEqual(a.solidCount(), n);
    try testing.expectEqualSlices(bool, a.solid, b.solid);
}

test "VTK export lists every point for each array" {
    var l = try lbm.Lattice.create(gpa, 5, 4);
    defer l.destroy(gpa);
    l.setObstacleCell(2, 2, true);
    var aw: std.Io.Writer.Allocating = .init(gpa);
    defer aw.deinit();
    try fieldio.writeVtk(&aw.writer, &l);
    const text = aw.written();
    try testing.expect(std.mem.startsWith(u8, text, "# vtk DataFile Version 3.0\n"));
    try testing.expect(std.mem.indexOf(u8, text, "DIMENSIONS 5 4 1\n") != null);
    try testing.expect(std.mem.indexOf(u8, text, "POINT_DATA 20\n") != null);
    var lines: usize = 0;
    var it = std.mem.splitScalar(u8, text, '\n');
    while (it.next()) |line| {
        if (line.len > 0) lines += 1;
    }
    // 8 header lines, 4 array headers (+3 lookup tables), 4 x 20 values.
    try testing.expectEqual(@as(usize, 8 + 7 + 80), lines);
}

test "PGM parser accepts P2, P5 and 16-bit rasters with comments" {
    {
        var img = try fieldio.parsePgm(gpa, "P2\n# comment\n3 2\n# another\n9\n0 9 4\n9 0 1\n");
        defer img.deinit(gpa);
        try testing.expectEqualSlices(u16, &.{ 0, 9, 4, 9, 0, 1 }, img.pixels);
    }
    {
        var img = try fieldio.parsePgm(gpa, "P5 2 2 255\n\x00\xff\x10\x80");
        defer img.deinit(gpa);
        try testing.expectEqualSlices(u16, &.{ 0, 255, 16, 128 }, img.pixels);
    }
    {
        var img = try fieldio.parsePgm(gpa, "P5\n2 1\n1000\n\x03\xe8\x00\x01");
        defer img.deinit(gpa);
        try testing.expectEqualSlices(u16, &.{ 1000, 1 }, img.pixels);
    }
}

test "PGM parser rejects hostile input" {
    const cases = [_]struct { bytes: []const u8, err: fieldio.PgmError }{
        .{ .bytes = "", .err = error.NotPgm },
        .{ .bytes = "P6\n1 1\n255\n\x00\x00\x00", .err = error.NotPgm },
        .{ .bytes = "P5x1 1 255\n\x00", .err = error.BadHeader },
        .{ .bytes = "P5\n0 4\n255\n", .err = error.BadDimensions },
        .{ .bytes = "P5\n99999 1\n255\n", .err = error.BadDimensions },
        .{ .bytes = "P5\n16384 16384\n255\n", .err = error.BadDimensions },
        .{ .bytes = "P5\n2 2\n0\n\x00\x00\x00\x00", .err = error.BadMaxval },
        .{ .bytes = "P5\n2 2\n70000\n", .err = error.BadMaxval },
        .{ .bytes = "P5\n2 2\n255\n\x00\x00\x00", .err = error.TruncatedData },
        .{ .bytes = "P5\n2 1\n100\n\x00\xc8", .err = error.BadPixel },
        .{ .bytes = "P5\n1 1\n1000\n\xff\xff", .err = error.BadPixel },
        .{ .bytes = "P2\n2 2\n10\n1 2 3", .err = error.TruncatedData },
        .{ .bytes = "P2\n2 1\n10\n1 11", .err = error.BadPixel },
        .{ .bytes = "P2\n2 1\n10\n1 x", .err = error.BadPixel },
        .{ .bytes = "P5\n2 2\n255", .err = error.BadHeader },
    };
    for (cases) |c| try testing.expectError(c.err, fieldio.parsePgm(gpa, c.bytes));
}

test "PGM parser survives random mutations" {
    const seed_img = "P2\n# mask\n4 3\n255\n0 255 0 255\n255 0 255 0\n12 34 56 78\n";
    var rng = std.Random.DefaultPrng.init(3);
    const r = rng.random();
    var buf: [128]u8 = undefined;
    for (0..5000) |_| {
        const len = r.intRangeAtMost(usize, 0, seed_img.len);
        @memcpy(buf[0..len], seed_img[0..len]);
        const flips = r.intRangeAtMost(usize, 0, 4);
        for (0..flips) |_| {
            if (len == 0) break;
            buf[r.uintLessThan(usize, len)] = r.int(u8);
        }
        if (fieldio.parsePgm(gpa, buf[0..len])) |img_const| {
            var img = img_const;
            defer img.deinit(gpa);
            try testing.expectEqual(img.width * img.height, img.pixels.len);
            for (img.pixels) |p| try testing.expect(p <= img.maxval);
        } else |_| {}
    }
}

test "mask import resamples with nearest neighbour and honours invert" {
    var img = try fieldio.parsePgm(gpa, "P2\n2 2\n255\n0 255\n255 0\n");
    defer img.deinit(gpa);
    var l = try lbm.Lattice.create(gpa, 4, 4);
    defer l.destroy(gpa);
    try testing.expectEqual(@as(usize, 8), fieldio.applyMask(&l, img, false));
    try testing.expect(l.solid[0] and l.solid[1] and !l.solid[2] and l.solid[15]);
    try testing.expectEqual(@as(usize, 8), fieldio.applyMask(&l, img, true));
    try testing.expect(!l.solid[0] and l.solid[2]);
}

// ---- command line ----------------------------------------------------------------

test "CLI options parse and reject bad input" {
    var bad: []const u8 = "";
    const o = try cli.parseRunOptions(&.{ "--width", "64", "--height", "32", "--collision", "trt", "--tau", "0.8", "--cylinder", "10,16,4", "--profile", "parabolic", "--json", "--every", "10", "--field", "vorticity" }, &bad);
    try testing.expectEqual(@as(usize, 64), o.width);
    try testing.expectEqual(lbm.Collision.trt, o.config.collision);
    try testing.expectEqual(lbm.toQ(0.8), o.tau_q16);
    try testing.expectEqual(@as(usize, 1), o.cylinder_count);
    try testing.expectEqual(@as(i32, 4), o.cylinders[0].r);
    try testing.expectEqual(lbm.InletProfile.parabolic, o.config.inlet_profile);
    try testing.expectEqual(fieldio.Field.vorticity, o.field);
    try testing.expect(o.json);

    try testing.expectError(error.UnknownOption, cli.parseRunOptions(&.{"--nope"}, &bad));
    try testing.expectEqualStrings("--nope", bad);
    try testing.expectError(error.MissingValue, cli.parseRunOptions(&.{"--steps"}, &bad));
    try testing.expectError(error.BadValue, cli.parseRunOptions(&.{ "--steps", "-4" }, &bad));
    try testing.expectError(error.BadValue, cli.parseRunOptions(&.{ "--tau", "nan" }, &bad));
    try testing.expectError(error.BadValue, cli.parseRunOptions(&.{ "--tau", "1e300" }, &bad));
    try testing.expectError(error.BadValue, cli.parseRunOptions(&.{ "--cylinder", "1,2" }, &bad));
    try testing.expectError(error.BadValue, cli.parseRunOptions(&.{ "--cylinder", "1,2,3,4" }, &bad));
    try testing.expectError(error.BadValue, cli.parseRunOptions(&.{ "--collision", "mrt" }, &bad));
    try testing.expectError(error.BadValue, cli.parseRunOptions(&.{ "--field", "pressure" }, &bad));
    try testing.expectError(error.BadValue, cli.parseRunOptions(&.{ "--width", "99999" }, &bad));
}

fn runCli(args: []const []const u8, out: *std.Io.Writer.Allocating, err: *std.Io.Writer.Allocating) !u8 {
    return cli.dispatch(testing.io, gpa, &out.writer, &err.writer, args);
}

test "run --json emits telemetry that parses and matches the library" {
    var out: std.Io.Writer.Allocating = .init(gpa);
    defer out.deinit();
    var err: std.Io.Writer.Allocating = .init(gpa);
    defer err.deinit();
    const code = try runCli(&.{ "run", "--width", "48", "--height", "24", "--steps", "60", "--json", "--every", "20", "--collision", "trt" }, &out, &err);
    try testing.expectEqual(@as(u8, 0), code);
    var lines = std.mem.splitScalar(u8, std.mem.trimEnd(u8, out.written(), "\n"), '\n');
    var count: usize = 0;
    var last_digest: u64 = 0;
    while (lines.next()) |line| {
        const parsed = try std.json.parseFromSlice(std.json.Value, gpa, line, .{});
        defer parsed.deinit();
        const obj = parsed.value.object;
        const tele = if (obj.get("telemetry")) |t| t.object else obj;
        try testing.expect(tele.get("stable").?.bool);
        last_digest = @intCast(tele.get("digest32").?.integer);
        count += 1;
    }
    try testing.expectEqual(@as(usize, 4), count); // 3 stream lines + summary

    var l = try cli.setupLattice(gpa, .{ .width = 48, .height = 24, .config = .{ .collision = .trt } }, null);
    defer l.destroy(gpa);
    for (0..60) |_| l.step();
    try testing.expectEqual(@as(u64, l.digest32()), last_digest);
}

test "bench --json reports both operators; bad commands return 2" {
    var out: std.Io.Writer.Allocating = .init(gpa);
    defer out.deinit();
    var err: std.Io.Writer.Allocating = .init(gpa);
    defer err.deinit();
    try testing.expectEqual(@as(u8, 0), try runCli(&.{ "bench", "--width", "32", "--height", "16", "--steps", "20", "--json" }, &out, &err));
    const parsed = try std.json.parseFromSlice(std.json.Value, gpa, out.written(), .{});
    defer parsed.deinit();
    const results = parsed.value.object.get("results").?.array;
    try testing.expectEqual(@as(usize, 2), results.items.len);
    try testing.expectEqualStrings("trt", results.items[1].object.get("collision").?.string);

    try testing.expectEqual(@as(u8, 2), try runCli(&.{"frobnicate"}, &out, &err));
    try testing.expectEqual(@as(u8, 2), try runCli(&.{ "run", "--width", "2" }, &out, &err));
    try testing.expectEqual(@as(u8, 2), try runCli(&.{ "run", "--inflow", "0.9" }, &out, &err));
    try testing.expectEqual(@as(u8, 2), try runCli(&.{ "run", "--tau", "0.5" }, &out, &err));
    try testing.expectEqual(@as(u8, 0), try runCli(&.{"version"}, &out, &err));
}
