// Copyright (c) 2026 jk9576. All rights reserved.
// Halite: Deterministic Lattice Boltzmann Fluid / Gas Automaton (D2Q9)
// Licensed under the MIT License.
//
// Field export (binary PGM, legacy VTK) and obstacle import from PGM masks.

const std = @import("std");
const lbm = @import("lbm.zig");

pub const Field = enum {
    speed,
    vorticity,
    density,
    solid,

    pub fn parse(s: []const u8) ?Field {
        return std.meta.stringToEnum(Field, s);
    }
};

/// Integer sample of a field at a cell. Speed is |u| and density is rho, both
/// Q16.16; vorticity is twice the curl in Q16.16; solid is 1 or 0.
pub fn sample(l: *const lbm.Lattice, field: Field, x: usize, y: usize) i64 {
    const idx = y * l.width + x;
    if (field == .solid) return @intFromBool(l.solid[idx]);
    if (l.solid[idx]) return 0;
    return switch (field) {
        .speed => blk: {
            const u = l.velocity(idx);
            // Clamp so a corrupted state cannot overflow the square.
            const ux = std.math.clamp(u[0], -(1 << 30), 1 << 30);
            const uy = std.math.clamp(u[1], -(1 << 30), 1 << 30);
            const sq: u64 = @intCast(ux * ux + uy * uy);
            break :blk @intCast(std.math.sqrt(sq));
        },
        .density => lbm.momentsOf(l.f[idx]).rho,
        .vorticity => l.vorticity2(x, y),
        .solid => unreachable,
    };
}

/// Writes a binary 8-bit PGM (P5). Fluid values are scaled linearly from the
/// field's minimum to maximum over fluid cells; solid cells are 0 except in
/// the `solid` field, where solid is 0 and fluid is 255.
pub fn writePgm(w: *std.Io.Writer, l: *const lbm.Lattice, field: Field) std.Io.Writer.Error!void {
    var lo: i64 = std.math.maxInt(i64);
    var hi: i64 = std.math.minInt(i64);
    if (field != .solid) {
        for (0..l.height) |y| {
            for (0..l.width) |x| {
                if (l.solid[y * l.width + x]) continue;
                const v = sample(l, field, x, y);
                lo = @min(lo, v);
                hi = @max(hi, v);
            }
        }
    }
    try w.print("P5\n# halite field={s} min={d} max={d} (Q16.16)\n{d} {d}\n255\n", .{
        @tagName(field),
        if (lo > hi) 0 else lo,
        if (lo > hi) 0 else hi,
        l.width,
        l.height,
    });
    // Image row r is lattice row y = r, matching the WebAssembly RGBA buffer.
    for (0..l.height) |row| {
        for (0..l.width) |x| {
            const idx = row * l.width + x;
            const byte: u8 = if (field == .solid)
                (if (l.solid[idx]) 0 else 255)
            else if (l.solid[idx] or hi <= lo)
                0
            else
                @intCast(@divTrunc((sample(l, field, x, row) - lo) * 255, hi - lo));
            try w.writeByte(byte);
        }
    }
}

/// Writes a legacy ASCII VTK STRUCTURED_POINTS file with density, velocity,
/// vorticity and the solid mask as point data.
pub fn writeVtk(w: *std.Io.Writer, l: *const lbm.Lattice) std.Io.Writer.Error!void {
    const n = l.cellCount();
    try w.print("# vtk DataFile Version 3.0\nhalite D2Q9 step {d}\nASCII\nDATASET STRUCTURED_POINTS\n", .{l.steps});
    try w.print("DIMENSIONS {d} {d} 1\nORIGIN 0 0 0\nSPACING 1 1 1\nPOINT_DATA {d}\n", .{ l.width, l.height, n });
    try w.writeAll("SCALARS density double 1\nLOOKUP_TABLE default\n");
    for (0..n) |idx| {
        const rho: i64 = if (l.solid[idx]) lbm.Q else lbm.momentsOf(l.f[idx]).rho;
        try w.print("{d:.6}\n", .{lbm.fromQ(rho)});
    }
    try w.writeAll("VECTORS velocity double\n");
    for (0..n) |idx| {
        const u: [2]i64 = if (l.solid[idx]) .{ 0, 0 } else l.velocity(idx);
        try w.print("{d:.6} {d:.6} 0\n", .{ lbm.fromQ(u[0]), lbm.fromQ(u[1]) });
    }
    try w.writeAll("SCALARS vorticity double 1\nLOOKUP_TABLE default\n");
    for (0..l.height) |y| {
        for (0..l.width) |x| {
            const v: i64 = if (l.solid[y * l.width + x]) 0 else l.vorticity2(x, y);
            try w.print("{d:.6}\n", .{lbm.fromQ(v) / 2.0});
        }
    }
    try w.writeAll("SCALARS solid int 1\nLOOKUP_TABLE default\n");
    for (0..n) |idx| try w.print("{d}\n", .{@intFromBool(l.solid[idx])});
}

pub const Image = struct {
    width: usize,
    height: usize,
    maxval: u16,
    /// Row-major, top row first.
    pixels: []u16,

    pub fn deinit(self: *Image, gpa: std.mem.Allocator) void {
        gpa.free(self.pixels);
        self.* = undefined;
    }
};

pub const PgmError = error{
    NotPgm,
    BadHeader,
    BadDimensions,
    BadMaxval,
    TruncatedData,
    BadPixel,
    OutOfMemory,
};

pub const MAX_IMAGE_DIM: usize = 16384;
pub const MAX_IMAGE_PIXELS: usize = 1 << 26;

const Cursor = struct {
    bytes: []const u8,
    pos: usize = 0,

    fn skipSpaceAndComments(c: *Cursor) void {
        while (c.pos < c.bytes.len) {
            const b = c.bytes[c.pos];
            if (b == '#') {
                while (c.pos < c.bytes.len and c.bytes[c.pos] != '\n' and c.bytes[c.pos] != '\r') c.pos += 1;
            } else if (std.ascii.isWhitespace(b)) {
                c.pos += 1;
            } else break;
        }
    }

    fn number(c: *Cursor, limit: u64) PgmError!u64 {
        c.skipSpaceAndComments();
        const start = c.pos;
        var v: u64 = 0;
        while (c.pos < c.bytes.len and std.ascii.isDigit(c.bytes[c.pos])) : (c.pos += 1) {
            v = v * 10 + (c.bytes[c.pos] - '0');
            if (v > limit) return error.BadHeader;
        }
        if (c.pos == start) return error.BadHeader;
        return v;
    }
};

/// Parses a plain (P2) or binary (P5) PGM, 8 or 16 bit. Rejects oversized,
/// truncated or out-of-range input instead of guessing.
pub fn parsePgm(gpa: std.mem.Allocator, bytes: []const u8) PgmError!Image {
    if (bytes.len < 2 or bytes[0] != 'P' or (bytes[1] != '2' and bytes[1] != '5')) return error.NotPgm;
    const binary = bytes[1] == '5';
    var c = Cursor{ .bytes = bytes, .pos = 2 };
    if (c.pos >= bytes.len or !(std.ascii.isWhitespace(bytes[c.pos]) or bytes[c.pos] == '#')) return error.BadHeader;
    const width = c.number(MAX_IMAGE_DIM) catch return error.BadDimensions;
    const height = c.number(MAX_IMAGE_DIM) catch return error.BadDimensions;
    if (width == 0 or height == 0 or width * height > MAX_IMAGE_PIXELS) return error.BadDimensions;
    const maxval = c.number(65535) catch return error.BadMaxval;
    if (maxval == 0) return error.BadMaxval;

    const n: usize = @intCast(width * height);
    const pixels = try gpa.alloc(u16, n);
    errdefer gpa.free(pixels);

    if (binary) {
        // Exactly one whitespace byte separates the header from the raster.
        if (c.pos >= bytes.len or !std.ascii.isWhitespace(bytes[c.pos])) return error.BadHeader;
        c.pos += 1;
        const bpp: usize = if (maxval > 255) 2 else 1;
        if (bytes.len - c.pos < n * bpp) return error.TruncatedData;
        for (pixels, 0..) |*p, i| {
            const v: u16 = if (bpp == 2)
                std.mem.readInt(u16, bytes[c.pos + 2 * i ..][0..2], .big)
            else
                bytes[c.pos + i];
            if (v > maxval) return error.BadPixel;
            p.* = v;
        }
    } else {
        for (pixels) |*p| {
            c.skipSpaceAndComments();
            if (c.pos >= bytes.len) return error.TruncatedData;
            const v = c.number(maxval) catch return error.BadPixel;
            p.* = @intCast(v);
        }
    }
    return .{ .width = @intCast(width), .height = @intCast(height), .maxval = @intCast(maxval), .pixels = pixels };
}

/// Replaces the lattice obstacles with an image mask, resampled with nearest
/// neighbour to the lattice size. Pixels darker than half of maxval are solid
/// (lighter ones with `invert`). Image row r maps to lattice row y = r, the
/// same orientation as the WebAssembly RGBA buffer.
/// Returns the number of solid cells.
pub fn applyMask(l: *lbm.Lattice, img: Image, invert: bool) usize {
    var count: usize = 0;
    for (0..l.height) |y| {
        const iy = y * img.height / l.height;
        for (0..l.width) |x| {
            const ix = x * img.width / l.width;
            const v = img.pixels[iy * img.width + ix];
            const dark = @as(u32, v) * 2 < img.maxval;
            const s = dark != invert;
            l.setObstacleCell(x, y, s);
            count += @intFromBool(s);
        }
    }
    return count;
}
