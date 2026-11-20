// Copyright (c) 2026 jk9576. All rights reserved.
// Halite: Deterministic Lattice Boltzmann Fluid / Gas Automaton (D2Q9)
// Licensed under the MIT License.

const std = @import("std");
pub const lbm = @import("lbm.zig");
pub const fieldio = @import("fieldio.zig");
pub const validate = @import("validate.zig");
pub const cli = @import("cli.zig");

pub fn main(init: std.process.Init) !u8 {
    const io = init.io;
    const gpa = init.gpa;
    const argv = try init.minimal.args.toSlice(init.arena.allocator());

    var out_buf: [16 * 1024]u8 = undefined;
    var out = std.Io.File.stdout().writerStreaming(io, &out_buf);
    var err_buf: [1024]u8 = undefined;
    var err = std.Io.File.stderr().writerStreaming(io, &err_buf);

    const args = try init.arena.allocator().alloc([]const u8, argv.len -| 1);
    for (args, 0..) |*a, i| a.* = argv[i + 1];

    const code = try cli.dispatch(io, gpa, &out.interface, &err.interface, args);
    try out.interface.flush();
    try err.interface.flush();
    return code;
}

test {
    _ = @import("tests.zig");
}
