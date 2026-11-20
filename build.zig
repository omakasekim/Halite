const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    // Native CLI executable
    const exe = b.addExecutable(.{
        .name = "halite",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/main.zig"),
            .target = target,
            .optimize = optimize,
        }),
    });
    b.installArtifact(exe);

    // Freestanding WebAssembly artifact
    const wasm_target = b.resolveTargetQuery(.{
        .cpu_arch = .wasm32,
        .os_tag = .freestanding,
    });
    const wasm = b.addExecutable(.{
        .name = "halite",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/wasm.zig"),
            .target = wasm_target,
            .optimize = .ReleaseSmall,
        }),
    });
    wasm.entry = .disabled;
    wasm.root_module.export_symbol_names = &.{
        "halite_init",
        "halite_reset",
        "halite_step",
        "halite_render",
        "halite_set_obstacle",
        "halite_get_rgba_ptr",
        "halite_get_width",
        "halite_get_height",
        // 0.2.0
        "halite_get_version",
        "halite_step_n",
        "halite_set_collision",
        "halite_set_trt_magic",
        "halite_set_tau",
        "halite_get_tau",
        "halite_set_boundary",
        "halite_set_inflow",
        "halite_set_outlet_density",
        "halite_set_body_force",
        "halite_set_circle",
        "halite_clear_obstacles",
        "halite_get_mask_ptr",
        "halite_apply_mask",
        "halite_set_render_mode",
        "halite_get_force_x",
        "halite_get_force_y",
        "halite_get_mean_density",
        "halite_get_max_speed",
        "halite_is_stable",
        "halite_get_step_count",
        "halite_get_digest",
    };
    b.installArtifact(wasm);

    // Test runner
    const unit_tests = b.addTest(.{
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/main.zig"),
            .target = target,
            .optimize = optimize,
        }),
    });
    const run_unit_tests = b.addRunArtifact(unit_tests);
    const test_step = b.step("test", "Run unit tests");
    test_step.dependOn(&run_unit_tests.step);
}
