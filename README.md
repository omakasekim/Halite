# Halite

Deterministic Lattice Boltzmann Fluid / Gas Automaton

A bit-packed 2D hydrodynamic simulation engine implementing the D2Q9 Lattice Boltzmann Method (LBM) with zero floating-point nondeterminism.

## Architecture

- Fixed-Point Determinism: Velocity vectors and particle distributions are discretized into fixed-point integer arithmetic (Q16.16), guaranteeing identical, bit-for-bit simulation frames across x86_64, AArch64, and WebAssembly targets.
- D2Q9 Lattice Scheme: Evaluates nine velocity directions and weights per grid cell with bounce-back boundary conditions for complex obstacle geometries.
- Direct Canvas Export: Exports raw memory slices directly into an RGBA byte buffer compatible with HTML5 Canvas `putImageData()` with minimal JavaScript glue.
- Freestanding WASM: Compiles to a freestanding WebAssembly binary (`halite.wasm`, about 11 KB in ReleaseSmall) with no imports and no C runtime.

## Solver (0.2.0)

- Collision operators: BGK and TRT (two relaxation times). TRT takes a magic parameter Lambda = (tau+ - 1/2)(tau- - 1/2); the default 3/16 places halfway bounce-back walls exactly for Poiseuille flow.
- Equilibrium: incompressible He-Luo form (rho0 = 1, u = j). Weights are exact multiples of 1/36 and every product is computed in 64-bit integers with one rounding (nearest, ties away from zero), which is odd-symmetric, so mirror-symmetric setups stay exactly mirror symmetric.
- Exact conservation: after rounding, the collision puts the residual back into the cell so that its mass and momentum (plus the body force) are exact. Mass is conserved to the last bit in periodic and closed domains.
- Boundaries: x is `channel` (Zou-He velocity inlet at x = 0, Zou-He pressure outlet at x = width - 1), `periodic`, or `closed` (no-slip walls); y is `wall` (halfway bounce-back) or `periodic`. The inlet profile is uniform or a parabola peaking at the inflow speed.
- Forcing: a uniform body force with Guo's scheme, split into even and odd parts for TRT.
- Forces on obstacles: momentum exchange on every fluid-solid link gives drag and lift each step. Walls are tracked separately. In integer arithmetic these account for every unit of momentum (`P(t+1) - P(t) = N F - F_obstacle - F_wall` holds exactly).
- Hardening: settings are range-checked (tau 0.505 to 16, |inflow| <= 0.3, |force| <= 1/64, outlet density 0.5 to 2, magic in (0, 4]). Intermediates are 64-bit with velocity clamping and saturating stores, so even a corrupted or diverging state cannot overflow. The step then reports itself unstable.

The 0.1.0 solver sealed the right edge with a wall and divided the momentum by the density before using it in an incompressible equilibrium. That broke momentum conservation. On the default workbench setup it overflowed 32-bit integers within about 50 steps. The default is now the open Zou-He channel, and the test suite runs that setup in a Debug build, which traps on overflow.

## Command Line

```bash
zig build -Doptimize=ReleaseFast

# 50-step demo on the default lattice; prints the state digest
./zig-out/bin/halite

# Simulate and print a summary; export speed as PGM and all fields as VTK
./zig-out/bin/halite run --collision trt --steps 3000 --pgm speed.pgm --field speed --vtk flow.vtk

# JSON telemetry: one line every 100 steps, then a summary object
./zig-out/bin/halite run --x-boundary periodic --force-x 0.00002 --collision trt --steps 200 --json --every 100

# Confined cylinder with a parabolic inlet; drag and lift coefficients in the JSON telemetry
./zig-out/bin/halite run --width 200 --height 61 --profile parabolic --inflow 0.06 --tau 0.8 --cylinder 50,30,6 --steps 4000 --json

# Obstacles from a PGM mask (dark pixels are solid, resampled to the lattice)
printf 'P2\n8 4\n255\n255 255 255 255 255 255 255 255\n255 255 0 0 255 255 255 255\n255 255 0 0 255 255 255 255\n255 255 255 255 255 255 255 255\n' > block.pgm
./zig-out/bin/halite run --obstacle block.pgm --steps 500

# Throughput of both collision operators in MLUPS
./zig-out/bin/halite bench --steps 500 --json

# Analytic and integer-exact validation (exit status 1 on failure)
./zig-out/bin/halite validate
```

`halite help` lists every option. Telemetry fields are `step`, `fluid_cells`, `mass`, `mean_density`, `momentum_x`, `momentum_y`, `max_speed`, `drag`, `lift`, `wall_force_x`, `wall_force_y`, `drag_coefficient`, `lift_coefficient`, `stable`, `digest` (FNV-1a 64 of the full state) and `digest32`. Coefficients use `C = 2 F / (U^2 D)`. Here D is the number of lattice rows that contain an obstacle, and U is the inflow speed, or two thirds of it for a parabolic inlet. They are `null` when U or D is zero.

PGM export writes an 8-bit binary image scaled from the field's minimum to its maximum over fluid cells; the header comment records both in Q16.16. Image row r is lattice row y = r, the same orientation as the WebAssembly RGBA buffer. Mask import accepts plain (P2) and binary (P5) PGM, 8 or 16 bit, and rejects truncated data, pixels above maxval, and images larger than 16384 per side or 2^26 pixels.

`halite validate` reports, on the current build:

| Check | Result |
| --- | --- |
| Poiseuille, TRT, tau 2 and tau 0.8 (H = 24) | L2 error 1.5e-3 and 1.4e-3 against the analytic parabola |
| Poiseuille, BGK, tau 2 | L2 error 3.5e-2 (bounce-back wall slip); TRT must be at least 4x better |
| Periodic box, BGK and TRT, 2000 steps | mass and momentum drift exactly 0 |
| Momentum-exchange balance, 1500 steps | residual exactly 0 |
| Cylinder on the centreline, parabolic inlet | lift exactly 0, no asymmetric cell |
| Zou-He channel | column mass flux within 0.06 % of the inlet, outlet density 1 |

On the development machine (aarch64 macOS, ReleaseFast) `halite bench` measures about 30 MLUPS for both operators on the 160 x 80 lattice.

## Exported WebAssembly API

The module imports nothing (`WebAssembly.instantiate(bytes, { env: {} })`). The lattice is fixed at 160 x 80 and the RGBA buffer lives in the exported `memory`.

```c
void halite_init(void);
void halite_reset(void);       // restarts the flow and the default cylinder; keeps settings
void halite_step(void);
void halite_render(void);
void halite_set_obstacle(uint32_t x, uint32_t y, uint32_t is_solid);
uint8_t* halite_get_rgba_ptr(void);
uint32_t halite_get_width(void);
uint32_t halite_get_height(void);
```

Added in 0.2.0 (setters return 1 if accepted and 0 if rejected, leaving the setting unchanged; Q16.16 means value * 65536):

```c
uint32_t halite_get_version(void);                      // 0x000200
uint32_t halite_step_n(uint32_t count);                 // up to 100000 steps; 1 if stable
uint32_t halite_set_collision(uint32_t op);             // 0 BGK, 1 TRT
uint32_t halite_set_trt_magic(int32_t magic_q16);
uint32_t halite_set_tau(int32_t tau_q16);
int32_t  halite_get_tau(void);
uint32_t halite_set_boundary(uint32_t x_mode, uint32_t y_mode); // x: 0 channel, 1 periodic, 2 closed; y: 0 wall, 1 periodic
uint32_t halite_set_inflow(int32_t ux_q16, uint32_t profile);   // profile: 0 uniform, 1 parabolic
uint32_t halite_set_outlet_density(int32_t rho_q16);
uint32_t halite_set_body_force(int32_t fx_q16, int32_t fy_q16);
void     halite_set_circle(int32_t cx, int32_t cy, int32_t r, uint32_t is_solid);
void     halite_clear_obstacles(void);
uint8_t* halite_get_mask_ptr(void);                     // 160*80 bytes, nonzero = solid
uint32_t halite_apply_mask(void);                       // returns the solid cell count
uint32_t halite_set_render_mode(uint32_t mode);         // 0 speed, 1 vorticity, 2 density
double   halite_get_force_x(void);                      // obstacle drag, last step
double   halite_get_force_y(void);                      // obstacle lift, last step
double   halite_get_mean_density(void);
double   halite_get_max_speed(void);
uint32_t halite_is_stable(void);
uint32_t halite_get_step_count(void);
uint32_t halite_get_digest(void);                       // equals digest32 from `halite run --json`
```

## Building and Testing

```bash
# Run unit, property and fuzz-style tests
zig build test

# Compile native binary and freestanding WASM
zig build -Doptimize=ReleaseFast

# Cross-target check: WASM and native reach identical digests; replays the workbench sequence
node scripts/wasm_parity.mjs zig-out/bin/halite.wasm zig-out/bin/halite
```

## License

MIT License - Copyright (c) 2026 Jeong Hoon Kim
