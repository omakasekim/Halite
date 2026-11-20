# Changelog

All notable changes to this project will be documented in this file.

## [0.2.0] - 2026-11-20

### Fixed
- The 0.1.0 solver overflowed 32-bit integers within about 50 steps on the default workbench setup. It sealed the right edge with a bounce-back wall, divided momentum by density before an incompressible equilibrium and so did not conserve momentum. The default is now an open channel, and the collision conserves mass and momentum exactly.

### Changed
- The solver core moved from a fixed 160 x 80 `Grid` to a `Lattice` over caller-owned buffers of 3 to 8192 cells per side and at most 2^24 cells. `Grid` remains as the static 160 x 80 wrapper used by WebAssembly.
- Equilibria use exact 1/36 weights and 64-bit intermediates with sign-symmetric round-to-nearest. The legacy truncating `qMul`/`qDiv` helpers are kept.
- `halite_reset` keeps solver settings (tau, collision, boundaries, forcing) and restores the flow and the default cylinder. Re-opened obstacle cells start as fluid at rest.
- The no-argument demo now prints the state digest it describes.
- The WebAssembly binary grows from 2.0 KB to about 11 KB. It still has no imports.

### Added
- TRT collision operator with a configurable magic parameter (default 3/16), alongside BGK.
- Zou-He velocity inlet (uniform or parabolic profile) and Zou-He pressure outlet. Periodic and closed x edges; wall or periodic y edges.
- Guo body forcing, split into even and odd parts for TRT.
- Momentum-exchange drag and lift on obstacles, with wall forces reported separately. Drag and lift coefficients appear in the telemetry.
- CLI subcommands: `run` (JSON summary, `--every` JSON lines, PGM and VTK field export, PGM obstacle import, cylinders), `bench` (MLUPS for BGK and TRT, text or JSON) and `validate` (Poiseuille against the analytic profile, exact conservation, exact momentum balance, exact mirror symmetry, Zou-He mass flux).
- Hardened PGM parser for P2/P5, 8 and 16 bit, with size limits and nearest-neighbour mask resampling.
- WebAssembly exports for collision, tau, TRT magic, boundaries, inflow, outlet density, body force, circles, bulk mask upload, obstacle clearing, render modes (speed, vorticity, density), forces, mean density, maximum speed, stability, step count, state digest and version. All 0.1.0 exports keep their names, signatures and meaning.
- Tests grow from 2 to 31: unit, property (rounding, equilibrium moments), integer-exact physics invariants, analytic Poiseuille, fuzz-style hostile configurations and PGM mutations, and CLI JSON round trips.
- `scripts/wasm_parity.mjs` checks that WebAssembly and native builds reach identical digests and replays the live workbench call sequence. CI runs it together with `halite validate`.

### CI
- Pinned `mlugg/setup-zig` to v2.2.1 (`d1434d08867e3ee9daa34448df10607b98908d29`); the previous pin did not exist.

## 0.1.0 - 2026-11-08

### Added
- Initial release of Halite deterministic 2D Lattice Boltzmann Method (D2Q9) fluid automaton.
- Fixed-point Q16.16 arithmetic guaranteeing bit-for-bit identical simulation results across targets.
- Direct RGBA frame buffer memory export for HTML5 Canvas putImageData rendering.
- Freestanding WebAssembly compilation target yielding a 2.0 KB binary.
