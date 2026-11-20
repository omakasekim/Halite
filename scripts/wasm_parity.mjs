// Cross-target determinism check: the freestanding WebAssembly build and the
// native CLI must reach bit-identical lattice states.
//
// usage: node scripts/wasm_parity.mjs zig-out/bin/halite.wasm zig-out/bin/halite
import { readFileSync } from "node:fs";
import { execFileSync } from "node:child_process";

const [wasmPath, cliPath] = process.argv.slice(2);
if (!wasmPath || !cliPath) {
  console.error("usage: node scripts/wasm_parity.mjs <halite.wasm> <halite>");
  process.exit(2);
}

const bytes = readFileSync(wasmPath);
const module = new WebAssembly.Module(bytes);
if (WebAssembly.Module.imports(module).length !== 0) {
  console.error("FAIL: the module must not import anything");
  process.exit(1);
}

function fresh() {
  return new WebAssembly.Instance(module, { env: {} }).exports;
}

function cliDigest(args) {
  const out = execFileSync(cliPath, ["run", "--json", ...args], { encoding: "utf8" });
  return JSON.parse(out).telemetry.digest32 >>> 0;
}

let failures = 0;
function check(name, ok, detail) {
  console.log(`${ok ? "PASS" : "FAIL"} ${name}: ${detail}`);
  if (!ok) failures++;
}

const Q = 65536;
const cases = [
  { name: "default bgk channel", steps: 300, setup: () => {}, args: [] },
  {
    name: "trt parabolic channel",
    steps: 250,
    setup: (e) => {
      e.halite_set_collision(1);
      e.halite_set_tau(Math.round(0.8 * Q));
      e.halite_set_inflow(Math.round(0.06 * Q), 1);
    },
    args: ["--collision", "trt", "--tau", "0.8", "--inflow", "0.06", "--profile", "parabolic"],
  },
  {
    name: "periodic forced trt",
    steps: 200,
    setup: (e) => {
      e.halite_set_collision(1);
      e.halite_set_boundary(1, 0);
      e.halite_set_body_force(7, 0);
    },
    args: ["--collision", "trt", "--x-boundary", "periodic", "--force-x", String(7 / Q)],
  },
];

for (const c of cases) {
  const e = fresh();
  e.halite_init();
  c.setup(e);
  e.halite_reset();
  for (let i = 0; i < c.steps; i++) e.halite_step();
  const wasm = e.halite_get_digest() >>> 0;
  const native = cliDigest(["--steps", String(c.steps), ...c.args]);
  check(c.name, wasm === native, `wasm ${wasm.toString(16)} native ${native.toString(16)}`);
}

// The live workbench sequence: init, clear, paint a cylinder cell by cell,
// two steps per frame, render, read RGBA through exported memory.
{
  const e = fresh();
  e.halite_init();
  const W = e.halite_get_width();
  const H = e.halite_get_height();
  for (let y = 0; y < H; y++) for (let x = 0; x < W; x++) e.halite_set_obstacle(x, y, 0);
  for (let y = 32; y <= 48; y++)
    for (let x = 27; x <= 43; x++)
      if ((x - 35) ** 2 + (y - 40) ** 2 <= 64) e.halite_set_obstacle(x, y, 1);
  let stable = true;
  for (let frame = 0; frame < 1500; frame++) {
    e.halite_step();
    e.halite_step();
    stable = stable && e.halite_is_stable() === 1;
  }
  e.halite_render();
  const rgba = new Uint8ClampedArray(e.memory.buffer, e.halite_get_rgba_ptr(), W * H * 4);
  const distinct = new Set();
  for (let i = 0; i < rgba.length; i += 4) distinct.add(rgba[i] | (rgba[i + 1] << 8) | (rgba[i + 2] << 16));
  check("workbench 3000 steps stable", stable, `max speed ${e.halite_get_max_speed().toFixed(4)}, mean density ${e.halite_get_mean_density().toFixed(4)}`);
  check("workbench frame has structure", distinct.size > 20 && rgba[3] === 255, `${distinct.size} distinct colours`);
  check("workbench drag positive", e.halite_get_force_x() > 0, `drag ${e.halite_get_force_x().toFixed(5)} lift ${e.halite_get_force_y().toFixed(5)}`);
  check("step counter", e.halite_get_step_count() === 3000, `${e.halite_get_step_count()}`);
}

// Rejected settings leave the configuration unchanged.
{
  const e = fresh();
  e.halite_init();
  const ok =
    e.halite_set_tau(100) === 0 &&
    e.halite_get_tau() === 78643 &&
    e.halite_set_collision(7) === 0 &&
    e.halite_set_boundary(3, 0) === 0 &&
    e.halite_set_inflow(Q, 0) === 0 &&
    e.halite_set_body_force(5000, 0) === 0 &&
    e.halite_set_render_mode(9) === 0 &&
    e.halite_get_version() === 0x000200;
  check("hostile settings rejected", ok, `tau ${e.halite_get_tau()}`);
}

// Mask upload and render modes.
{
  const e = fresh();
  e.halite_init();
  const W = e.halite_get_width();
  const H = e.halite_get_height();
  const mask = new Uint8Array(e.memory.buffer, e.halite_get_mask_ptr(), W * H);
  mask.fill(0);
  for (let y = 30; y < 50; y++) mask[y * W + 60] = 255;
  const solid = e.halite_apply_mask();
  let okModes = true;
  for (const mode of [0, 1, 2]) {
    okModes = okModes && e.halite_set_render_mode(mode) === 1;
    e.halite_step_n(20);
    e.halite_render();
  }
  check("mask upload", solid === 20, `${solid} solid cells`);
  check("render modes", okModes && e.halite_is_stable() === 1, "speed, vorticity, density");
}

if (failures > 0) {
  console.error(`${failures} parity check(s) failed`);
  process.exit(1);
}
console.log("all parity checks passed");
