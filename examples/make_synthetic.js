// Writes the two example spectra in examples/data/. Both are computed, so the answer is known.
//   Z = j w L + R0 + sum_k R_k / (1 + (j w tau_k)^alpha_k) + W (j w)^-0.5
//   R0 = 0.200 ohm; arcs (R, tau, alpha): contact 0.050 / 2e-4 s / 0.95, 0.080 / 2e-3 s / 0.85, 0.170 / 2e-2 s / 0.90;  W = 0.04
//   -> R_HFI 0.200, R_contact 0.050, R_ct 0.250, width 0.300 ohm
// synthetic_three_arcs.csv                 10 kHz - 10 mHz, no inductance, no noise
// synthetic_inductance_noise.csv           100 kHz - 10 mHz, L = 2e-7 H, 0.3 % noise (fixed seed, so the file is always the same)
'use strict';
const fs = require('fs'), path = require('path');
const R0 = 0.2, arcs = [[0.05, 2e-4, 0.95], [0.08, 2e-3, 0.85], [0.17, 2e-2, 0.9]], W = 0.04, p = 0.5;
function z(f, L) {
  const w = 2 * Math.PI * f; let re = R0, im = w * L;
  for (const [R, tau, a] of arcs) { const m = Math.pow(w * tau, a), c = 1 + m * Math.cos(a * Math.PI / 2), s = m * Math.sin(a * Math.PI / 2), d = c * c + s * s; re += R * c / d; im -= R * s / d; }
  const m = W * Math.pow(w, -p); return [re + m * Math.cos(p * Math.PI / 2), im - m * Math.sin(p * Math.PI / 2)];
}
let seed = 20261113; const rnd = () => { seed = (seed * 16807) % 2147483647; return seed / 2147483647; };
const gauss = () => Math.sqrt(-2 * Math.log(rnd() || 1e-12)) * Math.cos(2 * Math.PI * rnd());
function write(name, fHi, n, L, noise) {
  let out = 'frequency (Hz),Re(Z) (ohm),Im(Z) (ohm)\n';
  for (let i = 0; i < n; i++) { const f = Math.pow(10, Math.log10(fHi) - i * (Math.log10(fHi) + 2) / (n - 1)), [re, im] = z(f, L), m = Math.hypot(re, im);
    out += `${f.toPrecision(7)},${(re + noise * m * gauss()).toPrecision(8)},${(im + noise * m * gauss()).toPrecision(8)}\n`; }
  fs.writeFileSync(path.join(__dirname, 'data', name), out);
}
write('synthetic_three_arcs.csv', 1e4, 61, 0, 0);
write('synthetic_inductance_noise.csv', 1e5, 71, 2e-7, 0.003);
