// Runs the command on the example spectra and compares every number with the stored table (relative tolerance 1e-4),
// then checks the fitted values against the known answer of the synthetic spectra.
'use strict';
const {execFileSync} = require('child_process'), fs = require('fs'), path = require('path');
const root = path.join(__dirname, '..'), dir = path.join(root, 'examples', 'data');
const data = fs.readdirSync(dir).filter(f => f.endsWith('.csv')).sort().map(f => path.join('examples', 'data', f));
let bad = 0;
{ const got = execFileSync(process.execPath, [path.join('bin', 'cli.js'), ...data], {cwd: root, encoding: 'utf8'});
  const want = fs.readFileSync(path.join(root, 'examples', 'expected', 'ecm_3_semicircles.csv'), 'utf8');
  const tok = s => s.split(/[\n,]+/).filter(Boolean), a = tok(got), b = tok(want);
  let ok = a.length === b.length;
  for (let i = 0; ok && i < a.length; i++) { const x = Number(a[i]), y = Number(b[i]);
    ok = isFinite(x) && isFinite(y) ? Math.abs(x - y) <= 1e-4 * Math.max(Math.abs(x), Math.abs(y), 1e-12) : a[i] === b[i]; }
  console.log(`${ok ? 'ok  ' : 'FAIL'} stored table: ${data.length} spectra`); if (!ok) bad++; }
const EIS = require('../src/eis.js');
// exact data: every value within 0.1 %.  Noisy data with inductance: R_HFI and the width within 5 % (the split is not checked, see README)
for (const [file, tol, split] of [['synthetic_three_arcs.csv', 1e-3, true], ['synthetic_inductance_noise.csv', 0.05, false]]) {
  const d = EIS.parseEIS(fs.readFileSync(path.join(dir, file), 'utf8'), file)[0], f = EIS.fitSpectrum(d, {n: 3, tail: 'auto'});
  const close = (v, t) => Math.abs(v - t) <= tol * t, ok = close(f.R_HFI, 0.2) && close(f.R_tot, 0.3) && (!split || (close(f.R_contact, 0.05) && close(f.R_ct, 0.25)));
  console.log(`${ok ? 'ok  ' : 'FAIL'} ${file}: R_HFI ${f.R_HFI.toFixed(4)} (0.2000), R_contact ${f.R_contact.toFixed(4)} (0.0500), R_ct ${f.R_ct.toFixed(4)} (0.2500), width ${f.R_tot.toFixed(4)} (0.3000)`);
  if (!ok) bad++;
}
process.exit(bad ? 1 : 0);
