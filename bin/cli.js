#!/usr/bin/env node
// Command line: equivalent-circuit fit of impedance spectra. Prints a CSV table.
//   node bin/cli.js [--semicircles 3] [--inductive fit|skip] file1 [file2 ...]
'use strict';
const fs = require('fs'), path = require('path');
const EIS = require('../src/eis.js');

function usage(code) {
  console.error(`Usage:
  node bin/cli.js [--semicircles 2|3] [--inductive fit|skip] FILE...
                    R_HFI, R_contact, R_ct and semicircle width of every spectrum
Options
  --semicircles N   number of depressed semicircles in the circuit (default 3)
  --inductive fit   points with -Z'' < 0 at high frequency are fitted, with an inductor in the model (default)
  --inductive skip  these points are left out and the model has no inductor
Files: text with frequency (Hz), Z' and Z'' (or -Z'') columns; a header line is optional.`);
  process.exit(code);
}
const args = process.argv.slice(2);
if (!args.length || args[0] === '-h' || args[0] === '--help') usage(args.length ? 0 : 1);
const opt = {semicircles: 3, inductive: 'fit'}, files = [];
while (args.length) { const a = args.shift();
  if (a.startsWith('--')) { const k = a.slice(2); if (!(k in opt) || !args.length) usage(1);
    if (k === 'inductive') { opt[k] = args.shift(); if (opt[k] !== 'fit' && opt[k] !== 'skip') usage(1); }
    else { opt[k] = Number(args.shift()); if (!isFinite(opt[k])) usage(1); } }
  else files.push(a); }
if (!files.length) usage(1);

const n6 = v => v == null || !isFinite(v) ? '' : String(+(+v).toPrecision(6));
const csv = rows => rows.map(r => r.map(c => /[",\n]/.test(String(c)) ? '"' + String(c).replace(/"/g, '""') + '"' : c).join(',')).join('\n') + '\n';
const rows = [['spectrum', 'R_HFI (ohm)', 'R_contact (ohm)', 'R_ct (ohm)', 'semicircle width = R_contact + R_ct (ohm)', 'semicircles', 'inductive points', 'diffusion tail', 'tail exponent',
  "fit error (%), points with -Z'' >= 0", 'fit error (%), all points']];
for (const file of files) for (const d of EIS.parseEIS(fs.readFileSync(file, 'utf8'), path.basename(file))) {
  const f = EIS.fitSpectrum(d, {n: opt.semicircles === 2 ? 2 : 3, tail: 'auto', inductive: opt.inductive});
  rows.push(f ? [d.name, n6(f.R_HFI), n6(f.R_contact), n6(f.R_ct), n6(f.R_tot), f.n, f.nSkip ? f.nSkip + ' left out' : f.L ? 'fitted with an inductor' : 'none', f.tail ? 'yes' : 'no', f.tail ? n6(f.model.p) : '', f.err.toFixed(2), f.errAll.toFixed(2)]
    : [d.name, '', '', '', '', opt.semicircles, '', '', '', 'failed', '']);
}
process.stdout.write(csv(rows));
