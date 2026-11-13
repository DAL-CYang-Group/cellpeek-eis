# cellpeek-eis

Equivalent-circuit fitting of electrochemical impedance spectra (EIS) of Li-ion cells. For every spectrum it returns the
series resistance, the contact resistance and the charge-transfer resistance.

This is the analysis code of the "EIS ECM" tab of CellPeek (https://dal-cyang-group.github.io/cellpeek/), without the
web page. It has no dependencies and runs with Node.js 18 or newer.

## Quick start

```bash
node bin/cli.js examples/data/*.csv
```

```bash
npm test
```

The first command prints a CSV table. `npm test` fits the example spectra, compares the result with the table in
`examples/expected/`, and checks it against the known answer of the examples.

## Input files

Text files with three columns: frequency (Hz), Z′ (Ω) and Z″ (Ω). Comma, semicolon, tab or spaces as separator; a header
line is optional. Without a header the third column is taken as Im(Z) when it is mostly negative and as −Im(Z) when it is
mostly positive. Files with a header of BioLogic (`.mpt`) or Gamry (`.DTA`) exports are read as well; several sweeps in
one file are split.

## Equivalent circuit

```
Z(ω) = jωL + R_HFI + Σ_k R_k / (1 + (jωτ_k)^α_k) + W (jω)^(−p)        k = 1 … 3 (or 2)
```

| Element | Meaning |
|---|---|
| `L` | inductance of cables and cell hardware; only in the model when the spectrum has points with −Z″ < 0 |
| `R_HFI` | series resistance: high-frequency intercept with the real axis |
| `R_k, τ_k, α_k` | depressed semicircles (a resistor parallel to a constant-phase element), 0.5 ≤ α ≤ 1 |
| `W, p` | diffusion tail; only when −Z″ rises again at low frequency. p is fitted between 0.3 and 0.9 (0.5 = Warburg) |

Reported quantities:

| Quantity | Definition |
|---|---|
| `R_HFI` | as above |
| `R_contact` | diameter of the semicircle with the shortest time constant |
| `R_ct` | summed diameters of the other semicircle(s) |
| width | `R_contact + R_ct`, the whole span of the semicircles |

One model is used for all spectra of a series: three semicircles by default, two with `--semicircles 2`.

Points with −Z″ < 0 at the high-frequency end (inductive) can be treated in two ways:

| Option | What happens |
|---|---|
| `--inductive fit` (default) | all points are fitted, with the inductor `L` in the model |
| `--inductive skip` | these points are left out of the fit and the model has no inductor |

Measured spectra often show a real part that also changes with frequency in this region, which an ideal inductor cannot
follow. Leaving the points out then gives a smaller fit error, and somewhat different `R_HFI` and `R_contact`. Use the same
option for all spectra that are compared, and state it.

## Fit

Complex non-linear least squares (Levenberg–Marquardt) on the real and imaginary parts together, each point weighted by
1/|Z|, started from several initial guesses; the best result is refined. No random numbers are used, so the same file
always gives the same result.

Two errors are given, both as root-mean-square relative deviation between fit and data: over the points with −Z″ ≥ 0,
and over all points. They differ when the spectrum has inductive points at high frequency.

## How far to trust the split into R_contact and R_ct

When the semicircles overlap, `R_contact` and `R_ct` are not determined well separately, while their sum is. The second
example shows this: with 0.3 % noise the fit returns R_contact 0.043 and R_ct 0.257 Ω (true values 0.050 and 0.250 Ω),
but a width of 0.2999 Ω (true value 0.300 Ω). A model with too few semicircles moves the split further: with
`--semicircles 2` the same data give R_contact 0.097 and R_ct 0.211 Ω.

## Examples

Both spectra are computed (`examples/make_synthetic.js`) from the circuit above with R_HFI 0.200 Ω, three semicircles of
0.050, 0.080 and 0.170 Ω, and a Warburg tail, so R_contact = 0.050, R_ct = 0.250 and width = 0.300 Ω.

| File | Spectrum | Fitted R_HFI / R_contact / R_ct / width (Ω) |
|---|---|---|
| `synthetic_three_arcs.csv` | 10 kHz – 10 mHz, no noise | 0.2000 / 0.0500 / 0.2500 / 0.3000 |
| `synthetic_inductance_noise.csv` | 100 kHz – 10 mHz, inductance 0.2 µH, 0.3 % noise | 0.1996 / 0.0428 / 0.2571 / 0.2999 |

`examples/expected/ecm_3_semicircles.csv` is the full output for both.

## Use from JavaScript

```js
const EIS = require('./src/eis.js');
const [d] = EIS.parseEIS(text, 'name.csv');                 // {name, f[], zr[], zi[]}, high to low frequency
const fit = EIS.fitSpectrum(d, {n: 3, tail: 'auto', inductive: 'fit'});   // {R_HFI, R_contact, R_ct, R_tot, err, errAll, nSkip, model, …}
```

In a browser the file defines `window.CellPeekEIS`.

## Version

1.0.0 — taken from CellPeek commit `e0aa4a6` (2026-11-13), fitting method version 3.

## Licence

MIT, see `LICENSE`. References for the circuit and the definition of R_ct: `NOTICE.md`.
