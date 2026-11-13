// cellpeek-eis — equivalent-circuit fitting of impedance spectra (the analysis code of the EIS ECM tab of CellPeek).
// Model (series):  [L] + R_HFI + ZARC_1 (contact) + ZARC_2 [+ ZARC_3] [+ diffusion tail]
//   ZARC_k = R_k / (1 + (jωτ_k)^α_k)      a depressed semicircle (R parallel to a constant-phase element), diameter R_k
//   tail   = W (jω)^(-p)                  p = 0.5 is a Warburg line (45°); p is fitted between 0.3 and 0.9
// Output:  R_HFI = series resistance (high-frequency intercept with the real axis)
//          R_ct  = R_2 (+ R_3): summed diameters of the semicircle(s) after the first one.
// The first (highest-frequency) semicircle is the contact impedance (active particles / current collector) and is not part of R_ct.
// Fit: complex non-linear least squares (Levenberg–Marquardt), modulus weighting, several starting points.
'use strict';
(function () {
const TWO_PI = 2 * Math.PI;

// ---------------- reading files ----------------
// Returns [{name, f[], zr[], zi[]}] sorted from high to low frequency; zi = Im(Z) (negative for capacitive behaviour).
function parseEIS(text, name) {
  const lines = String(text).replace(/\r/g, '').split('\n');
  const split = ln => { const t = ln.trim(); if (t.includes('\t')) return t.split('\t'); if (t.includes(';')) return t.split(';'); const c = t.split(','); return c.length >= 3 ? c : t.split(/\s+/); };
  const num = s => { s = String(s).trim(); if (s === '') return NaN; if (/^-?\d+,\d+([eE][-+]?\d+)?$/.test(s)) s = s.replace(',', '.'); const v = Number(s); return Number.isFinite(v) ? v : NaN; };
  const isNumRow = c => c.length >= 3 && c.filter(x => !isNaN(num(x))).length >= Math.min(3, c.length);
  // header: the last non-numeric line before the first run of numeric rows
  let start = -1;
  for (let i = 0; i < lines.length; i++) { const c = split(lines[i]); if (isNumRow(c) && i + 2 < lines.length && isNumRow(split(lines[i + 1])) && isNumRow(split(lines[i + 2]))) { start = i; break; } }
  if (start < 0) throw new Error('no table of numbers found');
  let hdr = null;
  for (let i = start - 1; i >= Math.max(0, start - 3); i--) { const c = split(lines[i]).map(x => x.trim().toLowerCase()); if (c.length >= 3 && c.some(x => /freq|^f$|^f[\/ (]|hz/.test(x))) { hdr = c; break; } }
  let cf = 0, cr = 1, ci = 2, neg = null;        // neg: the column holds −Im(Z)
  if (hdr) {
    const find = re => hdr.findIndex(x => re.test(x));
    const a = find(/freq|^f$|^f[\/ (]/), b = find(/re\(z\)|zreal|z_?re|^z'(?!')|z′(?!′)|^re[ (\/]|real/), c = find(/im\(z\)|zimag|z_?im|z''|z″|^-?im[ (\/]|imag/);
    if (a >= 0 && b >= 0 && c >= 0) { cf = a; cr = b; ci = c; neg = /^-|−|-im|-z''/.test(hdr[c]) || /-im\(z\)/.test(hdr[c]); }
  }
  const rows = [];
  for (let i = start; i < lines.length; i++) { const c = split(lines[i]); if (c.length <= Math.max(cf, cr, ci)) continue; const f = num(c[cf]), r = num(c[cr]), m = num(c[ci]); if (f > 0 && isFinite(r) && isFinite(m)) rows.push([f, r, m]); }
  if (rows.length < 6) throw new Error('fewer than 6 data points');
  if (neg === null) { const s = rows.map(r => r[2]).sort((a, b) => a - b); neg = s[s.length >> 1] > 0; }   // no header: a mostly positive column is −Im(Z)
  // several sweeps in one file: a new sweep starts when the frequency jumps back
  const dir = Math.sign(Math.log(rows[Math.min(5, rows.length - 1)][0] / rows[0][0])) || -1, sweeps = [[]];
  for (let i = 0; i < rows.length; i++) { if (i && Math.sign(Math.log(rows[i][0] / rows[i - 1][0])) === -dir && Math.abs(Math.log(rows[i][0] / rows[i - 1][0])) > 2) sweeps.push([]); sweeps[sweeps.length - 1].push(rows[i]); }
  const base = String(name || 'spectrum').replace(/\.[A-Za-z0-9]{1,5}$/, '');
  return sweeps.filter(s => s.length >= 6).map((s, k, all) => { s.sort((a, b) => b[0] - a[0]);
    return {name: all.length > 1 ? `${base} #${k + 1}` : base, f: s.map(r => r[0]), zr: s.map(r => r[1]), zi: s.map(r => neg ? -r[2] : r[2])}; });
}

// ---------------- model ----------------
const sig = u => 1 / (1 + Math.exp(-u)), logit = v => Math.log(v / (1 - v));
const A_LO = 0.5, A_HI = 1, P_LO = 0.3, P_HI = 0.9;
// cfg = {n: 2|3, tail: bool, L: bool};  q = transformed parameter vector
function unpack(q, cfg) {
  let i = 0; const o = {R0: Math.exp(q[i++]), L: cfg.L ? Math.exp(q[i++]) : 0, arcs: [], W: 0, p: 0.5};
  for (let k = 0; k < cfg.n; k++) o.arcs.push({R: Math.exp(q[i++]), tau: Math.exp(q[i++]), a: A_LO + (A_HI - A_LO) * sig(q[i++])});
  if (cfg.tail) { o.W = Math.exp(q[i++]); o.p = P_LO + (P_HI - P_LO) * sig(q[i++]); }
  return o;
}
function pack(o, cfg) {
  const q = [Math.log(o.R0)]; if (cfg.L) q.push(Math.log(o.L));
  for (const a of o.arcs) q.push(Math.log(a.R), Math.log(a.tau), logit((Math.min(A_HI - 1e-3, Math.max(A_LO + 1e-3, a.a)) - A_LO) / (A_HI - A_LO)));
  if (cfg.tail) q.push(Math.log(o.W), logit((o.p - P_LO) / (P_HI - P_LO)));
  return q;
}
function zarc(R, tau, a, w) { const m = Math.pow(w * tau, a), c = 1 + m * Math.cos(a * Math.PI / 2), s = m * Math.sin(a * Math.PI / 2), d = c * c + s * s; return [R * c / d, -R * s / d]; }
function zModel(o, f) {                          // -> [re, im]
  const w = TWO_PI * f; let re = o.R0, im = w * o.L;
  for (const a of o.arcs) { const z = zarc(a.R, a.tau, a.a, w); re += z[0]; im += z[1]; }
  if (o.W) { const m = o.W * Math.pow(w, -o.p); re += m * Math.cos(o.p * Math.PI / 2); im -= m * Math.sin(o.p * Math.PI / 2); }
  return [re, im];
}

// ---------------- Levenberg–Marquardt ----------------
function solve(A, b) {
  const n = b.length, M = A.map((r, i) => r.concat(b[i]));
  for (let c = 0; c < n; c++) { let p = c; for (let r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[p][c])) p = r;
    if (Math.abs(M[p][c]) < 1e-300) return null; [M[c], M[p]] = [M[p], M[c]];
    for (let r = c + 1; r < n; r++) { const k = M[r][c] / M[c][c]; if (k) for (let j = c; j <= n; j++) M[r][j] -= k * M[c][j]; } }
  const x = new Array(n).fill(0);
  for (let r = n - 1; r >= 0; r--) { let s = M[r][n]; for (let j = r + 1; j < n; j++) s -= M[r][j] * x[j]; x[r] = s / M[r][r]; }
  return x.every(Number.isFinite) ? x : null;
}
function lm(res, q0, maxIt) {
  let q = q0.slice(), r = res(q), cost = r.reduce((s, v) => s + v * v, 0), lam = 1e-2; const n = q.length;
  for (let it = 0; it < maxIt; it++) {
    const J = []; for (let j = 0; j < n; j++) { const h = 1e-5 * Math.max(1, Math.abs(q[j])), q2 = q.slice(); q2[j] += h; const r2 = res(q2); J.push(r2.map((v, i) => (v - r[i]) / h)); }
    const A = [], g = new Array(n).fill(0);
    for (let a = 0; a < n; a++) { A.push(new Array(n).fill(0)); for (let i = 0; i < r.length; i++) g[a] += J[a][i] * r[i]; }
    for (let a = 0; a < n; a++) for (let b = a; b < n; b++) { let s = 0; for (let i = 0; i < r.length; i++) s += J[a][i] * J[b][i]; A[a][b] = A[b][a] = s; }
    let ok = false;
    for (let t = 0; t < 12; t++) {
      const B = A.map((row, i) => row.map((v, j) => i === j ? v * (1 + lam) + 1e-12 : v)), d = solve(B, g.map(v => -v));
      if (!d) { lam *= 10; continue; }
      const q2 = q.map((v, i) => v + Math.max(-3, Math.min(3, d[i]))), r2 = res(q2), c2 = r2.reduce((s, v) => s + v * v, 0);
      if (isFinite(c2) && c2 < cost) { const rel = (cost - c2) / (cost || 1); q = q2; r = r2; cost = c2; lam = Math.max(lam / 4, 1e-10); ok = true; if (rel < 1e-10) it = maxIt; break; }
      lam *= 6;
    }
    if (!ok) break;
  }
  return {q, cost};
}

// ---------------- shape analysis: tail, starting values ----------------
function shape(d) {
  const n = d.f.length, y = d.zi.map(v => -v);
  const ys = y.map((v, i) => { const a = [y[Math.max(0, i - 1)], v, y[Math.min(n - 1, i + 1)]].sort((p, q) => p - q); return a[1]; });   // 3-point median
  const ymax = Math.max(...ys), iInd = ys.findIndex(v => v > 0);                  // points before iInd are inductive (−Z″ < 0)
  // tail: at the low-frequency end −Z″ rises again after a minimum that follows a semicircle apex
  let run = n - 1; while (run > 0 && ys[run] >= ys[run - 1] - 1e-4 * ymax) run--;
  let apex = -1; for (let i = Math.max(1, iInd); i < run; i++) if (ys[i] >= ys[i - 1] && ys[i] >= ys[i + 1] && (apex < 0 || ys[i] > ys[apex])) apex = i;
  const tail = (n - 1 - run) >= 3 && apex >= 0 && ys[n - 1] - ys[run] > 0.06 * ymax && ys[run] < 0.92 * ys[apex];
  return {y, ys, ymax, iInd: Math.max(0, iInd), tail, iMin: tail ? run : n - 1, inductive: y.filter(v => v < 0).length >= 2};
}
function starts(d, sh, cfg) {
  const n = d.f.length, i0 = sh.iInd, i1 = cfg.tail ? sh.iMin : n - 1;
  const zmin = Math.min(...d.zr.slice(i0)), R0s = [Math.max(zmin * 0.97, 1e-9), Math.max(zmin - 0.5 * Math.max(0, sh.y[i0]), zmin * 0.5, 1e-9)];
  const Rtot = Math.max((cfg.tail ? d.zr[i1] : d.zr[n - 1] + 0.6 * Math.max(0, sh.y[n - 1])) - zmin, 1e-6 * Math.abs(zmin) + 1e-9);
  const lf = k => Math.log(d.f[Math.round(i0 + (i1 - i0) * k)]), tauAt = k => 1 / (TWO_PI * Math.exp(lf(k)));
  const peaks = []; for (let i = i0 + 1; i < i1; i++) if (sh.ys[i] > sh.ys[i - 1] && sh.ys[i] >= sh.ys[i + 1]) peaks.push(1 / (TWO_PI * d.f[i]));
  const N = cfg.n, tauSets = N === 2 ? [[tauAt(.25), tauAt(.7)], [tauAt(.1), tauAt(.55)], [tauAt(.4), tauAt(.85)]] : [[tauAt(.15), tauAt(.5), tauAt(.82)], [tauAt(.08), tauAt(.35), tauAt(.7)], [tauAt(.3), tauAt(.6), tauAt(.9)]];
  if (peaks.length >= 1) { const p = peaks.slice().sort((a, b) => a - b), s = []; for (let k = 0; k < N; k++) s.push(p[Math.round(k * (p.length - 1) / Math.max(1, N - 1))]);
    for (let k = 1; k < N; k++) if (s[k] <= s[k - 1] * 2) s[k] = s[k - 1] * 6; if (p.length === 1) { s[N - 1] = p[0]; for (let k = N - 2; k >= 0; k--) s[k] = s[k + 1] / 12; } tauSets.push(s); }
  const splits = N === 2 ? [[.15, .85], [.4, .6], [.7, .3]] : [[.1, .3, .6], [.2, .4, .4], [.1, .6, .3], [.34, .33, .33]];
  const wl = TWO_PI * d.f[n - 1], W0 = Math.max(sh.y[n - 1] - sh.y[sh.iMin], 0.05 * sh.ymax) * Math.sqrt(wl) / Math.SQRT1_2, L0 = sh.inductive ? Math.max(-sh.y[0], 1e-12) / (TWO_PI * d.f[0]) : 0;
  const out = [];
  for (const ts of tauSets) for (const sp of splits) out.push(pack({R0: R0s[out.length % 2], L: L0 || 1e-9, arcs: ts.map((t, k) => ({R: Rtot * sp[k], tau: t, a: 0.85})), W: W0, p: 0.5}, cfg));
  return out;
}

// ---------------- fitting one spectrum ----------------
function fitCfg(d, sh, cfg) {
  const n = d.f.length, mod = d.f.map((_, i) => Math.hypot(d.zr[i], d.zi[i]) || 1);
  const tLo = Math.log(0.05 / (TWO_PI * d.f[0])), tHi = Math.log((cfg.tail ? 1 : 8) / (TWO_PI * d.f[n - 1]));   // time constants stay inside the measured window
  const res = q => { const o = unpack(q, cfg), r = new Array(2 * n);
    for (let i = 0; i < n; i++) { const z = zModel(o, d.f[i]); r[2 * i] = (z[0] - d.zr[i]) / mod[i]; r[2 * i + 1] = (z[1] - d.zi[i]) / mod[i]; }
    for (const a of o.arcs) { const lt = Math.log(a.tau); r.push(lt < tLo ? 3 * (tLo - lt) : lt > tHi ? 3 * (lt - tHi) : 0); }
    return r; };
  let best = null;
  for (const q0 of starts(d, sh, cfg)) { const r = lm(res, q0, 40); if (isFinite(r.cost) && (!best || r.cost < best.cost)) best = r; }
  if (!best) return null;
  best = lm(res, best.q, 300);
  const o = unpack(best.q, cfg); o.arcs.sort((a, b) => a.tau - b.tau);          // arc 1 = highest frequency = contact
  let ss = 0; for (let i = 0; i < n; i++) { const z = zModel(o, d.f[i]); ss += ((z[0] - d.zr[i]) ** 2 + (z[1] - d.zi[i]) ** 2) / (mod[i] * mod[i]); }
  const k = best.q.length, N2 = 2 * n;
  return {cfg, o, cost: ss, err: 100 * Math.sqrt(ss / n), bic: N2 * Math.log(ss / N2 + 1e-300) + k * Math.log(N2)};
}
// opt = {n: 'auto'|2|3, tail: 'auto'|'on'|'off'}
// Inductive points at the high-frequency end (−Z″ < 0: cables, cell hardware), two ways, chosen with opt.inductive:
//   'fit'  (default) — all measured points are used and an inductor is in the model
//   'skip'           — these points are left out of the fit and the model has no inductor (nSkip = number of points left out).
//                      Smaller fit error, because there Z′ changes with frequency too, which an ideal inductor cannot follow; R_HFI and R_contact change with it.
const FITV = 3;                                  // version of the fitting method (stored fits of an older version are redone)
function fitSpectrum(d0, opt = {}) {
  const k0 = opt.inductive === 'skip' ? shape(d0).iInd : 0, nSkip = k0 > 0 && d0.f.length - k0 >= 10 ? k0 : 0;
  const d = nSkip ? {f: d0.f.slice(nSkip), zr: d0.zr.slice(nSkip), zi: d0.zi.slice(nSkip)} : d0;
  const sh = shape(d), tail = opt.tail === 'on' ? true : opt.tail === 'off' ? false : sh.tail, L = nSkip ? false : sh.inductive;
  const run = n => fitCfg(d, sh, {n, tail, L});
  let fit, note = '';
  if (opt.n === 2 || opt.n === 3) fit = run(opt.n);
  else {                                         // auto: three semicircles only when the data clearly ask for a third one
    const f2 = run(2), f3 = run(3);
    const sum3 = f3 ? f3.o.arcs.reduce((s, a) => s + a.R, 0) : 0;
    const ok3 = f3 && f2 && f3.bic < f2.bic - 10 && f3.err < 0.8 * f2.err && f3.o.arcs.every(a => a.R > 0.03 * sum3) && f3.o.arcs.every((a, i, A) => !i || a.tau > 3 * A[i - 1].tau);
    fit = ok3 ? f3 : f2 || f3;
  }
  if (!fit) return null;
  const A = fit.o.arcs, Rct = A.slice(1).reduce((s, a) => s + a.R, 0);
  if (L) note = 'inductive points at high frequency: an inductor is included';
  // error shown and used for the "fits poorly" mark: the points with −Z″ ≥ 0 only. The inductive points are fitted too, but there Z′ also
  // changes with frequency, which an ideal inductor cannot follow; their misfit says nothing about the semicircles. errAll = all points.
  let ss = 0, nb = 0; for (let i = sh.iInd; i < d.f.length; i++) { const z = zModel(fit.o, d.f[i]), m2 = d.zr[i] * d.zr[i] + d.zi[i] * d.zi[i] || 1; ss += ((z[0] - d.zr[i]) ** 2 + (z[1] - d.zi[i]) ** 2) / m2; nb++; }
  const errBody = nb ? 100 * Math.sqrt(ss / nb) : fit.err;
  return {v: FITV, n: A.length, tail, tailAuto: sh.tail, L, nSkip, R_HFI: fit.o.R0, R_contact: A[0].R, R_ct: Rct, R_tot: A[0].R + Rct, err: errBody, errAll: fit.err, model: fit.o, note};
}

const API = {parseEIS, fitSpectrum, zModel, zarc, shape, FITV};
if (typeof module !== 'undefined' && module.exports) module.exports = API;
else if (typeof window !== 'undefined') window.CellPeekEIS = API;
})();
