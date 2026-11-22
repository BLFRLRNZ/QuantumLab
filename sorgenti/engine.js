// MOTORE — Livello 0. Unico dato in ingresso: Z (numero di protoni).
// Leggi impostate: equazione di Schrodinger per gli elettroni (Kohn-Sham, unita atomiche),
// attrazione di Coulomb del nucleo, repulsione di Coulomb tra elettroni (Hartree),
// principio di Pauli (max 2l+1 elettroni per sottolivello e spin),
// scambio e correlazione LSDA: scambio esatto del gas di elettroni + correlazione PW92
// (formula adattata a simulazioni Monte Carlo quantistiche del gas di elettroni, non a esperimenti).
// Non contiene: tabelle di configurazioni, regola di Madelung, valori misurati, relativita.

var HA2EV = 27.211386245988;
var LETTERS = ['s', 'p', 'd', 'f'];

function makeGrid(Z, dx, rmax) {
  dx = dx || 0.005; rmax = rmax || 60;
  var xmin = -8.0, xmax = Math.log(Z * rmax);
  var mesh = Math.ceil((xmax - xmin) / dx);
  var r = new Float64Array(mesh + 1), r2 = new Float64Array(mesh + 1), sqr = new Float64Array(mesh + 1), w = new Float64Array(mesh + 1);
  for (var i = 0; i <= mesh; i++) {
    r[i] = Math.exp(xmin + i * dx) / Z; r2[i] = r[i] * r[i]; sqr[i] = Math.sqrt(r[i]);
    w[i] = 4 * Math.PI * r[i] * r2[i] * dx; // elemento di volume 4 pi r^2 dr con dr = r dx
  }
  w[0] *= 0.5; w[mesh] *= 0.5;
  return { Z: Z, dx: dx, mesh: mesh, r: r, r2: r2, sqr: sqr, w: w };
}

// Equazione radiale di Schrodinger (Numerov su griglia logaritmica, metodo shooting con correzione perturbativa).
function solveRadial(g, n, l, V, eGuess, Zn, noRetry) {
  var sol = solveRadial1(g, n, l, V, eGuess, Zn);
  if (sol || noRetry) return sol;
  // ripartenze da energie diverse (utile per gli orbitali f, che vivono in una doppia buca di potenziale)
  var tries = [-0.02, -0.08, -0.2, -0.5, -1, -2, -4, -8];
  for (var t = 0; t < tries.length; t++) {
    if (eGuess !== undefined && Math.abs(tries[t] - eGuess) < 1e-3) continue;
    sol = solveRadial1(g, n, l, V, tries[t], Zn);
    if (sol) return sol;
  }
  return null;
}
function solveRadial1(g, n, l, V, eGuess, Zn) {
  var mesh = g.mesh, r = g.r, r2 = g.r2, sqr = g.sqr, dx = g.dx;
  var ddx12 = dx * dx / 12, sqlhf = (l + 0.5) * (l + 0.5), x2l2 = 2 * l + 2;
  var nodes = n - l - 1;
  var eup = V[mesh], elw = 1e30;
  for (var i = 0; i <= mesh; i++) { var t = V[i] + sqlhf / (2 * r2[i]); if (t < elw) elw = t; }
  if (eup - elw < 1e-12) return null;
  var e = (eGuess !== undefined && eGuess > elw && eGuess < eup) ? eGuess : 0.5 * (elw + eup);
  var f = new Float64Array(mesh + 1), y = new Float64Array(mesh + 1);
  for (var iter = 0; iter < 300; iter++) {
    var icl = -1;
    f[0] = ddx12 * (sqlhf + 2 * r2[0] * (V[0] - e));
    for (i = 1; i <= mesh; i++) {
      f[i] = ddx12 * (sqlhf + 2 * r2[i] * (V[i] - e));
      if (f[i] === 0) f[i] = 1e-20;
      if ((f[i] > 0) !== (f[i - 1] > 0)) icl = i;
    }
    if (icl < 0 || icl >= mesh - 2) { eup = e; e = 0.5 * (eup + elw); if (eup - elw < 1e-12) return null; continue; }
    // punto di partenza dell'integrazione verso l'interno: dove la funzione e' ormai trascurabile
    var acc = 0, istart = mesh - 1;
    for (i = icl; i < mesh - 1; i++) { var gg = f[i] / ddx12; if (gg > 0) acc += Math.sqrt(gg) * dx; if (acc > 45 && i > icl + 2) { istart = i; break; } }
    for (i = 0; i <= mesh; i++) f[i] = 1 - f[i];
    for (i = 0; i <= mesh; i++) y[i] = 0;
    y[0] = Math.pow(r[0], l + 1) * (1 - 2 * Zn * r[0] / x2l2) / sqr[0];
    y[1] = Math.pow(r[1], l + 1) * (1 - 2 * Zn * r[1] / x2l2) / sqr[1];
    var ncross = 0;
    for (i = 1; i <= icl - 1; i++) {
      y[i + 1] = ((12 - 10 * f[i]) * y[i] - f[i - 1] * y[i - 1]) / f[i + 1];
      if ((y[i] > 0) !== (y[i + 1] > 0) && y[i] !== 0) ncross++;
    }
    var fac = y[icl];
    if (ncross !== nodes) {
      if (ncross > nodes) eup = e; else elw = e;
      e = 0.5 * (eup + elw);
      if (eup - elw < 1e-12) return null;
      continue;
    }
    y[istart + 1] = 0; y[istart] = dx;
    y[istart - 1] = (12 - 10 * f[istart]) * y[istart] / f[istart - 1];
    for (i = istart - 1; i >= icl + 1; i--) {
      y[i - 1] = ((12 - 10 * f[i]) * y[i] - f[i + 1] * y[i + 1]) / f[i - 1];
      if (Math.abs(y[i - 1]) > 1e100) { for (var k = i - 1; k <= istart; k++) y[k] *= 1e-100; }
    }
    fac = fac / y[icl];
    for (i = icl; i <= istart; i++) y[i] *= fac;
    var norm = 0;
    for (i = 0; i <= mesh; i++) norm += y[i] * y[i] * r2[i] * dx;
    norm = Math.sqrt(norm);
    for (i = 0; i <= mesh; i++) y[i] /= norm;
    i = icl;
    var ycusp = (y[i - 1] * f[i - 1] + f[i + 1] * y[i + 1] + 10 * f[i] * y[i]) / 12;
    var dfcusp = f[i] * (y[i] / ycusp - 1);
    var de = 0.5 * dfcusp / ddx12 * ycusp * ycusp * dx;
    if (de > 0) elw = e;
    if (de < 0) eup = e;
    var tol = 1e-10 * Math.max(1, Math.abs(e));
    if (Math.abs(de) < tol || eup - elw < tol) {
      // accettato solo se e' davvero uno stato legato e la derivata si raccorda
      if (Math.abs(dfcusp) > 1e-6 || e > V[mesh] - 1e-5) return null;
      return { e: e, y: y };
    }
    e = Math.max(Math.min(e + de, eup), elw);
  }
  return null;
}

// Scambio-correlazione LSDA (per unita di volume) e potenziali per spin.
var CX = Math.pow(6 / Math.PI, 1 / 3), FPP0 = 1.709921, D43 = Math.pow(2, 4 / 3) - 2;
function Gpw(rs, A, a1, b1, b2, b3, b4) {
  var srs = Math.sqrt(rs);
  var q0 = -2 * A * (1 + a1 * rs);
  var q1 = 2 * A * (b1 * srs + b2 * rs + b3 * rs * srs + b4 * rs * rs);
  var q1p = A * (b1 / srs + 2 * b2 + 3 * b3 * srs + 4 * b4 * rs);
  var lg = Math.log(1 + 1 / q1);
  return [q0 * lg, -2 * A * a1 * lg - q0 * q1p / (q1 * q1 + q1)];
}
function xc(nu, nd, out) {
  var n = nu + nd;
  if (n < 1e-30) { out[0] = 0; out[1] = 0; out[2] = 0; return; }
  if (nu < 0) nu = 0; if (nd < 0) nd = 0;
  var cu = Math.cbrt(nu), cd = Math.cbrt(nd);
  var ex = -0.75 * CX * (nu * cu + nd * cd);
  var vxu = -CX * cu, vxd = -CX * cd;
  var rs = Math.cbrt(3 / (4 * Math.PI * n));
  var z = (nu - nd) / n; if (z > 1) z = 1; if (z < -1) z = -1;
  var c0 = Gpw(rs, 0.031091, 0.21370, 7.5957, 3.5876, 1.6382, 0.49294);
  var c1 = Gpw(rs, 0.015545, 0.20548, 14.1189, 6.1977, 3.3662, 0.62517);
  var ca = Gpw(rs, 0.016887, 0.11125, 10.357, 3.6231, 0.88026, 0.49671);
  var ac = -ca[0], dac = -ca[1];
  var p1 = Math.cbrt(1 + z), m1 = Math.cbrt(1 - z);
  var fz = ((1 + z) * p1 + (1 - z) * m1 - 2) / D43;
  var dfz = (4 / 3) * (p1 - m1) / D43;
  var z3 = z * z * z, z4 = z3 * z;
  var ec = c0[0] + ac * fz / FPP0 * (1 - z4) + (c1[0] - c0[0]) * fz * z4;
  var decr = c0[1] * (1 - fz * z4) + c1[1] * fz * z4 + dac * fz / FPP0 * (1 - z4);
  var decz = 4 * z3 * fz * (c1[0] - c0[0] - ac / FPP0) + dfz * (z4 * (c1[0] - c0[0]) + (1 - z4) * ac / FPP0);
  var base = ec - rs / 3 * decr;
  out[0] = ex + n * ec;
  out[1] = vxu + base - (z - 1) * decz;
  out[2] = vxd + base - (z + 1) * decz;
}

// Potenziale di partenza: schermo statistico di Thomas-Fermi (approssimazione di Tietz) per Z-1 elettroni,
// piu' la carica netta +1 vista da lontano (legge di Gauss). Solo da Z.
function startPotential(g) {
  var Z = g.Z, b = 0.8853 * Math.pow(Z, -1 / 3), V = new Float64Array(g.mesh + 1);
  for (var i = 0; i <= g.mesh; i++) { var x = g.r[i] / b, phi = 1 / Math.pow(1 + 0.53625 * x, 2); V[i] = -(1 + (Z - 1) * phi) / g.r[i]; }
  return V;
}

function orbKey(o) { return o.n + ',' + o.l + ',' + o.s; }
function energiesOf(res) { var m = {}; (res.orbs || []).forEach(function (o) { m[orbKey(o)] = o.e; }); return m; }

// Calcolo autoconsistente per una configurazione fissata (occupazioni intere per (n,l,spin)).
function scf(g, occs, Vstart, opts) {
  opts = opts || {};
  var mesh = g.mesh, Z = g.Z, r = g.r, w = g.w;
  var N = 0; for (var k = 0; k < occs.length; k++) N += occs[k].occ;
  if (N === 0) return { E: 0, V: [Vstart[0].slice(), Vstart[1].slice()], orbs: [], it: 0, ok: true };
  var Vin = [Float64Array.from(Vstart[0]), Float64Array.from(Vstart[1])];
  var nu = new Float64Array(mesh + 1), nd = new Float64Array(mesh + 1);
  var eprev = Object.assign({}, opts.eg || {}), Eold = 0, out = [0, 0, 0];
  var VH = new Float64Array(mesh + 1), Vout = [new Float64Array(mesh + 1), new Float64Array(mesh + 1)];
  var ex = new Float64Array(mesh + 1);
  var alpha = opts.alpha || 0.35, maxit = opts.maxit || 200;
  var Fold = null, Vold = null, good = null, backs = 0;
  for (var it = 1; it <= maxit; it++) {
    nu.fill(0); nd.fill(0);
    var sumeps = 0, orbs = [], failed = false;
    for (k = 0; k < occs.length; k++) {
      var o = occs[k]; if (o.occ <= 0) continue;
      var key = orbKey(o);
      var sol = solveRadial(g, o.n, o.l, Vin[o.s], eprev[key], Z);
      if (!sol) { failed = true; break; }
      eprev[key] = sol.e;
      sumeps += o.occ * sol.e;
      var dens = o.s === 0 ? nu : nd, y = sol.y;
      for (var i = 0; i <= mesh; i++) dens[i] += o.occ * y[i] * y[i] / (4 * Math.PI * r[i]);
      orbs.push({ n: o.n, l: o.l, s: o.s, occ: o.occ, e: sol.e, y: sol.y });
    }
    if (failed) {
      // un orbitale occupato non e' legato in questo passo intermedio: si torna indietro e si procede piu' piano
      if (!good || backs > 8) return { ok: false };
      backs++;
      for (var s0 = 0; s0 < 2; s0++) for (i = 0; i <= mesh; i++) Vin[s0][i] = good[s0][i] + 0.25 * (Vin[s0][i] - good[s0][i]);
      alpha *= 0.6; Fold = null; Vold = null;
      continue;
    }
    good = [Float64Array.from(Vin[0]), Float64Array.from(Vin[1])];
    // Hartree: V_H(r) = Q(r)/r + integrale esterno
    var q = 0, qin = new Float64Array(mesh + 1);
    var prev = 0;
    for (i = 0; i <= mesh; i++) {
      var cur = 4 * Math.PI * (nu[i] + nd[i]) * r[i] * r[i] * r[i];
      if (i === 0) q = cur / 3; else q += 0.5 * (cur + prev) * g.dx;
      qin[i] = q; prev = cur;
    }
    var o2 = 0; prev = 0;
    for (i = mesh; i >= 0; i--) {
      var cur2 = 4 * Math.PI * (nu[i] + nd[i]) * r[i] * r[i];
      if (i < mesh) o2 += 0.5 * (cur2 + prev) * g.dx;
      prev = cur2;
      VH[i] = qin[i] / r[i] + o2;
    }
    var Exc = 0, EH = 0, Eext = 0, vinn = 0;
    for (i = 0; i <= mesh; i++) {
      xc(nu[i], nd[i], out);
      var ntot = nu[i] + nd[i];
      Exc += out[0] * w[i];
      EH += 0.5 * VH[i] * ntot * w[i];
      Eext += (-Z / r[i]) * ntot * w[i];
      vinn += (Vin[0][i] * nu[i] + Vin[1][i] * nd[i]) * w[i];
      Vout[0][i] = -Z / r[i] + VH[i] + out[1];
      Vout[1][i] = -Z / r[i] + VH[i] + out[2];
    }
    var Ts = sumeps - vinn;
    var E = Ts + Eext + EH + Exc;
    // residuo
    var res = 0;
    for (var s = 0; s < 2; s++) for (i = 0; i <= mesh; i++) { var d = Vout[s][i] - Vin[s][i]; res += d * d * (nu[i] + nd[i] + 1e-12) * w[i]; }
    res = Math.sqrt(res);
    if (it > 1 && Math.abs(E - Eold) < 1e-9 && res < 1e-6) {
      return { ok: true, E: E, V: Vin, orbs: orbs, it: it, nu: nu, nd: nd, Ts: Ts, Eext: Eext, EH: EH, Exc: Exc };
    }
    Eold = E;
    // miscelazione di Anderson (memoria 1)
    var F = [new Float64Array(mesh + 1), new Float64Array(mesh + 1)];
    for (s = 0; s < 2; s++) for (i = 0; i <= mesh; i++) F[s][i] = Vout[s][i] - Vin[s][i];
    var beta = 0;
    if (Fold) {
      var num = 0, den = 0;
      for (s = 0; s < 2; s++) for (i = 0; i <= mesh; i++) {
        var wt = (nu[i] + nd[i] + 1e-12) * w[i];
        var dF = F[s][i] - Fold[s][i];
        num += F[s][i] * dF * wt; den += dF * dF * wt;
      }
      if (den > 0) beta = num / den;
      if (beta > 0.9) beta = 0.9; if (beta < -0.9) beta = -0.9;
    }
    var Vnew = [new Float64Array(mesh + 1), new Float64Array(mesh + 1)];
    for (s = 0; s < 2; s++) for (i = 0; i <= mesh; i++) {
      var vb = Vin[s][i], fb = F[s][i];
      if (Fold) { vb = (1 - beta) * Vin[s][i] + beta * Vold[s][i]; fb = (1 - beta) * F[s][i] + beta * Fold[s][i]; }
      Vnew[s][i] = vb + alpha * fb;
    }
    Fold = F; Vold = Vin; Vin = Vnew;
  }
  return { ok: false };
}

// Elenco dei possibili sottolivelli (n,l): solo combinatoria, nessun ordine imposto.
function shells(maxn) {
  var list = [];
  for (var n = 1; n <= maxn; n++) for (var l = 0; l < Math.min(n, 4); l++) list.push([n, l]);
  return list;
}

function occsFromMap(map) {
  var occs = [];
  for (var key in map) { if (map[key] > 0) { var p = key.split(','); occs.push({ n: +p[0], l: +p[1], s: +p[2], occ: map[key] }); } }
  occs.sort(function (a, b) { return a.n - b.n || a.l - b.l || a.s - b.s; });
  return occs;
}

// Livelli di energia dei sottolivelli vuoti, nel potenziale attuale
function levels(g, V, maxn) {
  var lv = [];
  var sh = shells(maxn);
  for (var k = 0; k < sh.length; k++) for (var s = 0; s < 2; s++) {
    var sol = solveRadial(g, sh[k][0], sh[k][1], V[s], undefined, g.Z, true);
    if (sol && sol.e < 0) lv.push({ n: sh[k][0], l: sh[k][1], s: s, e: sol.e });
  }
  return lv;
}

// Riempimento iniziale: N elettroni nei livelli piu' bassi del potenziale di partenza (Pauli: 2l+1 per spin).
function initialFill(g, V, N, maxn) {
  var lv = levels(g, [V, V], maxn);
  lv.sort(function (a, b) { return a.e - b.e || a.s - b.s; });
  var map = {}, left = N;
  for (var k = 0; k < lv.length && left > 0; k++) {
    var cap = 2 * lv[k].l + 1, put = Math.min(cap, left);
    map[lv[k].n + ',' + lv[k].l + ',' + lv[k].s] = put; left -= put;
  }
  return map;
}

// Riempimento autoconsistente ("aufbau dinamico"): a ogni passo gli N elettroni occupano i livelli
// piu' bassi del potenziale attuale (Pauli: 2l+1 per spin), con occupazioni smorzate. Nessuna tabella.
function aufbauSCF(g, N, maxn) {
  var mesh = g.mesh, r = g.r, w = g.w, Z = g.Z;
  var V0 = startPotential(g);
  var Vin = [Float64Array.from(V0), Float64Array.from(V0)];
  var sh = shells(maxn), occ = {}, out = [0, 0, 0];
  var nu = new Float64Array(mesh + 1), nd = new Float64Array(mesh + 1), VH = new Float64Array(mesh + 1);
  var lastLv = null, good = null;
  for (var it = 0; it < 70; it++) {
    var lv = [], ys = {};
    for (var k = 0; k < sh.length; k++) for (var s = 0; s < 2; s++) {
      var key = sh[k][0] + ',' + sh[k][1] + ',' + s;
      var sol = solveRadial(g, sh[k][0], sh[k][1], Vin[s], undefined, Z, true);
      if (sol && sol.e < 0) { lv.push({ key: key, l: sh[k][1], s: s, e: sol.e }); ys[key] = sol.y; }
    }
    lv.sort(function (a, b) { return a.e - b.e || a.s - b.s; });
    var target = {}, left = N;
    for (k = 0; k < lv.length && left > 0; k++) { var put = Math.min(2 * lv[k].l + 1, left); target[lv[k].key] = put; left -= put; }
    if (left > 0) { if (!good) return null; for (s = 0; s < 2; s++) Vin[s] = Float64Array.from(good[s]); continue; }
    good = [Float64Array.from(Vin[0]), Float64Array.from(Vin[1])];
    var keys = {}; for (var a in occ) keys[a] = 1; for (a in target) keys[a] = 1;
    var beta = it === 0 ? 1 : 0.35;
    for (a in keys) { var v = (1 - beta) * (occ[a] || 0) + beta * (target[a] || 0); if (!ys[a]) v = 0; occ[a] = v; }
    // rinormalizza a N elettroni
    var tot = 0; for (a in occ) tot += occ[a];
    for (a in occ) occ[a] *= N / tot;
    nu.fill(0); nd.fill(0);
    for (a in occ) { if (occ[a] <= 0) continue; var y = ys[a], dens = a.slice(-1) === '0' ? nu : nd;
      for (var i = 0; i <= mesh; i++) dens[i] += occ[a] * y[i] * y[i] / (4 * Math.PI * r[i]); }
    var q = 0, prev = 0, qin = new Float64Array(mesh + 1);
    for (i = 0; i <= mesh; i++) { var cur = 4 * Math.PI * (nu[i] + nd[i]) * r[i] * r[i] * r[i]; if (i === 0) q = cur / 3; else q += 0.5 * (cur + prev) * g.dx; qin[i] = q; prev = cur; }
    var o2 = 0; prev = 0;
    for (i = mesh; i >= 0; i--) { var c2 = 4 * Math.PI * (nu[i] + nd[i]) * r[i] * r[i]; if (i < mesh) o2 += 0.5 * (c2 + prev) * g.dx; prev = c2; VH[i] = qin[i] / r[i] + o2; }
    for (i = 0; i <= mesh; i++) {
      xc(nu[i], nd[i], out);
      Vin[0][i] += 0.3 * ((-Z / r[i] + VH[i] + out[1]) - Vin[0][i]);
      Vin[1][i] += 0.3 * ((-Z / r[i] + VH[i] + out[2]) - Vin[1][i]);
    }
    lastLv = lv;
    var dmax = 0; for (a in keys) dmax = Math.max(dmax, Math.abs((target[a] || 0) - occ[a]));
    if (it > 12 && dmax < 0.02) break;
  }
  // configurazione intera: riempimento dei livelli finali
  var map = {}, left2 = N;
  for (k = 0; k < lastLv.length && left2 > 0; k++) { var p2 = Math.min(2 * lastLv[k].l + 1, left2); map[lastLv[k].key] = p2; left2 -= p2; }
  return { map: map, V: good, levels: lastLv, occ: occ };
}

// Partenze intere vicine al riempimento frazionario: arrotondamento e spostamenti di 1-3 elettroni
// dai livelli occupati piu' alti a quelli liberi piu' bassi.
function candidateStarts(au, N) {
  var lv = au.levels, cap = {}, e = {};
  lv.forEach(function (x) { cap[x.key] = 2 * x.l + 1; e[x.key] = x.e; });
  var seen = {}, out = [];
  function add(m) { var k = JSON.stringify(Object.keys(m).filter(function (q) { return m[q] > 0; }).sort().map(function (q) { return q + ':' + m[q]; })); if (!seen[k]) { seen[k] = 1; out.push(m); } }
  // arrotondamento delle occupazioni frazionarie (resti piu' grandi)
  var fl = {}, tot = 0, rem = [];
  for (var k in au.occ) { var f = Math.floor(au.occ[k] + 1e-9); fl[k] = f; tot += f; rem.push([k, au.occ[k] - f]); }
  rem.sort(function (a, b) { return b[1] - a[1]; });
  for (var i = 0; tot < N && i < rem.length; i++) { if (fl[rem[i][0]] < cap[rem[i][0]]) { fl[rem[i][0]]++; tot++; } }
  if (tot === N) add(fl);
  [au.map, fl].forEach(function (base) {
    for (var m = 1; m <= 3; m++) {
      var mm = Object.assign({}, base);
      var occd = Object.keys(mm).filter(function (q) { return mm[q] > 0 && e[q] !== undefined; }).sort(function (a, b) { return e[b] - e[a]; });
      var free = lv.filter(function (x) { return (mm[x.key] || 0) < cap[x.key]; }).map(function (x) { return x.key; });
      var moved = 0, oi = 0;
      while (moved < m && oi < occd.length) {
        var src = occd[oi]; if (mm[src] <= 0) { oi++; continue; }
        var dst = null; for (var j = 0; j < free.length; j++) { var fk = free[j]; if (fk !== src && (mm[fk] || 0) < cap[fk] && e[fk] > e[src] - 1e-9) { dst = fk; break; } }
        if (!dst) break;
        mm[src]--; mm[dst] = (mm[dst] || 0) + 1; moved++;
      }
      if (moved === m) add(mm);
    }
  });
  add(au.map);
  return out;
}

// Ricerca della configurazione di minima energia: spostamenti di un elettrone alla volta.
function searchConfig(g, startMap, Vstart, maxn, log, pre) {
  var best = pre || scf(g, occsFromMap(startMap), Vstart);
  var tries = 1;
  if (!best.ok) return { ok: false, tries: tries };
  best.map = startMap;
  for (var round = 0; round < 12; round++) {
    var occ = best.orbs;
    var homo = -1e9; for (var k = 0; k < occ.length; k++) if (occ[k].e > homo) homo = occ[k].e;
    var lv = levels(g, best.V, maxn);
    var eps = {}; for (k = 0; k < lv.length; k++) eps[lv[k].n + ',' + lv[k].l + ',' + lv[k].s] = lv[k].e;
    for (k = 0; k < occ.length; k++) eps[orbKey(occ[k])] = occ[k].e;
    var sources = occ.filter(function (o) { return o.e > homo - 0.6; });
    var targets = [];
    for (var key in eps) {
      var p = key.split(','), l = +p[1], cap = 2 * l + 1, have = best.map[key] || 0;
      if (have < cap && eps[key] < homo + 0.35) targets.push({ key: key, e: eps[key] });
    }
    var moves = [];
    for (var a = 0; a < sources.length; a++) for (var b = 0; b < targets.length; b++) {
      var sk = orbKey(sources[a]); if (sk === targets[b].key) continue;
      moves.push({ from: sk, to: targets[b].key, est: targets[b].e - sources[a].e });
    }
    moves.sort(function (x, y) { return x.est - y.est; });
    var improved = null;
    var lim = Math.min(moves.length, 10);
    for (var m = 0; m < lim; m++) {
      var nm = Object.assign({}, best.map);
      nm[moves[m].from]--; nm[moves[m].to] = (nm[moves[m].to] || 0) + 1;
      var res = scf(g, occsFromMap(nm), best.V, { eg: energiesOf(best) });
      tries++;
      if (res.ok && res.E < best.E - 1e-7 && (!improved || res.E < improved.E)) { res.map = nm; improved = res; }
    }
    if (!improved) break;
    if (log) log('  migliora: ' + configString(improved.map) + '  E=' + improved.E.toFixed(6));
    best = improved;
  }
  best.tries = tries;
  return best;
}

var SUP = { '0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴', '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹' };
function shellTotals(map) {
  var t = {};
  for (var key in map) { if (!map[key]) continue; var p = key.split(','); var k2 = p[0] + LETTERS[+p[1]]; t[k2] = (t[k2] || 0) + map[key]; }
  return t;
}
function configString(map, plain) {
  var t = shellTotals(map), keys = Object.keys(t);
  keys.sort(function (a, b) { var na = parseInt(a), nb = parseInt(b); return na - nb || LETTERS.indexOf(a.slice(-1)) - LETTERS.indexOf(b.slice(-1)); });
  return keys.map(function (k) { var v = String(t[k]); return k + (plain ? v : v.split('').map(function (c) { return SUP[c]; }).join('')); }).join(' ');
}

// Calcolo completo di un atomo: solo Z in ingresso.
function computeAtom(Z, opts) {
  opts = opts || {};
  var t0 = Date.now();
  var g = makeGrid(Z, opts.dx, opts.rmax);
  var maxn = Z <= 2 ? 3 : Z <= 10 ? 4 : Z <= 36 ? 5 : Z <= 54 ? 6 : 7;
  var au = aufbauSCF(g, Z, maxn);
  var neu = { ok: false }, nstart = 0;
  if (au) {
    var starts = candidateStarts(au, Z);
    // si valuta ogni partenza e si cerca a partire dalla migliore
    var bestStart = null;
    for (var j = 0; j < starts.length; j++) {
      var t = scf(g, occsFromMap(starts[j]), au.V); nstart++;
      if (t.ok && (!bestStart || t.E < bestStart.E)) { t.map = starts[j]; bestStart = t; }
    }
    if (bestStart) neu = searchConfig(g, bestStart.map, bestStart.V, maxn, opts.log, bestStart);
  }
  if (!neu.ok) {
    var V0 = startPotential(g);
    neu = searchConfig(g, initialFill(g, V0, Z, maxn), [V0, V0], maxn, opts.log);
  }
  if (!neu.ok) return { Z: Z, ok: false };
  // Catione: si toglie un elettrone dai livelli piu' esterni, poi si lascia riassestare.
  var homo = -1e9; neu.orbs.forEach(function (o) { if (o.e > homo) homo = o.e; });
  var cand = neu.orbs.filter(function (o) { return o.e > homo - 0.6; });
  var bestCat = null, tries = neu.tries + nstart;
  if (Z === 1) { bestCat = { ok: true, E: 0, map: {} }; }
  else {
    for (var k = 0; k < cand.length; k++) {
      var nm = Object.assign({}, neu.map); nm[orbKey(cand[k])]--;
      var c = scf(g, occsFromMap(nm), neu.V, { eg: energiesOf(neu) }); tries++;
      if (c.ok && (!bestCat || c.E < bestCat.E)) { c.map = nm; bestCat = c; }
    }
    if (bestCat) { var cs = searchConfig(g, bestCat.map, bestCat.V, maxn); tries += cs.tries; if (cs.ok && cs.E <= bestCat.E) bestCat = cs; }
  }
  if (!bestCat) return { Z: Z, ok: false };
  var IE = (bestCat.E - neu.E) * HA2EV;
  var Nu = 0, Nd = 0; for (var key in neu.map) { if (key.slice(-1) === '0') Nu += neu.map[key]; else Nd += neu.map[key]; }
  // densita' radiale per il grafico: 4 pi r^2 n(r), campionata
  var dens = [];
  var step = Math.max(1, Math.floor(g.mesh / 220));
  for (var i = 0; i <= g.mesh; i += step) { if (g.r[i] > 12) break; if (g.r[i] < 1e-3) continue; dens.push([+g.r[i].toPrecision(5), +(4 * Math.PI * g.r2[i] * (neu.nu[i] + neu.nd[i])).toPrecision(5)]); }
  var orbs = neu.orbs.map(function (o) {
    var rm = 0; for (var i2 = 0; i2 <= g.mesh; i2++) rm += o.y[i2] * o.y[i2] * g.r2[i2] * g.r[i2] * g.dx;
    return { n: o.n, l: o.l, s: o.s, occ: o.occ, e: o.e, r: rm };
  });
  return {
    Z: Z, ok: true, E: neu.E, Ecat: bestCat.E, IE: IE,
    map: neu.map, config: configString(neu.map), configPlain: configString(neu.map, true),
    catConfig: configString(bestCat.map), unpaired: Math.abs(Nu - Nd),
    orbs: orbs, dens: dens, tries: tries, ms: Date.now() - t0
  };
}

if (typeof module !== 'undefined') module.exports = { computeAtom: computeAtom, makeGrid: makeGrid, scf: scf, solveRadial: solveRadial, xc: xc, startPotential: startPotential, occsFromMap: occsFromMap, configString: configString, shellTotals: shellTotals, HA2EV: HA2EV };
