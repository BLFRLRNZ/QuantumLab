// MOTORE — Livello 0, versione 0.4. Unico dato in ingresso: Z (numero di protoni).
// Leggi impostate: equazione di Schrodinger per gli elettroni (Kohn-Sham, unita' atomiche),
// attrazione di Coulomb del nucleo, repulsione di Coulomb tra elettroni (Hartree),
// principio di Pauli (max 2l+1 elettroni per sottolivello e spin).
// Scambio e correlazione: PBE (gradiente generalizzato, costruito da vincoli fisici esatti) sopra la base
// LSDA (scambio esatto del gas di elettroni + correlazione PW92, tarata su simulazioni quantistiche, non su misure).
// Relativita': equazione radiale scalare-relativistica (Koelling-Harmon): massa che cresce con la velocita'
// e termine di Darwin, senza accoppiamento spin-orbita. Unica costante in piu': c = 137,036 in unita' atomiche.
// Livelli singoli (0.3): nei gusci aperti repulsione tra elettroni e spin-orbita risolte esattamente
// (interazione di configurazioni) con gli orbitali calcolati; si confronta il livello piu' basso.
// Scambio esatto (0.4): ibrido PBE0, un quarto di scambio esatto al posto di quello PBE, valutato sugli
// orbitali di ogni configurazione, con la media vera sulle disposizioni degli elettroni nei gusci aperti.
// Non contiene: tabelle di configurazioni, regola di Madelung, valori misurati.

var HA2EV = 27.211386245988;
var CLIGHT = 137.035999084;
var CFG = { rel: true, xc: 'pbe', lv: 'full', hyb: true };
// derivate del potenziale (prima e seconda rispetto a r), servono per i termini relativistici
// Potenziale di riferimento per i termini relativistici: nucleo + repulsione media degli elettroni (Hartree).
// Lo scambio-correlazione resta nell'equazione ma non nella massa relativistica: cosi' il rumore numerico
// del gradiente (PBE) non entra nei termini relativistici. Vicino al nucleo, dove la relativita' conta,
// il potenziale e' dominato da -Z/r e la differenza e' trascurabile.
function relRef(g) {
  if (g.relCache && g.relCache.src === g.VHref) return g.relCache;
  var m = g.mesh, V = new Float64Array(m + 1);
  for (var i = 0; i <= m; i++) V[i] = -g.Z / g.r[i] + (g.VHref ? g.VHref[i] : 0);
  var d = potDerivs(g, V);
  g.relCache = { src: g.VHref, Vp: d[0], Vpp: d[1], Vs: d[2] };
  return g.relCache;
}
function potDerivs(g, V) {
  // la parte coulombiana del nucleo -Z/r si deriva esattamente; il resto (elettroni) con un fit
  // quadratico locale su una finestra (Savitzky-Golay), che non amplifica il rumore numerico
  var m = g.mesh, r = g.r, dx = g.dx, Z = g.Z, Vp = new Float64Array(m + 1), Vpp = new Float64Array(m + 1), Vs = new Float64Array(m + 1);
  var S = new Float64Array(m + 1), K = 10, i, j;
  for (i = 0; i <= m; i++) S[i] = V[i] + Z / r[i];
  var sj2 = 0, sj4 = 0; for (j = -K; j <= K; j++) { sj2 += j * j; sj4 += j * j * j * j; }
  var mj2 = sj2 / (2 * K + 1), den2 = sj4 - 2 * mj2 * sj2 + mj2 * mj2 * (2 * K + 1);
  for (i = 0; i <= m; i++) {
    var c = Math.min(Math.max(i, K), m - K), c0 = 0, c1 = 0, c2 = 0;
    for (j = -K; j <= K; j++) { var sv = S[c + j]; c0 += sv; c1 += j * sv; c2 += (j * j - mj2) * sv; }
    c1 /= sj2 * dx; c2 /= den2 * dx * dx; c0 = c0 / (2 * K + 1) - c2 * mj2 * dx * dx;
    var t = (i - c) * dx, d1 = c1 + 2 * c2 * t, d2 = 2 * c2, ri = r[i];
    Vs[i] = -Z / ri + c0 + c1 * t + c2 * t * t;
    Vp[i] = Z / (ri * ri) + d1 / ri;
    Vpp[i] = -2 * Z / (ri * ri * ri) + (d2 - d1) / (ri * ri);
  }
  return [Vp, Vpp, Vs];
}
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
  var rel = CFG.rel, A2 = 1 / (CLIGHT * CLIGHT), Vp = null, Vpp = null, Vsm = null, lfac = l * (l + 1);
  if (rel) { var dd = relRef(g); Vp = dd.Vp; Vpp = dd.Vpp; Vsm = dd.Vs; }
  function gfun(i, e) {
    if (!rel) return sqlhf + 2 * r2[i] * (V[i] - e);
    // P = sqrt(M r) y ; M = 1 + (E - V)/(2c^2) : massa relativistica
    var M = 1 + 0.5 * A2 * (e - Vsm[i]), Mp = -0.5 * A2 * Vp[i], Mpp = -0.5 * A2 * Vpp[i];
    var a = Mp / M, ap = Mpp / M - a * a;
    return sqlhf + 2 * r2[i] * M * (V[i] - e) - a * r[i] + r2[i] * (0.25 * a * a - 0.5 * ap);
  }
  for (var iter = 0; iter < 300; iter++) {
    var icl = -1;
    f[0] = ddx12 * gfun(0, e);
    for (i = 1; i <= mesh; i++) {
      f[i] = ddx12 * gfun(i, e);
      if (f[i] === 0) f[i] = 1e-20;
      if ((f[i] > 0) !== (f[i - 1] > 0)) icl = i;
    }
    if (icl < 0 || icl >= mesh - 2) { eup = e; e = 0.5 * (eup + elw); if (eup - elw < 1e-12) return null; continue; }
    // punto di partenza dell'integrazione verso l'interno: dove la funzione e' ormai trascurabile
    var acc = 0, istart = mesh - 1;
    for (i = icl; i < mesh - 1; i++) { var gg = f[i] / ddx12; if (gg > 0) acc += Math.sqrt(gg) * dx; if (acc > 45 && i > icl + 2) { istart = i; break; } }
    for (i = 0; i <= mesh; i++) f[i] = 1 - f[i];
    for (i = 0; i <= mesh; i++) y[i] = 0;
    if (!rel) {
      y[0] = Math.pow(r[0], l + 1) * (1 - 2 * Zn * r[0] / x2l2) / sqr[0];
      y[1] = Math.pow(r[1], l + 1) * (1 - 2 * Zn * r[1] / x2l2) / sqr[1];
    } else {
      // vicino al nucleo la componente grande va come P ~ r^gamma (gamma = l+1 lontano dal regime relativistico)
      var M0 = 1 + 0.5 * A2 * (e - Vsm[0]), M1 = 1 + 0.5 * A2 * (e - Vsm[1]);
      var gam = Math.sqrt(Math.max(lfac + 1 - A2 * Zn * Zn, 0.01)), pexp = M0 > 3 ? gam : l + 1;
      // con il primo termine della serie (cuspide): P ~ r^p (1 - Z r/(l+1))
      y[0] = Math.pow(r[0], pexp) * (1 - Zn * r[0] / (l + 1)) / Math.sqrt(M0 * r[0]);
      y[1] = Math.pow(r[1], pexp) * (1 - Zn * r[1] / (l + 1)) / Math.sqrt(M1 * r[1]);
    }
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
      if (rel) {
        // componente grande P e piccola Q, normalizzate insieme; y restituito come sqrt((P^2+Q^2)/r)
        var P = new Float64Array(mesh + 1), Q = new Float64Array(mesh + 1), nn = 0;
        for (i = 0; i <= mesh; i++) P[i] = Math.sqrt((1 + 0.5 * A2 * (e - Vsm[i])) * r[i]) * y[i];
        for (i = 0; i <= mesh; i++) {
          var dPx = i === 0 ? (P[1] - P[0]) / dx : i === mesh ? (P[mesh] - P[mesh - 1]) / dx : (P[i + 1] - P[i - 1]) / (2 * dx);
          var Mi = 1 + 0.5 * A2 * (e - Vsm[i]);
          Q[i] = (dPx / r[i] - P[i] / r[i]) / (2 * Mi * CLIGHT);
          nn += (P[i] * P[i] + Q[i] * Q[i]) * r[i] * dx;
        }
        var yy = new Float64Array(mesh + 1), pp = new Float64Array(mesh + 1), np = 0;
        for (i = 0; i <= mesh; i++) { yy[i] = Math.sqrt((P[i] * P[i] + Q[i] * Q[i]) / (nn * r[i])); np += P[i] * P[i] * r[i] * dx; }
        np = Math.sqrt(np);
        for (i = 0; i <= mesh; i++) pp[i] = P[i] / np;
        return { e: e, y: yy, p: pp };
      }
      var p0 = new Float64Array(mesh + 1);
      for (i = 0; i <= mesh; i++) p0[i] = y[i] * sqr[i];
      return { e: e, y: y, p: p0 };
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

// PBE (Perdew-Burke-Ernzerhof 1996): approssimazione a gradiente generalizzato, costruita solo da vincoli
// fisici esatti (limite di Lieb-Oxford, limiti del gas di elettroni), senza parametri tarati su misure.
var KAPPA = 0.804, MUX = 0.2195149727645171, BETA = 0.06672455060314922, GAMC = (1 - Math.log(2)) / (Math.PI * Math.PI);
var CXU = 0.75 * Math.cbrt(3 / Math.PI);
function exUnpol(n, gabs) {
  if (n < 1e-30) return 0;
  var exu = -CXU * n * Math.cbrt(n);
  var kf = Math.cbrt(3 * Math.PI * Math.PI * n), sg = gabs / (2 * kf * n);
  return exu * (1 + KAPPA - KAPPA / (1 + MUX * sg * sg / KAPPA));
}
function ecUnif(rs, z) {
  var c0 = Gpw(rs, 0.031091, 0.21370, 7.5957, 3.5876, 1.6382, 0.49294)[0];
  var c1 = Gpw(rs, 0.015545, 0.20548, 14.1189, 6.1977, 3.3662, 0.62517)[0];
  var ac = -Gpw(rs, 0.016887, 0.11125, 10.357, 3.6231, 0.88026, 0.49671)[0];
  var fz = (Math.pow(1 + z, 4 / 3) + Math.pow(1 - z, 4 / 3) - 2) / D43, z4 = z * z * z * z;
  return c0 + ac * fz / FPP0 * (1 - z4) + (c1 - c0) * fz * z4;
}
// densita' di energia di scambio-correlazione PBE per unita' di volume; gu, gd = dn/dr per spin
function ePBE(nu, nd, gu, gd) {
  if (nu < 0) nu = 0; if (nd < 0) nd = 0;
  var n = nu + nd; if (n < 1e-30) return 0;
  var ex = 0.5 * exUnpol(2 * nu, 2 * Math.abs(gu)) + 0.5 * exUnpol(2 * nd, 2 * Math.abs(gd));
  var rs = Math.cbrt(3 / (4 * Math.PI * n)), z = (nu - nd) / n; if (z > 1) z = 1; if (z < -1) z = -1;
  var ecu = ecUnif(rs, z);
  var phi = 0.5 * (Math.cbrt((1 + z) * (1 + z)) + Math.cbrt((1 - z) * (1 - z))), phi3 = phi * phi * phi;
  var kf = Math.cbrt(3 * Math.PI * Math.PI * n), ks = Math.sqrt(4 * kf / Math.PI);
  var t = Math.abs(gu + gd) / (2 * phi * ks * n), t2 = t * t;
  var den = Math.exp(-ecu / (GAMC * phi3)) - 1;
  var A = den > 1e-300 ? BETA / GAMC / den : 1e300, At2 = A * t2;
  var H = GAMC * phi3 * Math.log(1 + BETA / GAMC * t2 * (1 + At2) / (1 + At2 + At2 * At2));
  return ex + n * (ecu + H);
}
// PBE con derivate analitiche: restituisce [e, de/dn_up, de/dn_dn, de/dg_up, de/dg_dn]
function exPart(n, gs, out) { // scambio per densita' non polarizzata n e gradiente con segno gs
  if (n < 1e-30) { out[0] = 0; out[1] = 0; out[2] = 0; return; }
  var n13 = Math.cbrt(n), kf = Math.cbrt(3 * Math.PI * Math.PI) * n13;
  var exu = -CXU * n * n13, q = 1 / (2 * kf * n), sg = Math.abs(gs) * q, s2 = sg * sg;
  var den = 1 + MUX * s2 / KAPPA, F = 1 + KAPPA - KAPPA / den;
  var dFds2 = MUX / (den * den);
  out[0] = exu * F;
  out[1] = -CXU * (4 / 3) * n13 * F + exu * dFds2 * (-(8 / 3) * s2 / n);   // ds2/dn = -(8/3) s2/n
  out[2] = exu * dFds2 * 2 * gs * q * q;                                     // ds2/dg = 2 g q^2
}
var EXA = [0, 0, 0], EXB = [0, 0, 0];
function pbeFull(nu, nd, gu, gd, out) {
  if (nu < 0) nu = 0; if (nd < 0) nd = 0;
  var n = nu + nd;
  if (n < 1e-30) { out[0] = out[1] = out[2] = out[3] = out[4] = 0; return; }
  exPart(2 * nu, 2 * gu, EXA); exPart(2 * nd, 2 * gd, EXB);
  var e = 0.5 * (EXA[0] + EXB[0]), dnu = EXA[1], dnd = EXB[1], dgu = EXA[2], dgd = EXB[2];
  // correlazione
  var rs = Math.cbrt(3 / (4 * Math.PI * n)), z = (nu - nd) / n;
  if (z > 1 - 1e-6) z = 1 - 1e-6; if (z < -1 + 1e-6) z = -1 + 1e-6;
  var c0 = Gpw(rs, 0.031091, 0.21370, 7.5957, 3.5876, 1.6382, 0.49294);
  var c1 = Gpw(rs, 0.015545, 0.20548, 14.1189, 6.1977, 3.3662, 0.62517);
  var ca = Gpw(rs, 0.016887, 0.11125, 10.357, 3.6231, 0.88026, 0.49671);
  var ac = -ca[0], dac = -ca[1];
  var p1 = Math.cbrt(1 + z), m1 = Math.cbrt(1 - z);
  var fz = ((1 + z) * p1 + (1 - z) * m1 - 2) / D43, dfz = (4 / 3) * (p1 - m1) / D43;
  var z3 = z * z * z, z4 = z3 * z;
  var ec = c0[0] + ac * fz / FPP0 * (1 - z4) + (c1[0] - c0[0]) * fz * z4;
  var decr = c0[1] * (1 - fz * z4) + c1[1] * fz * z4 + dac * fz / FPP0 * (1 - z4);
  var decz = 4 * z3 * fz * (c1[0] - c0[0] - ac / FPP0) + dfz * (z4 * (c1[0] - c0[0]) + (1 - z4) * ac / FPP0);
  var decn = -decr * rs / (3 * n);
  var phi = 0.5 * (p1 * p1 + m1 * m1), dphi = (1 / 3) * (1 / p1 - 1 / m1), phi2 = phi * phi, phi3 = phi2 * phi;
  var g = gu + gd, kf = Math.cbrt(3 * Math.PI * Math.PI * n);
  var t2 = g * g * Math.PI / (16 * phi2 * kf * n * n);
  var bg = BETA / GAMC, gp3 = GAMC * phi3;
  var E = Math.exp(-ec / gp3), Em1 = E - 1;
  var A = Em1 > 1e-300 ? bg / Em1 : 1e300;
  var dAde = bg * E / gp3 / (Em1 * Em1);                 // dA/d eps
  var dAdphi = -bg * E * 3 * ec / (GAMC * phi3 * phi) / (Em1 * Em1);
  var u = t2, N1 = 1 + A * u, D1 = 1 + A * u + A * A * u * u;
  var Q = bg * u * N1 / D1;
  var dQdu = bg * (N1 * D1 + u * A * D1 - u * N1 * (A + 2 * A * A * u)) / (D1 * D1);
  var dQdA = bg * u * (u * D1 - N1 * (u + 2 * A * u * u)) / (D1 * D1);
  var L = Math.log(1 + Q), H = gp3 * L;
  var dHdQ = gp3 / (1 + Q), dHdu = dHdQ * dQdu, dHdA = dHdQ * dQdA;
  var dHdn = dHdA * dAde * decn + dHdu * (-(7 / 3) * u / n);
  var dHdz = (3 * GAMC * phi2 * L + dHdA * dAdphi + dHdu * (-2 * u / phi)) * dphi + dHdA * dAde * decz;
  var dHdg = g !== 0 ? dHdu * 2 * u / g : 0;
  e += n * (ec + H);
  var dEdn = ec + H + n * (decn + dHdn), dEdz = n * (decz + dHdz);
  dnu += dEdn + (1 - z) / n * dEdz;
  dnd += dEdn - (1 + z) / n * dEdz;
  dgu += n * dHdg; dgd += n * dHdg;
  out[0] = e; out[1] = dnu; out[2] = dnd; out[3] = dgu; out[4] = dgd;
}

// energia e potenziali di scambio-correlazione su tutta la griglia
function xcArrays(g, nu, nd, Ed, Vu, Vd) {
  var mesh = g.mesh, r = g.r, dx = g.dx, out = [0, 0, 0], i;
  if (CFG.xc !== 'pbe') {
    for (i = 0; i <= mesh; i++) { xc(nu[i], nd[i], out); Ed[i] = out[0]; Vu[i] = out[1]; Vd[i] = out[2]; }
    return;
  }
  var gu = new Float64Array(mesh + 1), gd = new Float64Array(mesh + 1), Fu = new Float64Array(mesh + 1), Fd = new Float64Array(mesh + 1);
  for (i = 0; i <= mesh; i++) {
    var a = i === 0 ? 0 : i - 1, b = i === mesh ? mesh : i + 1, span = (b - a) * dx * r[i];
    gu[i] = (nu[b] - nu[a]) / span; gd[i] = (nd[b] - nd[a]) / span;
  }
  var P5 = [0, 0, 0, 0, 0];
  for (i = 0; i <= mesh; i++) {
    var n = nu[i] + nd[i];
    if (n < 1e-14) { Ed[i] = 0; Vu[i] = 0; Vd[i] = 0; Fu[i] = 0; Fd[i] = 0; continue; }
    pbeFull(nu[i], nd[i], gu[i], gd[i], P5);
    Ed[i] = P5[0]; Vu[i] = P5[1]; Vd[i] = P5[2];
    Fu[i] = r[i] * r[i] * P5[3]; Fd[i] = r[i] * r[i] * P5[4];
  }
  // termine di gradiente: - (1/r^2) d/dr (r^2 de/dg)
  for (i = 0; i <= mesh; i++) {
    var a2 = i === 0 ? 0 : i - 1, b2 = i === mesh ? mesh : i + 1, sp = (b2 - a2) * dx;
    var r3 = r[i] * r[i] * r[i];
    Vu[i] -= (Fu[b2] - Fu[a2]) / sp / r3;
    Vd[i] -= (Fd[b2] - Fd[a2]) / sp / r3;
  }
  // i primissimi punti (derivate unilaterali) si raccordano con la forma A + B/r, che e' il comportamento
  // del potenziale PBE vicino al nucleo
  [Vu, Vd].forEach(function (Vx) {
    var B = (Vx[3] - Vx[4]) / (1 / r[3] - 1 / r[4]), A = Vx[3] - B / r[3];
    for (var k = 0; k < 3; k++) Vx[k] = A + B / r[k];
  });
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
  var XE = new Float64Array(mesh + 1), XU = new Float64Array(mesh + 1), XD = new Float64Array(mesh + 1);
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
      orbs.push({ n: o.n, l: o.l, s: o.s, occ: o.occ, e: sol.e, y: sol.y, p: sol.p });
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
    g.VHref = Float64Array.from(VH);
    var Exc = 0, EH = 0, Eext = 0, vinn = 0;
    xcArrays(g, nu, nd, XE, XU, XD);
    for (i = 0; i <= mesh; i++) {
      var ntot = nu[i] + nd[i];
      Exc += XE[i] * w[i];
      EH += 0.5 * VH[i] * ntot * w[i];
      Eext += (-Z / r[i]) * ntot * w[i];
      vinn += (Vin[0][i] * nu[i] + Vin[1][i] * nd[i]) * w[i];
      Vout[0][i] = -Z / r[i] + VH[i] + XU[i];
      Vout[1][i] = -Z / r[i] + VH[i] + XD[i];
    }
    var Ts = sumeps - vinn;
    var E = Ts + Eext + EH + Exc;
    // residuo
    var res = 0;
    for (var s = 0; s < 2; s++) for (i = 0; i <= mesh; i++) { var d = Vout[s][i] - Vin[s][i]; res += d * d * (nu[i] + nd[i] + 1e-12) * w[i]; }
    res = Math.sqrt(res);
    if (it > 1 && Math.abs(E - Eold) < 1e-9 && res < 1e-6) {
      if (POOL) poolAdd(occs, E, orbs, g.VHref);
      return { ok: true, E: E, V: Vin, orbs: orbs, it: it, nu: nu, nd: nd, Ts: Ts, Eext: Eext, EH: EH, Exc: Exc, VH: g.VHref };
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
    g.VHref = Float64Array.from(VH);
    var XE = new Float64Array(mesh + 1), XU = new Float64Array(mesh + 1), XD = new Float64Array(mesh + 1);
    xcArrays(g, nu, nd, XE, XU, XD);
    for (i = 0; i <= mesh; i++) {
      Vin[0][i] += 0.3 * ((-Z / r[i] + VH[i] + XU[i]) - Vin[0][i]);
      Vin[1][i] += 0.3 * ((-Z / r[i] + VH[i] + XD[i]) - Vin[1][i]);
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

// ---------------------------------------------------------------------------------------------
// LIVELLI SINGOLI: MULTIPLETTI E SPIN-ORBITA (dalla versione 0.3)
// Il calcolo a simmetria sferica da' l'energia media delle disposizioni con un dato numero di elettroni
// con spin su e giu' in ogni guscio. La busta pero' registra il singolo livello piu' basso. Dentro i gusci
// aperti si risolve quindi il problema esatto a molti elettroni (interazione di configurazioni) con:
//  - la repulsione tra gli elettroni dei gussci aperti (integrali di Slater dagli orbitali del motore);
//  - l'accoppiamento spin-orbita, con la costante zeta = < (1/(2 M^2 c^2)) (1/r) dV/dr > (stesso potenziale).
// Correzione = livello piu' basso - media delle disposizioni con gli stessi spin del calcolo sferico.
// Nessun parametro esterno: tutto viene dagli orbitali calcolati.

// raccolta delle disposizioni calcolate (per scegliere alla fine quella col livello piu' basso)
var POOL = null;
function poolAdd(occs, E, orbs, VH) {
  var map = {}; occs.forEach(function (o) { if (o.occ > 0) map[orbKey(o)] = o.occ; });
  var key = JSON.stringify(Object.keys(map).sort().map(function (k) { return k + ':' + map[k]; }));
  for (var i = 0; i < POOL.length; i++) if (POOL[i].key === key) { if (E < POOL[i].E) { POOL[i].E = E; POOL[i].orbs = orbs; POOL[i].VH = VH; } return; }
  POOL.push({ key: key, map: map, E: E, orbs: orbs, VH: VH });
  POOL.sort(function (a, b) { return a.E - b.E; });
  if (POOL.length > 10) POOL.length = 10;
}

var LF = [0]; for (var q_ = 1; q_ < 80; q_++) LF[q_] = LF[q_ - 1] + Math.log(q_);
function w3j(j1, j2, j3, m1, m2, m3) { // simboli 3j di Wigner (argomenti interi)
  if (m1 + m2 + m3 !== 0) return 0;
  if (j3 < Math.abs(j1 - j2) || j3 > j1 + j2) return 0;
  if (Math.abs(m1) > j1 || Math.abs(m2) > j2 || Math.abs(m3) > j3) return 0;
  var pre = 0.5 * (LF[j1 + j2 - j3] + LF[j1 - j2 + j3] + LF[-j1 + j2 + j3] - LF[j1 + j2 + j3 + 1]
    + LF[j1 + m1] + LF[j1 - m1] + LF[j2 + m2] + LF[j2 - m2] + LF[j3 + m3] + LF[j3 - m3]);
  var tmin = Math.max(0, j2 - j3 - m1, j1 - j3 + m2), tmax = Math.min(j1 + j2 - j3, j1 - m1, j2 + m2), sum = 0;
  for (var t = tmin; t <= tmax; t++) {
    var den = LF[t] + LF[j3 - j2 + t + m1] + LF[j3 - j1 + t - m2] + LF[j1 + j2 - j3 - t] + LF[j1 - t - m1] + LF[j2 - t + m2];
    sum += ((t % 2) ? -1 : 1) * Math.exp(pre - den);
  }
  return (((j1 - j2 - m3) % 2 + 2) % 2 ? -1 : 1) * sum;
}
// coefficienti angolari di Condon-Shortley c^k(l m, l' m')
function ckc(k, l, m, lp, mp) {
  return ((Math.abs(m) % 2) ? -1 : 1) * Math.sqrt((2 * l + 1) * (2 * lp + 1)) * w3j(l, k, lp, 0, 0, 0) * w3j(l, k, lp, -m, m - mp, mp);
}
// integrale radiale di Slater R^k(ab;cd) = int int Pa(1)Pc(1) r<^k/r>^(k+1) Pb(2)Pd(2)
function rPowers(g, k) {
  if (!g.pw) g.pw = {};
  if (!g.pw[k]) { var m = g.mesh, a = new Float64Array(m + 1), b = new Float64Array(m + 1); for (var i = 0; i <= m; i++) { a[i] = Math.pow(g.r[i], k); b[i] = 1 / (a[i] * g.r[i]); } g.pw[k] = [a, b]; }
  return g.pw[k];
}
function slaterR(g, Pa, Pb, Pc, Pd, k) {
  var m = g.mesh, r = g.r, dx = g.dx, i, inner = new Float64Array(m + 1), outer = new Float64Array(m + 1), acc = 0, prev = 0;
  var pw = rPowers(g, k), rk = pw[0], rmk = pw[1];
  for (i = 0; i <= m; i++) { var f = Pb[i] * Pd[i] * rk[i] * r[i]; if (i) acc += 0.5 * (f + prev) * dx; prev = f; inner[i] = acc; }
  acc = 0; prev = 0;
  for (i = m; i >= 0; i--) { var f2 = Pb[i] * Pd[i] * rmk[i] * r[i]; if (i < m) acc += 0.5 * (f2 + prev) * dx; prev = f2; outer[i] = acc; }
  var R = 0;
  for (i = 0; i <= m; i++) { var Y = inner[i] * rmk[i] + rk[i] * outer[i]; R += Pa[i] * Pc[i] * Y * r[i] * dx * (i === 0 || i === m ? 0.5 : 1); }
  return R;
}
// autovalore piu' basso di una matrice simmetrica sparsa (Lanczos con riortogonalizzazione completa)
function lowestEig(n, rows, cols, vals, diag) {
  if (n === 1) return diag[0];
  function mul(v, out) { for (var i = 0; i < n; i++) out[i] = diag[i] * v[i]; for (var k = 0; k < rows.length; k++) { out[rows[k]] += vals[k] * v[cols[k]]; out[cols[k]] += vals[k] * v[rows[k]]; } }
  var mmax = Math.min(n, 150), Qv = [], alpha = [], beta = [], v = new Float64Array(n), w = new Float64Array(n), i, j;
  for (i = 0; i < n; i++) v[i] = 1 / Math.sqrt(n) * (1 + 0.37 * Math.sin(1.7 * i + 0.3));
  var nv = 0; for (i = 0; i < n; i++) nv += v[i] * v[i]; nv = Math.sqrt(nv); for (i = 0; i < n; i++) v[i] /= nv;
  var last = 1e30, ev = 0;
  for (j = 0; j < mmax; j++) {
    Qv.push(Float64Array.from(v)); mul(v, w);
    var a = 0; for (i = 0; i < n; i++) a += w[i] * v[i]; alpha.push(a);
    // ortogonalizzazione ripetuta due volte ("due volte bastano") contro tutti i vettori precedenti
    for (var pass = 0; pass < 2; pass++)
      for (var q = 0; q < Qv.length; q++) { var d = 0, Qq = Qv[q]; for (i = 0; i < n; i++) d += w[i] * Qq[i]; for (i = 0; i < n; i++) w[i] -= d * Qq[i]; }
    var b = 0; for (i = 0; i < n; i++) b += w[i] * w[i]; b = Math.sqrt(b);
    ev = tridiagLowest(alpha, beta);
    var scale = 0; for (i = 0; i < alpha.length; i++) scale = Math.max(scale, Math.abs(alpha[i]));
    if (b < 1e-9 * Math.max(1, scale) || (j > 8 && Math.abs(ev - last) < 1e-11)) break;
    last = ev; beta.push(b); for (i = 0; i < n; i++) v[i] = w[i] / b;
  }
  return ev;
}
function tridiagLowest(a, b) { // bisezione di Sturm sul tridiagonale
  var n = a.length, lo = 1e30, hi = -1e30, i;
  for (i = 0; i < n; i++) { var rad = (i > 0 ? Math.abs(b[i - 1]) : 0) + (i < n - 1 && b[i] !== undefined ? Math.abs(b[i]) : 0); lo = Math.min(lo, a[i] - rad); hi = Math.max(hi, a[i] + rad); }
  for (var it = 0; it < 200; it++) {
    var mid = 0.5 * (lo + hi), cnt = 0, d = 1;
    for (i = 0; i < n; i++) { d = a[i] - mid - (i > 0 ? b[i - 1] * b[i - 1] / d : 0); if (d === 0) d = 1e-300; if (d < 0) cnt++; }
    if (cnt >= 1) hi = mid; else lo = mid;
    if (hi - lo < 1e-13 * Math.max(1, Math.abs(mid))) break;
  }
  return 0.5 * (lo + hi);
}
function popc(x) { x = x - ((x >>> 1) & 0x55555555); x = (x & 0x33333333) + ((x >>> 2) & 0x33333333); return (((x + (x >>> 4)) & 0x0F0F0F0F) * 0x01010101) >>> 24; }

// Correzione per una disposizione: { dE (Hartree), size, shells }
// shells: [{l, P, zeta, nUp, nDn}] gusci aperti; Iext: integrali radiali (per i test si possono passare a mano)
function multipletCore(shells, Rk) {
  // spin-orbitali: per ogni guscio m = -l..l, spin +1/2 (indice pari) e -1/2 (dispari)
  var so = [];
  shells.forEach(function (sh, si) { for (var m = -sh.l; m <= sh.l; m++) { so.push({ sh: si, l: sh.l, m: m, s: 1 }); so.push({ sh: si, l: sh.l, m: m, s: -1 }); } });
  var NS = so.length;
  if (NS > 31) return null;
  // tabella degli integrali a due elettroni <pq|rs>
  var T = new Float64Array(NS * NS * NS * NS);
  for (var p = 0; p < NS; p++) for (var q = 0; q < NS; q++) for (var r = 0; r < NS; r++) {
    if (so[p].s !== so[r].s) continue;
    for (var t = 0; t < NS; t++) {
      if (so[q].s !== so[t].s || so[p].m + so[q].m !== so[r].m + so[t].m) continue;
      var A = so[p], B = so[q], C = so[r], D = so[t], val = 0;
      var kmin = Math.max(Math.abs(A.l - C.l), Math.abs(B.l - D.l)), kmax = Math.min(A.l + C.l, B.l + D.l);
      for (var k = kmin; k <= kmax; k++) {
        if ((A.l + C.l + k) % 2 || (B.l + D.l + k) % 2) continue;
        var c1 = ckc(k, A.l, A.m, C.l, C.m); if (!c1) continue;
        var c2 = ckc(k, D.l, D.m, B.l, B.m); if (!c2) continue;
        val += c1 * c2 * Rk(A.sh, B.sh, C.sh, D.sh, k);
      }
      T[((p * NS + q) * NS + r) * NS + t] = val;
    }
  }
  function V2(p, q, r, t) { return T[((p * NS + q) * NS + r) * NS + t]; }
  // spin-orbita a un corpo nello stesso guscio
  function hso(a, b) {
    var A = so[a], B = so[b]; if (A.sh !== B.sh) return 0;
    var z = shells[A.sh].zeta, l = A.l; if (!z || !l) return 0;
    if (a === b) return z * A.m * A.s * 0.5;
    if (A.m === B.m + 1 && A.s === -1 && B.s === 1) return z * 0.5 * Math.sqrt(l * (l + 1) - B.m * (B.m + 1));
    if (A.m === B.m - 1 && A.s === 1 && B.s === -1) return z * 0.5 * Math.sqrt(l * (l + 1) - B.m * (B.m - 1));
    return 0;
  }
  // determinanti: combinazioni per guscio
  function combos(si) {
    var idx = []; for (var i = 0; i < NS; i++) if (so[i].sh === si) idx.push(i);
    var N = shells[si].nUp + shells[si].nDn, out = [];
    (function rec(start, left, bits) { if (!left) { out.push(bits); return; } for (var i = start; i <= idx.length - left; i++) rec(i + 1, left - 1, bits | (1 << idx[i])); })(0, N, 0);
    return out;
  }
  var lists = shells.map(function (_, si) { return combos(si); });
  var total = lists.reduce(function (a, l) { return a * l.length; }, 1);
  if (total > 400000) return null;
  var Ntot = shells.reduce(function (a, sh) { return a + sh.nUp + sh.nDn; }, 0), target2 = Ntot % 2;
  function m2(bits) { var mm = 0; for (var i = 0; i < NS; i++) if (bits & (1 << i)) mm += 2 * so[i].m + so[i].s; return mm; }
  function spinsOk(bits) { for (var si = 0; si < shells.length; si++) { var u = 0, d = 0; for (var i = 0; i < NS; i++) if ((bits & (1 << i)) && so[i].sh === si) { if (so[i].s > 0) u++; else d++; } if (u !== shells[si].nUp || d !== shells[si].nDn) return false; } return true; }
  function diagE(bits) {
    var occ = []; for (var i = 0; i < NS; i++) if (bits & (1 << i)) occ.push(i);
    var e = 0; for (var a = 0; a < occ.length; a++) { e += hso(occ[a], occ[a]); for (var b = a + 1; b < occ.length; b++) e += V2(occ[a], occ[b], occ[a], occ[b]) - V2(occ[a], occ[b], occ[b], occ[a]); }
    return e;
  }
  var dets = [], refSum = 0, refN = 0;
  (function prod(si, bits) {
    if (si === lists.length) {
      if (spinsOk(bits)) { refSum += diagE(bits); refN++; }
      if (Math.abs(m2(bits)) === target2 && m2(bits) >= 0) dets.push(bits);
      return;
    }
    var L = lists[si]; for (var i = 0; i < L.length; i++) prod(si + 1, bits | L[i]);
  })(0, 0);
  var n = dets.length; if (!n || !refN) return null;
  if (n > 6000) return null;
  var diag = new Float64Array(n), rows = [], cols = [], vals = [];
  for (var i = 0; i < n; i++) diag[i] = diagE(dets[i]);
  function below(bits, pos) { return popc(bits & ((1 << pos) - 1)); }
  for (i = 0; i < n; i++) for (var j = i + 1; j < n; j++) {
    var I = dets[i], J = dets[j], x = I ^ J, deg = popc(x);
    if (deg !== 2 && deg !== 4) continue;
    var rem = I & x, add = J & x, h;
    if (deg === 2) {
      var ii = 31 - Math.clz32(rem), aa = 31 - Math.clz32(add);
      var sgn = ((below(I, ii) + below(I ^ (1 << ii), aa)) % 2) ? -1 : 1;
      h = hso(aa, ii);
      var occ = I & ~(1 << ii);
      for (var kk = 0; kk < NS; kk++) if (occ & (1 << kk)) h += V2(aa, kk, ii, kk) - V2(aa, kk, kk, ii);
      h *= sgn;
    } else {
      var i1 = 31 - Math.clz32(rem & -rem), i2 = 31 - Math.clz32(rem & ~(rem & -rem));
      var a1 = 31 - Math.clz32(add & -add), a2 = 31 - Math.clz32(add & ~(add & -add));
      // segno di a+_a1 a+_a2 a_i2 a_i1 |I>
      var bts = I, sg = 0;
      sg += below(bts, i1); bts ^= (1 << i1);
      sg += below(bts, i2); bts ^= (1 << i2);
      sg += below(bts, a2); bts ^= (1 << a2);
      sg += below(bts, a1); bts ^= (1 << a1);
      h = ((sg % 2) ? -1 : 1) * (V2(a1, a2, i1, i2) - V2(a1, a2, i2, i1));
    }
    if (h !== 0 && Math.abs(h) > 1e-14) { rows.push(i); cols.push(j); vals.push(h); }
  }
  var low = lowestEig(n, rows, cols, vals, diag);
  return { dE: low - refSum / refN, size: n, low: low, ref: refSum / refN };
}

// Gusci aperti di una disposizione calcolata e loro correzione
function multipletCorrection(g, entry) {
  var tot = {}, ups = {}, dns = {}, k;
  for (k in entry.map) { var p = k.split(','), nl = p[0] + ',' + p[1]; tot[nl] = (tot[nl] || 0) + entry.map[k]; if (p[2] === '0') ups[nl] = entry.map[k]; else dns[nl] = entry.map[k]; }
  var shells = [];
  for (k in tot) {
    var q = k.split(','), n = +q[0], l = +q[1];
    if (tot[k] > 0 && tot[k] < 2 * (2 * l + 1)) {
      var o = null; entry.orbs.forEach(function (x) { if (x.n === n && x.l === l && (!o || x.occ > o.occ)) o = x; });
      if (!o || !o.p) return { dE: 0, skipped: true };
      shells.push({ n: n, l: l, nUp: ups[k] || 0, nDn: dns[k] || 0, P: o.p, e: o.e });
    }
  }
  if (!shells.length) return { dE: 0, shells: 0 };
  // costante spin-orbita: zeta = int P^2 (1/(2 M^2 c^2)) (1/r) dV/dr dr, con V = nucleo + Hartree
  var m = g.mesh, r = g.r, dx = g.dx, Vr = new Float64Array(m + 1), i;
  for (i = 0; i <= m; i++) Vr[i] = -g.Z / r[i] + (entry.VH ? entry.VH[i] : 0);
  var dd = potDerivs(g, Vr), Vp = dd[0], A2 = 1 / (CLIGHT * CLIGHT);
  shells.forEach(function (sh) {
    if (!sh.l) { sh.zeta = 0; return; }
    var z = 0;
    for (i = 0; i <= m; i++) { var M = CFG.rel ? 1 + 0.5 * A2 * (sh.e - Vr[i]) : 1; z += sh.P[i] * sh.P[i] * 0.5 * A2 / (M * M) * Vp[i] / r[i] * r[i] * dx; }
    sh.zeta = z;
  });
  var cache = {};
  function Rk(a, b, c, d, kk) {
    if (CFG.lv === 'so' && kk > 0) return 0; // variante: solo spin-orbita, senza la struttura dei multipletti
    var key = a + ',' + b + ',' + c + ',' + d + ',' + kk;
    if (cache[key] === undefined) cache[key] = slaterR(g, shells[a].P, shells[b].P, shells[c].P, shells[d].P, kk);
    return cache[key];
  }
  var res = multipletCore(shells, Rk);
  if (!res) return { dE: 0, skipped: true };
  res.zetas = shells.map(function (sh) { return [sh.n + LETTERS[sh.l], sh.zeta]; });
  return res;
}

// ---------------------------------------------------------------------------------------------
// SCAMBIO ESATTO: IBRIDO PBE0 (dalla versione 0.4)
// Nello scambio esatto ogni elettrone non interagisce con se stesso. PBE0 sostituisce un quarto dello
// scambio PBE con quello esatto: E = E_PBE + 1/4 (Ex_esatto - Ex_PBE). Il quarto viene da un argomento
// di teoria delle perturbazioni (Perdew, Ernzerhof, Burke 1996), non da un adattamento ai dati.
// Qui la correzione si valuta sugli orbitali PBE di ogni configurazione (non autoconsistente).
var ANGX = {};
function angX(la, lb, k) { // somma su m, m' di c^k(la m, lb m')^2
  var key = la + ',' + lb + ',' + k;
  if (ANGX[key] === undefined) { var sm = 0; for (var ma = -la; ma <= la; ma++) for (var mb = -lb; mb <= lb; mb++) { var c = ckc(k, la, ma, lb, mb); sm += c * c; } ANGX[key] = sm; }
  return ANGX[key];
}
function exactExchange(g, orbs) {
  var Ex = 0;
  for (var sp = 0; sp < 2; sp++) {
    var list = orbs.filter(function (o) { return o.s === sp && o.occ > 0; });
    for (var a = 0; a < list.length; a++) for (var b = a; b < list.length; b++) {
      var A = list[a], B = list[b], fa = A.occ / (2 * A.l + 1), fb = B.occ / (2 * B.l + 1), sum = 0;
      for (var k = Math.abs(A.l - B.l); k <= A.l + B.l; k += 2) { var co = angX(A.l, B.l, k); if (co) sum += co * slaterR(g, A.p, B.p, B.p, A.p, k); }
      Ex += -0.5 * fa * fb * sum * (a === b ? 1 : 2);
    }
    // gusci aperti: media vera sulle disposizioni degli elettroni (non occupazioni frazionarie).
    // Termine 1/2 * N(N-L)/(L^2 (L-1)) * somma_{m != m'} (J_mm' - K_mm'), L = 2l+1
    list.forEach(function (o) {
      var L = 2 * o.l + 1, N = o.occ; if (N <= 0 || N >= L || L < 2) return;
      var tJ = 0, tK = 0;
      for (var k = 0; k <= 2 * o.l; k += 2) {
        var Fk = slaterR(g, o.p, o.p, o.p, o.p, k), diag2 = 0;
        for (var m = -o.l; m <= o.l; m++) { var c = ckc(k, o.l, m, o.l, m); diag2 += c * c; }
        tJ += (k === 0 ? L * L - L : -diag2) * Fk;
        tK += (angX(o.l, o.l, k) - diag2) * Fk;
      }
      Ex += 0.5 * N * (N - L) / (L * L * (L - 1)) * (tJ - tK);
    });
  }
  return Ex;
}
function pbeExchange(g, orbs) {
  var m = g.mesh, r = g.r, dx = g.dx, nu = new Float64Array(m + 1), nd = new Float64Array(m + 1), i;
  orbs.forEach(function (o) { var d = o.s === 0 ? nu : nd; for (var j = 0; j <= m; j++) d[j] += o.occ * o.y[j] * o.y[j] / (4 * Math.PI * r[j]); });
  var Ex = 0, A = [0, 0, 0], B = [0, 0, 0];
  for (i = 0; i <= m; i++) {
    var lo = i === 0 ? 0 : i - 1, hi = i === m ? m : i + 1, span = (hi - lo) * dx * r[i];
    exPart(2 * nu[i], 2 * (nu[hi] - nu[lo]) / span, A); exPart(2 * nd[i], 2 * (nd[hi] - nd[lo]) / span, B);
    Ex += 0.5 * (A[0] + B[0]) * g.w[i];
  }
  return Ex;
}
function hybridCorrection(g, entry) {
  if (!CFG.hyb || !entry.orbs || !entry.orbs.length || !entry.orbs[0].p) return 0;
  return 0.25 * (exactExchange(g, entry.orbs) - pbeExchange(g, entry.orbs));
}

// Tra le disposizioni raccolte sceglie quella con il livello piu' basso (energia media + correzione)
function pickLowestLevel(g, pool, fallback) {
  var best = null;
  var E0 = pool.length ? pool[0].E : fallback.E;
  // per ogni configurazione si usa la disposizione degli spin che il calcolo sferico trova piu' bassa
  // (di solito quella con piu' elettroni paralleli): e' il riferimento in cui l'approssimazione e' affidabile
  var seen = {};
  pool.forEach(function (en) {
    if (en.E > E0 + 0.08) return;
    var ck = configString(en.map, true);
    if (seen[ck]) return; seen[ck] = 1;
    var c = multipletCorrection(g, en);
    var h = hybridCorrection(g, en);
    var lev = en.E + c.dE + h;
    if (!best || lev < best.level) best = { level: lev, entry: en, corr: c, hyb: h };
  });
  return best;
}

// Calcolo completo di un atomo: solo Z in ingresso.
function computeAtom(Z, opts) {
  opts = opts || {};
  CFG.rel = opts.rel !== undefined ? !!opts.rel : true;
  CFG.xc = opts.xc || 'pbe';
  var t0 = Date.now();
  var g = makeGrid(Z, opts.dx, opts.rmax);
  var maxn = Z <= 2 ? 3 : Z <= 10 ? 4 : Z <= 36 ? 5 : Z <= 54 ? 6 : 7;
  var useLevels = opts.levels !== false;
  CFG.lv = opts.levels === 'so' ? 'so' : 'full';
  CFG.hyb = opts.hybrid !== undefined ? !!opts.hybrid : true;
  POOL = useLevels ? [] : null;
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
  if (!neu.ok) { POOL = null; return { Z: Z, ok: false }; }
  var poolN = POOL; if (useLevels) POOL = [];
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
  var poolC = POOL; POOL = null;
  if (!bestCat) return { Z: Z, ok: false };
  // livello piu' basso: energia della disposizione + correzione per multipletti e spin-orbita
  var lvN = { level: neu.E, entry: { map: neu.map, orbs: neu.orbs }, corr: { dE: 0 } }, lvC = { level: bestCat.E, entry: { map: bestCat.map }, corr: { dE: 0 } };
  if (useLevels) {
    var pn = pickLowestLevel(g, poolN, neu); if (pn) lvN = pn;
    if (Z > 1) { var pc = pickLowestLevel(g, poolC, bestCat); if (pc) lvC = pc; }
  }
  var IE = (lvC.level - lvN.level) * HA2EV;
  var chosen = lvN.entry, cmap = chosen.map;
  var Nu = 0, Nd = 0; for (var key in cmap) { if (key.slice(-1) === '0') Nu += cmap[key]; else Nd += cmap[key]; }
  // densita' radiale per il grafico: 4 pi r^2 n(r), campionata, dagli orbitali della disposizione scelta
  var ntot = new Float64Array(g.mesh + 1);
  chosen.orbs.forEach(function (o) { for (var i3 = 0; i3 <= g.mesh; i3++) ntot[i3] += o.occ * o.y[i3] * o.y[i3] / (4 * Math.PI * g.r[i3]); });
  var dens = [];
  var step = Math.max(1, Math.floor(g.mesh / 220));
  for (var i = 0; i <= g.mesh; i += step) { if (g.r[i] > 12) break; if (g.r[i] < 1e-3) continue; dens.push([+g.r[i].toPrecision(5), +(4 * Math.PI * g.r2[i] * ntot[i]).toPrecision(5)]); }
  var orbs = chosen.orbs.map(function (o) {
    var rm = 0; for (var i2 = 0; i2 <= g.mesh; i2++) rm += o.y[i2] * o.y[i2] * g.r2[i2] * g.r[i2] * g.dx;
    return { n: o.n, l: o.l, s: o.s, occ: o.occ, e: o.e, r: rm };
  });
  return {
    Z: Z, ok: true, E: lvN.level, Ecat: lvC.level, IE: IE,
    map: cmap, config: configString(cmap), configPlain: configString(cmap, true),
    catConfig: configString(lvC.entry.map), unpaired: Math.abs(Nu - Nd),
    corrN: lvN.corr.dE * HA2EV, corrC: lvC.corr.dE * HA2EV,
    zetas: (lvN.corr.zetas || []).map(function (z) { return [z[0], z[1] * HA2EV]; }),
    orbs: orbs, dens: dens, tries: tries, ms: Date.now() - t0
  };
}

if (typeof module !== 'undefined') module.exports = { CFG: CFG, computeAtom: computeAtom, makeGrid: makeGrid, scf: scf, solveRadial: solveRadial, xc: xc, startPotential: startPotential, occsFromMap: occsFromMap, configString: configString, shellTotals: shellTotals, HA2EV: HA2EV };
