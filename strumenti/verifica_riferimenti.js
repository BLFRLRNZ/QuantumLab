const E = require('../sorgenti/engine.js');
// Riferimenti LDA/LSD (NIST, Kotochigova et al.): H -0.478671 (LSD), He -2.834836, Ne -128.233481, Ar -525.946195
const cases = [
  [1, { '1,0,0': 1 }, -0.478671],
  [2, { '1,0,0': 1, '1,0,1': 1 }, -2.834836],
  [10, { '1,0,0': 1, '1,0,1': 1, '2,0,0': 1, '2,0,1': 1, '2,1,0': 3, '2,1,1': 3 }, -128.233481],
  [18, { '1,0,0': 1, '1,0,1': 1, '2,0,0': 1, '2,0,1': 1, '2,1,0': 3, '2,1,1': 3, '3,0,0': 1, '3,0,1': 1, '3,1,0': 3, '3,1,1': 3 }, -525.946195],
];
for (const [Z, map, ref] of cases) {
  for (const dx of [0.012, 0.008, 0.005]) {
    const g = E.makeGrid(Z, dx, 60);
    const V0 = E.startPotential(g);
    const t0 = Date.now();
    const r = E.scf(g, E.occsFromMap(map), [V0, V0]);
    console.log(Z, 'dx', dx, 'mesh', g.mesh, r.ok, r.E && r.E.toFixed(6), 'ref', ref, 'diff', r.E && (r.E - ref).toExponential(2), 'it', r.it, (Date.now() - t0) + 'ms',
      r.orbs && r.orbs.filter(o=>o.s===0).map(o => o.n + 'spdf'[o.l] + ':' + o.e.toFixed(5)).join(' '));
  }
}
