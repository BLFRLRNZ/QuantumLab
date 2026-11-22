const E = require('../sorgenti/engine.js');
const { BUSTA } = require('../sorgenti/busta.js');
const CORES = { He: '1s2', Ne: '1s2 2s2 2p6', Ar: '1s2 2s2 2p6 3s2 3p6' };
CORES.Kr = CORES.Ar + ' 3d10 4s2 4p6'; CORES.Xe = CORES.Kr + ' 4d10 5s2 5p6'; CORES.Rn = CORES.Xe + ' 4f14 5d10 6s2 6p6'; CORES.Hg = CORES.Xe + ' 4f14 5d10 6s2';
function expand(c){ return c.replace(/\((\w+)\)/, (m,k)=>CORES[k]).trim().split(/\s+/).reduce((t,tok)=>{const m=tok.match(/^(\d)([spdf])(\d*)$/); t[m[1]+m[2]]=(t[m[1]+m[2]]||0)+(m[3]?+m[3]:1); return t;},{}); }
function same(a,b){ const k=new Set([...Object.keys(a),...Object.keys(b)]); for(const x of k) if((a[x]||0)!==(b[x]||0)) return false; return true; }
const [z1,z2]=[+process.argv[2]||1,+process.argv[3]||36];
let ok=0,tot=0,err=[],T=0;
for (let Z=z1; Z<=z2; Z++){
  const r=E.computeAtom(Z); T+=r.ms;
  const b=BUSTA[Z-1];
  if(!r.ok){console.log(Z,'FAIL');continue;}
  const match=same(E.shellTotals(r.map),expand(b[2]));
  tot++; if(match) ok++;
  err.push(r.IE-b[3]);
  console.log(String(Z).padStart(3),b[1].padEnd(3),(match?'OK ':'-- ')+r.configPlain.padEnd(38),'mis:',b[2].padEnd(22),'IE',r.IE.toFixed(3),'mis',b[3].toFixed(3),'d',(r.IE-b[3]).toFixed(2),'spaiati',r.unpaired,'prove',r.tries,r.ms+'ms');
}
const mae=err.reduce((a,b)=>a+Math.abs(b),0)/err.length;
console.log('config ok',ok+'/'+tot,'MAE IE',mae.toFixed(3),'eV','tempo tot',(T/1000).toFixed(1)+'s');
