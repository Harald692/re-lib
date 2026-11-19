"use strict";
/* =====================================================================
   re.lib – WÄRMETAUSCHER-BERECHNUNG
   Rechenkern (HX) + Rechenkarte (hxCard) für den Bereich „Berechnungen“
   Stoffwerte: wtdata.js (CoolProp-Tabellen)
   ===================================================================== */
const HX=(()=>{
const G0=9.80665,lg=Math.log10,ln=Math.log,sq=Math.sqrt,PI=Math.PI,cl=(x,a,b)=>Math.min(Math.max(x,a),b);

/* ---------- Stoffwerte ---------- */
function ip(rows,t,logCols){t=cl(t,rows[0][0],rows[rows.length-1][0]);let i=0;while(i<rows.length-2&&t>rows[i+1][0])i++;
  const a=rows[i],b=rows[i+1],f=(t-a[0])/(b[0]-a[0]);
  return a.map((v,k)=>{const w=b[k];if(v==null)return w;if(w==null)return v;return logCols.includes(k)&&v>0&&w>0?Math.exp(ln(v)+f*(ln(w)-ln(v))):v+f*(w-v)})}
const f1=x=>(Math.round(x*10)/10).toString().replace(".",",");
function props(fl,t,W){
  if(fl.custom){const c=fl.custom;return{rho:c.rho,cp:c.cp,mu:c.mu,k:c.k,Pr:c.cp*c.mu/c.k}}
  const rows=HXDATA.liq[fl.name].t;
  if(W&&(t<rows[0][0]-0.5||t>rows[rows.length-1][0]+0.5))W.add(fl.name+": "+f1(t)+" °C liegt ausserhalb des Stoffdatenbereichs ("+rows[0][0]+" … "+rows[rows.length-1][0]+" °C) – Werte am Rand genommen.");
  const r=ip(rows,t,[3]);return{rho:r[1],cp:r[2],mu:r[3],k:r[4],Pr:r[2]*r[3]/r[4]}}
function sat(name,t,W){const d=HXDATA.ref[name],rows=d.t;
  if(W&&(t<rows[0][0]||t>rows[rows.length-1][0]))W.add(name+": Sättigungstemperatur "+f1(t)+" °C ausserhalb der Tabelle ("+rows[0][0]+" … "+rows[rows.length-1][0]+" °C).");
  const r=ip(rows,t,[1,7,8]);const S={t,p:r[1],rl:r[2],rv:r[3],r:r[4],cpl:r[5],cpv:r[6],mul:r[7],muv:r[8],kl:r[9],kv:r[10],sig:r[11]||0.01,pr:r[1]/d.pc,M:d.M,tc:d.tc};
  S.liq={rho:S.rl,cp:S.cpl,mu:S.mul,k:S.kl,Pr:S.cpl*S.mul/S.kl};S.vap={rho:S.rv,cp:S.cpv,mu:S.muv,k:S.kv,Pr:S.cpv*S.muv/S.kv};return S}

/* ---------- Feuchte Luft (p = 1,01325 bar) ---------- */
const PA=101325;
const pws=t=>t>=0?611.2*Math.exp(17.62*t/(243.12+t)):611.2*Math.exp(22.46*t/(272.62+t));
const xAir=(t,phi)=>0.622*phi*pws(t)/(PA-phi*pws(t));
const hAir=(t,x)=>1006*t+x*(2501e3+1860*t);
const hSat=t=>hAir(t,xAir(t,1));
const bis=(f,a,b,n=60)=>{let fa=f(a);for(let i=0;i<n;i++){const m=(a+b)/2,fm=f(m);if((fm>0)===(fa>0)){a=m;fa=fm}else b=m}return(a+b)/2};
const tOfHsat=h=>bis(t=>hSat(t)-h,-60,100);
const tDew=x=>bis(t=>xAir(t,1)-x,-60,100);
const phiOf=(t,x)=>x*PA/((0.622+x)*pws(t));

/* ---------- Wirkungsgrad ε(NTU) ---------- */
const E={
  counter:(N,C)=>C<1e-6?1-Math.exp(-N):Math.abs(1-C)<1e-6?N/(1+N):(1-Math.exp(-N*(1-C)))/(1-C*Math.exp(-N*(1-C))),
  parallel:(N,C)=>(1-Math.exp(-N*(1+C)))/(1+C),
  crossUU:(N,C)=>C<1e-6?1-Math.exp(-N):1-Math.exp(Math.pow(N,0.22)/C*(Math.exp(-C*Math.pow(N,0.78))-1)),
  crossMaxMixed:(N,C)=>C<1e-6?1-Math.exp(-N):(1/C)*(1-Math.exp(-C*(1-Math.exp(-N)))),
  crossMinMixed:(N,C)=>C<1e-6?1-Math.exp(-N):1-Math.exp(-(1/C)*(1-Math.exp(-C*N))),
  shell12:(N,C)=>{const s=sq(1+C*C),e=Math.exp(-N*s);return 2/(1+C+s*(1+e)/(1-e))},
  multi:(n,ep,C)=>{if(n<=1)return ep;if(Math.abs(1-C)<1e-6)return n*ep/(1+(n-1)*ep);const a=Math.pow((1-ep*C)/(1-ep),n);return(a-1)/(a-C)}
};

/* ---------- Korrelationen: einphasig ---------- */
function nuLam(Re,Pr,dL){const z=Re*Pr*dL,n2=1.615*Math.cbrt(z),n3=Math.pow(2/(1+22*Pr),1/6)*sq(z);return Math.cbrt(3.66**3+0.7**3+(n2-0.7)**3+n3**3)}
function nuTurb(Re,Pr,dL){const x=Math.pow(1.8*lg(Re)-1.5,-2);return(x/8)*Re*Pr/(1+12.7*sq(x/8)*(Math.pow(Pr,2/3)-1))*(1+Math.pow(dL,2/3))}
function nuTube(Re,Pr,dL){if(Re<=2300)return nuLam(Re,Pr,dL);if(Re>=1e4)return nuTurb(Re,Pr,dL);const y=(Re-2300)/7700;return(1-y)*nuLam(2300,Pr,dL)+y*nuTurb(1e4,Pr,dL)}
function lamTube(Re){const pt=R=>Math.pow(0.79*ln(R)-1.64,-2);if(Re<2300)return 64/Re;if(Re<4000){const y=(Re-2300)/1700;return(1-y)*64/2300+y*pt(4000)}return pt(Re)}
// Rohr, einphasig: m gesamt, n parallele Rohre, Länge L je Strömungsweg
function tubeSP(P,m,d,L,n){const A=PI*d*d/4,w=m/(P.rho*A*n),Re=P.rho*w*d/P.mu,Nu=nuTube(Re,P.Pr,d/L);
  return{h:Nu*P.k/d,Re,w,Nu,dp:lamTube(Re)*L/d*P.rho*w*w/2,q2:P.rho*w*w/2}}
// Ringspalt (Doppelrohr), Wärmeübergang am Innenrohr (VDI-Wärmeatlas G2)
function annSP(P,m,Di,da,L){const a=da/Di,Dh=Di-da,A=PI*(Di*Di-da*da)/4,w=m/(P.rho*A),Re=P.rho*w*Dh/P.mu,z=Re*P.Pr*Dh/L;
  const lam=R=>Math.cbrt((3.66+1.2*Math.pow(a,-0.8))**3+(1.615*(1+0.14*Math.pow(a,-0.5))*Math.cbrt(R*P.Pr*Dh/L))**3);
  const tur=R=>nuTurb(R,P.Pr,Dh/L)*0.86*Math.pow(a,-0.16);
  const Nu=Re<=2300?lam(Re):Re>=1e4?tur(Re):(1-(Re-2300)/7700)*lam(2300)+(Re-2300)/7700*tur(1e4);
  return{h:Nu*P.k/Dh,Re,w,Nu,dp:lamTube(Re)*L/Dh*P.rho*w*w/2,q2:P.rho*w*w/2}}
// Mantelraum nach Kern
function kern(P,m,g){const De=g.layout==="sq"?4*(g.Pt*g.Pt-PI*g.do*g.do/4)/(PI*g.do):4*(sq(3)/4*g.Pt*g.Pt-PI*g.do*g.do/8)/(PI*g.do/2);
  const As=g.Ds*(g.Pt-g.do)*g.B/g.Pt,Gs=m/As,Re=Math.max(Gs*De/P.mu,1),Nu=0.36*Math.pow(Re,0.55)*Math.cbrt(P.Pr),Nb=Math.max(1,Math.round(g.L/g.B)-1);
  const f=Math.exp(0.576-0.19*ln(Re));return{h:Nu*P.k/De,Re,w:Gs/P.rho,Nu,De,dp:f*Gs*Gs*g.Ds*(Nb+1)/(2*P.rho*De)*g.shells,Nb}}
// Plattenwärmetauscher nach Martin (VDI-Wärmeatlas)
function xiMartin(Re,phi){const xi0=Re<2000?64/Re:Math.pow(1.8*lg(Re)-1.5,-2),xi1=Re<2000?597/Re+3.85:39/Math.pow(Re,0.289),c=Math.cos(phi);
  return 1/Math.pow(c/sq(0.18*Math.tan(phi)+0.36*Math.sin(phi)+xi0/c)+(1-c)/sq(3.8*xi1),2)}
function pheSP(P,m,n,g){const Dh=2*g.b/g.Phi,w=m/(P.rho*n*g.b*g.B),Re=Math.max(P.rho*w*Dh/P.mu,1),xi=xiMartin(Re,g.phi);
  const Nu=0.122*Math.cbrt(P.Pr)*Math.pow(xi*Re*Re*Math.sin(2*g.phi),0.374),wp=m/(P.rho*PI*g.Dp*g.Dp/4);
  return{h:Nu*P.k/Dh,Re,w,Nu,dp:xi*g.Lp/Dh*P.rho*w*w/2+1.5*P.rho*wp*wp/2,wp}}

/* ---------- Korrelationen: zweiphasig ---------- */
const xs=(a,b,n=12)=>{const r=[];for(let i=0;i<n;i++)r.push(a+(b-a)*(i+0.5)/n);return r};
const avg=(a,b,fn)=>{const v=xs(Math.min(a,b),Math.max(a,b));return v.reduce((s,x)=>s+fn(x),0)/v.length};
// Kondensation im Rohr: Shah (1979), bei kleinen Massenströmen mind. Schwerkraftkondensation nach Chato
function condTube(S,G,D,xa,xb,q){const Re=G*D/S.mul,hlo=0.023*Math.pow(Re,0.8)*Math.pow(S.liq.Pr,0.4)*S.kl/D;
  const C=0.555*Math.pow(S.rl*(S.rl-S.rv)*G0*S.kl**3*S.r/(S.mul*D),0.25),hch=Math.pow(C,4/3)/Math.cbrt(Math.max(q||10000,100));
  return avg(xa,xb,x=>Math.max(hlo*(Math.pow(1-x,0.8)+3.8*Math.pow(x,0.76)*Math.pow(1-x,0.04)/Math.pow(S.pr,0.38)),hch*(x>0.05?1:x/0.05+0.2)))}
// Strömungssieden im Rohr (Gungor-Winterton 1987), Dryout ab x ≈ 0,85 vereinfacht
function evapTube(S,G,D,xa,xb,q){const Bo=q/(G*S.r);
  return avg(xa,xb,x=>{const xe=Math.min(x,0.85),Re=Math.max(G*(1-xe)*D/S.mul,50),hl=0.023*Math.pow(Re,0.8)*Math.pow(S.liq.Pr,0.4)*S.kl/D;
    let h=Math.max(hl*(1+3000*Math.pow(Bo,0.86)+1.12*Math.pow(xe/(1-xe),0.75)*Math.pow(S.rl/S.rv,0.41)),pool(S,q));
    if(x>0.85){const hv=tubeSP(S.vap,G*PI*D*D/4,D,1e9,1).h;h=h+(hv-h)*(x-0.85)/0.15}return h})}
// Zweiphasen-Druckverlust Rohr (Müller-Steinhagen/Heck)
function dpTP(S,G,D,L,xa,xb){const fl=lamTube(G*D/S.mul),fg=lamTube(G*D/S.muv),A=fl*G*G/(2*D*S.rl),B=fg*G*G/(2*D*S.rv);
  return L*avg(xa,xb,x=>{const Gx=A+2*(B-A)*x;return Gx*Math.cbrt(1-x)+B*x*x*x})}
// Platte: Verdampfung / Kondensation nach Longo et al. (2015, gelötete Platten-WT)
//  Verdampfung: max(Blasensieden 0,58·Φ·α_Cooper ; erzwungene Konvektion 0,122·Φ·λ/Dh·Re_eq^0,8·Pr^1/3)
//  Kondensation: Re_eq < 1600 Schwerkraft (Nusselt, 0,943·Φ), sonst 1,875·Φ·λ/Dh·Re_eq^0,445·Pr^1/3
const geq=(S,G,x)=>G*((1-x)+x*sq(S.rl/S.rv));
function pheTP(S,G,Dh,xa,xb,q,cond,g){const Phi=g.Phi,Lp=g.Lp;return avg(xa,xb,x=>{const Re=geq(S,G,x)*Dh/S.mul,fc=S.kl/Dh*Math.cbrt(S.liq.Pr);
  if(!cond)return Math.max(0.58*Phi*pool(S,q),0.122*Phi*fc*Math.pow(Re,0.8));
  if(Re>=1600)return 1.875*Phi*fc*Math.pow(Re,0.445);
  const C=0.943*Phi*Math.pow(S.rl*(S.rl-S.rv)*G0*S.r*S.kl**3/(S.mul*Lp),0.25);return Math.pow(C,4/3)/Math.cbrt(Math.max(q,100))})}
function dpPheTP(S,G,g,xa,xb){const Dh=2*g.b/g.Phi;return avg(xa,xb,x=>{const rh=1/(x/S.rv+(1-x)/S.rl),mu=1/(x/S.muv+(1-x)/S.mul),Re=G*Dh/mu;return xiMartin(Re,g.phi)*g.Lp/Dh*G*G/(2*rh)})}
// Mantelseitige Kondensation am waagrechten Rohrbündel (Nusselt, Reihenkorrektur Kern)
function condBundle(S,dO,Nv,q){const c=0.725*Math.pow(S.rl*(S.rl-S.rv)*G0*S.r*S.kl**3/(S.mul*dO),0.25)*Math.pow(Math.max(Nv,1),-1/6);
  return Math.pow(c,4/3)/Math.cbrt(Math.max(q,100))}         // h = c·ΔT^-1/4 mit ΔT = q/h
// Behältersieden (Cooper) – überfluteter Verdampfer
const pool=(S,q)=>55*Math.pow(S.pr,0.12)*Math.pow(-lg(S.pr),-0.55)*Math.pow(S.M,-0.5)*Math.pow(Math.max(q,100),0.67);

/* ---------- Bauarten (Geometrie-Modelle) ---------- */
// Jedes Modell: Aref (Bezugsfläche), Rw (Wandwiderstand bezogen auf Aref), side(loc,st) → {h,A,eta,dp,Re,w}, eps(NTU,Cr,cHot,cCold)
function mPlate(c){const g=c.geo,nh=Math.ceil((g.np-1)/2),nc=Math.floor((g.np-1)/2),A=Math.max(g.np-2,1)*g.B*g.Lp*g.Phi,Dh=2*g.b/g.Phi;
  return{A,Aref:A,refLabel:"wärmeübertragende Plattenfläche",Rw:g.s/g.kw,loc:r=>r,info:{nh,nc,Dh},
    side(loc,st){const n=loc==="hot"?nh:nc;
      if(st.kind==="sp"){const r=pheSP(st.P,st.m,n,g);return{...r,A,eta:1}}
      const G=st.m/(n*g.b*g.B),S=st.S,h=pheTP(S,G,Dh,st.xa,st.xb,st.q||5000,st.kind==="cond",g);
      return{h,A,eta:1,G,Re:G*Dh/S.mul,dp:dpPheTP(S,G,g,st.xa,st.xb),w:G/S.rl}},
    eps:(N,C)=>g.flow==="parallel"?E.parallel(N,C):E.counter(N,C),series:g.flow==="parallel"?-1:1}}
function mShell(c){const g=c.geo,sh=g.shells,Ao=PI*g.do*g.L*g.Nt*sh,Ai=PI*g.di*g.L*g.Nt*sh,nTp=g.Nt/g.npp,Lp=g.L*g.npp*sh,Nv=Math.max(1,Math.round(sq(g.Nt)));
  const CL=g.layout==="sq"?1:0.87,CTP=g.npp===1?0.93:g.npp===2?0.9:0.85,NtMax=Math.floor(0.785*(CTP/CL)*g.Ds*g.Ds/(g.Pt*g.Pt));
  const tubeLoc=c.hotLoc==="tube"?"hot":"cold";
  return{A:Ao,Aref:Ao,refLabel:"Rohraussenfläche",Rw:g.do*ln(g.do/g.di)/(2*g.kw),loc:r=>r===tubeLoc?"tube":"shell",info:{NtMax,Nv},
    side(loc,st){
      if(loc==="tube"){
        if(st.kind==="sp"){const r=tubeSP(st.P,st.m,g.di,Lp,nTp);r.dp+=4*g.npp*sh*r.q2;return{...r,A:Ai,eta:1}}
        const G=st.m/(nTp*PI*g.di*g.di/4),S=st.S;
        const h=st.kind==="cond"?condTube(S,G,g.di,st.xa,st.xb,st.q):evapTube(S,G,g.di,st.xa,st.xb,st.q||10000);
        return{h,A:Ai,eta:1,G,Re:G*g.di/S.mul,dp:dpTP(S,G,g.di,Lp,st.xa,st.xb),w:G/S.rl}}
      if(st.kind==="sp"){const r=kern(st.P,st.m,{...g,shells:sh});return{...r,A:Ao,eta:1}}
      const S=st.S,h=st.kind==="cond"?condBundle(S,g.do,Nv,st.q||15000):pool(S,st.q||15000);return{h,A:Ao,eta:1,Re:NaN,dp:0,w:NaN}},
    eps(N,C){if(g.npp===1)return E.counter(N,C);if(sh===2)return E.multi(2,E.shell12(N/2,C),C);return E.shell12(N,C)},series:1}}
function mDpipe(c){const g=c.geo,L=g.Lel*g.n,Ao=PI*g.da*L,Ai=PI*g.di*L,inLoc=c.hotLoc==="in"?"hot":"cold";
  return{A:Ao,Aref:Ao,refLabel:"Aussenfläche Innenrohr",Rw:g.da*ln(g.da/g.di)/(2*g.kw),loc:r=>r===inLoc?"in":"ann",info:{},
    side(loc,st){
      if(loc==="in"){if(st.kind==="sp"){const r=tubeSP(st.P,st.m,g.di,L,1);r.dp+=(g.n-1)*1.2*r.q2;return{...r,A:Ai,eta:1}}
        const G=st.m/(PI*g.di*g.di/4),S=st.S,h=st.kind==="cond"?condTube(S,G,g.di,st.xa,st.xb,st.q):evapTube(S,G,g.di,st.xa,st.xb,st.q||10000);
        return{h,A:Ai,eta:1,G,Re:G*g.di/S.mul,dp:dpTP(S,G,g.di,L,st.xa,st.xb),w:G/S.rl}}
      const r=annSP(st.P,st.m,g.Di,g.da,L);r.dp+=(g.n-1)*2.5*r.q2;return{...r,A:Ao,eta:1}},
    eps:(N,C)=>g.flow==="parallel"?E.parallel(N,C):E.counter(N,C),series:g.flow==="parallel"?-1:1}}
// Lamellenrohr (Fin-and-Tube): Luftseite Wang/Chi/Chang (2000) für glatte Lamellen, Faktoren für andere Lamellenformen
const FINF={plain:[1,1],wavy:[1.25,1.45],louver:[1.55,1.7],slit:[1.65,1.9]};
function mFin(c){const g=c.geo,Dc=g.do+2*g.delta,di=g.do-2*g.tw,N=g.N,H=g.Nt*g.Pt,Afr=H*g.L,Nf=g.L/g.Fp;
  const Af=2*Nf*(H*N*g.Pl-g.Nt*N*PI*Dc*Dc/4),At=PI*Dc*g.L*g.Nt*N*(1-g.delta/g.Fp),Ao=Af+At,Ai=PI*di*g.L*g.Nt*N;
  const a=g.Pt-Dc,cd=2*(sq((g.Pt/2)**2+g.Pl**2)-Dc),sig=(g.stag?Math.min(a,cd):a)/g.Pt*(g.Fp-g.delta)/g.Fp,Amin=sig*Afr,Dh=4*Amin*N*g.Pl/Ao;
  const ff=FINF[g.fin]||FINF.plain,fi=g.stag?[1,1]:[0.85,0.8],mi=g.inner==="micro";
  const Lc=g.Nt*N*g.L/g.circ,bends=Math.max(g.Nt*N/g.circ-1,0),airLoc=c.airSide;
  // Lamellenwirkungsgrad nach Schmidt
  const etaFin=h=>{const r=Dc/2,Xm=g.Pt/2,XL=g.stag?sq((g.Pt/2)**2+g.Pl**2)/2:g.Pl/2,Rr=g.stag?1.27*Xm/r*sq(Math.max(XL/Xm-0.3,0.01)):1.28*Xm/r*sq(Math.max(XL/Xm-0.2,0.01));
    const ph=(Rr-1)*(1+0.35*ln(Rr)),m=sq(2*h/(g.kf*g.delta)),z=m*r*ph;return z<1e-6?1:Math.tanh(z)/z};
  return{A:Ao,Aref:Ao,refLabel:"luftseitige Aussenfläche",Rw:Ao*ln(g.do/di)/(2*PI*g.kt*g.L*g.Nt*N),loc:r=>r===airLoc?"air":"tube",
    info:{Afr,Amin,Ao,Ai,Af,Dh,Dc,sig,Lc},
    side(loc,st){
      if(loc==="air"){const P=st.P,Gc=st.m/Amin,Re=Gc*Dc/P.mu,lr=ln(Re);let j,f;
        const r1=g.Fp/Dc,r2=g.Fp/Dh,r3=g.Fp/g.Pt;
        if(N===1){const P1=1.9-0.23*lr,P2=-0.236+0.126*lr;j=0.108*Math.pow(Re,-0.29)*Math.pow(g.Pt/g.Pl,P1)*Math.pow(r1,-1.084)*Math.pow(r2,-0.786)*Math.pow(r3,P2)}
        else{const P3=-0.361-0.042*N/lr+0.158*ln(N*Math.pow(r1,0.41)),P4=-1.224-0.076*Math.pow(g.Pl/Dh,1.42)/lr,P5=-0.083+0.058*N/lr,P6=-5.735+1.21*ln(Re/N);
          j=0.086*Math.pow(Re,P3)*Math.pow(N,P4)*Math.pow(r1,P5)*Math.pow(r2,P6)*Math.pow(r3,-0.93)}
        const F1=-0.764+0.739*g.Pt/g.Pl+0.177*r1-0.00758/N,F2=-15.689+64.021/lr,F3=1.696-15.695/lr;f=0.0267*Math.pow(Re,F1)*Math.pow(g.Pt/g.Pl,F2)*Math.pow(r1,F3);
        j*=ff[0]*fi[0];f*=ff[1]*fi[1];const h=j*Gc*P.cp/Math.pow(P.Pr,2/3),ef=etaFin(h),eta=1-Af/Ao*(1-ef);
        return{h,A:Ao,eta,etaFin:ef,Re,w:st.m/(P.rho*Afr),wmax:Gc/P.rho,dp:f*Gc*Gc/(2*P.rho)*Ao/Amin,j,f}}
      const fh=mi?(st.kind==="sp"?1.5:1.8):1,fd=mi?1.3:1;
      if(st.kind==="sp"){const r=tubeSP(st.P,st.m,di,Lc,g.circ);r.h*=fh;r.dp=r.dp*fd+bends*1.0*r.q2;return{...r,A:Ai,eta:1}}
      const G=st.m/(g.circ*PI*di*di/4),S=st.S,h=(st.kind==="cond"?condTube(S,G,di,st.xa,st.xb,st.q):evapTube(S,G,di,st.xa,st.xb,st.q||8000))*fh;
      return{h,A:Ai,eta:1,G,Re:G*di/S.mul,dp:dpTP(S,G,di,Lc,st.xa,st.xb)*fd,w:G/S.rl}},
    eps(Nt_,C,cAir,cTube){const mixedMax=cTube>=cAir,ep=E[mixedMax?"crossMaxMixed":"crossMinMixed"](Nt_/N,C);return E.multi(N,ep,C)},airPar:true}}
// Microchannel: Luftseite Kiemenlamelle (Kim & Bullard 2002), Flachrohr mit Mehrkanal-Profil
function mMchx(c){const g=c.geo,Nt=Math.max(1,Math.floor(g.H/g.Tp)),Fl=g.Tp-g.Dt,Ng=Nt,sl=g.slabs;
  const Af=Ng*(g.W/g.Fp)*2*Fl*g.Td*sl,At=Nt*2*(g.Td+g.Dt)*g.W*sl,Ao=Af+At,Amin=Ng*g.W*Fl*(1-g.delta/g.Fp),Afr=g.W*g.H;
  const Ap=g.pw*g.ph,Pp=2*(g.pw+g.ph),dh=4*Ap/Pp,Ai=Nt*g.Np*Pp*g.W*sl,nT=Math.max(1,Nt/g.passes),Lp=g.W*g.passes*sl,Ll=g.Ll||0.85*Fl,tw=Math.max((g.Dt-g.ph)/2,1e-4),airLoc=c.airSide;
  return{A:Ao,Aref:Ao,refLabel:"luftseitige Aussenfläche",Rw:tw/g.kt*Ao/Ai,loc:r=>r===airLoc?"air":"tube",info:{Nt,Fl,Afr,Amin,Ao,Ai,Af,dh},
    side(loc,st){
      if(loc==="air"){const P=st.P,Gc=st.m/Amin,Re=Gc*g.Lp/P.mu,t=g.theta/90,L=g.Lp;
        const j=Math.pow(Re,-0.487)*Math.pow(t,0.257)*Math.pow(g.Fp/L,-0.13)*Math.pow(Fl/L,-0.29)*Math.pow(g.Td/L,-0.235)*Math.pow(Ll/L,0.68)*Math.pow(g.Tp/L,-0.279)*Math.pow(g.delta/L,-0.05);
        const f=Math.pow(Re,-0.781)*Math.pow(t,0.444)*Math.pow(g.Fp/L,-1.682)*Math.pow(Fl/L,-1.22)*Math.pow(g.Td/L,0.818)*Math.pow(Ll/L,1.97);
        const h=j*Gc*P.cp/Math.pow(P.Pr,2/3),m=sq(2*h/(g.kf*g.delta)),z=m*Fl/2,ef=Math.tanh(z)/z,eta=1-Af/Ao*(1-ef);
        return{h,A:Ao,eta,etaFin:ef,Re,w:st.m/(P.rho*Afr),wmax:Gc/P.rho,dp:f*Gc*Gc/(2*P.rho)*(Ao/sl)/Amin*sl,j,f}}
      if(st.kind==="sp"){const A1=Ap,w=st.m/(st.P.rho*A1*nT*g.Np),Re=st.P.rho*w*dh/st.P.mu,Nu=nuTube(Re,st.P.Pr,dh/Lp);
        return{h:Nu*st.P.k/dh,Re,w,A:Ai,eta:1,dp:lamTube(Re)*Lp/dh*st.P.rho*w*w/2+(g.passes*sl-1)*1.5*st.P.rho*w*w/2}}
      const G=st.m/(nT*g.Np*Ap),S=st.S,h=st.kind==="cond"?condTube(S,G,dh,st.xa,st.xb,st.q):evapTube(S,G,dh,st.xa,st.xb,st.q||8000);
      return{h,A:Ai,eta:1,G,Re:G*dh/S.mul,dp:dpTP(S,G,dh,Lp,st.xa,st.xb),w:G/S.rl}},
    eps(N_,C,cAir,cTube){const one=(n,c)=>E[cTube>=cAir?"crossMaxMixed":"crossMinMixed"](n,c);return sl>1?E.multi(sl,one(N_/sl,C),C):one(N_,C)},airPar:true}}
// k-Wert vorgegeben
function mKman(c){const k=c.kman;return{A:k.A,Aref:k.A,refLabel:"Wärmeübertragungsfläche",Rw:0,loc:r=>r,info:{},Ufix:k.U,
  side:()=>({h:NaN,A:k.A,eta:1,dp:NaN,Re:NaN,w:NaN}),
  eps:(N,C)=>({counter:E.counter,parallel:E.parallel,cross:E.crossUU,shell:E.shell12})[k.arr](N,C),series:k.arr==="parallel"?-1:1}}
const MODELS={plate:mPlate,shell:mShell,dpipe:mDpipe,fin:mFin,mchx:mMchx,kman:mKman};

/* ---------- Wärmedurchgang ---------- */
// liefert U (bezogen auf Aref) und Teilwiderstände (m²K/W bezogen auf Aref)
function Ucalc(M,c,sH,sC){if(M.Ufix)return{U:M.Ufix,Rh:0,Rc:0,Rw:0,Rf:1/M.Ufix};
  const A=M.Aref,Rh=A/(sH.eta*sH.h*sH.A),Rc=A/(sC.eta*sC.h*sC.A),Rf=c.Rf.hot*A/sH.A+c.Rf.cold*A/sC.A;
  return{U:1/(Rh+Rc+M.Rw+Rf),Rh,Rc,Rw:M.Rw,Rf}}
const cpAir=1006;

/* ---------- Lösung ohne Phasenwechsel ---------- */
function solveSP(c,M,W){
  const H=c.hot,C=c.cold;let Tho=H.tin-0.3*(H.tin-C.tin),Tco=C.tin+0.3*(H.tin-C.tin),out;
  for(let it=0;it<12;it++){
    const Ph=props(H.fl,(H.tin+Tho)/2,it>9?W:null),Pc=props(C.fl,(C.tin+Tco)/2,it>9?W:null);
    const sH=M.side(M.loc("hot"),{kind:"sp",P:Ph,m:H.m}),sC=M.side(M.loc("cold"),{kind:"sp",P:Pc,m:C.m});
    const u=Ucalc(M,c,sH,sC),Ch=H.m*Ph.cp,Cc=C.m*Pc.cp,Cmin=Math.min(Ch,Cc),Cr=Cmin/Math.max(Ch,Cc),NTU=u.U*M.Aref/Cmin;
    const cAir=M.airPar?(c.airSide==="hot"?Ch:Cc):0,cTube=M.airPar?(c.airSide==="hot"?Cc:Ch):0;
    const eps=M.eps(NTU,Cr,cAir,cTube),Q=eps*Cmin*(H.tin-C.tin),nTho=H.tin-Q/Ch,nTco=C.tin+Q/Cc;
    out={Q,Tho:nTho,Tco:nTco,eps,NTU,Cr,U:u.U,R:u,sH,sC,Ph,Pc,Ch,Cc};
    if(Math.abs(nTho-Tho)<1e-4&&Math.abs(nTco-Tco)<1e-4)break;Tho=nTho;Tco=nTco}
  // Luftkühler: Entfeuchtung prüfen (vollständig nasser Kühler, Enthalpie-Methode nach Braun/Threlkeld)
  if(M.airPar&&c.airSide==="hot"&&c.air&&c.air.rh>0){const xin=xAir(H.tin,c.air.rh),td=tDew(xin),u=out.R;
    const Ts=C.tin+(out.Tho-C.tin)*(u.Rc+u.Rw+u.Rf)/(u.Rh+u.Rc+u.Rw+u.Rf);
    if(Ts<td){const hin=hAir(H.tin,xin),sH=out.sH,sC=out.sC,Pc=out.Pc;let Two=out.Tco,wet;
      for(let it=0;it<30;it++){const cs=Math.abs(Two-C.tin)>0.05?(hSat(Two)-hSat(C.tin))/(Two-C.tin):(hSat(C.tin+0.05)-hSat(C.tin-0.05))/0.1;
        const Ro=M.Aref/(sH.eta*sH.h*sH.A),Ri=u.Rc+u.Rw+u.Rf,Ustar=1/(cpAir*Ro+cs*Ri);
        const Ca=H.m,Cw=C.m*Pc.cp/cs,Cmin=Math.min(Ca,Cw),Cr=Cmin/Math.max(Ca,Cw),NTU=Ustar*M.Aref/Cmin;
        const eps=M.eps(NTU,Cr,Ca,Cw),Q=eps*Cmin*(hin-hSat(C.tin)),n=C.tin+Q/(C.m*Pc.cp);wet={Q,Two:n,eps,NTU,cs,Ustar};if(Math.abs(n-Two)<1e-4)break;Two=n}
      if(wet.Q>out.Q){const hout=hin-wet.Q/H.m,NTUo=sH.eta*sH.h*sH.A/(H.m*cpAir),hse=hin-(hin-hout)/(1-Math.exp(-NTUo)),Tse=tOfHsat(hse),
        Tao=Tse+(H.tin-Tse)*Math.exp(-NTUo);let xo=(hout-1006*Tao)/(2501e3+1860*Tao);xo=Math.min(xo,xin);
        out={...out,Q:wet.Q,Tho:Tao,Tco:wet.Two,eps:wet.eps,NTU:wet.NTU,wet:true,air:{xin,xout:xo,phiOut:Math.min(phiOf(Tao,xo),1),hin,hout,Qs:H.m*(cpAir+1860*xin)*(H.tin-Tao),cond:H.m*(xin-xo),td}};}}
    if(!out.air)out.air={xin,xout:xin,phiOut:phiOf(out.Tho,xin),td,Qs:out.Q,cond:0}}
  return out}

/* ---------- Lösung mit Phasenwechsel (Zonenmodell) ---------- */
function solveTP(c,M,W){
  const cond=c.phase==="cond",rR=cond?"hot":"cold",sR=cond?"cold":"hot",R=c.ref,S=sat(R.name,R.tsat,W),T=R.tsat,sec=c[sR],rl=M.loc(rR),sl=M.loc(sR);
  let zones;
  if(cond){const tg=Math.max(R.tgas,T),sc=Math.max(R.sc,0);zones=[{k:"dsh",n:"Enthitzung",sp:S.vap,Ta:tg,Tb:T,dh:S.cpv*(tg-T)},{k:"tp",n:"Kondensation",xa:1,xb:0,dh:S.r},{k:"sc",n:"Unterkühlung",sp:S.liq,Ta:T,Tb:T-sc,dh:S.cpl*sc}]}
  else{const xin=cl(S.cpl*(R.tliq-T)/S.r,0,0.95),shh=R.flooded?0:Math.max(R.sh,0);zones=[{k:"tp",n:"Verdampfung",xa:xin,xb:1,dh:(1-xin)*S.r},{k:"sh",n:"Überhitzung",sp:S.vap,Ta:T,Tb:T+shh,dh:S.cpv*shh}];R.xin=xin}
  zones=zones.filter(z=>z.dh>1e-9);
  const air=M.airPar&&sl==="air",wetCheck=air&&!cond&&c.air&&c.air.rh>0,xin=wetCheck?xAir(sec.tin,c.air.rh):0,td=wetCheck?tDew(xin):-99;
  let Tsm=sec.tin+(cond?3:-3),q={tp:cond?15000:8000},res=null;
  for(let it=0;it<6;it++){
    const Ps=props(sec.fl,Tsm,it===5?W:null),sS=M.side(sl,{kind:"sp",P:Ps,m:sec.m}),Cs=sec.m*Ps.cp;
    const ev=mr=>{const Z=[];let sum=0,Tsec=sec.tin;const order=M.series===-1?zones:[...zones].reverse();
      const det={};for(const z of zones){const st=z.sp?{kind:"sp",P:z.sp,m:mr}:{kind:cond?"cond":"evap",S,m:mr,xa:z.xa,xb:z.xb,q:q.tp};
        const sr=M.side(rl,st),u=cond?Ucalc(M,c,sr,sS):Ucalc(M,c,sS,sr);det[z.k]={sr,u,Q:mr*z.dh}}
      if(air){for(const z of zones){const d=det[z.k],Q=d.Q;let f;
          if(!z.sp){const Uw=d.u;let wet=false,fd;const Ts=T+(sec.tin-T)*(cond?(Uw.Rh+Uw.Rw+Uw.Rf):(Uw.Rc+Uw.Rw+Uw.Rf))/(Uw.Rh+Uw.Rc+Uw.Rw+Uw.Rf);
            const epsD=1-Math.exp(-d.u.U*M.Aref/Cs);fd=Q/(epsD*Cs*Math.abs(sec.tin-T));f=fd;
            if(wetCheck&&Ts<td){const cs=(hSat(td)-hSat(T))/Math.max(td-T,0.1),Ro=M.Aref/(sS.eta*sS.h*sS.A),Ri=1/d.u.U-Ro,Us=1/(cpAir*Ro+cs*Ri),
              hin=hAir(sec.tin,xin),epsW=1-Math.exp(-Us*M.Aref/sec.m),fw=Q/(epsW*sec.m*(hin-hSat(T)));if(fw<fd){f=fw;wet=true}}
            d.wet=wet}
          else{const Cr_=mr*z.sp.cp,dT=Math.abs(z.Ta-sec.tin);if(Q>=Cr_*dT*0.999){f=Infinity}else f=bis(ff=>{const Ca=ff*Cs,Cmin=Math.min(Ca,Cr_),Cr=Cmin/Math.max(Ca,Cr_);
              return E.crossUU(d.u.U*M.Aref*ff/Cmin,Cr)*Cmin*dT-Q},1e-9,50)}
          d.f=f;sum+=f}}
      else{for(const z of order){const d=det[z.k],Q=d.Q,Tn=Tsec+(cond?1:-1)*Q/Cs;
          const ra=M.series===-1?z.Ta??T:z.Tb??T,rb=M.series===-1?z.Tb??T:z.Ta??T;      // Kältemitteltemp. an Ende „Sekundär-Eintritt“ / „-Austritt“ der Zone
          const d1=cond?ra-Tsec:Tsec-ra,d2=cond?rb-Tn:Tn-rb;let f;
          if(d1<=0.01||d2<=0.01)f=Infinity;else{const lm=Math.abs(d1-d2)<1e-6?d1:(d1-d2)/ln(d1/d2);f=Q/(d.u.U*lm)/M.Aref}
          d.f=f;sum+=f;d.Tsa=Tsec;d.Tsb=Tn;Tsec=Tn}}
      return{sum,det}};
    // Kältemittel-Massenstrom so bestimmen, dass die Zonen genau die vorhandene Fläche füllen
    let hi=1e-5,eh=ev(hi);while(eh.sum<1&&hi<1e4){hi*=2;eh=ev(hi)}
    let lo=hi/2;for(let i=0;i<70;i++){const m=(lo+hi)/2,e=ev(m);if(e.sum<1)lo=m;else hi=m}
    const mr=lo,e=ev(mr),Q=zones.reduce((s,z)=>s+mr*z.dh,0);
    const Tso=sec.tin+(cond?1:-1)*Q/Cs;Tsm=(sec.tin+Tso)/2;
    const tz=zones.find(z=>!z.sp);if(tz){const d=e.det[tz.k],qn=mr*tz.dh/(Math.max(d.f,1e-6)*d.sr.A);if(isFinite(qn)&&qn>0)q.tp=0.5*q.tp+0.5*qn}
    const Ubar=zones.reduce((s,z)=>s+e.det[z.k].u.U*e.det[z.k].f,0)/Math.max(zones.reduce((s,z)=>s+e.det[z.k].f,0),1e-9);
    res={Q,mr,S,zones:zones.map(z=>({k:z.k,n:z.n,Q:mr*z.dh,f:e.det[z.k].f,U:e.det[z.k].u.U,h:e.det[z.k].sr.h,dp:e.det[z.k].sr.dp*e.det[z.k].f,wet:e.det[z.k].wet})),
      sS,Ps,Cs,Tso,U:Ubar,sumf:e.sum,det:e.det,q:q.tp}}
  // Luftzustand nach Verdampfer
  if(air&&!cond){const hin=hAir(sec.tin,c.air?xAir(sec.tin,c.air.rh||0):0),x0=c.air&&c.air.rh?xin:0;let hout=hin,xo=0,To=0,F=0;
    const NTUo=res.sS.eta*res.sS.h*res.sS.A/(sec.m*cpAir);
    for(const z of res.zones){const qz=z.Q/(z.f*sec.m),f=z.f;let Tz,xz;
      if(z.wet){const hz=hin-qz,hse=hin-qz/(1-Math.exp(-NTUo)),Tse=tOfHsat(hse);Tz=Tse+(sec.tin-Tse)*Math.exp(-NTUo);xz=Math.min((hz-1006*Tz)/(2501e3+1860*Tz),x0)}
      else{Tz=sec.tin-qz/(cpAir+1860*x0);xz=x0}To+=f*Tz;xo+=f*xz;F+=f}
    To/=F;xo/=F;res.air={xin:x0,xout:xo,phiOut:Math.min(phiOf(To,xo),1),td,Qs:sec.m*(cpAir+1860*x0)*(sec.tin-To),cond:sec.m*(x0-xo)};res.Tso=To}
  const dpR=res.zones.reduce((s,z)=>s+z.dp,0),sR_={h:res.zones.find(z=>!z.sp&&z.k==="tp")?.h,dp:dpR,Re:NaN,w:NaN};
  return{...res,Q:res.Q,cond,rR,sR,dpRef:dpR}}

/* ---------- Hauptfunktion ---------- */
function solve(c){const W=new Set(),M=MODELS[c.type](c);
  if(c.phase==="none"){const r=solveSP(c,M,W);return finish(c,M,W,{mode:"sp",...r})}
  const r=solveTP(c,M,W);return finish(c,M,W,{mode:"tp",...r})}
function finish(c,M,W,r){
  const o={M,warn:W,mode:r.mode,Q:r.Q,Aref:M.Aref,refLabel:M.refLabel,info:M.info};
  if(r.mode==="sp"){Object.assign(o,{Thi:c.hot.tin,Tho:r.Tho,Tci:c.cold.tin,Tco:r.Tco,mh:c.hot.m,mc:c.cold.m,U:r.U,eps:r.eps,NTU:r.NTU,hot:r.sH,cold:r.sC,R:r.R,air:r.air,wet:!!r.wet,Ph:r.Ph,Pc:r.Pc})}
  else{const cond=r.cond,sec=cond?c.cold:c.hot,T=c.ref.tsat;
    const refIn=cond?Math.max(c.ref.tgas,T):T,refOut=cond?T-Math.max(c.ref.sc,0):T+(c.ref.flooded?0:Math.max(c.ref.sh,0));
    const tpz=r.zones.find(z=>z.k==="tp");
    const refSide={h:tpz?tpz.h:NaN,dp:r.dpRef,Re:NaN,w:NaN,tp:true};
    Object.assign(o,{zones:r.zones,mref:r.mr,S:r.S,U:r.U,air:r.air,wet:r.zones.some(z=>z.wet),sumf:r.sumf});
    if(cond)Object.assign(o,{Thi:refIn,Tho:refOut,Tci:sec.tin,Tco:r.Tso,mh:r.mr,mc:sec.m,hot:refSide,cold:r.sS,Pc:r.Ps});
    else Object.assign(o,{Thi:sec.tin,Tho:r.Tso,Tci:refIn,Tco:refOut,mh:sec.m,mc:r.mr,hot:r.sS,cold:refSide,Ph:r.Ps});
    o.eps=NaN;o.NTU=NaN;
    {const S2=sat(c.ref.name,T+0.5,null),S1=sat(c.ref.name,T-0.5,null),dTs=r.dpRef/(S2.p-S1.p);o.dTsat=dTs;
      if(dTs>1.5)W.add("Kältemittel-Druckverlust "+f1(r.dpRef/1000)+" kPa entspricht ≈ "+f1(dTs)+" K Sättigungstemperatur-Abfall – mehr Kreisläufe / Pässe bzw. grösserer Querschnitt empfohlen.")}
    if(r.sumf<0.98||!isFinite(r.sumf))W.add("Zonenmodell konvergiert nicht sauber – Eingaben prüfen (Temperaturkreuzung?).")}
  // mittlere log. Temperaturdifferenz (wirksam) = Q/(U·A)
  o.lmtd=o.Q/(o.U*o.Aref);
  const d1=o.Thi-o.Tco,d2=o.Tho-o.Tci;o.lmtdCounter=Math.abs(d1-d2)<1e-6?d1:(d1-d2)/ln(d1/d2);
  if(o.Tho<o.Tci-0.05||o.Tco>o.Thi+0.05)W.add("Temperaturkreuzung: Austrittstemperatur über/unter der Eintrittstemperatur der anderen Seite.");
  validity(c,M,o,W);return o}
function validity(c,M,o,W){
  const chk=(s,lab)=>{if(!s||!isFinite(s.Re))return;
    if(c.type==="plate"&&s.Re<10)W.add(lab+": Re < 10 – Korrelation nach Martin ausserhalb des Gültigkeitsbereichs.");
    if(c.type==="shell"&&M.loc(lab==="Warme Seite"?"hot":"cold")==="shell"&&s.Re<2000)W.add(lab+" (Mantel): Re < 2000 – Kern-Verfahren hier ungenau (Bell-Delaware empfohlen).")};
  chk(o.hot,"Warme Seite");chk(o.cold,"Kalte Seite");
  if(c.type==="fin"){const s=M.loc("hot")==="air"?o.hot:o.cold;if(s&&(s.Re<300||s.Re>20000))W.add("Luftseite: Re(Dc) = "+Math.round(s.Re)+" ausserhalb 300 … 20 000 (Wang-Korrelation).");
    if(s&&s.w>4.5)W.add("Luftseite: Anströmgeschwindigkeit > 4,5 m/s – Tropfenmitriss / Geräusch prüfen.")}
  if(c.type==="mchx"){const s=M.loc("hot")==="air"?o.hot:o.cold;if(s&&(s.Re<100||s.Re>600))W.add("Luftseite: Re(Lp) = "+Math.round(s.Re)+" ausserhalb 100 … 600 (Kim & Bullard).")}
  if(c.type==="shell"&&c.geo.Nt>M.info.NtMax)W.add("Rohranzahl "+c.geo.Nt+" passt vermutlich nicht in den Mantel (max. ca. "+M.info.NtMax+" Rohre bei dieser Teilung).");
}

/* ---------- Vergleichsrechnung ---------- */
// target: {k:"Q"|"Tho"|"Tco", v}, knob: {k:"size"|"mh"|"mc"|"tsat"}, sizeKey: Geometriegrösse, int: ganzzahlig
function compare(c,target,knob,sizeKey,isInt){
  const get=r=>target.k==="Q"?r.Q:target.k==="Tho"?r.Tho:r.Tco;
  const setK=(cc,x)=>{const n=JSON.parse(JSON.stringify(cc));if(knob.k==="size"){if(n.type==="kman")n.kman.A=x;else n.geo[sizeKey]=x}
    else if(knob.k==="mh")n.hot.m=x;else if(knob.k==="mc")n.cold.m=x;else if(knob.k==="tsat")n.ref.tsat=x;return n};
  const x0=knob.k==="size"?(c.type==="kman"?c.kman.A:c.geo[sizeKey]):knob.k==="mh"?c.hot.m:knob.k==="mc"?c.cold.m:c.ref.tsat;
  const lin=knob.k==="tsat",xmin=knob.k==="size"&&isInt?(sizeKey==="np"?3:1):-Infinity;
  const F=x=>{try{const v=get(solve(setK(c,x)));return isFinite(v)?v-target.v:NaN}catch(e){return NaN}};
  const f0=F(x0);if(!isFinite(f0))return{ok:false,msg:"Ausgangsrechnung ungültig."};
  let a=NaN,b=NaN,fa,fb;
  if(Math.abs(f0)<1e-9*Math.max(1,Math.abs(target.v))){a=b=x0}
  else for(const dir of[1,-1]){let prev=x0,fp=f0;
    for(let i=1;i<=45;i++){const x=lin?x0+dir*i*2:x0*Math.pow(1.3,dir*i);if(!lin&&x<=0||x<xmin-1e-9)break;const fx=F(x);if(!isFinite(fx))continue;
      if((fx>0)!==(fp>0)){a=prev;fa=fp;b=x;fb=fx;break}prev=x;fp=fx}
    if(isFinite(a))break}
  if(!isFinite(a))return{ok:false,msg:"Wunschwert mit dieser Stellgrösse nicht erreichbar (geprüfter Bereich ausgeschöpft)."};
  if(a!==b)for(let i=0;i<60;i++){const m=lin||a*b<=0?(a+b)/2:Math.sign(a)*sq(a*b),fm=F(m);if(!isFinite(fm))break;if((fm>0)===(fa>0)){a=m;fa=fm}else{b=m;fb=fm}}
  let x=a===b?a:(Math.abs(fa)<Math.abs(fb)?a:b);
  if(isInt&&knob.k==="size"){const good=v=>{const d=F(v);return target.k==="Q"?d>=-1e-6*Math.abs(target.v):target.k==="Tho"?d<=1e-6:d>=-1e-6};
    const dn=Math.max(Math.floor(x+1e-9),xmin),up=Math.max(Math.ceil(x-1e-9),xmin);x=good(dn)?dn:up}
  const cc=setK(c,x);return{ok:true,x,c:cc,r:solve(cc)}}

return{solve,compare,props,sat,xAir,hAir,tDew,MODELS,E};
})();

/* =====================================================================
   RECHENKARTE (Oberfläche)
   ===================================================================== */
// Bauarten, Varianten, Voreinstellungen
const HX_TYPES=[["plate","Plattenwärmetauscher"],["shell","Rohrbündelwärmetauscher"],["dpipe","Doppelrohr (Rohr-in-Rohr)"],["fin","Lamellenrohr (Fin-and-Tube)"],["mchx","Microchannel (Flachrohr, Kiemenlamelle)"],["kman","Freie Eingabe (k-Wert und Fläche vorgegeben)"]];
const HX_FINPRE=[["7 mm (Teilung 21 × 12,7 mm)",7e-3,21e-3,12.7e-3,0.35e-3],["3/8″ 9,52 mm (25,4 × 22 mm)",9.52e-3,25.4e-3,22e-3,0.35e-3],["1/2″ 12,7 mm (31,75 × 27,5 mm)",12.7e-3,31.75e-3,27.5e-3,0.4e-3],["5/8″ 15,88 mm (38,1 × 33 mm)",15.88e-3,38.1e-3,33e-3,0.4e-3]];
const HX_WALL=[["Edelstahl 1.4404 / 1.4301",15],["Kupfer",390],["Aluminium",205],["Titan",21],["Stahl",50],["Kupfer-Nickel CuNi10Fe",45],["Hastelloy C-276",10],["Messing",110]];

function hxCard(){
  const c=el("div",{class:"card wide"});c.dataset.s="wärmetauscher wärmeübertrager platten rohrbündel doppelrohr lamellen fin tube microchannel verdampfer verflüssiger kondensator luftkühler rückkühler k-wert lmtd ntu";
  c.dataset.pdf="Wärmetauscher-Berechnung";
  c.append(el("div",{class:"grp"},"Wärmeübertragung"),el("h3",{},"Wärmetauscher – Auslegung & Vergleichsrechnung"));
  const fig=el("div",{class:"fig"});c.append(fig);
  const fx=el("div",{class:"f"});c.append(fx);
  const F={};let busy=false;
  // --- Bausteine ---
  const flowQ={"m³/h":["v",1/3600],"l/min":["v",1/60000],"l/s":["v",1e-3],"l/h":["v",1/3.6e6],"m³/s":["v",1],"kg/h":["m",1/3600],"kg/s":["m",1],"t/h":["m",1/3.6]};
  const row=(par,lab,...ch)=>{const r=el("div",{class:"row"},el("label",{},lab),...ch);par.append(r);return r};
  const num=(par,id,lab,q,u,v,opt={})=>{const inp=el("input",{inputmode:"decimal",value:String(v).replace(".",",")});let us,get,last;
    if(q==="flow"){us=el("select");Object.keys(flowQ).forEach(k=>us.append(el("option",{value:k},k)));us.value=u}else if(q)us=unitSel(q,u);else us=el("span",{},u||"");
    const r=row(par,lab,inp,us);last=us.value;
    const rd=()=>parse(inp.value.replace(/[’']/g,""));
    if(q==="flow")get=()=>({x:rd()*flowQ[us.value][1],t:flowQ[us.value][0]});else get=()=>q?toSI(q,us.value,rd()):rd();
    if(q&&q!=="flow")us.addEventListener("change",()=>{const x=rd();if(isFinite(x))inp.value=fmt(fromSI(q,us.value,toSI(q,last,x))).replace(/[’']/g,"");last=us.value;calc()});
    if(q==="flow")us.addEventListener("change",calc);
    inp.addEventListener("input",calc);
    const o={row:r,inp,us,get,set:x=>{inp.value=q&&q!=="flow"?fmt(fromSI(q,us.value,x)).replace(/[’']/g,""):fmt(x).replace(/[’']/g,"");last=us.value},opt};F[id]=o;return o};
  const sel=(par,id,lab,opts,def)=>{const s=el("select");opts.forEach(([v,t])=>s.append(el("option",{value:v},t)));if(def!=null)s.value=def;const r=row(par,lab,s);s.style.gridColumn="2 / 4";s.addEventListener("change",()=>{onSel(id);calc()});
    const o={row:r,s,get:()=>s.value,set:v=>{s.value=v}};F[id]=o;return o};
  const blk=(par,title)=>{const b=el("div");if(title)b.append(el("div",{class:"sub"},title));par.append(b);return b};
  const outv=(par,lab,q,u)=>{const v=el("div",{class:"v"},"–"),us=q?unitSel(q,u):el("span",{},u||"");const r=el("div",{class:"row out"},el("label",{},lab),v,us);par.append(r);let si=NaN;
    if(q)us.addEventListener("change",()=>o.set(si));const f4=x=>{if(!isFinite(x))return "–";const a=Math.abs(x);return x.toLocaleString("de-CH",a>=1000?{maximumFractionDigits:0}:{maximumSignificantDigits:4})};
    const o={row:r,us,set:x=>{si=x;v.textContent=q?f4(fromSI(q,us.value,x)):(typeof x==="string"?x:f4(x))}};return o};

  // --- Allgemein ---
  const gen=blk(c,"");
  sel(gen,"type","Bauart",HX_TYPES,"plate");
  sel(gen,"phase","Betriebsart",[["none","ohne Phasenwechsel"],["cond","Verflüssiger (KM kondensiert)"],["evap","Verdampfer (KM verdampft)"]],"none");
  sel(gen,"app","Anwendung",[["cool","Luftkühler (Kaltwasser/Sole)"],["heat","Lufterhitzer / Rückkühler"],["evap","Verdampfer (DX)"],["cond","Verflüssiger (luftgekühlt)"]],"cool");
  sel(gen,"kmode","k-Wert",[["geo","aus Geometrie berechnen"],["man","vorgeben (k-Wert und Fläche)"]],"geo");
  // --- Varianten je Bauart ---
  const vp=blk(c,"Variante");
  sel(vp,"pl_kind","Bauform",[["brazed","gelötet (Cu / Ni)"],["gasket","gedichtet (verschraubt, zerlegbar)"],["welded","verschweisst"]],"brazed");
  sel(vp,"pl_ang","Prägung (Winkel zur Strömung)",[["60","H – 60° (hohe Leistung, hoher Δp)"],["45","M – gemischt, ≈ 45°"],["30","L – 30° (geringer Δp)"]],"60");
  sel(vp,"pl_flow","Strömungsführung",[["counter","Gegenstrom"],["parallel","Gleichstrom"]],"counter");
  sel(vp,"sh_pass","Durchgänge (Mantel – Rohr)",[["1-1","1 – 1 (reiner Gegenstrom)"],["1-2","1 – 2"],["1-4","1 – 4"],["U","U-Rohr-Bündel (1 – 2)"],["2-4","2 – 4 (zwei Mäntel in Reihe)"]],"1-2");
  sel(vp,"sh_lay","Rohranordnung",[["tri","Dreieck 30°"],["rtri","Dreieck gedreht 60°"],["sq","Quadrat 90°"],["rsq","Quadrat gedreht 45°"]],"tri");
  sel(vp,"sh_hot","Warmes Medium",[["tube","in den Rohren"],["shell","im Mantel"]],"shell");
  sel(vp,"sh_evap","Verdampferart",[["dx","Trockenexpansion (KM in Rohren)"],["flooded","überflutet (KM im Mantel)"]],"dx");
  sel(vp,"dp_flow","Strömungsführung",[["counter","Gegenstrom"],["parallel","Gleichstrom"]],"counter");
  sel(vp,"dp_hot","Warmes Medium",[["in","im Innenrohr"],["ann","im Ringspalt"]],"in");
  sel(vp,"fi_fin","Lamellenform",[["plain","glatt"],["wavy","gewellt (Herringbone)"],["louver","Kiemen (Louver)"],["slit","geschlitzt (Slit)"]],"wavy");
  sel(vp,"fi_arr","Rohranordnung",[["stag","versetzt"],["inl","fluchtend"]],"stag");
  sel(vp,"fi_in","Rohr innen",[["smooth","glatt"],["micro","innenberippt (Microfin)"]],"smooth");
  sel(vp,"fi_mat","Werkstoff Rohr / Lamelle",[["cual","Kupfer / Aluminium"],["cucu","Kupfer / Kupfer"],["alal","Aluminium / Aluminium"],["ssal","Edelstahl / Aluminium"],["ssss","Edelstahl / Edelstahl"]],"cual");
  sel(vp,"mc_slab","Ausführung",[["1","1 Block (Slab)"],["2","2 Blöcke hintereinander"]],"1");
  sel(vp,"km_arr","Strömungsform",[["counter","Gegenstrom"],["parallel","Gleichstrom"],["cross","Kreuzstrom"],["shell","Rohrbündel 1 – 2"]],"counter");

  // --- Medien ---
  const mc=el("div",{class:"cols"}),mH=el("div"),mC=el("div");mc.append(mH,mC);c.append(mc);
  const fluids=[...Object.keys(HXDATA.liq).map(n=>[n,n]),["x","Eigene Stoffwerte …"]];
  const side=(par,p,title,tin,fl,flow,fu)=>{par.append(el("div",{class:"sub"},title));
    sel(par,p+"fl","Medium",fluids,fl);num(par,p+"tin","Eintrittstemperatur","temp","°C",tin);num(par,p+"flow","Durchfluss","flow",fu,flow);
    num(par,p+"rh","Rel. Luftfeuchte Eintritt","pct","%",50);
    num(par,p+"rho","Dichte ρ","dens","kg/m³",1000);num(par,p+"cp","Spez. Wärmekapazität cp","cp","kJ/(kg·K)",4.19);num(par,p+"mu","Dyn. Viskosität μ","visc","mPa·s",1);num(par,p+"k","Wärmeleitfähigkeit λ","","W/(m·K)",0.6);
    num(par,p+"rf","Verschmutzungsfaktor","","m²·K/W",0)};
  side(mH,"h_","Warme Seite",60,"Wasser",5,"m³/h");side(mC,"c_","Kalte Seite",20,"Ethylenglykol (MEG) 30 %",6,"m³/h");
  // Kältemittel
  const rb=blk(c,"Kältemittel");const rc=el("div",{class:"cols"}),r1=el("div"),r2=el("div");rc.append(r1,r2);rb.append(rc);
  sel(r1,"ref","Kältemittel",Object.keys(HXDATA.ref).map(n=>[n,n]),"R134a");
  num(r1,"tsat","Sättigungstemperatur","temp","°C",40);
  num(r2,"tgas","Heissgastemperatur Eintritt","temp","°C",70);num(r2,"sc","Unterkühlung","dtemp","K",3);
  num(r2,"tliq","Flüssigkeitstemp. vor Expansionsventil","temp","°C",35);num(r2,"sh","Überhitzung","dtemp","K",5);num(r1,"rrf","Verschmutzungsfaktor KM-Seite","","m²·K/W",0);

  // --- Geometrie ---
  const gb=blk(c,"Geometrie");const gc=el("div",{class:"cols"}),g1=el("div"),g2=el("div");gc.append(g1,g2);gb.append(gc);
  const G={plate:[],shell:[],dpipe:[],fin:[],mchx:[],kman:[]};
  const gn=(t,par,id,lab,q,u,v)=>{const o=num(par,id,lab,q,u,v);G[t].push(o.row);return o};
  const gs=(t,par,id,lab,opts,def)=>{const o=sel(par,id,lab,opts,def);G[t].push(o.row);return o};
  // Platte
  gn("plate",g1,"pl_np","Plattenanzahl","num","–",40);gn("plate",g1,"pl_B","Kanalbreite (Plattenbreite innen)","len","mm",115);gn("plate",g1,"pl_L","Plattenlänge (Anschlussmitte)","len","mm",470);
  gn("plate",g1,"pl_b","Prägetiefe / Plattenabstand","len","mm",2.2);gn("plate",g2,"pl_s","Plattendicke","len","mm",0.4);gn("plate",g2,"pl_phi","Flächenvergrösserung Φ","num","–",1.17);
  gn("plate",g2,"pl_dp","Anschlussdurchmesser","len","mm",30);gs("plate",g2,"pl_w","Plattenwerkstoff",HX_WALL.map(([n,k])=>[k,n]),"15");
  // Rohrbündel
  gn("shell",g1,"sh_Ds","Mantel-Innendurchmesser","len","mm",305);gn("shell",g1,"sh_Nt","Rohranzahl (Bohrungen je Mantel)","num","–",120);gn("shell",g1,"sh_do","Rohr-Aussendurchmesser","len","mm",19.05);
  gn("shell",g1,"sh_t","Rohr-Wandstärke","len","mm",1.24);gn("shell",g2,"sh_L","Rohrlänge (je Mantel)","len","m",2.5);gn("shell",g2,"sh_Pt","Rohrteilung","len","mm",23.8);
  gn("shell",g2,"sh_B","Umlenkblech-Abstand","len","mm",150);gs("shell",g2,"sh_w","Rohrwerkstoff",HX_WALL.map(([n,k])=>[k,n]),"15");
  // Doppelrohr
  gn("dpipe",g1,"dp_da","Innenrohr Aussen-Ø","len","mm",33.7);gn("dpipe",g1,"dp_t","Innenrohr Wandstärke","len","mm",2);gn("dpipe",g1,"dp_Di","Aussenrohr Innen-Ø","len","mm",54.5);
  gn("dpipe",g2,"dp_L","Länge je Element","len","m",3);gn("dpipe",g2,"dp_n","Anzahl Elemente in Reihe","num","–",4);gs("dpipe",g2,"dp_w","Rohrwerkstoff",HX_WALL.map(([n,k])=>[k,n]),"15");
  // Lamellenrohr
  const fpre=el("select");fpre.append(el("option",{value:""},"Rohr / Teilung wählen …"));HX_FINPRE.forEach((p,i)=>fpre.append(el("option",{value:i},p[0])));
  const fpr=el("div",{class:"preset"},el("span"),fpre);g1.append(fpr);G.fin.push(fpr);
  gn("fin",g1,"fi_do","Rohr-Aussendurchmesser","len","mm",9.52);gn("fin",g1,"fi_t","Rohr-Wandstärke","len","mm",0.35);gn("fin",g1,"fi_Pt","Rohrteilung quer (in Reihe)","len","mm",25.4);
  gn("fin",g1,"fi_Pl","Rohrteilung längs (Reihenabstand)","len","mm",22);gn("fin",g1,"fi_N","Rohrreihen (in Luftrichtung)","num","–",4);gn("fin",g1,"fi_Nt","Rohre je Reihe","num","–",24);
  gn("fin",g2,"fi_L","Berippte Länge","len","mm",1200);gn("fin",g2,"fi_Fp","Lamellenabstand (Teilung)","len","mm",2.5);gn("fin",g2,"fi_d","Lamellendicke","len","mm",0.12);
  gn("fin",g2,"fi_c","Anzahl Kreisläufe (Einspritzungen)","num","–",8);
  fpre.addEventListener("change",()=>{if(fpre.value==="")return;const p=HX_FINPRE[+fpre.value];F.fi_do.set(p[1]);F.fi_Pt.set(p[2]);F.fi_Pl.set(p[3]);F.fi_t.set(p[4]);fpre.value="";calc()});
  // Microchannel
  gn("mchx",g1,"mc_W","Blockbreite (Rohrlänge)","len","mm",800);gn("mchx",g1,"mc_H","Blockhöhe","len","mm",600);gn("mchx",g1,"mc_Td","Blocktiefe (Rohrtiefe)","len","mm",20);
  gn("mchx",g1,"mc_Dt","Rohrhöhe (Flachrohr aussen)","len","mm",1.8);gn("mchx",g1,"mc_Tp","Rohrteilung","len","mm",10);gn("mchx",g1,"mc_Np","Kanäle je Rohr","num","–",16);
  gn("mchx",g1,"mc_pw","Kanalbreite","len","mm",1.0);gn("mchx",g1,"mc_ph","Kanalhöhe","len","mm",1.2);
  gn("mchx",g2,"mc_Fp","Lamellenabstand (Teilung)","len","mm",1.3);gn("mchx",g2,"mc_d","Lamellendicke","len","mm",0.1);gn("mchx",g2,"mc_Lp","Kiemenabstand","len","mm",1.2);
  gn("mchx",g2,"mc_th","Kiemenwinkel","angle","°",27);gn("mchx",g2,"mc_pa","Anzahl Kältemittel-Pässe","num","–",3);
  // freie Eingabe
  gn("kman",g1,"km_U","k-Wert","uval","W/(m²·K)",2500);gn("kman",g1,"km_A","Fläche","area","m²",5);

  // --- Vergleichsrechnung ---
  const vb=blk(c,"Vergleichsrechnung (Wunschwert)");
  sel(vb,"cmp_t","Wunschwert",[["off","keine Vergleichsrechnung"],["Q","Leistung"],["Tho","Austrittstemperatur warme Seite"],["Tco","Austrittstemperatur kalte Seite"]],"off");
  num(vb,"cmp_q","Wunsch-Leistung","power","kW",100);num(vb,"cmp_T","Wunsch-Temperatur","temp","°C",30);
  sel(vb,"cmp_k","Angepasst wird",[["size","Baugrösse"],["mh","Durchfluss warme Seite"],["mc","Durchfluss kalte Seite"],["tsat","Sättigungstemperatur Kältemittel"]],"size");


  // --- Auswahlassistent (Fin-and-Tube / Microchannel) ---
  const ab=blk(c,"Auswahlassistent: passende Bauart & Ausführung vorschlagen");
  const abi=el("div",{class:"hint",style:"margin:0 0 6px"},"Nutzt die Luft-, Medium- bzw. Kältemitteldaten von oben und rechnet viele Ausführungen durch: Rohrdurchmesser, Lamellenform und -teilung, Rohrreihen bzw. Microchannel-Tiefe, Lamellenteilung und Blöcke.");ab.append(abi);
  sel(ab,"aw_typ","Vorauswahl Bauart",[["fin","Lamellenrohr (Fin-and-Tube)"],["mchx","Microchannel"],["both","beide vergleichen"]],"both");
  sel(ab,"aw_mode","Vorgabe",[["dim","Einbaumass vorgegeben (Breite × Höhe)"],["cap","Leistung / Temperatur vorgegeben – Grösse frei"]],"dim");
  sel(ab,"aw_tgt","Ziel",[["Q","Leistung"],["T","Luft-Austrittstemperatur"],["max","grösstmögliche Leistung im Einbaumass"]],"Q");
  num(ab,"aw_q","Gewünschte Leistung","power","kW",30);num(ab,"aw_T","Gewünschte Luft-Austrittstemperatur","temp","°C",14);
  const awc=el("div",{class:"cols"}),aw1=el("div"),aw2=el("div");awc.append(aw1,aw2);ab.append(awc);
  num(aw1,"aw_W","Breite (Stirnfläche, Rohrlänge)","len","mm",1200);num(aw1,"aw_H","Höhe (Stirnfläche)","len","mm",600);
  num(aw2,"aw_D","Max. Bautiefe (0 = beliebig)","len","mm",0);num(aw2,"aw_v","Max. Anströmgeschwindigkeit","vel","m/s",2.5);
  const awBtn=el("button",{type:"button",class:"pdfbtn",style:"margin:6px 0"},"Vorschläge berechnen");ab.append(awBtn);
  const awBox=el("div");ab.append(awBox);
  // --- Ergebnisse ---
  c.append(el("div",{class:"sub"},"Ergebnis"));const oc=el("div",{class:"cols"}),o1=el("div"),o2=el("div");oc.append(o1,o2);c.append(oc);
  const O={};const ov=(par,id,lab,q,u)=>{O[id]=outv(par,lab,q,u);return O[id]};
  ov(o1,"Q","Wärmeleistung Q","power","kW");ov(o1,"Tho","Austritt warme Seite","temp","°C");ov(o1,"Tco","Austritt kalte Seite","temp","°C");
  ov(o1,"mh","Massenstrom warme Seite","mflow","kg/h");ov(o1,"mc","Massenstrom kalte Seite","mflow","kg/h");
  ov(o1,"U","k-Wert (Wärmedurchgang)","uval","W/(m²·K)");ov(o1,"A","Bezugsfläche","area","m²");ov(o1,"lm","Mittlere Temperaturdifferenz","dtemp","K");ov(o1,"eps","Wirkungsgrad ε / NTU","","–");
  ov(o2,"ah","α warme Seite","uval","W/(m²·K)");ov(o2,"ac","α kalte Seite","uval","W/(m²·K)");ov(o2,"wh","Geschwindigkeit warme Seite","vel","m/s");ov(o2,"wc","Geschwindigkeit kalte Seite","vel","m/s");
  ov(o2,"Reh","Reynoldszahl warme / kalte Seite","","–");ov(o2,"dph","Druckverlust warme Seite","press","kPa");ov(o2,"dpc","Druckverlust kalte Seite","press","kPa");
  ov(o2,"etaf","Lamellenwirkungsgrad","pct","%");ov(o2,"wf","Luft: Anströmgeschwindigkeit","vel","m/s");ov(o2,"air","Luft Austritt: Feuchte / Kondensat","","");ov(o2,"shr","Sensible Leistung (SHR)","power","kW");
  ov(o2,"zon","Zonen (Flächenanteil)","","");
  [O.air,O.zon].forEach(o=>{const v=o.row.querySelector(".v");v.style.gridColumn="2 / 4";v.style.textAlign="right";o.row.lastChild.remove()});ov(o2,"mref","Kältemittel-Massenstrom","mflow","kg/h");ov(o2,"psat","Sättigungsdruck","pabs","bar (abs)");
  const gi=el("div",{class:"hint"});c.append(gi);
  const cmpBox=el("div");c.append(cmpBox);const varBox=el("div");c.append(varBox);const mfBox=el("div");c.append(mfBox);
  const note=el("div",{class:"note"}),hint=el("div",{class:"hint"});c.append(note,hint);

  // --- Sichtbarkeit ---
  const vis=(o,on)=>{(o.row||o).classList.toggle("hide",!on)};
  // Startwerte je Bauart / Betriebsart (werden nur beim Umschalten gesetzt)
  const sv_=(id,v,u)=>{const o=F[id];if(u&&o.us&&o.us.tagName==="SELECT")o.us.value=u;o.inp.value=String(v).replace(".",",")};
  const DEF={
    plate:{none:{h:["Wasser",60,5],c:["Ethylenglykol (MEG) 30 %",20,6]},cond:{c:["Wasser",30,5],r:["R134a",40]},evap:{h:["Wasser",12,4],r:["R134a",2]}},
    shell:{none:{h:["Wasser",80,36],c:["Wasser",20,36]},cond:{c:["Wasser",30,25],r:["R134a",40]},evap:{h:["Wasser",12,25],r:["R134a",3]}},
    dpipe:{none:{h:["Wasser",80,1.8],c:["Wasser",15,2.2]},cond:{c:["Wasser",30,1.5],r:["R134a",40]},evap:{h:["Wasser",12,1.5],r:["R134a",3]}},
    fin:{cool:{h:[null,27,8000],c:["Ethylenglykol (MEG) 30 %",6,2.5]},heat:{h:["Wasser",70,1.5],c:[null,10,8000]},evap:{h:[null,27,8000],r:["R410A",5]},cond:{c:[null,35,8000],r:["R410A",50]}},
    mchx:{cool:{h:[null,27,4000],c:["Ethylenglykol (MEG) 30 %",6,1.5]},heat:{h:["Wasser",70,1],c:[null,10,4000]},evap:{h:[null,27,4000],r:["R32",5]},cond:{c:[null,35,4000],r:["R32",50]}},
    kman:{none:{h:["Wasser",60,5],c:["Wasser",20,5]},cond:{c:["Wasser",30,5],r:["R134a",40]},evap:{h:["Wasser",12,5],r:["R134a",2]}}};
  function onSel(id){if(!["type","phase","app","sh_evap"].includes(id))return;const st=refresh(),d=DEF[st.t][st.coil?st.app:st.ph];if(!d)return;
    for(const [k,p] of [["h","h_"],["c","c_"]])if(d[k]){const[fl,t,q]=d[k];if(fl)F[p+"fl"].set(fl);sv_(p+"tin",t,"°C");sv_(p+"flow",q,"m³/h")}
    if(st.t==="fin")sv_("fi_c",{cool:8,heat:8,evap:24,cond:12}[st.app]);if(st.t==="mchx")sv_("mc_pa",{cool:1,heat:1,evap:1,cond:3}[st.app]);
    O.dph.us.value=st.coil&&st.airH?"Pa":"kPa";O.dpc.us.value=st.coil&&st.airC?"Pa":"kPa";
    if(d.r){F.ref.set(d.r[0]);sv_("tsat",d.r[1],"°C");if(st.ph==="cond"){sv_("tgas",d.r[1]+30,"°C");sv_("sc",3,"K")}else{sv_("tliq",40,"°C");sv_("sh",5,"K")}}
    refresh()}
  function refresh(){
    const t=F.type.get(),coil=t==="fin"||t==="mchx";let ph=F.phase.get();
    if(coil){const a=F.app.get();ph=a==="evap"?"evap":a==="cond"?"cond":"none";F.phase.set(ph)}
    if(t==="kman"&&F.kmode.get()!=="man")F.kmode.set("man");if(t!=="kman")F.kmode.set("geo");
    vis(F.phase,!coil);vis(F.app,coil);vis(F.kmode,false);
    const show={pl_kind:t==="plate",pl_ang:t==="plate",pl_flow:t==="plate",sh_pass:t==="shell",sh_lay:t==="shell",sh_hot:t==="shell"&&ph==="none",sh_evap:t==="shell"&&ph==="evap",
      dp_flow:t==="dpipe",dp_hot:t==="dpipe"&&ph==="none",fi_fin:t==="fin",fi_arr:t==="fin",fi_in:t==="fin",fi_mat:t==="fin",mc_slab:t==="mchx",km_arr:t==="kman"};
    for(const k in show)vis(F[k],show[k]);
    Object.keys(G).forEach(k=>G[k].forEach(r=>r.classList.toggle("hide",k!==t)));
    const app=F.app.get(),airH=coil&&(app==="cool"||app==="evap"),airC=coil&&(app==="heat"||app==="cond");
    if(coil){const AIR="Luft (1,013 bar)";if(airH){F.h_fl.set(AIR);if(F.c_fl.get()===AIR)F.c_fl.set("Ethylenglykol (MEG) 30 %")}if(airC){F.c_fl.set(AIR);if(F.h_fl.get()===AIR)F.h_fl.set("Ethylenglykol (MEG) 30 %")}}
    // Medienseiten
    ["h_","c_"].forEach(p=>{const isRef=(p==="h_"&&ph==="cond")||(p==="c_"&&ph==="evap"),isAir=coil&&((p==="h_"&&airH)||(p==="c_"&&airC)),cust=F[p+"fl"].get()==="x";
      [p+"fl",p+"tin",p+"flow",p+"rf"].forEach(k=>vis(F[k],!isRef));vis(F[p+"fl"],!isRef&&!isAir);
      [p+"rho",p+"cp",p+"mu",p+"k"].forEach(k=>vis(F[k],!isRef&&cust&&!isAir));
      vis(F[p+"rh"],!isRef&&isAir&&p==="h_");});
    mH.querySelector(".sub").textContent=ph==="cond"?"Warme Seite: Kältemittel (siehe unten)":coil&&airH?"Warme Seite: Luft":"Warme Seite";
    mC.querySelector(".sub").textContent=ph==="evap"?"Kalte Seite: Kältemittel (siehe unten)":coil&&airC?"Kalte Seite: Luft":"Kalte Seite";
    vis(rb,ph!=="none");vis(F.tgas,ph==="cond");vis(F.sc,ph==="cond");vis(F.tliq,ph==="evap");vis(F.sh,ph==="evap"&&!(t==="shell"&&F.sh_evap.get()==="flooded"));
    F.tsat.row.querySelector("label").textContent=ph==="cond"?"Kondensationstemperatur":"Verdampfungstemperatur";
    // Vergleich
    const ct=F.cmp_t.get();vis(F.cmp_q,ct==="Q");vis(F.cmp_T,ct==="Tho"||ct==="Tco");vis(F.cmp_k,ct!=="off");
    const ko=F.cmp_k.s.options;ko[0].text="Baugrösse ("+({plate:"Plattenanzahl",shell:"Rohrlänge",dpipe:"Anzahl Elemente",fin:"berippte Länge",mchx:"Blockbreite",kman:"Fläche"})[t]+")";
    ko[1].disabled=ph==="cond";ko[2].disabled=ph==="evap";ko[3].disabled=ph==="none";if(F.cmp_k.s.selectedOptions[0].disabled)F.cmp_k.set("size");
    // Ergebnisse
    vis(O.etaf,coil);vis(O.wf,coil);vis(O.air,coil&&airH&&(app==="cool"||app==="evap"));vis(O.shr,coil&&(app==="cool"||app==="evap"));
    vis(O.zon,ph!=="none");vis(O.mref,ph!=="none");vis(O.psat,ph!=="none");vis(O.eps,ph==="none");
    [O.ah,O.ac,O.wh,O.wc,O.Reh,O.dph,O.dpc].forEach(o=>vis(o,t!=="kman"));
    vis(ab,coil);const am=F.aw_mode.get(),at=F.aw_tgt.get();const mo=F.aw_tgt.s.options[2];mo.disabled=am==="cap";if(am==="cap"&&at==="max")F.aw_tgt.set("Q");
    vis(F.aw_q,F.aw_tgt.get()==="Q");vis(F.aw_T,F.aw_tgt.get()==="T");vis(F.aw_W,am==="dim");F.aw_H.row.querySelector("label").textContent=am==="dim"?"Höhe (Stirnfläche)":"Höhe vorgeben (0 = automatisch)";
    vis(F.aw_v,am==="cap");
    return{t,ph,coil,app,airH,airC}}

  // --- Eingaben einsammeln ---
  function flowSide(p,P){const f=F[p+"flow"].get();return f.t==="m"?f.x:f.x*P.rho}
  function fluid(p){const n=F[p+"fl"].get();if(n!=="x")return{name:n};return{custom:{rho:F[p+"rho"].get(),cp:F[p+"cp"].get(),mu:F[p+"mu"].get(),k:parse(F[p+"k"].inp.value)}}}
  function cfg(st){const t=st.t,ph=st.ph,v=id=>F[id].get(),n=id=>parse(F[id].inp.value.replace(/[’']/g,""));
    const c={type:t,phase:ph,hot:{},cold:{},Rf:{hot:0,cold:0},geo:{}};
    for(const [p,role] of [["h_","hot"],["c_","cold"]]){const isRef=(role==="hot"&&ph==="cond")||(role==="cold"&&ph==="evap");
      c.Rf[role]=isRef?n("rrf")||0:n(p+"rf")||0;if(isRef)continue;const fl=fluid(p),tin=v(p+"tin")-273.15,P=HX.props(fl,tin);c[role]={fl,tin,m:flowSide(p,P)}}
    if(ph!=="none")c.ref={name:v("ref"),tsat:v("tsat")-273.15,tgas:v("tgas")-273.15,sc:v("sc"),tliq:v("tliq")-273.15,sh:v("sh"),flooded:t==="shell"&&v("sh_evap")==="flooded"};
    if(st.coil){c.airSide=st.airH?"hot":"cold";c.air={rh:st.airH?v("h_rh"):0}}
    const L=x=>v(x);
    if(t==="plate")c.geo={np:Math.max(3,Math.round(n("pl_np"))),B:L("pl_B"),Lp:L("pl_L"),b:L("pl_b"),s:L("pl_s"),Phi:n("pl_phi"),Dp:L("pl_dp"),kw:+v("pl_w"),phi:+v("pl_ang")*Math.PI/180,flow:v("pl_flow")};
    if(t==="shell"){const ps=v("sh_pass"),lay=v("sh_lay");c.geo={Ds:L("sh_Ds"),Nt:Math.round(n("sh_Nt")),do:L("sh_do"),di:L("sh_do")-2*L("sh_t"),L:L("sh_L"),Pt:L("sh_Pt"),B:L("sh_B"),kw:+v("sh_w"),
      layout:lay==="sq"||lay==="rsq"?"sq":"tri",npp:ps==="1-1"?1:ps==="1-4"?4:2,shells:ps==="2-4"?2:1};
      c.hotLoc=ph==="cond"?"shell":ph==="evap"?(c.ref.flooded?"tube":"shell"):v("sh_hot")}
    if(t==="dpipe"){c.geo={da:L("dp_da"),di:L("dp_da")-2*L("dp_t"),Di:L("dp_Di"),Lel:L("dp_L"),n:Math.max(1,Math.round(n("dp_n"))),kw:+v("dp_w"),flow:v("dp_flow")};c.hotLoc=ph==="cond"?"in":ph==="evap"?"ann":v("dp_hot")}
    if(t==="fin"){const m=v("fi_mat"),kt={cual:390,cucu:390,alal:205,ssal:15,ssss:15}[m],kf={cual:205,cucu:390,alal:205,ssal:205,ssss:15}[m];
      c.geo={do:L("fi_do"),tw:L("fi_t"),Pt:L("fi_Pt"),Pl:L("fi_Pl"),N:Math.max(1,Math.round(n("fi_N"))),Nt:Math.max(1,Math.round(n("fi_Nt"))),L:L("fi_L"),Fp:L("fi_Fp"),delta:L("fi_d"),
        circ:Math.max(1,Math.round(n("fi_c"))),fin:v("fi_fin"),stag:v("fi_arr")==="stag",inner:v("fi_in"),kt,kf}}
    if(t==="mchx")c.geo={W:L("mc_W"),H:L("mc_H"),Td:L("mc_Td"),Dt:L("mc_Dt"),Tp:L("mc_Tp"),Np:Math.max(1,Math.round(n("mc_Np"))),pw:L("mc_pw"),ph:L("mc_ph"),Fp:L("mc_Fp"),delta:L("mc_d"),
      Lp:L("mc_Lp"),theta:v("mc_th")*180/Math.PI,passes:Math.max(1,Math.round(n("mc_pa"))),slabs:+v("mc_slab"),kt:160,kf:200};
    if(t==="kman")c.kman={U:v("km_U"),A:v("km_A"),arr:v("km_arr")};
    return c}
  const SIZE={plate:["np",1],shell:["L",0],dpipe:["n",1],fin:["L",0],mchx:["W",0],kman:["A",0]};

  // --- Formel / Skizze je Bauart ---
  const FXT={plate:"Q = k·A·ΔTm   ·   ε-NTU-Methode\nα Plattenkanal: Martin (VDI-Wärmeatlas), ξ und Nu aus Prägewinkel φ\nVerdampfung / Verflüssigung: Longo et al. (2015)",
    shell:"Q = k·A·F·ΔTm,log   ·   ε-NTU für 1–2n und 2–4 Durchgänge\nα Rohr: Gnielinski / VDI G1 · α Mantel: Kern-Verfahren\nKM in Rohren: Shah/Chato (Kondensation), Gungor-Winterton/Cooper (Verdampfung) · Mantel: Nusselt-Bündel / Cooper",
    dpipe:"Q = k·A·ΔTm,log   ·   1/k = 1/αa + da·ln(da/di)/(2λ) + da/(di·αi)\nα Rohr: Gnielinski / VDI G1 · α Ringspalt: Gnielinski mit Ringspaltfaktor (VDI G2)\nKältemittel im Innenrohr: Shah/Chato bzw. Gungor-Winterton/Cooper",
    fin:"Q = ε·Cmin·(t_w,ein − t_k,ein)   ·   Gegen-Kreuzstrom über n Rohrreihen\nLuftseite: Wang/Chi/Chang (2000) + Faktoren je Lamellenform · Rippenwirkungsgrad: Schmidt\nRohr: Gnielinski · KM: Shah/Chato bzw. Gungor-Winterton/Cooper · Entfeuchtung: Enthalpie-Methode",
    mchx:"Q = ε·Cmin·ΔTmax   ·   Kreuzstrom (2 Blöcke: Gegen-Kreuzstrom)\nLuftseite Kiemenlamelle: Kim & Bullard (2002) · Kanal: Gnielinski / laminar\nKM: Shah/Chato bzw. Gungor-Winterton/Cooper, Δp Müller-Steinhagen/Heck",
    kman:"Q = k·A·ΔTm   ·   ε-NTU-Methode für die gewählte Strömungsform"};
  const SKZ={plate:`<rect x="70" y="10" width="120" height="90" rx="6" class="lf"/>`+[0,1,2,3,4,5,6].map(i=>`<line x1="${85+i*15}" y1="16" x2="${85+i*15}" y2="94" class="ln" style="stroke-width:1"/>`).join("")+`<circle cx="84" cy="22" r="5" class="lf"/><circle cx="176" cy="22" r="5" class="lf"/><circle cx="84" cy="88" r="5" class="lf"/><circle cx="176" cy="88" r="5" class="lf"/>`,
    shell:`<rect x="40" y="30" width="180" height="50" rx="20" class="lf"/><line x1="40" y1="45" x2="220" y2="45" class="ln" style="stroke-width:1"/><line x1="40" y1="55" x2="220" y2="55" class="ln" style="stroke-width:1"/><line x1="40" y1="65" x2="220" y2="65" class="ln" style="stroke-width:1"/>`+[90,130,170].map((x,i)=>`<line x1="${x}" y1="${i%2?52:30}" x2="${x}" y2="${i%2?80:58}" class="ac"/>`).join("")+`<rect x="70" y="18" width="10" height="12" class="lf"/><rect x="180" y="80" width="10" height="12" class="lf"/>`,
    dpipe:`<rect x="30" y="35" width="200" height="40" class="lf"/><rect x="20" y="48" width="220" height="14" class="lf" style="fill:var(--bg)"/>`,
    fin:Array.from({length:14},(_,i)=>`<line x1="${60+i*10}" y1="15" x2="${60+i*10}" y2="95" class="ln" style="stroke-width:1"/>`).join("")+[25,45,65,85].map(y=>`<circle cx="0" cy="0" r="0"/><line x1="50" y1="${y}" x2="200" y2="${y}" class="ac"/>`).join(""),
    mchx:Array.from({length:6},(_,i)=>`<rect x="40" y="${15+i*16}" width="180" height="4" class="lf"/>`).join("")+Array.from({length:5},(_,i)=>`<path d="M40,${19+i*16} `+Array.from({length:18},(_,k)=>`L${45+k*10},${(k%2?19:31)+i*16}`).join(" ")+`" class="ln" style="stroke-width:.8"/>`).join("")+`<rect x="28" y="10" width="12" height="90" class="lf"/><rect x="220" y="10" width="12" height="90" class="lf"/>`,
    kman:`<rect x="80" y="25" width="100" height="60" rx="8" class="lf"/><text x="130" y="60">k · A</text>`};
  const arrows=`<line x1="5" y1="40" x2="40" y2="40" class="ac" marker-end="url(#ah)"/><line x1="220" y1="72" x2="255" y2="72" class="ac" marker-end="url(#ah)"/>`;

  // --- Rechnen ---
  const fmtS=(x,d=1)=>isFinite(x)?x.toLocaleString("de-CH",{maximumFractionDigits:d,minimumFractionDigits:d}):"–";
  let lastT="";
  function calc(){if(busy)return;busy=true;try{run()}catch(e){note.textContent="Berechnung nicht möglich: Eingaben prüfen."+(e&&e.message?" ("+e.message+")":"");console.error(e)}busy=false}
  function run(){note.textContent="";hint.textContent="";gi.textContent="";cmpBox.innerHTML="";varBox.innerHTML="";mfBox.innerHTML="";
    const st=refresh();if(st.t!==lastT){lastT=st.t;fig.innerHTML=sv(SKZ[st.t]+arrows);}
    fx.textContent=FXT[st.t];c.dataset.pdf="Wärmetauscher – "+HX_TYPES.find(x=>x[0]===st.t)[1].replace(/ \(.*$/,"");
    const c0=cfg(st);
    for(const r of ["hot","cold"])if(c0[r]&&c0[r].m!==undefined&&!(c0[r].m>0))throw new Error("Durchfluss "+(r==="hot"?"warm":"kalt")+" fehlt");
    const r=HX.solve(c0);show(r,st,c0);
    const ct=F.cmp_t.get();
    if(ct!=="off"){const tv=ct==="Q"?F.cmp_q.get():F.cmp_T.get()-273.15,kk=F.cmp_k.get(),[sk,isInt]=SIZE[st.t];
      const cr=HX.compare(c0,{k:ct,v:tv},{k:kk},sk,isInt);compareTable(r,cr,c0,st,ct,tv,kk,sk);
      if(st.t==="fin"||st.t==="mchx")variantTable(r,c0,st,ct,tv)}
    makerBox(st);
    const W=[...r.warn];if(W.length)note.textContent=W.join(" ");
    hint.textContent="Richtwerte – Auslegung durch Fachperson bzw. Hersteller-Software bestätigen. Stoffwerte: CoolProp (Glykole in Massen-%). Kältemittel: Zonenmodell mit Sättigungsstoffwerten, Temperaturgleit zeotroper Gemische nicht berücksichtigt."}
  function show(r,st,c0){
    O.ah.row.querySelector("label").textContent="α warme Seite";O.ac.row.querySelector("label").textContent="α kalte Seite";
    O.Q.set(r.Q);O.Tho.set(r.Tho+273.15);O.Tco.set(r.Tco+273.15);O.mh.set(r.mh);O.mc.set(r.mc);O.U.set(r.U);O.A.set(r.Aref);O.lm.set(r.lmtd);
    O.eps.set(isFinite(r.eps)?fmtS(r.eps*100,1)+" %  /  "+fmtS(r.NTU,2):"–");
    const hs=r.hot||{},cs=r.cold||{};O.ah.set(hs.h*(hs.eta||1));O.ac.set(cs.h*(cs.eta||1));O.wh.set(hs.w);O.wc.set(cs.w);
    O.Reh.set((isFinite(hs.Re)?Math.round(hs.Re).toLocaleString("de-CH"):"–")+"  /  "+(isFinite(cs.Re)?Math.round(cs.Re).toLocaleString("de-CH"):"–"));
    O.dph.set(hs.dp);O.dpc.set(cs.dp);
    if(st.coil){const a=st.airH?hs:cs;O.etaf.set(a.etaFin);O.wf.set(a.w);
      if(r.air){O.air.set("φ = "+fmtS(r.air.phiOut*100,0)+" %, x = "+fmtS(r.air.xout*1000,1)+" g/kg, Kondensat "+fmtS(r.air.cond*3600,1)+" kg/h");O.shr.set(r.air.Qs);
        O.shr.row.querySelector("label").textContent="Sensible Leistung (SHR "+fmtS(Math.min(r.air.Qs/r.Q,1),2)+")"}}
    if(r.mode==="tp"){O.zon.set(r.zones.map(z=>z.n+" "+fmtS(z.f*100,0)+" %").join(" · "));O.mref.set(st.ph==="cond"?r.mh:r.mc);O.psat.set(r.S.p);
      O[st.ph==="cond"?"ah":"ac"].row.querySelector("label").textContent="α "+(st.ph==="cond"?"warme":"kalte")+" Seite (Zweiphasen)";
      O[st.ph==="cond"?"wh":"wc"].set(NaN);
      O[st.ph==="cond"?"ah":"ac"].set(r[st.ph==="cond"?"hot":"cold"].h)}
    // Geometrie-Info
    const i=r.info,mm=x=>fmtS(x*1000,1),m2=x=>fmtS(x,2);
    const gl={plate:()=>`Kanäle warm/kalt: ${i.nh}/${i.nc} · hydraul. Durchmesser ${mm(i.Dh)} mm · Fläche ${m2(r.Aref)} m²`,
      shell:()=>`Rohrinnen-Ø ${mm(c0.geo.di)} mm · Rohre je Durchgang ${fmtS(c0.geo.Nt/c0.geo.npp,0)} · Umlenkbleche ≈ ${Math.max(1,Math.round(c0.geo.L/c0.geo.B)-1)} je Mantel · max. Rohranzahl im Mantel ≈ ${i.NtMax}`,
      dpipe:()=>`Gesamtlänge ${fmtS(c0.geo.Lel*c0.geo.n,1)} m · Ringspalt-Dh ${mm(c0.geo.Di-c0.geo.da)} mm`,
      fin:()=>`Stirnfläche ${m2(i.Afr)} m² · Aussenfläche ${m2(i.Ao)} m² (Lamellen ${fmtS(i.Af/i.Ao*100,0)} %) · Innenfläche ${m2(i.Ai)} m² · Verhältnis ${fmtS(i.Ao/i.Ai,1)} · Rohrlänge je Kreislauf ${fmtS(i.Lc,1)} m`,
      mchx:()=>`${i.Nt} Flachrohre · Lamellenhöhe ${mm(i.Fl)} mm · Kanal-Dh ${mm(i.dh)} mm · Stirnfläche ${m2(i.Afr)} m² · Aussenfläche ${m2(i.Ao)} m² · Innenfläche ${m2(i.Ai)} m²`,kman:()=>""}[st.t]();
    gi.textContent=(gl?gl+" · ":"")+"k-Wert bezogen auf "+r.refLabel+(r.wet?" · Kühler arbeitet nass (Entfeuchtung).":"");
  }
  function compareTable(r,cr,c0,st,ct,tv,kk,sk){
    cmpBox.append(el("div",{class:"sub"},"Vergleich: Auslegung ↔ Wunschwert"));
    if(!cr.ok){cmpBox.append(el("div",{class:"note"},cr.msg));return}
    const v=cr.r,unitK={plate:["Plattenanzahl","–",1],shell:["Rohrlänge je Mantel","m",1],dpipe:["Anzahl Elemente","–",1],fin:["Berippte Länge","mm",1000],mchx:["Blockbreite","mm",1000],kman:["Fläche","m²",1]}[st.t];
    const x0=kk==="size"?(st.t==="kman"?c0.kman.A:c0.geo[sk]):kk==="mh"?c0.hot.m:kk==="mc"?c0.cold.m:c0.ref.tsat,x1=cr.x;
    const rows=[];const pct=(a,b)=>isFinite(a)&&isFinite(b)&&Math.abs(a)>1e-12?((b-a)/Math.abs(a)*100):NaN;
    const add=(lab,a,b,u,d=1,rel=true)=>rows.push([lab,fmtS(a,d),fmtS(b,d),u,rel?(isFinite(pct(a,b))?(pct(a,b)>=0?"+":"")+fmtS(pct(a,b),1)+" %":"–"):(isFinite(b-a)?(b-a>=0?"+":"")+fmtS(b-a,d)+" "+u:"–")]);
    if(kk==="size")add(unitK[0],x0*unitK[2],x1*unitK[2],unitK[1],unitK[1]==="–"?0:(unitK[1]==="m"?2:0));
    if(kk==="mh")add("Durchfluss warme Seite",x0*3600,x1*3600,"kg/h",0);
    if(kk==="mc")add("Durchfluss kalte Seite",x0*3600,x1*3600,"kg/h",0);
    if(kk==="tsat")add("Sättigungstemperatur",x0,x1,"°C",1,false);
    add("Wärmeleistung",r.Q/1000,v.Q/1000,"kW",2);
    add("Austritt warme Seite",r.Tho,v.Tho,"°C",1,false);add("Austritt kalte Seite",r.Tco,v.Tco,"°C",1,false);
    add("Bezugsfläche",r.Aref,v.Aref,"m²",2);add("k-Wert",r.U,v.U,"W/(m²·K)",0);add("Mittl. Temperaturdifferenz",r.lmtd,v.lmtd,"K",2,false);
    if(st.t!=="kman"){const pa=k=>st.coil&&(k==="h"?st.airH:st.airC);
      add("Druckverlust warme Seite",r.hot.dp/(pa("h")?1:1000),v.hot.dp/(pa("h")?1:1000),pa("h")?"Pa":"kPa",pa("h")?0:2);
      add("Druckverlust kalte Seite",r.cold.dp/(pa("c")?1:1000),v.cold.dp/(pa("c")?1:1000),pa("c")?"Pa":"kPa",pa("c")?0:2)}
    if(r.mode==="tp")add("Kältemittel-Massenstrom",(st.ph==="cond"?r.mh:r.mc)*3600,(st.ph==="cond"?v.mh:v.mc)*3600,"kg/h",1);
    const t=el("table",{class:"cmp"});t.append(el("tr",{},...["Grösse","Auslegung","Vergleich","Einheit","Abweichung"].map(h=>el("th",{},h))));
    rows.forEach(rw=>t.append(el("tr",{},...rw.map(x=>el("td",{},x)))));cmpBox.append(el("div",{class:"tw"},t));
    // Reserve
    let msg="";const tl=ct==="Q"?"Wunsch-Leistung "+fmtS(tv/1000,2)+" kW":"Wunsch-Austrittstemperatur "+fmtS(tv,1)+" °C";
    if(kk==="size"){const res=(x0/x1-1)*100;msg=tl+": "+(res>=0?"Auslegung hat "+fmtS(res,1)+" % Reserve (Baugrösse).":"Auslegung ist "+fmtS(-res,1)+" % zu klein – Baugrösse erhöhen.")}
    else msg=tl+" wird mit der angepassten Grösse erreicht.";
    cmpBox.append(el("div",{class:"pdfmsg hint",style:"color:var(--fg);font-weight:600"},msg))}


  /* --- Variantenvergleich Fin-and-Tube / Microchannel ---
     Gleiche Stirnhöhe, gleiche Reihen-/Blockzahl, gleiche Betriebsdaten; je Ausführung wird die
     Baugrösse (berippte Länge bzw. Blockbreite) gesucht, die den Wunschwert erreicht. */
  function variantTable(r,c0,st,ct,tv){
    const g0=c0.geo,vars=[],FN={plain:"glatt",wavy:"gewellt",louver:"Kiemen",slit:"geschlitzt"};
    if(st.t==="fin"){const H=g0.Nt*g0.Pt;
      HX_FINPRE.forEach(p=>Object.keys(FN).forEach(fin=>{const g=Object.assign({},g0,{do:p[1],Pt:p[2],Pl:p[3],tw:p[4],fin,Nt:Math.max(1,Math.round(H/p[2]))});
        vars.push({name:fmtS(p[1]*1000,p[1]*1000%1?2:0)+" mm "+FN[fin]+" ("+g.Nt+" R.)",geo:g,cur:Math.abs(p[1]-g0.do)<1e-5&&fin===g0.fin})}))}
    else{[16,20,25.4,32].forEach(td=>[1.0,1.2,1.4,1.7].forEach(fp=>{const Td=td/1000,Fp=fp/1000,g=Object.assign({},g0,{Td,Fp,Np:Math.max(1,Math.round(g0.Np*Td/g0.Td))});
        vars.push({name:"T "+fmtS(td,td%1?1:0)+" · "+g.Np+" K. · Lam. "+fmtS(fp,1)+" mm",geo:g,cur:Math.abs(Td-g0.Td)<1e-5&&Math.abs(Fp-g0.Fp)<1e-5})}))}
    if(!vars.some(v=>v.cur))vars.push({name:"aktuelle Eingabe",geo:Object.assign({},g0),cur:true});
    const sk=st.t==="fin"?"L":"W",out=[];
    for(const v of vars){const cc=JSON.parse(JSON.stringify(c0));cc.geo=v.geo;let r1,cr;
      try{r1=HX.solve(cc);cr=HX.compare(cc,{k:ct,v:tv},{k:"size"},sk,false)}catch(e){continue}
      const air=st.airH?"hot":"cold",inn=st.airH?"cold":"hot",rr=cr&&cr.ok?cr.r:null;
      out.push({v,Q:r1.Q,x:cr&&cr.ok?cr.x:NaN,dpa:rr?rr[air].dp:NaN,dpi:rr?(rr.mode==="tp"?rr.dpRef:rr[inn].dp):NaN})}
    if(!out.length)return;
    const ok=out.filter(o=>isFinite(o.x)).sort((a,b)=>a.x-b.x);const best=ok[0];
    out.sort((a,b)=>(isFinite(a.x)?a.x:1e9)-(isFinite(b.x)?b.x:1e9));
    const lab=st.t==="fin"?"Berippte Länge":"Blockbreite";
    varBox.append(el("div",{class:"sub"},"Variantenvergleich: Rohrdimension & Lamellen ("+(st.t==="fin"?"gleiche Stirnhöhe und Reihenzahl":"gleiche Blockhöhe und Blockzahl")+")"));
    const t=el("table",{class:"cmp"});t.append(el("tr",{},...["Ausführung","Q akt. Grösse","Erf. "+(st.t==="fin"?"Länge":"Breite"),"Δp Luft","Δp innen"].map(h=>el("th",{},h))));
    out.forEach(o=>t.append(el("tr",{style:o===best?"font-weight:700":o.v.cur?"font-style:italic":""},
      el("td",{},(o===best?"▶ ":"")+o.v.name+(o.v.cur?" ◀":"")),el("td",{},fmtS(o.Q/1000,2)+" kW"),
      el("td",{},isFinite(o.x)?fmtS(o.x*1000,0)+" mm":"nicht erreichbar"),el("td",{},isFinite(o.dpa)?fmtS(o.dpa,0)+" Pa":"–"),el("td",{},isFinite(o.dpi)?fmtS(o.dpi/1000,1)+" kPa":"–"))));
    varBox.append(el("div",{class:"tw"},t));varBox.append(el("div",{class:"hint",style:"margin-top:2px"},st.t==="fin"?"R. = Rohre je Reihe (gleiche Stirnhöhe) · ▶ kompakteste Ausführung · ◀ aktuelle Eingabe":"T = Rohrtiefe mm · K. = Kanäle je Rohr · Lam. = Lamellenteilung · ▶ kompakteste Ausführung · ◀ aktuelle Eingabe"));
    const cur=out.find(o=>o.v.cur);
    varBox.append(el("div",{class:"pdfmsg hint",style:"color:var(--fg);font-weight:600"},best?"Kompakteste Ausführung: "+best.v.name+" mit "+fmtS(best.x*1000,0)+" mm "+lab.toLowerCase()+
      (cur&&isFinite(cur.x)?" (aktuelle Ausführung: "+fmtS(cur.x*1000,0)+" mm, "+(cur.x>best.x?fmtS((cur.x/best.x-1)*100,0)+" % grösser)":"bereits die kompakteste)"):"")+
      ". Druckverluste bei der jeweils erforderlichen Grösse; höherer Luft-Δp bedeutet mehr Ventilatorleistung.":"Wunschwert mit keiner Ausführung erreichbar – Betriebsdaten oder Blockhöhe anpassen."))}

  /* --- Hersteller-Vorschläge (Auswahl, ohne Gewähr) --- */
  const MAKERS={plate:"Kelvion (DE), Danfoss (DK), SWEP (SE), Alfa Laval (SE) – gelötete und gedichtete Plattenwärmetauscher, jeweils mit eigener Auslegungssoftware.",
    shell:"Kelvion (DE), Alfa Laval (SE), Danfoss (DK) – Rohrbündel und Rohrbündelverdampfer/-verflüssiger, auch als Sonderanfertigung.",
    dpipe:"Kelvion (DE), Danfoss (DK, Rohr-in-Rohr HE) – Doppelrohr- bzw. Rohr-in-Rohr-Wärmetauscher.",
    fin:"Güntner (DE), Kelvion (DE), LU-VE (IT) – Lamellenrohr-Wärmetauscher (Coils) nach Mass, verschiedene Rohrdurchmesser, Teilungen und Lamellenformen; Kelvion z. B. Cu-Rohr 7,94 mm (25 × 21,65 mm) und 12 mm (32 × 27,71 mm).",
    mchx:"Kaltra (DE/AT) – Microchannel-Wärmetauscher nach Mass mit vielen Flachrohr- und Lamellengeometrien, eigene Auslegungssoftware MCHEselect; Danfoss (DK) – MCHE-Standardreihe und kundenspezifisch.",
    kman:"Hersteller je nach Bauart (siehe andere Bauarten)."};
  function makerBox(st){mfBox.append(el("div",{class:"sub"},"Hersteller (Auswahl, ohne Gewähr)"),el("div",{class:"pdfmsg hint",style:"color:var(--fg)"},MAKERS[st.t]+" Für Angebot und verbindliche Auslegung die Geometrie und Betriebsdaten aus diesem Bericht an den Hersteller senden."))}


  /* ===== Auswahlassistent: Kandidaten erzeugen, rechnen, bewerten ===== */
  const FNn={plain:"glatt",wavy:"gewellt",louver:"Kiemen",slit:"geschlitzt"};
  function awCircuits(c,g,Qest){ // Kreisläufe so wählen, dass die Innenseite sinnvoll durchströmt wird
    const di=g.do-2*g.tw,Ai=PI_*di*di/4,tubes=g.Nt*g.N;
    if(c.ref){const S=HX.sat(c.ref.name,c.ref.tsat),mr=Math.max(Qest,100)/(S.r*(c.phase==="cond"?1.15:0.9)),G=c.phase==="cond"?250:150;return Math.max(1,Math.min(tubes,Math.round(mr/(G*Ai))))}
    const side=c.airSide==="hot"?"cold":"hot",P=HX.props(c[side].fl,c[side].tin);return Math.max(1,Math.min(tubes,Math.round(c[side].m/(P.rho*1.0*Ai))))}
  const PI_=Math.PI;
  function awCandidates(c0,st,typ,W,H,Dmax){const out=[],evap=st.app==="evap";
    if(typ!=="mchx"){const kt=c0.type==="fin"?c0.geo.kt:390,kf=c0.type==="fin"?c0.geo.kf:205,dl=c0.type==="fin"?c0.geo.delta:0.12e-3;
      for(const p of HX_FINPRE)for(const fin of Object.keys(FNn))for(const fp of(evap?[2.1,2.5,3.2]:[1.8,2.5,3.2]))for(let N=1;N<=8;N++){
        const D=N*p[3];if(Dmax>0&&D>Dmax+1e-9)break;const Nt=Math.max(1,Math.floor(H/p[2]+1e-9));
        out.push({t:"fin",D,face:[W,Nt*p[2]],name:fmtS(p[1]*1000,p[1]*1000%1?2:0)+" mm "+FNn[fin]+" · "+N+" R. · Lam. "+fmtS(fp,1),
          geo:{do:p[1],tw:p[4],Pt:p[2],Pl:p[3],N,Nt,L:W,Fp:fp/1000,delta:dl,circ:1,fin,stag:true,inner:c0.ref?"micro":"smooth",kt,kf}})}}
    if(typ!=="fin"){const g0=c0.type==="mchx"?c0.geo:{Dt:1.8e-3,Tp:10e-3,Np:16,Td:0.02,pw:1e-3,ph:1.2e-3,delta:0.1e-3,Lp:1.2e-3,theta:27,kt:160,kf:200};
      const pas=c0.type==="mchx"?c0.geo.passes:({cool:1,heat:1,evap:1,cond:3}[st.app]||1);
      for(const td of[16,20,25.4,32])for(const fp of[1.0,1.2,1.4,1.7])for(const sl of[1,2]){const D=sl*td/1000+(sl-1)*0.01;if(Dmax>0&&D>Dmax+1e-9)continue;const Np=Math.max(1,Math.round(g0.Np/g0.Td*td/1000));
        out.push({t:"mchx",D,face:[W,H],name:"T "+fmtS(td,td%1?1:0)+" · "+Np+" K. · Lam. "+fmtS(fp,1)+" · "+sl+" Block",
          geo:Object.assign({},g0,{W,H,Td:td/1000,Fp:fp/1000,Np,slabs:sl,passes:pas})})}}
    return out}
  function awEval(c0,st,k){const cc=JSON.parse(JSON.stringify(c0));cc.type=k.t;cc.geo=JSON.parse(JSON.stringify(k.geo));
    try{if(k.t==="fin"){cc.geo.circ=awCircuits(cc,cc.geo,k.Qest||20000);let r=HX.solve(cc);cc.geo.circ=awCircuits(cc,cc.geo,r.Q);r=HX.solve(cc);return{cc,r}}
      return{cc,r:HX.solve(cc)}}catch(e){return null}}
  async function awRun(){
    const st=refresh();if(!st.coil)return;awBox.innerHTML="";awBtn.disabled=true;awBtn.textContent="Rechne Varianten …";await new Promise(r=>setTimeout(r,30));
    try{
      const c0=cfg(st),typ=F.aw_typ.get(),mode=F.aw_mode.get(),tgt=F.aw_tgt.get(),air=st.airH?"hot":"cold",inn=st.airH?"cold":"hot";
      c0.type=st.t;const Pa=HX.props(c0[air].fl,c0[air].tin),Vair=c0[air].m/Pa.rho;
      const tq=F.aw_q.get(),tT=F.aw_T.get()-273.15,Dmax=F.aw_D.get()||0;
      let W=F.aw_W.get(),H=F.aw_H.get();
      if(mode==="cap"){const A=Vair/Math.max(F.aw_v.get(),0.5);if(!(H>0))H=Math.sqrt(A/1.6);H=Math.ceil(H*100)/100;W=Math.ceil(A/H*100)/100}
      const meets=r=>tgt==="Q"?r.Q>=tq*0.995:tgt==="T"?(st.airH?r.Tho<=tT+0.05:r.Tco>=tT-0.05):true;
      const cand=awCandidates(c0,st,typ,W,H,Dmax);const res=[];let n=0;
      for(const k of cand){if(tgt==="Q")k.Qest=tq;const e=awEval(c0,st,k);if(e&&isFinite(e.r.Q)){const r=e.r;
          res.push({k,cc:e.cc,Q:r.Q,ta:st.airH?r.Tho:r.Tco,dpa:r[air].dp,dpi:isFinite(r.dpRef)?r.dpRef:r[inn].dp,vf:r[air].w,ok:meets(r),res:tgt==="Q"?r.Q/tq-1:tgt==="T"?Math.abs((st.airH?c0[air].tin-r.Tho:r.Tco-c0[air].tin)/(st.airH?c0[air].tin-tT:tT-c0[air].tin))-1:0})}
        if(++n%60===0)await new Promise(r=>setTimeout(r,0))}
      // je Grundausführung nur die kleinste Bautiefe, die das Ziel erreicht
      let list;
      if(tgt==="max")list=res.sort((a,b)=>b.Q-a.Q||a.dpa-b.dpa);
      else{const seen=new Set(),bk=x=>x.res<=0.25?0:1;list=res.filter(x=>x.ok).sort((a,b)=>bk(a)-bk(b)||(bk(a)?a.res-b.res:a.k.D-b.k.D||a.dpa-b.dpa)).filter(x=>{const key=x.k.t+x.k.name.replace(/ · \d+ R\./,"");if(seen.has(key))return false;seen.add(key);return true})}
      let top=list.slice(0,12);
      if(typ==="both"){const f=list.filter(x=>x.k.t==="fin").slice(0,6),m=list.filter(x=>x.k.t==="mchx").slice(0,6);
        const bk=x=>x.res<=0.25?0:1;top=[...f,...m].sort(tgt==="max"?(a,b)=>b.Q-a.Q||a.dpa-b.dpa:(a,b)=>bk(a)-bk(b)||(bk(a)?a.res-b.res:a.k.D-b.k.D||a.dpa-b.dpa))}
      const vf=Vair/(W*H);
      awBox.append(el("div",{class:"sub"},"Vorschläge"+(mode==="cap"?" (Stirnfläche "+fmtS(W*1000,0)+" × "+fmtS(H*1000,0)+" mm aus Luftmenge "+fmtS(Vair*3600,0)+" m³/h)":" für "+fmtS(W*1000,0)+" × "+fmtS(H*1000,0)+" mm")));
      if(!top.length){awBox.append(el("div",{class:"note"},"Kein Vorschlag erreicht das Ziel"+(Dmax>0?" innerhalb der max. Bautiefe":"")+". Stirnfläche vergrössern, mehr Bautiefe zulassen oder Betriebsdaten (Luftmenge, Temperaturen) prüfen. Grösste gefundene Leistung: "+fmtS(Math.max(...res.map(x=>x.Q))/1000,2)+" kW."));return}
      const t=el("table",{class:"cmp"});t.append(el("tr",{},...["Nr.","Bauart / Ausführung","Stirnmass B×H mm","Tiefe mm","Leistung kW",tgt==="max"?"Luft aus °C":"Reserve","Δp Luft Pa","Δp innen kPa"].map(h=>el("th",{},h))));
      top.forEach((x,i)=>t.append(el("tr",{style:i===0?"font-weight:700":""},el("td",{},(i===0?"▶ ":"")+(i+1)),el("td",{style:"text-align:left"},(x.k.t==="fin"?"FT ":"MC ")+x.k.name),
        el("td",{},fmtS(x.k.face[0]*1000,0)+"×"+fmtS(x.k.face[1]*1000,0)),el("td",{},fmtS(x.k.D*1000,0)),el("td",{},fmtS(x.Q/1000,2)),el("td",{},tgt==="max"?fmtS(x.ta,1):"+"+fmtS(Math.max(0,x.res)*100,0)+" %"),el("td",{},fmtS(x.dpa,0)),el("td",{},isFinite(x.dpi)?fmtS(x.dpi/1000,1):"–"))));
      awBox.append(el("div",{class:"tw"},t),el("div",{class:"hint",style:"margin-top:2px"},"FT = Fin-and-Tube (Rohr-Ø · Lamellenform · Rohrreihen · Lamellenteilung mm) · MC = Microchannel (Rohrtiefe mm · Kanäle · Lamellenteilung mm · Blöcke). "+(tgt==="max"?"Sortiert nach grösster Leistung.":"Zuerst Ausführungen mit 0–25 % Reserve (kleinste Bautiefe, dann geringster Luft-Druckverlust), danach stärker überdimensionierte.")));
      const best=top[0],bf=list.find(x=>x.k.t==="fin"),bm=list.find(x=>x.k.t==="mchx");
      let msg="Empfehlung: "+(best.k.t==="fin"?"Fin-and-Tube ":"Microchannel ")+best.k.name+", "+fmtS(best.k.D*1000,0)+" mm Bautiefe, "+fmtS(best.Q/1000,2)+" kW, Δp Luft "+fmtS(best.dpa,0)+" Pa.";
      if(typ==="both"&&bf&&bm)msg+=" Bester "+(best.k.t==="fin"?"Microchannel: "+bm.k.name+" ("+fmtS(bm.k.D*1000,0)+" mm, "+fmtS(bm.dpa,0)+" Pa)":"Fin-and-Tube: "+bf.k.name+" ("+fmtS(bf.k.D*1000,0)+" mm, "+fmtS(bf.dpa,0)+" Pa)")+".";
      if(st.airH&&(st.app==="cool"||st.app==="evap")&&vf>2.8)msg+=" Achtung: Anströmgeschwindigkeit "+fmtS(vf,1)+" m/s – bei nassem Kühler über ca. 2,5–2,8 m/s Tropfenabriss, Tropfenabscheider vorsehen.";
      if(typ==="both"&&bm&&bm.res>0.5)msg+=" Microchannel ist hier selbst in der kleinsten Ausführung "+fmtS(bm.res*100,0)+" % überdimensioniert – mit kleinerer Blockbreite bzw. -höhe möglich (Vorgabe Einbaumass nutzen).";
      if(st.ph!=="none")msg+=" Kältemittelseite: Anzahl Kreisläufe bzw. Pässe vom Hersteller optimieren lassen.";
      awBox.append(el("div",{class:"pdfmsg hint",style:"color:var(--fg);font-weight:600"},msg));
      // übernehmen
      const ps=el("select");top.forEach((x,i)=>ps.append(el("option",{value:i},(i+1)+": "+(x.k.t==="fin"?"FT ":"MC ")+x.k.name)));
      const pb=el("button",{type:"button",class:"updbtn",style:"color:var(--fg);border-color:var(--line)"},"In die Berechnung übernehmen");
      pb.addEventListener("click",()=>awApply(top[+ps.value]));
      ps.style.flex="1";const pw=el("div",{style:"display:flex;gap:8px;align-items:center;flex-wrap:wrap;margin-top:8px"},el("span",{style:"font-size:13.5px"},"Vorschlag übernehmen:"),ps,pb);awBox.append(pw);
    }finally{awBtn.disabled=false;awBtn.textContent="Vorschläge berechnen"}}
  function awApply(x){const g=x.cc.geo,setL=(id,v)=>F[id].set(v);
    F.type.set(x.k.t);refresh();
    if(x.k.t==="fin"){setL("fi_do",g.do);setL("fi_t",g.tw);setL("fi_Pt",g.Pt);setL("fi_Pl",g.Pl);F.fi_N.set(g.N);F.fi_Nt.set(g.Nt);setL("fi_L",g.L);setL("fi_Fp",g.Fp);F.fi_c.set(g.circ);F.fi_fin.set(g.fin);F.fi_arr.set("stag");F.fi_in.set(g.inner)}
    else{setL("mc_W",g.W);setL("mc_H",g.H);setL("mc_Td",g.Td);F.mc_Np.set(g.Np);setL("mc_Fp",g.Fp);F.mc_slab.set(String(g.slabs));F.mc_pa.set(g.passes)}
    lastT="";calc();toast("Vorschlag übernommen – Berechnung aktualisiert");fig.scrollIntoView({behavior:"smooth",block:"start"})}
  awBtn.addEventListener("click",awRun);

  onSel("type");setTimeout(calc,0);
  return c;
}
