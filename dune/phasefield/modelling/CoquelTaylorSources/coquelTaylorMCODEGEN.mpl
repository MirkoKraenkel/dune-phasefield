restart; with(codegen); outstring2 := "coquelTaylorMRho.cc"; outstring := "coquelTaylorMmaple.cc"; outstring3 := "coquelTaylorMSource.cc";beta := 0; with(plots); Digits := 20;
F0 := proc (rho) options operator, arrow; (b-c)*rho+c*rho*ln(rho)+d end proc;
F1 := proc (rho) options operator, arrow; (f-e)*rho+e*rho*ln(rho)+g end proc;

c := 1.5; b := log(2.); d := 0; e := 1; f := 0; g := .5;
G1 := D(F1); G0 := D(F0); p1 := proc (rho) options operator, arrow; -F1(rho)+rho*G1(rho) end proc; p0 := proc (rho) options operator, arrow; -F0(rho)+rho*G0(rho) end proc;
sols := fsolve([p1(rho) = p0(x), G1(rho) = G0(x)], {rho = 0 .. 4, x = 0 .. 4}); sol1 := solve(sols[1]); sol2 := solve(sols[2]); mwg := proc (rho) options operator, arrow; G1(sol1)*rho-p1(sol1) end proc;
nn := proc (phi) options operator, arrow; 6.*phi^5+(-1)*15.*phi^4+10.*phi^3 end proc; W := proc (phi) options operator, arrow; 2*A*(phi^4-2*phi^3+phi^2) end proc; dW := D(W); F := proc (rho, phi) options operator, arrow; .1*(W(phi)/delta+nn(phi)*F1(rho)+(1-nn(phi))*F0(rho)) end proc; Pressure := proc (rho, phi) options operator, arrow; nn(phi)*p1(rho)+(1-nn(phi))*p0(rho) end proc; P2 := proc (rho, phi) options operator, arrow; Pressure(rho, phi)-W(phi)/delta end proc;
solrho1 := proc (x) options operator, arrow; exp((Const-nn(x)*(f-b)-b)/(nn(x)*(e-c)+c)) end proc; Const := G1(sol1); solproc1 := makeproc(solrho1(x), x); evalRho1 := optimize(solproc1);


Fproc := makeproc(F(rho, phi), [rho, phi]);
helmholtz := optimize(Fproc); C(helmholtz, filename = outstring, ansi);
# 
fp3 := proc (rho, phi) options operator, arrow; diff(helmholtz(rho, phi), phi, phi, phi) end proc; fp1 := proc (rho, phi) options operator, arrow; diff(helmholtz(rho, phi), phi) end proc; Sourceproc := makeproc(fp1(x, y), [x, y]); Sproc := makeproc(fp1(rho, mid)+(1/24)*fp3(rho, mid)*(phi-old)^2, [rho, phi, old, mid]); SSproc := makeproc(Sproc(rho, phi, old, .5*(phi+old)), [rho, phi, old]); reactionSource := optimize(SSproc); C(reactionSource, filename = outstring, ansi);
dphiS := proc (rho, phi, old) options operator, arrow; diff(SSproc(rho, phi, old), phi) end proc; dphiSproc := makeproc(dphiS(rho, phi, old, mid), [rho, phi, old, mid]); dphireactionSource := optimize(dphiSproc); C(dphireactionSource, filename = outstring, ansi);
drhoS := proc (rho, phi, old) options operator, arrow; diff(SSproc(rho, phi, old), rho) end proc; drhoSproc := makeproc(drhoS(rho, phi, old, mid), [rho, phi, old, mid]); drhoreactionSource := optimize(drhoSproc); C(drhoreactionSource, filename = outstring, ansi);
fr3 := proc (rho, phi) options operator, arrow; diff(helmholtz(rho, phi), rho, rho, rho) end proc;
fr1 := proc (rho, phi) options operator, arrow; diff(helmholtz(rho, phi), rho) end proc;
CProc := makeproc(fr1(mid, phi)+(1/24)*fr3(mid, phi)*(rho-old)^2, [rho, phi, old, mid]); CCProc := makeproc(CProc(rho, phi, old, .5*(rho+old)), [rho, phi, old]);
chemicalPotential := optimize(CCProc); C(chemicalPotential, filename = outstring, ansi); chemproc := makeproc(fr1(rho, phi), [rho, phi]);
drhoPotential := proc (rho, phi, old, mid) options operator, arrow; diff(CCProc(rho, phi, old), rho) end proc; drhomuproc := makeproc(drhoPotential(rho, phi, old), [rho, phi, old]); drhochemicalPotential := optimize(drhomuproc); C(drhochemicalPotential, filename = outstring, ansi);
dphiPotential := proc (rho, phi, old, mid) options operator, arrow; diff(CCProc(rho, phi, old), phi) end proc; dphimuproc := makeproc(dphiPotential(rho, phi, old), [rho, phi, old]); dphichemicalPotential := optimize(dphimuproc); C(dphichemicalPotential, filename = outstring, ansi);
Pproc := makeproc(P2(rho, phi), [rho, phi]); pressure := optimize(Pproc); C(pressure, filename = outstring, ansi);
wavespeed := proc (rho, phi) options operator, arrow; sqrt(diff(Pressure(rho, phi), rho)) end proc; wproc := makeproc(simplify(wavespeed(rho, phi)), [rho, phi]); a := optimize(wproc); C(a, filename = outstring, ansi);
psol := proc (t, x, y) options operator, arrow; .5*cos(2*Pi*t)*cos(2*Pi*x)+.5 end proc; exactphi := makeproc(psol(t, x, y), [t, x, y]); exactphi := optimize(exactphi); C(exactphi, filename = outstring2, ansi);
v1sol := proc (t, x, y) options operator, arrow; cos(2*Pi*t)*sin(2*Pi*x) end proc; exactv1 := makeproc(v1sol(t, x, y), [t, x, y]); exactv1 := optimize(exactv1); C(exactv1, filename = outstring2, ansi);
v2sol := proc (t, x, y) options operator, arrow; cos(2*Pi*t)*sin(2*Pi*x) end proc; exactv2 := makeproc(v2sol(t, x, y), [t, x, y]); exactv2 := optimize(exactv2); C(exactv2, filename = outstring2, ansi);

rsol := proc (t, x, y) options operator, arrow; solproc1(psol(t, x, y)) end proc; .5*cos(2*Pi*t)*cos(2*Pi*x)+1.5; exactrho := makeproc(rsol(t, x, y), [t, x, y]); exactrho := optimize(exactrho); C(exactrho, filename = outstring2, ansi);
dtr := diff(rsol(t, x, y), t);
divrv := diff(rsol(t, x, y)*v1sol(t, x, y), x)+diff(rsol(t, x, y)*v2sol(t, x, y), y);
dxphi := diff(psol(t, x, y), x); dyphi := diff(psol(t, x, y), y);
exactsigma1 := makeproc(dxphi, [t, x, y]); exactsigma2 := makeproc(dyphi, [t, x, y]); C(exactsigma2, filename = outstring2, ansi); C(exactsigma1, filename = outstring2, ansi);
dtphi := diff(psol(t, x, y), t); lapphi := A*(diff(psol(t, x, y), x, x)+diff(psol(t, x, y), y, y));
phitrans := v1sol(t, x, y)*dxphi+v2sol(t, x, y)*dyphi;
tau := Sourceproc(rsol(t, x, y), psol(t, x, y)); exacttau := makeproc(tau-delta*lapphi, [t, x, y]); extacttau := optimize(exacttau); C(exacttau, filename = outstring2, ansi); exactmu := makeproc(chemproc(rsol(t, x, y), psol(t, x, y))+.5*(v1sol(t, x, y)^2+v2sol(t, x, y)^2), [t, x, y]); extactmu := optimize(exactmu); C(exactmu, filename = outstring2, ansi);
dxmu := diff(exactmu(t, x, y), x); dymu := diff(exactmu(t, x, y), y);
dtmom1 := rsol(t, x, y)*(diff(v1sol(t, x, y), t)); dtmom2 := rsol(t, x, y)*(diff(v2sol(t, x, y), t));
momtrans1 := diff(rsol(t, x, y)*v1sol(t, x, y)*v1sol(t, x, y), x)+diff(rsol(t, x, y)*v1sol(t, x, y)*v2sol(t, x, y), y);
momtrans2 := diff(rsol(t, x, y)*v1sol(t, x, y)*v2sol(t, x, y), x)+diff(rsol(t, x, y)*v2sol(t, x, y)*v2sol(t, x, y), y);
lapv1 := visc1*(diff(v1sol(t, x, y), x, x)+diff(v1sol(t, x, y), y, y));
lapv2 := visc1*(diff(v2sol(t, x, y), x, x)+diff(v2sol(t, x, y), y, y));
RHSrho := dtr+divrv;
RHSv1 := dtmom1+momtrans1-divrv*v1sol(t, x, y)+rsol(t, x, y)*dxmu-.5*(diff(v1sol(t, x, y)^2+v2sol(t, x, y)^2, x))*rsol(t, x, y)-dxphi*(tau-delta*lapphi)-lapv1;
RHSv2 := dtmom2+momtrans2-divrv*v2sol(t, x, y)+rsol(t, x, y)*dymu-.5*(diff(v1sol(t, x, y)^2+v2sol(t, x, y)^2, y))*rsol(t, x, y)-dyphi*(tau-delta*lapphi)-lapv2;
RHSphi := dtphi+phitrans+(tau-delta*lapphi)/rsol(t, x, y);
rhsRho := makeproc(RHSrho, [t, x, y]); rhsRho := optimize(rhsRho); C(rhsRho, filename = outstring3, filename = outstring3, ansi);
rhsV1 := makeproc(RHSv1, [t, x, y]); rhsV1 := optimize(rhsV1); C(rhsV1, filename = outstring3, ansi);
rhsV2 := makeproc(RHSv2, [t, x, y]); rhsV2 := optimize(rhsV2); C(rhsV2, filename = outstring3, ansi);
rhsPhi := makeproc(RHSphi, [t, x, y]); rhsPhi := optimize(rhsPhi); C(rhsPhi, filename = outstring3, ansi);


# 
# 
# 
