restart; with(codegen); outstring2 := "realRho.cc"; outstring := "realmaple.cc";
F11 := proc (rho) options operator, arrow; (b-c)*rho+c*rho*ln(rho)+d end proc;
F00 := proc (rho) options operator, arrow; (f-e)*rho+e*rho*ln(rho)+g end proc;
c := 28.9445110522; b := -13.9748379249; d := 20.9488105766; e := 9.50205321859; f := 6.09644617041; g := 0.649044609107e-1;
G1 := D(F11); G0 := D(F00); p1 := proc (rho) options operator, arrow; -F11(rho)+rho*G1(rho) end proc; p0 := proc (rho) options operator, arrow; -F00(rho)+rho*G0(rho) end proc;
sols := fsolve([p1(rho) = p0(x), G1(rho) = G0(x)], {rho = 0 .. 4, x = 0 .. 4}); sol1 := solve(sols[1]); sol2 := solve(sols[2]); mwg := proc (x) options operator, arrow; G1(sol1)*x-p1(sol1) end proc; F1 := F11-mwg; F0 := F00-mwg;
nn := proc (phi) options operator, arrow; 6.*phi^5+(-1)*15.*phi^4+10.*phi^3 end proc; W := proc (phi) options operator, arrow; 2*A*(phi^4-2*phi^3+phi^2) end proc; Pressure := proc (rho, phi) options operator, arrow; nn(phi)*p1(rho)+(1-nn(phi))*p0(rho) end proc; Potential := proc (rho, phi) options operator, arrow; nn(phi)*G1(rho)+(1-nn(phi))*G0(rho)-G0(sol0) end proc;
F := proc (rho, phi) options operator, arrow; W(phi)/delta+nn(phi)*F1(rho)+(1-nn(phi))*F0(rho) end proc;
solrho1 := proc (x) options operator, arrow; exp((Const-nn(x)*(b-f)-f)/(nn(x)*(c-e)+e)) end proc; Const := G1(sol1); solproc1 := makeproc(solrho1(x), x); evalRho := optimize(solproc1);
Fbulk := proc (rho, phi) options operator, arrow; beta*(nn(phi)*F1(rho)+(1-nn(phi))*F0(rho)) end proc;

Fproc := makeproc(F(rho, phi), [rho, phi]);
helmholtz := optimize(Fproc); C(helmholtz, filename = outstring, ansi);
P2 := proc (rho, phi) options operator, arrow; Pressure(rho, phi)-W(phi)/delta end proc;
# 
fp3 := proc (rho, phi) options operator, arrow; diff(helmholtz(rho, phi), phi, phi, phi) end proc; fp1 := proc (rho, phi) options operator, arrow; diff(helmholtz(rho, phi), phi) end proc; Sourceproc := makeproc(fp1(x, y), [x, y]); Sproc := makeproc(fp1(rho, mid)+(1/24)*fp3(rho, mid)*(phi-old), [rho, phi, old, mid]); SSproc := makeproc(Sproc(rho, phi, old, .5*(phi+old)), [rho, phi, old]); reactionSource := optimize(SSproc); C(reactionSource, filename = outstring, ansi);
dphiS := proc (rho, phi, old) options operator, arrow; diff(SSproc(rho, phi, old), phi) end proc; dphiSproc := makeproc(dphiS(rho, phi, old, mid), [rho, phi, old, mid]); dphireactionSource := optimize(dphiSproc); C(dphireactionSource, filename = outstring, ansi);
drhoS := proc (rho, phi, old) options operator, arrow; diff(SSproc(rho, phi, old), rho) end proc; drhoSproc := makeproc(drhoS(rho, phi, old, mid), [rho, phi, old, mid]); drhoreactionSource := optimize(drhoSproc); C(drhoreactionSource, filename = outstring, ansi);
fr3 := proc (rho, phi) options operator, arrow; diff(helmholtz(rho, phi), rho, rho, rho) end proc; fr1 := proc (rho, phi) options operator, arrow; diff(helmholtz(rho, phi), rho) end proc; CProc := makeproc(fr1(mid, phi)+(1/24)*fr3(mid, phi)*(rho-old), [rho, phi, old, mid]); CCProc := makeproc(CProc(rho, phi, old, .5*(rho+old)), [rho, phi, old]);
chemicalPotential := optimize(CCProc); C(chemicalPotential, filename = outstring, ansi);
drhoPotential := proc (rho, phi, old, mid) options operator, arrow; diff(CCProc(rho, phi, old), rho) end proc; drhomuproc := makeproc(drhoPotential(rho, phi, old), [rho, phi, old]); drhochemicalPotential := optimize(drhomuproc); C(drhochemicalPotential, filename = outstring, ansi);
dphiPotential := proc (rho, phi, old, mid) options operator, arrow; diff(CCProc(rho, phi, old), phi) end proc; dphimuproc := makeproc(dphiPotential(rho, phi, old), [rho, phi, old]); dphichemicalPotential := optimize(dphimuproc); C(dphichemicalPotential, filename = outstring, ansi);
Pproc := makeproc(P2(rho, phi), [rho, phi]); pressure := optimize(Pproc); C(pressure, filename = outstring, ansi);
wavespeed := proc (rho, phi) options operator, arrow; sqrt(diff(Pressure(rho, phi), rho)) end proc; wproc := makeproc(simplify(wavespeed(rho, phi)), [rho, phi]); a := optimize(wproc); C(a, filename = outstring, ansi);
A := 0.1e-1; alpha := 0; inttanh := int(sqrt(2*A*W(x)), x = 0. .. 1.);
sols := fsolve([p1(rho)-inttanh/(.25) = p0(x), G1(rho) = G0(x)], {rho = 0 .. 1, x = 0 .. 1});
sols := fsolve([p1(rho)+inttanh/(.25) = p0(x), G1(rho) = G0(x)], {rho = 0 .. 1, x = 0 .. 1}); bubblesol := solve(sols[1]); bubblesol2 := solve(sols[2]); ml := proc (a) options operator, arrow; bubblesol end proc; mwpliq := makeproc(ml(x), x);
mv := proc (a) options operator, arrow; bubblesol2 end proc; mwpvap := makeproc(mv(x), x); C(mwpliq, filename = outstring2, ansi); C(mwpvap, filename = outstring2, ansi); Const := G1(bubblesol);
solrho := proc (x) options operator, arrow; exp((Const-nn(x)*(b-f)-f)/(nn(x)*(c-e)+e)) end proc; solproc := makeproc(solrho(x), x); evalRho := optimize(solproc); C(evalRho, filename = outstring2, ansi);


NULL;

# 
