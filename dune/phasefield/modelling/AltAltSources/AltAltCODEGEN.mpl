restart; with(codegen); outstring2 := "balanced.cc"; outstring := "maple.cc";
c0 := 0; c1 := 0; a0 := .1; r0 := 2; r1 := .1292; a1 := 10;
F0 := proc (rho) options operator, arrow; a0*rho*(log(rho/r0)-1)+a0*r0+c0 end proc;
g0 := D(F0); p0 := proc (rho) options operator, arrow; -F0(rho)+rho*g0(rho) end proc; G0 := makeproc(g0(x), x);
F1 := proc (rho) options operator, arrow; beta*(a1*rho*(log(rho/r1)-1)+a1*r1+c1) end proc;

g1 := D(F1); p1 := proc (rho) options operator, arrow; -F1(rho)+rho*g1(rho) end proc; G1 := makeproc(g1(x), x);
nn := proc (phi) options operator, arrow; 6.*phi^5+(-1)*15.*phi^4+10.*phi^3 end proc;
W := proc (phi) options operator, arrow; 2*A*(phi^4-2*phi^3+phi^2) end proc;
dW := D(W); F := proc (rho, phi) options operator, arrow; rho*W(phi)/delta+nn(phi)*F1(rho)+(1-nn(phi))*F0(rho) end proc;
Pressure := proc (rho, phi) options operator, arrow; beta*(nn(phi)*p1(rho)+(1-nn(phi))*p0(rho)) end proc; Potential := proc (rho, phi) options operator, arrow; nn(phi)*G1(rho)+(1-nn(phi))*G0(rho) end proc; P2 := proc (rho, phi) options operator, arrow; Pressure(rho, phi) end proc;
solrho := proc (x) options operator, arrow; (r0-r1)*x+r1 end proc; Const := G1(sol1); solproc := makeproc(solrho(x), x); evalRho := optimize(solproc); C(evalRho, filename = outstring, ansi);
Fproc := makeproc(F(rho, phi), [rho, phi]); helmholtz := optimize(Fproc); C(helmholtz, filename = outstring, ansi);
fp3 := proc (rho, phi) options operator, arrow; diff(helmholtz(rho, phi), phi, phi, phi) end proc; fp1 := proc (rho, phi) options operator, arrow; diff(helmholtz(rho, phi), phi) end proc; Sproc := makeproc(fp1(rho, mid)+(1/24)*fp3(rho, mid)*(phi-old), [rho, phi, old, mid]); SSproc := makeproc(Sproc(rho, phi, old, .5*(phi+old)), [rho, phi, old]);
reactionSource := optimize(SSproc); C(reactionSource, filename = outstring, ansi);
dphiS := proc (rho, phi, old) options operator, arrow; diff(SSproc(rho, phi, old), phi) end proc; dphiSproc := makeproc(dphiS(rho, phi, old), [rho, phi, old]);
dphireactionSource := optimize(dphiSproc); C(dphireactionSource, filename = outstring, ansi);
fr3 := proc (rho, phi) options operator, arrow; diff(helmholtz(rho, phi), rho, rho, rho) end proc; fr1 := proc (rho, phi) options operator, arrow; diff(helmholtz(rho, phi), rho) end proc; CProc := makeproc(fr1(mid, phi)+(1/24)*fr3(mid, phi)*(rho-old), [rho, phi, old, mid]); CCProc := makeproc(CProc(rho, phi, old, .5*(rho+old)), [rho, phi, old]);
chemicalPotential := optimize(CCProc); C(chemicalPotential, filename = outstring, ansi);
drhoPotential := proc (rho, phi, old) options operator, arrow; diff(CCProc(rho, phi, old), rho) end proc; drhomuproc := makeproc(drhoPotential(rho, phi, old), [rho, phi, old]); drhochemicalPotential := optimize(drhomuproc); C(drhochemicalPotential, filename = outstring, ansi);
dphiPotential := proc (rho, phi, old) options operator, arrow; diff(CCProc(rho, phi, old), phi) end proc; dphimuproc := makeproc(dphiPotential(rho, phi, old), [rho, phi, old]); dphichemicalPotential := optimize(dphimuproc); C(dphichemicalPotential, filename = outstring, ansi);
Pproc := makeproc(P2(rho, phi), [rho, phi]); pressure := optimize(Pproc); C(pressure, filename = outstring, ansi);
wavespeed := proc (rho, phi) options operator, arrow; sqrt(diff(Pressure(rho, phi), rho)) end proc; wproc := makeproc(simplify(wavespeed(rho, phi)), [rho, phi]); a := optimize(wproc); C(a, filename = outstring, ansi);




