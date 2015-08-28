restart; with(codegen); outstring2 := "balanced.cc"; outstring := "maple.cc";
c0 := 0; c1 := 0; a0 := 5; r0 := 2.; r1 := .1292; a1 := 10;
F0 := proc (rho) options operator, arrow; a0*rho*(log(rho/r0)-1)+a0*r0+c0 end proc;
g0 := D(F0); p0 := proc (rho) options operator, arrow; -F0(rho)+rho*g0(rho) end proc; G0 := makeproc(g0(x), x);
F1 := proc (rho) options operator, arrow; a1*rho*(log(rho/r1)-1)+a1*r1+c1 end proc;
g1 := D(F1); p1 := proc (rho) options operator, arrow; -F1(rho)+rho*g1(rho) end proc; G1 := makeproc(g1(x), x);
nn := proc (phi) options operator, arrow; 6.*phi^5+(-1)*15.*phi^4+10.*phi^3 end proc;
W := proc (phi) options operator, arrow; 2*A*(phi^4+(2*alpha-2)*phi^3+(-3*alpha+1)*phi^2+alpha) end proc;
dW := D(W); F := proc (rho, phi) options operator, arrow; rho*W(phi)/delta+beta*(nn(phi)*F1(rho)+(1-nn(phi))*F0(rho)) end proc;
Pressure := proc (rho, phi) options operator, arrow; beta*(nn(phi)*p1(rho)+(1-nn(phi))*p0(rho)) end proc; Potential := proc (rho, phi) options operator, arrow; W(phi)/delta+beta*(nn(phi)*G1(rho)+(1-nn(phi))*G0(rho)) end proc; P2 := proc (rho, phi) options operator, arrow; Pressure(rho, phi) end proc;
solrho := proc (x) options operator, arrow; (r0-r1)*x+r1 end proc; Const := G1(sol1); solproc := makeproc(solrho(x), x); evalRho := optimize(solproc); C(evalRho, filename = outstring, ansi);
Fproc := makeproc(F(rho, phi), [rho, phi]); helmholtz := optimize(Fproc); C(helmholtz, filename = outstring, ansi);
S := proc (rho, phi) options operator, arrow; diff(F(rho, phi), phi) end proc; Sproc := makeproc(S(rho, phi), [rho, phi]); reactionSource := optimize(Sproc); C(reactionSource, filename = outstring, ansi);
dphiS := proc (rho, phi) options operator, arrow; diff(S(rho, phi), phi) end proc; dphiSproc := makeproc(dphiS(rho, phi), [rho, phi]); dphireactionSource := optimize(dphiSproc); C(dphireactionSource, filename = outstring, ansi);
CProc := makeproc(simplify(Potential(rho, phi)), [rho, phi]); chemicalPotential := optimize(CProc); C(chemicalPotential, filename = outstring, ansi);
drhoPotential := proc (rho, phi) options operator, arrow; diff(Potential(rho, phi), rho) end proc; drhomuproc := makeproc(simplify(drhoPotential(rho, phi)), [rho, phi]); drhochemicalPotential := optimize(drhomuproc); C(drhochemicalPotential, filename = outstring, ansi);
dphiPotential := proc (rho, phi) options operator, arrow; diff(Potential(rho, phi), phi) end proc; dphimuproc := makeproc(simplify(dphiPotential(rho, phi)), [rho, phi]); dphichemicalPotential := optimize(dphimuproc); C(dphichemicalPotential, filename = outstring, ansi);
Pproc := makeproc(P2(rho, phi), [rho, phi]); pressure := optimize(Pproc); C(pressure, filename = outstring, ansi);
wavespeed := proc (rho, phi) options operator, arrow; sqrt(diff(Pressure(rho, phi), rho)) end proc; wproc := makeproc(simplify(wavespeed(rho, phi)), [rho, phi]); a := optimize(wproc); C(a, filename = outstring, ansi);
solphi := proc (x) options operator, arrow; .5*tanh(x/delta)+.5 end proc; rs := proc (x) options operator, arrow; solrho(solphi(x)) end proc; rhosol := makeproc(rs(x), x); rhosol := optimize(rhosol); gr := proc (x) options operator, arrow; diff(rhosol(x), x) end proc; gradrho := makeproc(gr(x), x); gradrho := optimize(gradrho); C(rhosol, filename = outstring2, ansi); C(gradrho, filename = outstring2, ansi);

gp := D(solphi); gradphi := makeproc(gp(x), x); C(gradphi, filename = outstring2, ansi);
lp := D(gradphi); laplacephi := makeproc(lp(x), x);
tht := proc (x) options operator, arrow; reactionSource(rhosol(x), solphi(x))-delta*laplacephi(x) end proc; thetasol := makeproc(tht(x), x); thetasol := optimize(thetasol); C(thetasol, filename = outstring2, ansi);
RsP := proc (x) options operator, arrow; thetasol(rhosol(x), solphi(x))/rhosol(x) end proc; phiSource := makeproc(RsP(x), x); phiSource := optimize(phiSource); C(phiSource, filename = outstring2, ansi);
ms := proc (x) options operator, arrow; chemicalPotential(rhosol(x), solphi(x)) end proc;
musol := makeproc(ms(x), x); musol := optimize(musol); C(musol, filename = outstring2, ansi);
dsolmu := proc (x) options operator, arrow; diff(musol(x), x) end proc; gradmu := makeproc(dsolmu(x), x); gradmu := optimize(gradmu);
RsV := proc (x) options operator, arrow; -gradphi(x)*thetasol(x)+rhosol(x)*gradmu(x) end proc;
veloSource := makeproc(RsV(x), x); veloSource := optimize(veloSource); C(veloSource, filename = outstring2, ansi);


