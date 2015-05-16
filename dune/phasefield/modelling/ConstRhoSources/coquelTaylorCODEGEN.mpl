restart; with(codegen); outstring2 := "balanced.cc"; outstring := "maple.cc";
F11 := proc (rho) options operator, arrow; (b-c)*rho+c*rho*ln(rho)+d end proc; g1 := D(F11);
F00 := proc (rho) options operator, arrow; (f-e)*rho+e*rho*ln(rho)+g end proc; g0 := D(F00);
c := 1.5; b := log(2); d := 0; e := 1; f := 0; z1 := solve(g1(rho) = 0, rho); z2 := solve(g0(rho) = 0.*rho); gg := F1(z1)-e*z2*ln(z2)-(f-e)*z2; g := .5;
G1 := D(F11); G0 := D(F00); p1 := proc (rho) options operator, arrow; -F11(rho)+rho*G1(rho) end proc; p0 := proc (rho) options operator, arrow; -F00(rho)+rho*G0(rho) end proc;
sols := fsolve([p1(rho) = p0(x), G1(rho) = G0(x)], {rho = 0 .. 4, x = 0 .. 4}); sol1 := solve(sols[1]); sol2 := solve(sols[2]); mwg := proc (rho) options operator, arrow; G1(sol1)*rho-p1(sol1) end proc; F1 := F11-mwg; F0 := F00-mwg;
nn := proc (phi) options operator, arrow; 6.*phi^5+(-1)*15.*phi^4+10.*phi^3 end proc; W := proc (phi) options operator, arrow; 2*A*(phi^4+(2*alpha-2)*phi^3+(-3*alpha+1)*phi^2+alpha) end proc; dW := D(W); F := proc (rho, phi) options operator, arrow; W(phi)/delta+beta*(nn(phi)*F1(rho)+(1-nn(phi))*F0(rho)) end proc; Pressure := proc (rho, phi) options operator, arrow; beta*(nn(phi)*p1(rho)+(1-nn(phi))*p0(rho)) end proc; Potential := proc (rho, phi) options operator, arrow; beta*(nn(phi)*G1(rho)+(1-nn(phi))*G0(rho)) end proc; P2 := proc (rho, phi) options operator, arrow; Pressure(rho, phi)-W(phi)/delta end proc;
solrho := proc (x) options operator, arrow; exp((Const-nn(x)*(b-f)-f)/(nn(x)*(c-e)+e)) end proc; Const := G1(sol1); solproc := makeproc(solrho(x), x); evalRho := optimize(solproc); C(evalRho, filename = outstring, ansi);

Fproc := makeproc(F(rho, phi), [rho, phi]); helmholtz := optimize(Fproc); C(helmholtz, filename = outstring, ansi);
fp3 := proc (rho, phi) options operator, arrow; diff(helmholtz(rho, phi), phi, phi, phi) end proc; fp1 := proc (rho, phi) options operator, arrow; diff(helmholtz(rho, phi), phi) end proc; Sproc := makeproc(fp1(rho, phiMid)+(1/24)*fp3(rho, phiMid)*(phi-phiOld), [rho, phi, phiOld, phiMid]); reactionSource := optimize(Sproc); C(reactionSource, filename = outstring, ansi);
dphiS := proc (rho, phi) options operator, arrow; 0 end proc; dphiSproc := makeproc(dphiS(rho, phi), [rho, phi]); dphireactionSource := optimize(dphiSproc); C(dphireactionSource, filename = outstring, ansi);
fr3 := proc (rho, phi) options operator, arrow; diff(helmholtz(rho, phi), rho, rho, rho) end proc; fr1 := proc (rho, phi) options operator, arrow; diff(helmholtz(rho, phi), rho) end proc; CProc := makeproc(fr1(rhoMid, phi)+(1/24)*fr3(rhoMid, phi)*(rho-rhoOld), [rho, phi, rhoOld, phiMid]);
chemicalPotential := optimize(CProc); C(chemicalPotential, filename = outstring, ansi);
drhoPotential := proc (rho, phi) options operator, arrow; diff(Potential(rho, phi), rho) end proc; drhomuproc := makeproc(simplify(drhoPotential(rho, phi)), [rho, phi]); drhochemicalPotential := optimize(drhomuproc); C(drhochemicalPotential, filename = outstring, ansi);
dphiPotential := proc (rho, phi) options operator, arrow; diff(Potential(rho, phi), phi) end proc; dphimuproc := makeproc(simplify(dphiPotential(rho, phi)), [rho, phi]); dphichemicalPotential := optimize(dphimuproc); C(dphichemicalPotential, filename = outstring, ansi);
Pproc := makeproc(P2(rho, phi), [rho, phi]); pressure := optimize(Pproc); C(pressure, filename = outstring, ansi);
wavespeed := proc (rho, phi) options operator, arrow; sqrt(diff(Pressure(rho, phi), rho)) end proc; wproc := makeproc(simplify(wavespeed(rho, phi)), [rho, phi]); a := optimize(wproc); C(a, filename = outstring, ansi);
solphi := proc (x) options operator, arrow; .5*tanh(x/delta)+.5 end proc; rs := proc (x) options operator, arrow; solrho(solphi(x)) end proc; rhosol := makeproc(rs(x), x); rhosol := optimize(rhosol); gr := proc (x) options operator, arrow; diff(rhosol(x), x) end proc; gradrho := makeproc(gr(x), x); gradrho := optimize(gradrho); C(rhosol, filename = outstring2, ansi); C(gradrho, filename = outstring2, ansi);
gp := D(solphi); gradphi := makeproc(gp(x), x); C(gradphi, filename = outstring2, ansi);
lp := D(gradphi); laplacephi := makeproc(lp(x), x);
tht := proc (x) options operator, arrow; 0 end proc; thetasol := makeproc(tht(x), x); thetasol := optimize(thetasol); C(thetasol, filename = outstring2, ansi);
RsP := proc (x) options operator, arrow; thetasol(rhosol(x), solphi(x))/rhosol(x) end proc; phiSource := makeproc(RsP(x), x); phiSource := optimize(phiSource); C(phiSource, filename = outstring2, ansi);
ms := proc (x) options operator, arrow; chemicalPotential(rhosol(x), solphi(x)) end proc;
musol := makeproc(ms(x), x); musol := optimize(musol); C(musol, filename = outstring2, ansi);
dsolmu := proc (x) options operator, arrow; diff(musol(x), x) end proc; gradmu := makeproc(dsolmu(x), x); gradmu := optimize(gradmu);
RsV := proc (x) options operator, arrow; -gradphi(x)*thetasol(x)+rhosol(x)*gradmu(x) end proc;
veloSource := makeproc(RsV(x), x); veloSource := optimize(veloSource); C(veloSource, filename = outstring2, ansi);

fr3 := diff(helmholtz(rho, phi), rho, rho, rho); fr1 := diff(helmholtz(rho, phi), rho); Tfr := fr1+(1/24)*fr3*(r1-r2)^2;


# 
