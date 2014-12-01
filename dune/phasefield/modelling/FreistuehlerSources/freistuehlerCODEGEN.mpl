??print(); # input placeholder
restart; with(codegen);
F1 := proc (rho) options operator, arrow; b*rho*ln(rho) end proc; g1 := D(F1);
F0 := proc (rho) options operator, arrow; e*rho*ln(rho) end proc; g0 := D(F0); solve(g0(rho) = 0, rho);
b := 1; e := 2;
G1 := D(F1); G0 := D(F0); p1 := proc (rho) options operator, arrow; -F1(rho)+rho*G1(rho) end proc; p0 := proc (rho) options operator, arrow; -F0(rho)+rho*G0(rho) end proc;
sols := fsolve([p1(rho) = p0(x), G1(rho) = G0(x)], {rho = 0 .. 4, x = 0 .. 10}); outstring := "maple.cc"; sol1 := solve(sols[1]); sol2 := solve(sols[2]);
nu := proc (phi) options operator, arrow; phi end proc; W := proc (phi) options operator, arrow; phi*log(phi)+(1-phi)*log(1-phi)+(-1)*0*phi^2 end proc; dW := D(W); F := proc (rho, phi) options operator, arrow; theta*rho*W(phi)+nu(phi)*F1(rho)+(1-nu(phi))*F0(rho) end proc; Pressure := proc (rho, phi) options operator, arrow; nu(phi)*p1(rho)+(1-nu(phi))*p0(rho) end proc; Potential := proc (rho, phi) options operator, arrow; W(phi)+nu(phi)*G1(rho)+(1-nu(phi))*G0(rho) end proc;
solrho := proc (x) options operator, arrow; 0 end proc; Const := G1(sol1); solproc := makeproc(solrho(x), x); evalRho := optimize(solproc); C(evalRho, filename = outstring, ansi);
Fproc := makeproc(F(rho, phi), [rho, phi]); helmholtz := optimize(Fproc); C(helmholtz, filename = outstring, ansi);
S := proc (rho, phi) options operator, arrow; diff(F(rho, phi), phi) end proc; Sproc := makeproc(S(rho, phi), [rho, phi]); reactionSource := optimize(Sproc); C(reactionSource, filename = outstring, ansi);
dphiS := proc (rho, phi) options operator, arrow; diff(S(rho, phi), phi) end proc; dphiSproc := makeproc(dphiS(rho, phi), [rho, phi]); dphireactionSource := optimize(dphiSproc); C(dphireactionSource, filename = outstring, ansi);
CProc := makeproc(simplify(Potential(rho, phi)), [rho, phi]); chemicalPotential := optimize(CProc); C(chemicalPotential, filename = outstring, ansi);
drhoPotential := proc (rho, phi) options operator, arrow; diff(Potential(rho, phi), rho) end proc; drhomuproc := makeproc(simplify(drhoPotential(rho, phi)), [rho, phi]); drhochemicalPotential := optimize(drhomuproc); C(drhochemicalPotential, filename = outstring, ansi);
dphiPotential := proc (rho, phi) options operator, arrow; diff(Potential(rho, phi), phi) end proc; dphimuproc := makeproc(simplify(dphiPotential(rho, phi)), [rho, phi]); dphichemicalPotential := optimize(dphimuproc); C(dphichemicalPotential, filename = outstring, ansi);
Pproc := makeproc(simplify(Pressure(rho, phi)), [rho, phi]); pressure := optimize(Pproc); C(pressure, filename = outstring, ansi);
wavespeed := proc (rho, phi) options operator, arrow; sqrt(diff(Pressure(rho, phi), rho)) end proc; wproc := makeproc(simplify(wavespeed(rho, phi)), [rho, phi]); a := optimize(wproc); C(a, filename = outstring, ansi);


