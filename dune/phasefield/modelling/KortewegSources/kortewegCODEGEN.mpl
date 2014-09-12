restart; with(codegen); c := 4;
F1 := proc (rho) options operator, arrow; rho*(-3*rho+(8/3)*theta*log(rho/(3-rho))+c*theta*(1-log(theta))) end proc; theta := .9;
G1 := D(F1); p1 := proc (rho) options operator, arrow; -F1(rho)+rho*G1(rho) end proc;
pp := makeproc(p1(rho), rho); gg := makeproc(G1(rho), [rho]);

outstring := "maple.c";
F := proc (rho) options operator, arrow; F1(rho) end proc; Pressure := proc (rho) options operator, arrow; p1(rho) end proc; Potential := proc (rho) options operator, arrow; G1(rho) end proc; simplify(Potential(rho));

Fproc := makeproc(F(rho), [rho]); helmholtz := optimize(Fproc); C(helmholtz, filename = outstring, ansi);


CProc := makeproc(Potential(rho), [rho]); chemicalPotential := optimize(CProc); C(chemicalPotential, filename = outstring, ansi);
drhoPotential := proc (rho) options operator, arrow; diff(Potential(rho), rho) end proc; drhomuproc := makeproc(simplify(drhoPotential(rho)), [rho]); drhochemicalPotential := optimize(drhomuproc); C(drhochemicalPotential, filename = outstring, ansi);

Pproc := makeproc(Pressure(rho), [rho]); pressure := optimize(Pproc); C(pressure, filename = outstring, ansi);
wavespeed := proc (rho, phi) options operator, arrow; diff(Pressure(rho, phi), rho) end proc; wproc := makeproc(simplify(wavespeed(rho, phi)), [rho, phi]); a := optimize(wproc); C(a, filename = outstring, ansi);


# 
