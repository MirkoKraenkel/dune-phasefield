restart; with(codegen);
w := proc (rho) options operator, arrow; -3*rho+(8/3)*theta*log(rho/(3-rho))+c*theta*(1-log(theta)) end proc; g := proc (rho) options operator, arrow; diff(rho*w(rho), rho) end proc;
c := 4; theta := .9;
p := proc (rho) options operator, arrow; rho^2*(diff(w(rho), rho)) end proc; pp := makeproc(p(rho), rho); gg := makeproc(g(rho), [rho]);
sols := fsolve([p(x) = p(y), g(x) = g(y)], {x = 0 .. 1.5, y = 1.5 .. 3}); r1 := solve(sols[1]); r2 := solve(sols[2]);
F1 := proc (rho) options operator, arrow; a2*rho*ln(rho)+(b2-a2)*rho+k end proc; G1 := proc (rho) options operator, arrow; diff(F1(rho), rho) end proc; gg1 := makeproc(G1(rho), rho); p1 := proc (rho) options operator, arrow; -F1(rho)+rho*gg1(rho) end proc; pp1 := makeproc(p1(rho), rho);
F0 := proc (rho) options operator, arrow; a*rho*ln(rho)+(b-a)*rho end proc; G0 := D(F0); p0 := proc (rho) options operator, arrow; -F1(rho)+rho*gg1(rho) end proc; pp0 := makeproc(p1(rho), rho);
W := proc (rho) options operator, arrow; rho*w(rho) end proc; a := pp(r1)/r1; b := gg(r1)-a*ln(r1);
a2 := (pp(r2)+k)/r2; b2 := gg(r2)-a2*ln(r2); k := 11;

nn := proc (phi) options operator, arrow; 6.*phi^5+(-1)*15.*phi^4+10.*phi^3 end proc; W := proc (phi) options operator, arrow; 2.*(phi^2*(1-phi)^2) end proc; dW := D(W);
F := proc (rho, phi) options operator, arrow; W(phi)/delta+nn(phi)*F1(rho)+(1-nn(phi))*F0(rho) end proc;
Pressure := proc (rho, phi) options operator, arrow; nn(phi)*p1(rho)+(1-nn(phi))*p0(rho) end proc;
P2 := proc (rho, phi) options operator, arrow; Pressure(rho, phi)-W(phi)/delta end proc;
Potential := proc (rho, phi) options operator, arrow; nn(phi)*G1(rho)+(1-nn(phi))*G0(rho) end proc;

outstring := "maple.c";

solrho := proc (x) options operator, arrow; 0 end proc; solproc := makeproc(solrho(x), x); evalRho := optimize(solproc); C(evalRho, filename = outstring, ansi);
Fproc := makeproc(F(rho, phi), [rho, phi]); helmholtz := optimize(Fproc); C(helmholtz, filename = outstring, ansi);
S := proc (rho, phi) options operator, arrow; diff(F(rho, phi), phi) end proc; Sproc := makeproc(S(rho, phi), [rho, phi]); reactionSource := optimize(Sproc); C(reactionSource, filename = outstring, ansi);
dphiS := proc (rho, phi) options operator, arrow; diff(S(rho, phi), phi) end proc; dphiSproc := makeproc(dphiS(rho, phi), [rho, phi]); dphireactionSource := optimize(dphiSproc); C(dphireactionSource, filename = outstring, ansi);
CProc := makeproc(simplify(Potential(rho, phi)), [rho, phi]); chemicalPotential := optimize(CProc); C(chemicalPotential, filename = outstring, ansi);
drhoPotential := proc (rho, phi) options operator, arrow; diff(Potential(rho, phi), rho) end proc; drhomuproc := makeproc(simplify(drhoPotential(rho, phi)), [rho, phi]); drhochemicalPotential := optimize(drhomuproc); C(drhochemicalPotential, filename = outstring, ansi);
dphiPotential := proc (rho, phi) options operator, arrow; diff(Potential(rho, phi), phi) end proc; dphimuproc := makeproc(simplify(dphiPotential(rho, phi)), [rho, phi]); dphichemicalPotential := optimize(dphimuproc); C(dphichemicalPotential, filename = outstring, ansi);
Pproc := makeproc(P2(rho, phi), [rho, phi]); pressure := optimize(Pproc); C(pressure, filename = outstring, ansi);
wavespeed := proc (rho, phi) options operator, arrow; sqrt(diff(Pproc(rho, phi), rho)) end proc; wproc := makeproc(simplify(wavespeed(rho, phi)), [rho, phi]); a := optimize(wproc); C(a, filename = outstring, ansi);


# 
