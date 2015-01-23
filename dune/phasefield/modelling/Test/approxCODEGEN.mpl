restart; with(codegen); outstring := "maple.cc";
# 
w := proc (rho) options operator, arrow; -3*rho+(8/3)*theta*log(rho/(3-rho))+d*theta*(1-log(theta)) end proc; d := 0; theta := .7;
g := proc (rho) options operator, arrow; diff(rho*w(rho), rho) end proc; gg := makeproc(g(rho), [rho]); p := proc (rho) options operator, arrow; rho^2*(diff(w(rho), rho)) end proc; pp := makeproc(p(rho), rho);
# solve for maxwell points
sols := fsolve([p(x) = p(y), g(x) = g(y)], {x = 0 .. 1.5, y = 1.5 .. 3}); r1 := solve(sols[1]); r2 := solve(sols[2]);
# second derivative of the free Energy
ddw := proc (rho) options operator, arrow; diff(g(rho), rho) end proc; ddwproc := makeproc(ddw(rho), [rho]);
# 
psiV := r1*w(r1);
# Chemical Potential at first Maxwell Point
dpsiV := gg(r1);
# Second derivative at first Maxwell Point
ddpsiV := ddwproc(r1);
# Energy at second MaxwellPoint
psiL := r2*w(r2);
# Chemical Potential at second Maxwell Point
dpsiL := gg(r2);
# Second derivative at second Maxwell Point
ddpsiL := ddwproc(r2);
# 
fV := proc (rho) options operator, arrow; av*rho*log(rho)+bv-av+cv end proc;

# 
dfV := D(fV);
# 
ddfV := D(dfV);
# 
av := solve(ddfV(r1) = ddpsiV, av);
# 
bv := solve(dfV(r1) = dpsiV, bv);
# 
cv := solve(fV(r1) = psiV, cv);
pv := proc (rho) options operator, arrow; rho*dfV(rho)-fV(rho) end proc;
# 
fL := proc (rho) options operator, arrow; al*rho*log(rho)+(bl-al)*rho+cl end proc;
# 
dfL := D(fL);
# 
ddfL := D(dfL);
pL := proc (rho) options operator, arrow; rho*dfL(rho)-fL(rho) end proc;
# 
al := solve(ddfL(r2) = ddpsiL, al);
# 
bl := solve(dfL(r2) = dpsiL, bl);
cl := solve(fL(r2) = psiL, cl);
f1 := proc (rho) options operator, arrow; piecewise(rho <= r1, rho*w(rho), fV(rho)) end proc; f0 := proc (rho) options operator, arrow; piecewise(rho < r2, fL(rho), rho*w(rho)) end proc;
mm1 := proc (rho) options operator, arrow; piecewise(rho <= r1, gg(rho), dfV(rho)) end proc; mm0 := proc (rho) options operator, arrow; piecewise(rho < r2, dfL(rho), gg(rho)) end proc;
pp1 := proc (rho) options operator, arrow; piecewie(rho <= r1, pp(rho), pv(rho)) end proc; pp0 := proc (rho) options operator, arrow; piecewise(rho < r2, pL(rho), pp(rho)) end proc;
ddm1 := proc (rho) options operator, arrow; piecewise(rho < r1, ddwproc(rho), ddfV(rho)) end proc; ddm2 := proc (rho) options operator, arrow; piecewise(rho < r2, ddfV(rho), ddwproc(rho)) end proc;
nn := proc (phi) options operator, arrow; 6.*phi^5+(-1)*15.*phi^4+10.*phi^3 end proc; W := proc (phi) options operator, arrow; 2*A*(phi^4+(2*alpha-2)*phi^3+(-3*alpha+1)*phi^2+alpha) end proc; dW := D(W); F := proc (rho, phi) options operator, arrow; W(phi)/delta+nn(phi)*ff1+(1-nn(phi))*ff2 end proc; Pressure := proc (rho, phi) options operator, arrow; nn(phi)*ppp1+(1-nn(phi))*ppp0 end proc; Potential := proc (rho, phi) options operator, arrow; nn(phi)*m1+(1-nn(phi))*m0 end proc; P2 := proc (rho, phi) options operator, arrow; Pressure(rho, phi)-W(phi)/delta end proc;
Fproc := makeproc(F(rho, phi), [rho, phi]); helmholtz := optimize(Fproc); C(helmholtz, filename = outstring, ansi);
S := proc (rho, phi) options operator, arrow; diff(F(rho, phi), phi) end proc; Sproc := makeproc(S(rho, phi), [rho, phi]); reactionSource := optimize(Sproc); C(reactionSource, filename = outstring, ansi);
dphiS := proc (rho, phi) options operator, arrow; diff(S(rho, phi), phi) end proc; dphiSproc := makeproc(dphiS(rho, phi), [rho, phi]); dphireactionSource := optimize(dphiSproc); C(dphireactionSource, filename = outstring, ansi);
CProc := makeproc(Potential(rho, phi), [rho, phi]); chemicalPotential := optimize(CProc); C(chemicalPotential, filename = outstring, ansi);
dphiPotential := proc (rho, phi) options operator, arrow; diff(Potential(rho, phi), phi) end proc; dphimuproc := makeproc(dphiPotential(rho, phi), [rho, phi]); dphichemicalPotential := optimize(dphimuproc); C(dphichemicalPotential, filename = outstring, ansi);
drhoPotential := proc (rho, phi) options operator, arrow; nn(phi)*dmu1+(1-nn(phi))*dmu0 end proc; drhomuproc := makeproc(drhoPotential(rho, phi), [rho, phi]); drhochemicalPotential := optimize(drhomuproc); C(drhochemicalPotential, filename = outstring, ansi);
Pproc := makeproc(P2(rho, phi), [rho, phi]); pressure := optimize(Pproc); C(pressure, filename = outstring, ansi);
wavespeed := proc (rho, phi) options operator, arrow; sqrt(nn(phi)*dp1+(1-nn(phi))*dp0) end proc; wproc := makeproc(wavespeed(rho, phi), [rho, phi]); a := optimize(wproc); C(a, filename = outstring, ansi);
f1proc := makeproc(f1(rho), rho); psi1 := optimize(f1proc); psi1 := prep2trans(psi1); C(psi1, filename = outstring, ansi);
f2proc := makeproc(f0(rho), rho); psi2 := optimize(f2proc); psi2 := prep2trans(psi1); C(psi2, filename = outstring, ansi);
mm1proc := makeproc(mm1(rho), rho); mu1 := optimize(mm1proc); mu1 := prep2trans(mu1); C(mu1, filename = outstring, ansi);
mm2proc := makeproc(mm0(rho), rho); mu2 := optimize(mm2proc); mu2 := prep2trans(mu2); C(mu2, filename = outstring, ansi);
pp1proc := makeproc(mm1(rho), rho); p1 := optimize(pp1proc); p1 := prep2trans(p1); C(p1, filename = outstring, ansi);
pp2proc := makeproc(mm2(rho), rho); p2 := optimize(pp1proc); p2 := prep2trans(p2); C(p2, filename = outstring, ansi);
ddm1proc := makeproc(ddm1(rho), rho); dmu1 := optimize(ddm1proc); dmu1 := prep2trans(dmu1); C(dmu1, filename = outstring, ansi);
ddm2proc := makeproc(ddm2(rho), rho); dmu2 := optimize(ddm2proc); dmu2 := prep2trans(dmu2); C(dmu2, filename = outstring, ansi);


