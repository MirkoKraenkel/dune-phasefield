restart; with(codegen); with(plots); outstring2 := "balanced.cc"; outstring := "maple.cc"; with(CurveFitting);
FvdW := proc (rho) options operator, arrow; rho*(-rho+(8/27)*theta*log(rho/(1-rho))) end proc; theta := .85; alpha := 0;
gvdW := proc (rho) options operator, arrow; diff(FvdW(rho), rho) end proc; pp := proc (rho) options operator, arrow; -FvdW(rho)+rho*gvdW(rho) end proc; pvdW := makeproc(pp(rho), rho); GvdW := makeproc(gvdW(rho), rho);
sols := fsolve([pvdW(x) = pvdW(y), gvdW(x) = gvdW(y)], {x = 0 .. 3, y = .4 .. 1}); r1 := solve(sols[1]); r2 := solve(sols[2]); dpvdW := proc (rho) options operator, arrow; diff(pvdW(rho), rho) end proc; mwg := proc (rho) options operator, arrow; GvdW(r1)*rho-pvdW(r2) end proc; F2 := proc (rho) options operator, arrow; FvdW(rho)-mwg(rho) end proc; G2 := D(F2);
psols := solve(dpvdW(x), x); pp1 := psols[1]; pp2 := psols[2]; Gv2 := D(F2); p2 := proc (rho) options operator, arrow; -F2(rho)+rho*Gv2(rho) end proc; C2 := makeproc(diff(pvdW(rho), rho), [rho]);
F1 := proc (rho) options operator, arrow; i*rho/(rho-1)+c*rho*log(rho)+(b-c)*rho+d end proc; F1(x); g1 := D(F1); g1(r2); g1(pp2); p1 := proc (rho) options operator, arrow; -F1(rho)+rho*g1(rho) end proc; G1 := makeproc(g1(x), x); G1(r2); C1 := makeproc(diff(p1(rho), rho), [rho]);
sol3 := fsolve({G2(x) = 0}, x = .3 .. .4); r3 := solve(sol3[1]);
koeffs := fsolve({C1(pp1) = 0, G1(r2) = Gv2(r2), p1(pp2) = p2(pp2), p1(r2) = p2(r2)});
b := solve(koeffs[1]); c := solve(koeffs[2]); d := solve(koeffs[3]); i := solve(koeffs[4]);
F0 := proc (rho) options operator, arrow; h*rho^2+e*rho*log(rho)+(f-e)*rho+g end proc; g0 := D(F0); p0 := proc (rho) options operator, arrow; -F0(rho)+rho*g0(rho) end proc; G0 := makeproc(g0(x), x); C0 := makeproc(diff(p0(rho), rho), [rho]);
koeffs2 := solve({C0(pp2) = 0, G0(r1) = Gv2(r1), p0(pp1) = p2(pp1), p0(r1) = p2(r1)});
e := solve(koeffs2[1]); f := solve(koeffs2[2]); g := solve(koeffs2[3]); h := solve(koeffs2[4]);

nn := proc (phi) options operator, arrow; 6.*phi^5+(-1)*15.*phi^4+10.*phi^3 end proc; W := proc (phi) options operator, arrow; 2*phi^4+2*(2*alpha-2)*phi^3+2*(-3*alpha+1)*phi^2+2*alpha end proc; dW := D(W); F := proc (rho, phi) options operator, arrow; A*W(phi)/delta+nn(phi)*F1(rho)+(1-nn(phi))*F0(rho) end proc; Pressure := proc (rho, phi) options operator, arrow; nn(phi)*p1(rho)+(1-nn(phi))*p0(rho) end proc; Potential := proc (rho, phi) options operator, arrow; nn(phi)*G1(rho)+(1-nn(phi))*G0(rho) end proc; P2 := proc (rho, phi) options operator, arrow; Pressure(rho, phi)-W(phi)/delta end proc; rhosol := proc (phi) options operator, arrow; phi*r2+(1-phi)*r1 end proc;
Fproc := makeproc(F(rho, phi), [rho, phi]); helmholtz := optimize(Fproc); C(helmholtz, filename = outstring, ansi);
fp3 := proc (rho, phi) options operator, arrow; diff(helmholtz(rho, phi), phi, phi, phi) end proc; fp1 := proc (rho, phi) options operator, arrow; diff(helmholtz(rho, phi), phi) end proc; SourceProc := makeproc(fp1(x, y), [x, y]); Sproc := makeproc(fp1(rho, mid)+(1/24)*fp3(rho, mid)*(phi-old), [rho, phi, old, mid]); SSproc := makeproc(Sproc(rho, phi, old, .5*(phi+old)), [rho, phi, old]); reactionSource := optimize(SSproc); C(reactionSource, filename = outstring, ansi);

Fbulk := proc (rho, phi) options operator, arrow; beta*(nn(phi)*F1(rho)+(1-nn(phi))*F0(rho)) end proc;
dphiS := proc (rho, phi, old) options operator, arrow; diff(SSproc(rho, phi, old), phi) end proc; dphiSproc := makeproc(dphiS(rho, phi, old, mid), [rho, phi, old, mid]); dphireactionSource := optimize(dphiSproc); C(dphireactionSource, filename = outstring, ansi);
fr3 := proc (rho, phi) options operator, arrow; diff(helmholtz(rho, phi), rho, rho, rho) end proc; fr1 := proc (rho, phi) options operator, arrow; diff(helmholtz(rho, phi), rho) end proc; CProc := makeproc(fr1(mid, phi)+(1/24)*fr3(mid, phi)*(rho-old), [rho, phi, old, mid]); CCProc := makeproc(CProc(rho, phi, old, .5*(rho+old)), [rho, phi, old]);
chemicalPotential := optimize(CCProc); C(chemicalPotential, filename = outstring, ansi);
drhoPotential := proc (rho, phi, old, mid) options operator, arrow; diff(CCProc(rho, phi, old), rho) end proc; drhomuproc := makeproc(drhoPotential(rho, phi, old), [rho, phi, old]); drhochemicalPotential := optimize(drhomuproc); C(drhochemicalPotential, filename = outstring, ansi);
dphiPotential := proc (rho, phi, old, mid) options operator, arrow; diff(CCProc(rho, phi, old), phi) end proc; dphimuproc := makeproc(dphiPotential(rho, phi, old), [rho, phi, old]); dphichemicalPotential := optimize(dphimuproc); C(dphichemicalPotential, filename = outstring, ansi);
# 
Pproc := makeproc(P2(rho, phi), [rho, phi]); pressure := optimize(Pproc); C(pressure, filename = outstring, ansi);
wavespeed := proc (rho, phi) options operator, arrow; sqrt(diff(Pressure(rho, phi), rho)) end proc; wproc := makeproc(simplify(wavespeed(rho, phi)), [rho, phi]); a := optimize(wproc); C(a, filename = outstring, ansi);
A := 0.1e-2; delta := 0.5e-1; inttanh := 0; Re(evalf(int(sqrt(2*A*delta*Fproc(rhosol(x), x)), x = 0. .. 1.))); Re(evalf(int(sqrt(2*W(x)), x = 0. .. 1.)));
phisol := proc (x) options operator, arrow; .5*tanh(x/delta)+.5 end proc; lap := proc (x) options operator, arrow; diff(phisol(x), x, x) end proc;
sols := fsolve([p1(rho)+inttanh/(.25) = p0(x), G1(rho) = G0(x)], {rho = .2 .. 1, x = 0 .. 4}); bubblesol := solve(sols[1]); bubblesol2 := solve(sols[2]); ml := proc (a) options operator, arrow; bubblesol end proc; mwpliq := makeproc(ml(x), x);
mv := proc (a) options operator, arrow; bubblesol2 end proc; mwpvap := makeproc(mv(x), x); C(mwpliq, filename = outstring2, ansi); C(mwpvap, filename = outstring2, ansi); Const := G1(bubblesol); G0(bubblesol2);
numP := 20; CC := Array(1 .. numP); B := Array(1 .. numP);

for k to numP do CC[k] := (k-1)/(numP-1); B[k] := fsolve(Potential(x, (k-1)/numP) = G1(bubblesol), x = 0.1e-1 .. 1) end do; pair := proc (x, y) options operator, arrow; [x, y] end proc; rhoI := makeproc(Spline(CC, B, z), [z]);
evalRho1 := optimize(rhoI); evalRho2 := prep2trans(evalRho1); C(evalRho2, filename = outstring2, ansi);

# 

# 



# 
