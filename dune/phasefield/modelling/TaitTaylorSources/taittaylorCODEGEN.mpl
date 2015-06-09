restart; with(codegen); outstring2 := "0balanced.cc"; outstring := "maple.cc";
FvdW := proc (rho) options operator, arrow; rho*(-rho+(8/27)*theta*log(rho/(1-rho))) end proc; theta := .85;
gvdW := proc (rho) options operator, arrow; diff(FvdW(rho), rho) end proc; pp := proc (rho) options operator, arrow; -FvdW(rho)+rho*gvdW(rho) end proc; pvdW := makeproc(pp(rho), rho); GvdW := makeproc(gvdW(rho), rho);
sols := fsolve([pvdW(x) = pvdW(y), gvdW(x) = gvdW(y)], {x = 0 .. 3, y = .4 .. 1}); r1 := solve(sols[1]); r2 := solve(sols[2]); dpvdW := proc (rho) options operator, arrow; diff(pvdW(rho), rho) end proc; mwg := proc (rho) options operator, arrow; GvdW(r1)*rho-pvdW(r2) end proc; F2 := proc (rho) options operator, arrow; FvdW(rho)-mwg(rho) end proc; G2 := D(F2);
psols := solve(dpvdW(x), x); pp1 := psols[1]; pp2 := psols[2]; Gv2 := D(F2); p2 := proc (rho) options operator, arrow; -F2(rho)+rho*Gv2(rho) end proc;
sol3 := fsolve({G2(x) = 0}, x = .3 .. .4); r3 := solve(sol3[1]);
F1 := proc (rho) options operator, arrow; a*rho^7+b*rho+K end proc; F1(x); g1 := D(F1); g1(r2); g1(pp2); p1 := proc (rho) options operator, arrow; -F1(rho)+rho*g1(rho) end proc; G1 := makeproc(g1(x), x); G1(r2);

koeffs := fsolve({F1(r3) = F2(r3), G1(r2) = Gv2(r2), p1(r2) = p2(r2)});
K := solve(koeffs[1]); a := solve(koeffs[2]); b := solve(koeffs[3]);
F0 := proc (rho) options operator, arrow; av*rho*log(rho)+(bv-av)*rho+cv end proc; g0 := D(F0); p0 := proc (rho) options operator, arrow; -F0(rho)+rho*g0(rho) end proc; G0 := makeproc(g0(x), x);
koeffs2 := solve({F0(r3) = F2(r3), G0(r1) = Gv2(r1), p0(r1) = p2(r1)});
av := solve(koeffs2[1]); bv := solve(koeffs2[2]); cv := solve(koeffs2[3]);

nn := proc (phi) options operator, arrow; 6.*phi^5+(-1)*15.*phi^4+10.*phi^3 end proc; W := proc (phi) options operator, arrow; 2*A*(phi^4+(2*alpha-2)*phi^3+(-3*alpha+1)*phi^2+alpha) end proc; dW := D(W); F := proc (rho, phi) options operator, arrow; W(phi)/delta+beta*(nn(phi)*F1(rho)+(1-nn(phi))*F0(rho)) end proc; Pressure := proc (rho, phi) options operator, arrow; beta*(nn(phi)*p1(rho)+(1-nn(phi))*p0(rho)) end proc; Potential := proc (rho, phi) options operator, arrow; beta*(nn(phi)*G1(rho)+(1-nn(phi))*G0(rho)) end proc; P2 := proc (rho, phi) options operator, arrow; Pressure(rho, phi)-W(phi)/delta end proc;
solrho := proc (x) options operator, arrow; (r2-r1)*x+r1 end proc; Const := G1(sol1); solproc := makeproc(solrho(x), x); evalRho := optimize(solproc); C(evalRho, filename = outstring, ansi);
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

# 


