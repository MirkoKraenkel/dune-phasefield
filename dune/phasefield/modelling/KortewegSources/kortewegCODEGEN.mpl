restart; with(codegen); with(plots); outstring2 := "nskRho.cc"; outstring := "nskmaple.cc";
FvdW := proc (rho) options operator, arrow; rho*(-rho+(8/27)*theta*log(rho/(1-rho))) end proc; theta := .85; alpha := 0;
gvdW := proc (rho) options operator, arrow; diff(FvdW(rho), rho) end proc; pp := proc (rho) options operator, arrow; -FvdW(rho)+rho*gvdW(rho) end proc; pvdW := makeproc(pp(rho), rho); GvdW := makeproc(gvdW(rho), rho);
sols := fsolve([pvdW(x) = pvdW(y), gvdW(x) = gvdW(y)], {x = 0 .. 3, y = .4 .. 1}); r1 := solve(sols[1]); r2 := solve(sols[2]); dpvdW := proc (rho) options operator, arrow; diff(pvdW(rho), rho) end proc; mwg := proc (rho) options operator, arrow; GvdW(r1)*rho-pvdW(r2) end proc; F2 := proc (rho) options operator, arrow; FvdW(rho)-mwg(rho) end proc; G2 := D(F2);
psols := solve(dpvdW(x), x); pp1 := psols[1]; pp2 := psols[2]; Gv2 := D(F2); p2 := proc (rho) options operator, arrow; -F2(rho)+rho*Gv2(rho) end proc; C2 := makeproc(diff(pvdW(rho), rho), [rho]);


Fproc := makeproc(FvdW(rho), [rho]); helmholtz := optimize(Fproc); C(helmholtz, filename = outstring, ansi);


fr3 := proc (rho) options operator, arrow; diff(helmholtz(rho), rho, rho, rho) end proc; fr1 := proc (rho) options operator, arrow; diff(helmholtz(rho), rho) end proc; CProc := makeproc(fr1(mid)+(1/24)*fr3(mid)*(rho-old), [rho, old, mid]); CCProc := makeproc(CProc(rho, old, .5*(rho+old)), [rho, old]);
chemicalPotential := optimize(CCProc); C(chemicalPotential, filename = outstring, ansi);
drhoPotential := proc (rho, old, mid) options operator, arrow; diff(CCProc(rho, old), rho) end proc; drhomuproc := makeproc(drhoPotential(rho, old), [rho, old]); drhochemicalPotential := optimize(drhomuproc); C(drhochemicalPotential, filename = outstring, ansi);

Pproc := makeproc(pvdW(rho), [rho]); pressure := optimize(Pproc); C(pressure, filename = outstring, ansi);
wavespeed := proc (rho, phi) options operator, arrow; sqrt(max(diff(pressure(rho), rho), 0)) end proc; wproc := makeproc(simplify(wavespeed(rho, phi)), [rho, phi]); a := optimize(wproc); C(a, filename = outstring, ansi);
sols := fsolve([pvdW(rho) = pvdW(x), gvdW(rho) = gvdW(x)], {rho = .2 .. 1, x = 0 .. 4}); bubblesol := solve(sols[1]); bubblesol2 := solve(sols[2]); ml := proc (a) options operator, arrow; bubblesol end proc; mwpliq := makeproc(ml(x), x);
mv := proc (a) options operator, arrow; bubblesol2 end proc; mwpvap := makeproc(mv(x), x); C(mwpliq, filename = outstring2, ansi); C(mwpvap, filename = outstring2, ansi); Const := G1(bubblesol); G0(bubblesol2);


# 

# 



# 
