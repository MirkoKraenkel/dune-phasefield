restart; with(codegen); outstring2 := "stiffenedgasRho.cc"; outstring := "stiffenedgasmaple.cc";
F11 := proc (rho) options operator, arrow; (bl-cl)*rho+cl*rho*ln(rho)+dl end proc; g1 := D(F11);
F00 := proc (rho) options operator, arrow; (bv-cv)*rho+cv*rho*ln(rho)+dv end proc; g0 := D(F00);
G1 := D(F11); G0 := D(F00); p1 := proc (rho) options operator, arrow; -F11(rho)+rho*G1(rho) end proc; p0 := proc (rho) options operator, arrow; -F00(rho)+rho*G0(rho) end proc;
nn := proc (phi) options operator, arrow; 6.*phi^5+(-1)*15.*phi^4+10.*phi^3 end proc; W := proc (phi) options operator, arrow; 2*A*(phi^4+(2*alpha-2)*phi^3+(-3*alpha+1)*phi^2+alpha) end proc; dW := D(W); F := proc (rho, phi) options operator, arrow; W(phi)/delta+beta*(nn(phi)*F11(rho)+(1-nn(phi))*F00(rho)) end proc; Pressure := proc (rho, phi) options operator, arrow; beta*(nn(phi)*p1(rho)+(1-nn(phi))*p0(rho)) end proc; Potential := proc (rho, phi) options operator, arrow; beta*(nn(phi)*G1(rho)+(1-nn(phi))*G0(rho)) end proc; P2 := proc (rho, phi) options operator, arrow; Pressure(rho, phi)-W(phi)/delta end proc;
Fproc := makeproc(F(rho, phi), [rho, phi]); helmholtz := optimize(Fproc); C(helmholtz, filename = outstring, ansi);
fp3 := proc (rho, phi) options operator, arrow; diff(helmholtz(rho, phi), phi, phi, phi) end proc; fp1 := proc (rho, phi) options operator, arrow; diff(helmholtz(rho, phi), phi) end proc; Sproc := makeproc(fp1(rho, mid)+(1/24)*fp3(rho, mid)*(phi-old), [rho, phi, old, mid]); SSproc := makeproc(Sproc(rho, phi, old, .5*(phi+old)), [rho, phi, old]); reactionSource := optimize(SSproc); C(reactionSource, filename = outstring, ansi); diff(Sproc(rho, phi, old, mid), phi);
dphiS := proc (rho, phi, old) options operator, arrow; diff(SSproc(rho, phi, old), phi) end proc; dphiSproc := makeproc(dphiS(rho, phi, old, mid), [rho, phi, old, mid]); dphireactionSource := optimize(dphiSproc); C(dphireactionSource, filename = outstring, ansi);
drhoS := proc (rho, phi, old) options operator, arrow; diff(SSproc(rho, phi, old), rho) end proc; drhoSproc := makeproc(drhoS(rho, phi, old, mid), [rho, phi, old, mid]); drhoreactionSource := optimize(drhoSproc); C(drhoreactionSource, filename = outstring, ansi);
fr3 := proc (rho, phi) options operator, arrow; diff(helmholtz(rho, phi), rho, rho, rho) end proc; fr1 := proc (rho, phi) options operator, arrow; diff(helmholtz(rho, phi), rho) end proc; CProc := makeproc(fr1(mid, phi)+(1/24)*fr3(mid, phi)*(rho-old), [rho, phi, old, mid]); CCProc := makeproc(CProc(rho, phi, old, .5*(rho+old)), [rho, phi, old]);
chemicalPotential := optimize(CCProc); C(chemicalPotential, filename = outstring, ansi);
drhoPotential := proc (rho, phi, old, mid) options operator, arrow; diff(CCProc(rho, phi, old), rho) end proc; drhomuproc := makeproc(drhoPotential(rho, phi, old), [rho, phi, old]); drhochemicalPotential := optimize(drhomuproc); C(drhochemicalPotential, filename = outstring, ansi);
# 
dphiPotential := proc (rho, phi, old, mid) options operator, arrow; diff(CCProc(rho, phi, old), phi) end proc; dphimuproc := makeproc(dphiPotential(rho, phi, old), [rho, phi, old]); dphichemicalPotential := optimize(dphimuproc); C(dphichemicalPotential, filename = outstring, ansi);
Pproc := makeproc(P2(rho, phi), [rho, phi]); pressure := optimize(Pproc); C(pressure, filename = outstring, ansi);
wavespeed := proc (rho, phi) options operator, arrow; sqrt(diff(Pressure(rho, phi), rho)) end proc; wproc := makeproc(simplify(wavespeed(rho, phi)), [rho, phi]); a := optimize(wproc); C(a, filename = outstring, ansi);
ml := proc (a) options operator, arrow; 0 end proc; mwpliq := makeproc(ml(x), x); mv := proc (a) options operator, arrow; 2 end proc; mwpvap := makeproc(mv(x), x); mwpvap := makeproc(mv(x), x); C(mwpliq, filename = outstring2, ansi); C(mwpvap, filename = outstring2, ansi);

solrho := proc (x) options operator, arrow; exp((-nn(x)*(bv-bl)-bl)/(nn(x)*(cv-cl)+cl)) end proc; solproc := makeproc(solrho(x), x); evalRho := optimize(solproc); C(evalRho, filename = outstring2, ansi);




