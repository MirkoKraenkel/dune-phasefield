
inline double helmholtz ( double rho ,double phi ) const
{
  double t1;
  double t10;
  double t11;
  double t12;
  double t13;
  double t14;
  double t15;
  double t16;
  double t19;
  double t2;
  double t3;
  double t34;
  double t35;
  double t36;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t1*phi;
    t10 = t2*phi;
    t11 = 6.0*t10;
    t12 = 15.0*t2;
    t13 = 10.0*t3;
    t14 = t11-t12+t13;
    t15 = log(rho);
    t16 = rho*t15;
    t19 = 1.0-t11+t12-t13;
    t34 = exp((0.4158883084E1*t10-0.1039720771E2*t2+0.6931471806E1*t3)/(-0.3E1*
t10+0.75E1*t2-0.5E1*t3+0.15E1));
    t35 = log(t34);
    t36 = t34*t35;
    return(2.0*A_*(t2-2.0*t3+t1)/delta_+t14*(-rho+t16+0.5)+t19*(-0.8068528194*rho
+0.15E1*t16)-beta_*(t14*(-t34+t36+0.5)+t19*(-0.8068528194*t34+0.15E1*t36)
-0.6931471806*t34+0.15E1));
  }
}


inline double reactionSource ( double rho ,double phi ,double old ) const
{
  double t101;
  double t105;
  double t106;
  double t107;
  double t114;
  double t116;
  double t117;
  double t118;
  double t12;
  double t133;
  double t147;
  double t15;
  double t150;
  double t155;
  double t156;
  double t167;
  double t172;
  double t180;
  double t184;
  double t186;
  double t19;
  double t20;
  double t203;
  double t21;
  double t22;
  double t24;
  double t27;
  double t29;
  double t3;
  double t33;
  double t37;
  double t38;
  double t4;
  double t40;
  double t41;
  double t42;
  double t43;
  double t45;
  double t46;
  double t47;
  double t48;
  double t5;
  double t52;
  double t54;
  double t55;
  double t56;
  double t60;
  double t62;
  double t67;
  double t69;
  double t70;
  double t71;
  double t72;
  double t74;
  double t87;
  double t89;
  double t91;
  double t92;
  double t93;
  double t99;
  {
    t3 = 0.5*phi+0.5*old;
    t4 = t3*t3;
    t5 = t4*t3;
    t12 = 1/delta_;
    t15 = t4*t4;
    t19 = 30.0*t15-60.0*t5+30.0*t4;
    t20 = log(rho);
    t21 = rho*t20;
    t22 = -rho+t21+0.5;
    t24 = -t19;
    t27 = -0.8068528194*rho+0.15E1*t21;
    t29 = t15*t3;
    t33 = 0.4158883084E1*t29-0.1039720771E2*t15+0.6931471806E1*t5;
    t37 = -0.3E1*t29+0.75E1*t15-0.5E1*t5+0.15E1;
    t38 = 1/t37;
    t40 = exp(t33*t38);
    t41 = log(t40);
    t42 = t40*t41;
    t43 = -t40+t42+0.5;
    t45 = 6.0*t29;
    t46 = 15.0*t15;
    t47 = 10.0*t5;
    t48 = t45-t46+t47;
    t52 = 0.2079441542E2*t15-0.4158883084E2*t5+0.2079441542E2*t4;
    t54 = t37*t37;
    t55 = 1/t54;
    t56 = t33*t55;
    t60 = -0.15E2*t15+0.3E2*t5-0.15E2*t4;
    t62 = t52*t38-t56*t60;
    t67 = -0.8068528194*t40+0.15E1*t42;
    t69 = 1.0-t45+t46-t47;
    t70 = t62*t40;
    t71 = 0.6931471806*t70;
    t72 = t70*t41;
    t74 = t71+0.15E1*t72;
    t87 = 360.0*t4-0.18E3*phi-0.18E3*old+60.0;
    t89 = -t87;
    t91 = t62*t62;
    t92 = t91*t62;
    t93 = t92*t40;
    t99 = 0.8317766168E2*t5-0.1247664925E3*t4+0.2079441542E2*phi+0.2079441542E2
*old;
    t101 = t52*t55;
    t105 = 1/t54/t37;
    t106 = t33*t105;
    t107 = t60*t60;
    t114 = -0.6E2*t5+0.9E2*t4-0.15E2*phi-0.15E2*old;
    t116 = t99*t38-2.0*t101*t60+2.0*t106*t107-t56*t114;
    t117 = t116*t62;
    t118 = t117*t40;
    t133 = t54*t54;
    t147 = (0.249532985E3*t4-0.1247664925E3*phi-0.1247664925E3*old+
0.4158883084E2)*t38-3.0*t99*t55*t60+6.0*t52*t105*t107-3.0*t101*t114-6.0*t33/
t133*t107*t60+6.0*t106*t60*t114-t56*(-0.18E3*t4+0.9E2*phi+0.9E2*old-0.3E2);
    t150 = t19*t91;
    t155 = t147*t40;
    t156 = 0.6931471806*t155;
    t167 = t48*t116;
    t172 = t48*t92;
    t180 = -120.0*t5+180.0*t4-0.3E2*phi-0.3E2*old;
    t184 = t116*t40;
    t186 = t91*t40;
    t203 = -0.6931471806*t93-0.2079441542E1*t118+t48*t147*t42+3.0*t150*t42+3.0*
t150*t40+t69*(t156+0.6579441543E1*t118+0.3693147181E1*t93+0.15E1*t155*t41+
0.45E1*t117*t42+0.15E1*t93*t41)-t156+3.0*t167*t72+3.0*t167*t70+2.0*t172*t40+t89
*t67+3.0*t180*t74+t172*t42+3.0*t24*(0.6931471806*t184+0.2193147181E1*t186+
0.15E1*t184*t41+0.15E1*t186*t41)+t87*t43-3.0*t180*t62*t42+3.0*t19*t116*t42;
    return(2.0*A_*(4.0*t5-6.0*t4+0.1E1*phi+0.1E1*old)*t12+t19*t22+t24*t27-beta_*(
t19*t43+t48*t62*t42+t24*t67+t69*t74-t71)+(2.0*A_*(0.12E2*phi+0.12E2*old-12.0)*
t12+t87*t22+t89*t27-beta_*t203)*(phi-old)/24.0);
  }
}


inline double dphireactionSource ( double rho ,double phi ,double old ) const
{
  double t10;
  double t101;
  double t103;
  double t106;
  double t107;
  double t108;
  double t113;
  double t115;
  double t118;
  double t120;
  double t122;
  double t124;
  double t126;
  double t127;
  double t128;
  double t129;
  double t13;
  double t130;
  double t132;
  double t136;
  double t139;
  double t145;
  double t147;
  double t150;
  double t153;
  double t157;
  double t161;
  double t162;
  double t163;
  double t165;
  double t166;
  double t167;
  double t168;
  double t171;
  double t172;
  double t173;
  double t175;
  double t178;
  double t18;
  double t183;
  double t184;
  double t185;
  double t186;
  double t189;
  double t19;
  double t195;
  double t197;
  double t198;
  double t199;
  double t20;
  double t201;
  double t206;
  double t209;
  double t21;
  double t216;
  double t219;
  double t222;
  double t226;
  double t23;
  double t232;
  double t234;
  double t235;
  double t236;
  double t242;
  double t245;
  double t247;
  double t251;
  double t26;
  double t28;
  double t29;
  double t3;
  double t31;
  double t312;
  double t313;
  double t316;
  double t317;
  double t319;
  double t320;
  double t322;
  double t33;
  double t336;
  double t338;
  double t341;
  double t342;
  double t343;
  double t344;
  double t346;
  double t354;
  double t356;
  double t358;
  double t363;
  double t37;
  double t370;
  double t378;
  double t379;
  double t38;
  double t383;
  double t4;
  double t40;
  double t401;
  double t403;
  double t409;
  double t41;
  double t414;
  double t42;
  double t420;
  double t428;
  double t43;
  double t436;
  double t439;
  double t440;
  double t447;
  double t450;
  double t457;
  double t473;
  double t48;
  double t51;
  double t516;
  double t53;
  double t54;
  double t55;
  double t59;
  double t61;
  double t67;
  double t69;
  double t71;
  double t73;
  double t75;
  double t78;
  double t79;
  double t80;
  double t81;
  double t85;
  double t87;
  double t89;
  double t92;
  double t93;
  double t94;
  {
    t3 = 0.5*phi+0.5*old;
    t4 = t3*t3;
    t10 = 1/delta_;
    t13 = t4*t3;
    t18 = 0.6E2*t13-0.9E2*t4+0.15E2*phi+0.15E2*old;
    t19 = log(rho);
    t20 = rho*t19;
    t21 = -rho+t20+0.5;
    t23 = -t18;
    t26 = -0.8068528194*rho+0.15E1*t20;
    t28 = t4*t4;
    t29 = t28*t3;
    t31 = 0.1039720771E2*t28;
    t33 = 0.4158883084E1*t29-t31+0.6931471806E1*t13;
    t37 = -0.3E1*t29+0.75E1*t28-0.5E1*t13+0.15E1;
    t38 = 1/t37;
    t40 = exp(t33*t38);
    t41 = log(t40);
    t42 = t40*t41;
    t43 = -t40+t42+0.5;
    t48 = 30.0*t28-60.0*t13+30.0*t4;
    t51 = t31-0.2079441542E2*t13+0.1039720771E2*t4;
    t53 = t37*t37;
    t54 = 1/t53;
    t55 = t33*t54;
    t59 = -0.75E1*t28+0.15E2*t13-0.75E1*t4;
    t61 = t51*t38-t55*t59;
    t67 = 0.15E2*t28-0.3E2*t13+0.15E2*t4;
    t69 = 0.4158883084E2*t13;
    t71 = 0.2079441542E2*t28-t69+0.2079441542E2*t4;
    t73 = -t67;
    t75 = t71*t38-t55*t73;
    t78 = 6.0*t29;
    t79 = 15.0*t28;
    t80 = 10.0*t13;
    t81 = t78-t79+t80;
    t85 = t69-0.6238324626E2*t4+0.1039720771E2*phi+0.1039720771E2*old;
    t87 = t71*t54;
    t89 = t51*t54;
    t92 = 1/t53/t37;
    t93 = t33*t92;
    t94 = t73*t59;
    t101 = -0.3E2*t13+0.45E2*t4-0.75E1*phi-0.75E1*old;
    t103 = t85*t38-t87*t59-t89*t73+2.0*t93*t94-t55*t101;
    t106 = t81*t75;
    t107 = t61*t40;
    t108 = t107*t41;
    t113 = -0.8068528194*t40+0.15E1*t42;
    t115 = -t48;
    t118 = 0.6931471806*t107+0.15E1*t108;
    t120 = t75*t40;
    t122 = t120*t41;
    t124 = 0.6931471806*t120+0.15E1*t122;
    t126 = 1.0-t78+t79-t80;
    t127 = t103*t40;
    t128 = 0.6931471806*t127;
    t129 = t75*t61;
    t130 = t129*t40;
    t132 = t127*t41;
    t136 = t128+0.2193147181E1*t130+0.15E1*t132+0.15E1*t129*t42;
    t139 = t18*t43+t48*t61*t42+t67*t75*t42+t81*t103*t42+t106*t108+t106*t107+t23
*t113+t115*t118+t73*t124+t126*t136-t128-0.6931471806*t130;
    t145 = -0.18E3+0.18E3*phi+0.18E3*old;
    t147 = -t145;
    t150 = 0.1247664925E3*t4;
    t153 = 0.8317766168E2*t13-t150+0.2079441542E2*phi+0.2079441542E2*old;
    t157 = t73*t73;
    t161 = t153*t38-2.0*t87*t73+2.0*t93*t157-t55*t23;
    t162 = t161*t103;
    t163 = t162*t40;
    t165 = t75*t75;
    t166 = t165*t75;
    t167 = t166*t61;
    t168 = t167*t40;
    t171 = 0.1247664925E3*phi;
    t172 = 0.1247664925E3*old;
    t173 = 0.4158883084E2+0.249532985E3*t4-t171-t172;
    t175 = t153*t54;
    t178 = t71*t92;
    t183 = t53*t53;
    t184 = 1/t183;
    t185 = t33*t184;
    t186 = t157*t73;
    t189 = t73*t23;
    t195 = -0.18E3*t4+0.9E2*phi+0.9E2*old-0.3E2;
    t197 = t173*t38-3.0*t175*t73+6.0*t178*t157-3.0*t87*t23-6.0*t185*t186+6.0*
t93*t189-t55*t195;
    t198 = t197*t61;
    t199 = t198*t40;
    t201 = t48*t75;
    t206 = 0.2079441542E2+t150-0.6238324625E2*phi-0.6238324625E2*old;
    t209 = t85*t54;
    t216 = t51*t92;
    t219 = t157*t59;
    t222 = t73*t101;
    t226 = t23*t59;
    t232 = -0.15E2-0.9E2*t4+0.45E2*phi+0.45E2*old;
    t234 = t206*t38-t175*t59-2.0*t209*t73+4.0*t178*t94-2.0*t87*t101+2.0*t216*
t157-6.0*t185*t219+4.0*t93*t222-t89*t23+2.0*t93*t226-t55*t232;
    t235 = t234*t75;
    t236 = t235*t40;
    t242 = -120.0*t13+180.0*t4-0.3E2*phi-0.3E2*old;
    t245 = t81*t197;
    t247 = t81*t165;
    t251 = t18*t165;
    t312 = -3.0*t87*t232-6.0*t51*t184*t186+24.0*t33/t183/t37*t186*t59-18.0*t185
*t157*t101+6.0*t216*t189-18.0*t185*t189*t59+6.0*t93*t101*t23+6.0*t93*t73*t232-
t89*t195+2.0*t93*t195*t59-t55*(0.9E2-0.9E2*phi-0.9E2*old);
    t313 = (-0.1247664925E3+t171+t172)*t38-t173*t54*t59-3.0*t206*t54*t73+6.0*
t153*t92*t94-3.0*t175*t101+6.0*t85*t92*t157-18.0*t71*t184*t219+12.0*t178*t222
-3.0*t209*t23+6.0*t178*t226+t312;
    t316 = t165*t40;
    t317 = t316*t103;
    t319 = -0.2079441542E1*t163-0.6931471806*t168-0.6931471806*t199+6.0*t201*
t132-0.2079441542E1*t236+3.0*t242*t136+t245*t108+6.0*t247*t127+t245*t107+3.0*
t251*t42+t81*t313*t42-0.2079441542E1*t317;
    t320 = t81*t166;
    t322 = t81*t161;
    t336 = t67*t166;
    t338 = t81*t234;
    t341 = t197*t40;
    t342 = 0.6931471806*t341;
    t343 = t161*t75;
    t344 = t343*t40;
    t346 = t166*t40;
    t354 = t342+0.6579441543E1*t344+0.3693147181E1*t346+0.15E1*t341*t41+0.45E1*
t343*t42+0.15E1*t346*t41;
    t356 = t234*t40;
    t358 = t161*t61;
    t363 = t165*t61;
    t370 = t41*t103;
    t378 = t313*t40;
    t379 = 0.6931471806*t378;
    t383 = t343*t107;
    t401 = t379+0.2193147181E1*t199+0.6579441543E1*t236+0.6579441543E1*t163+
0.1107944154E2*t383+0.1107944154E2*t317+0.5193147181E1*t168+0.15E1*t378*t41+
0.15E1*t198*t42+0.45E1*t235*t42+0.45E1*t162*t42+0.45E1*t343*t108+0.45E1*t316*
t370+0.15E1*t167*t42;
    t403 = t320*t108+3.0*t322*t75*t108+3.0*t247*t132+3.0*t251*t40+3.0*t195*t124
+3.0*t322*t132+3.0*t320*t107+t336*t42+3.0*t338*t120+t73*t354+3.0*t115*(
0.6931471806*t356+0.2193147181E1*t358*t40+0.4386294362E1*t120*t103+
0.3693147181E1*t363*t40+0.15E1*t356*t41+0.15E1*t358*t42+0.3E1*t120*t370+0.15E1*
t363*t42)+t126*t401;
    t409 = -t242;
    t414 = t48*t161;
    t420 = 360.0*t4-0.18E3*phi-0.18E3*old+60.0;
    t428 = t409*t75;
    t436 = t48*t165;
    t439 = -3.0*t195*t75*t42+3.0*t409*t103*t42+t147*t113+3.0*t414*t107+t420*t61
*t42+3.0*t414*t108+6.0*t201*t127-0.2079441542E1*t383+3.0*t428*t107+3.0*t18*t161
*t42+3.0*t322*t127+6.0*t436*t107;
    t440 = t161*t40;
    t447 = 0.6931471806*t440+0.2193147181E1*t316+0.15E1*t440*t41+0.15E1*t316*
t41;
    t450 = t67*t161;
    t457 = -t420;
    t473 = 3.0*t23*t447+3.0*t450*t120+2.0*t336*t40-t379+t67*t197*t42+t457*t118+
t145*t43+3.0*t436*t108+3.0*t48*t234*t42+3.0*t450*t122+6.0*t322*t130+3.0*t428*
t108+3.0*t338*t122;
    t516 = -0.6931471806*t346-0.2079441542E1*t344+t245*t42+3.0*t436*t42+3.0*
t436*t40+t126*t354-t342+3.0*t322*t122+3.0*t322*t120+2.0*t320*t40+t457*t113+3.0*
t242*t124+t320*t42+3.0*t115*t447+t420*t43+3.0*t428*t42+3.0*t414*t42;
    return(2.0*A_*(0.6E1*t4-0.3E1*phi-0.3E1*old+0.1E1)*t10+t18*t21+t23*t26-beta_*
t139+(0.24E2*A_*t10+t145*t21+t147*t26-beta_*(t319+t403+t439+t473))*(phi-old)/24.0
+A_*(0.12E2*phi+0.12E2*old-12.0)*t10/12.0+t420*t21/24.0+t457*t26/24.0-beta_*t516/
24.0);
  }
}


inline double chemicalPotential ( double rho ,double phi ,double old ) const
{
  double t1;
  double t11;
  double t12;
  double t14;
  double t18;
  double t19;
  double t2;
  double t4;
  double t5;
  double t7;
  double t8;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t8 = t4-t5+t7;
    t11 = 0.5*rho+0.5*old;
    t12 = log(t11);
    t14 = 1.0-t4+t5-t7;
    t18 = t11*t11;
    t19 = 1/t18;
    return(t8*t12+t14*(0.6931471806+0.15E1*t12)+(-t8*t19-0.15E1*t14*t19)*(rho-
old)/24.0);
  }
}

inline double drhochemicalPotential ( double rho ,double phi ,double old ) const
{
  double t1;
  double t11;
  double t12;
  double t15;
  double t18;
  double t2;
  double t20;
  double t29;
  double t4;
  double t5;
  double t7;
  double t8;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t8 = t4-t5+t7;
    t11 = 0.5*rho+0.5*old;
    t12 = 1/t11;
    t15 = 1.0-t4+t5-t7;
    t18 = t11*t11;
    t20 = 1/t18/t11;
    t29 = 1/t18;
    return(0.5*t8*t12+0.75*t15*t12+(0.1E1*t8*t20+0.15E1*t15*t20)*(rho-old)/24.0
-t8*t29/24.0-0.625E-1*t15*t29);
  }
}


inline double dphichemicalPotential ( double rho ,double phi ,double old ) const
{
  double t1;
  double t10;
  double t11;
  double t13;
  double t17;
  double t18;
  double t2;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t7 = 30.0*t2-60.0*t1*phi+30.0*t1;
    t10 = 0.5*rho+0.5*old;
    t11 = log(t10);
    t13 = -t7;
    t17 = t10*t10;
    t18 = 1/t17;
    return(t7*t11+t13*(0.6931471806+0.15E1*t11)+(-t7*t18-0.15E1*t13*t18)*(rho-
old)/24.0);
  }
}


inline double pressure ( double rho ,double phi ) const
{
  double t1;
  double t13;
  double t2;
  double t4;
  double t5;
  double t6;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t6 = t1*phi;
    t7 = 10.0*t6;
    t13 = log(rho);
    return((t4-t5+t7)*(rho-0.5)+(1.0-t4+t5-t7)*(0.8068528194*rho-0.15E1*rho*t13
+rho*(0.6931471806+0.15E1*t13))-2.0*A_*(t2-2.0*t6+t1)/delta_);
  }
}


inline double a ( double rho ,double phi ) const
{
  double t1;
  double t2;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t9 = sqrt(-12.0*t2*phi+30.0*t2-20.0*t1*phi+6.0);
    return(0.5*t9);
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t15;
  double t16;
  double t17;
  double t18;
  double t19;
  double t21;
  double t22;
  double t26;
  double t3;
  double t4;
  double t41;
  double t43;
  double t44;
  double t7;
  {
    t3 = phi*phi;
    t4 = t3*t3;
    t7 = t3*phi;
    t15 = t4*phi;
    t16 = 6.0*t15;
    t17 = 15.0*t4;
    t18 = 10.0*t7;
    t19 = t16-t17+t18;
    t21 = log(rho);
    t22 = rho*t21;
    t26 = 1.0-t16+t17-t18;
    t41 = exp((-0.2889626582E2+0.1204277045E3*t15-0.3010692614E3*t4+
0.2007128409E3*t7)/(0.9502053219E1+0.116654747E3*t15-0.2916368674E3*t4+
0.1944245783E3*t7));
    t43 = log(t41);
    t44 = t41*t43;
    return(2.0/delta*A*(t4+(2.0*alpha-2.0)*t7+(-3.0*alpha+1.0)*t3+alpha)+t19*(
-0.2011952932E2*rho+0.289445110522E2*t22+0.2133795654E2)+t26*(0.193942126E2*rho
+0.950205321859E1*t22+0.4540504209)-beta*(t19*(0.2133795654E2-0.2011952932E2*
t41+0.289445110522E2*t44)+t26*(0.4540504209+0.193942126E2*t41+0.950205321859E1*
t44)));
  }
}

#include <math.h>
double reactionSource(double rho,double phi,double old)
{
  double t100;
  double t107;
  double t114;
  double t116;
  double t120;
  double t121;
  double t122;
  double t129;
  double t131;
  double t132;
  double t134;
  double t135;
  double t137;
  double t139;
  double t157;
  double t172;
  double t174;
  double t175;
  double t178;
  double t180;
  double t182;
  double t184;
  double t2;
  double t20;
  double t24;
  double t26;
  double t27;
  double t29;
  double t31;
  double t34;
  double t36;
  double t40;
  double t44;
  double t45;
  double t47;
  double t49;
  double t5;
  double t50;
  double t52;
  double t54;
  double t55;
  double t56;
  double t57;
  double t6;
  double t61;
  double t63;
  double t64;
  double t65;
  double t69;
  double t7;
  double t71;
  double t72;
  double t74;
  double t76;
  double t80;
  double t82;
  double t85;
  double t98;
  {
    t2 = 1/delta*A;
    t5 = 0.5*phi+0.5*old;
    t6 = t5*t5;
    t7 = t6*t5;
    t20 = t6*t6;
    t24 = 30.0*t20-60.0*t7+30.0*t6;
    t26 = log(rho);
    t27 = rho*t26;
    t29 = -0.2011952932E2*rho+0.289445110522E2*t27+0.2133795654E2;
    t31 = -t24;
    t34 = 0.193942126E2*rho+0.950205321859E1*t27+0.4540504209;
    t36 = t20*t5;
    t40 = -0.2889626582E2+0.1204277045E3*t36-0.3010692614E3*t20+0.2007128409E3*
t7;
    t44 = 0.9502053219E1+0.116654747E3*t36-0.2916368674E3*t20+0.1944245783E3*t7
;
    t45 = 1/t44;
    t47 = exp(t40*t45);
    t49 = log(t47);
    t50 = t47*t49;
    t52 = 0.2133795654E2-0.2011952932E2*t47+0.289445110522E2*t50;
    t54 = 6.0*t36;
    t55 = 15.0*t20;
    t56 = 10.0*t7;
    t57 = t54-t55+t56;
    t61 = 0.6021385225E3*t20-0.1204277046E4*t7+0.6021385227E3*t6;
    t63 = t44*t44;
    t64 = 1/t63;
    t65 = t40*t64;
    t69 = 0.583273735E3*t20-0.116654747E4*t7+0.5832737349E3*t6;
    t71 = t61*t45-t65*t69;
    t72 = t71*t47;
    t74 = t72*t49;
    t76 = 0.882498173E1*t72+0.289445110522E2*t74;
    t80 = 0.4540504209+0.193942126E2*t47+0.950205321859E1*t50;
    t82 = 1.0-t54+t55-t56;
    t85 = 0.2889626582E2*t72+0.950205321859E1*t74;
    t98 = 360.0*t6-0.18E3*phi-0.18E3*old+60.0;
    t100 = -t98;
    t107 = 120.0*t7-180.0*t6+0.3E2*phi+0.3E2*old;
    t114 = 0.240855409E4*t7-0.3612831138E4*t6+0.6021385225E3*phi+0.6021385225E3
*old;
    t116 = t61*t64;
    t120 = 1/t63/t44;
    t121 = t40*t120;
    t122 = t69*t69;
    t129 = 0.233309494E4*t7-0.349964241E4*t6+0.583273735E3*phi+0.583273735E3*
old;
    t131 = t114*t45-2.0*t116*t69+2.0*t121*t122-t65*t129;
    t132 = t131*t47;
    t134 = t71*t71;
    t135 = t134*t47;
    t137 = t132*t49;
    t139 = t135*t49;
    t157 = t63*t63;
    t172 = ((0.722566227E4*t6-0.3612831138E4*phi-0.3612831138E4*old+
0.1204277045E4)*t45-3.0*t114*t64*t69+6.0*t61*t120*t122-3.0*t116*t129-6.0*t40/
t157*t122*t69+6.0*t121*t69*t129-t65*(0.699928482E4*t6-0.349964241E4*phi
-0.349964241E4*old+0.116654747E4))*t47;
    t174 = t131*t71;
    t175 = t174*t47;
    t178 = t134*t71*t47;
    t180 = t172*t49;
    t182 = t174*t50;
    t184 = t178*t49;
    return(2.0*t2*(4.0*t7+3.0*(2.0*alpha-2.0)*t6+2.0*(-3.0*alpha+1.0)*t5)+t24*
t29+t31*t34-beta*(t24*t52+t57*t76+t31*t80+t82*t85)+(2.0*t2*(0.12E2*phi+0.12E2*
old+12.0*alpha-12.0)+t98*t29+t100*t34-beta*(t98*t52+3.0*t107*t76+3.0*t24*(
0.882498173E1*t132+0.3776949278E2*t135+0.289445110522E2*t137+0.289445110522E2*
t139)+t57*(0.882498173E1*t172+0.1133084783E3*t175+0.6671400383E2*t178+
0.289445110522E2*t180+0.8683353315E2*t182+0.289445110522E2*t184)+t100*t80-3.0*
t107*t85+3.0*t31*(0.2889626582E2*t132+0.3839831904E2*t135+0.950205321859E1*t137
+0.950205321859E1*t139)+t82*(0.2889626582E2*t172+0.1151949571E3*t175+
0.4790037226E2*t178+0.950205321859E1*t180+0.2850615966E2*t182+0.950205321859E1*
t184)))*(phi-old)/24.0);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi,double old,double mid)
{
  double t100;
  double t105;
  double t107;
  double t109;
  double t112;
  double t113;
  double t114;
  double t120;
  double t122;
  double t123;
  double t125;
  double t126;
  double t128;
  double t130;
  double t132;
  double t136;
  double t138;
  double t141;
  double t143;
  double t146;
  double t148;
  double t153;
  double t16;
  double t160;
  double t162;
  double t168;
  double t173;
  double t180;
  double t187;
  double t191;
  double t195;
  double t198;
  double t2;
  double t200;
  double t201;
  double t203;
  double t204;
  double t206;
  double t208;
  double t21;
  double t210;
  double t216;
  double t218;
  double t220;
  double t223;
  double t228;
  double t23;
  double t231;
  double t232;
  double t233;
  double t234;
  double t237;
  double t24;
  double t241;
  double t246;
  double t248;
  double t249;
  double t251;
  double t252;
  double t254;
  double t256;
  double t257;
  double t259;
  double t26;
  double t261;
  double t263;
  double t264;
  double t266;
  double t274;
  double t28;
  double t282;
  double t285;
  double t289;
  double t290;
  double t291;
  double t293;
  double t294;
  double t296;
  double t297;
  double t299;
  double t300;
  double t302;
  double t304;
  double t306;
  double t308;
  double t31;
  double t33;
  double t34;
  double t368;
  double t370;
  double t372;
  double t373;
  double t375;
  double t376;
  double t378;
  double t379;
  double t38;
  double t381;
  double t383;
  double t385;
  double t386;
  double t388;
  double t390;
  double t392;
  double t394;
  double t396;
  double t398;
  double t400;
  double t402;
  double t405;
  double t410;
  double t417;
  double t42;
  double t43;
  double t437;
  double t45;
  double t453;
  double t455;
  double t47;
  double t48;
  double t5;
  double t50;
  double t55;
  double t59;
  double t6;
  double t61;
  double t62;
  double t63;
  double t67;
  double t69;
  double t70;
  double t72;
  double t74;
  double t79;
  double t83;
  double t86;
  double t88;
  double t90;
  double t91;
  double t93;
  double t95;
  double t97;
  double t98;
  double t99;
  {
    t2 = 1/delta*A;
    t5 = 0.5*phi+0.5*old;
    t6 = t5*t5;
    t16 = t6*t5;
    t21 = 0.6E2*t16-0.9E2*t6+0.15E2*phi+0.15E2*old;
    t23 = log(rho);
    t24 = rho*t23;
    t26 = -0.2011952932E2*rho+0.289445110522E2*t24+0.2133795654E2;
    t28 = -t21;
    t31 = 0.193942126E2*rho+0.950205321859E1*t24+0.4540504209;
    t33 = t6*t6;
    t34 = t33*t5;
    t38 = -0.2889626582E2+0.1204277045E3*t34-0.3010692614E3*t33+0.2007128409E3*
t16;
    t42 = 0.9502053219E1+0.116654747E3*t34-0.2916368674E3*t33+0.1944245783E3*
t16;
    t43 = 1/t42;
    t45 = exp(t38*t43);
    t47 = log(t45);
    t48 = t45*t47;
    t50 = 0.2133795654E2-0.2011952932E2*t45+0.289445110522E2*t48;
    t55 = 30.0*t33-60.0*t16+30.0*t6;
    t59 = 0.3010692612E3*t33-0.6021385228E3*t16+0.3010692614E3*t6;
    t61 = t42*t42;
    t62 = 1/t61;
    t63 = t38*t62;
    t67 = 0.2916368675E3*t33-0.5832737348E3*t16+0.2916368674E3*t6;
    t69 = t59*t43-t63*t67;
    t70 = t69*t45;
    t72 = t70*t47;
    t74 = 0.882498173E1*t70+0.289445110522E2*t72;
    t79 = 0.15E2*t33-0.3E2*t16+0.15E2*t6;
    t83 = 0.6021385225E3*t33-0.1204277046E4*t16+0.6021385227E3*t6;
    t86 = 0.116654747E4*t16;
    t88 = 0.583273735E3*t33-t86+0.5832737349E3*t6;
    t90 = t83*t43-t63*t88;
    t91 = t90*t45;
    t93 = t91*t47;
    t95 = 0.882498173E1*t91+0.289445110522E2*t93;
    t97 = 6.0*t34;
    t98 = 15.0*t33;
    t99 = 10.0*t16;
    t100 = t97-t98+t99;
    t105 = 0.1204277045E4*t16-0.1806415569E4*t6+0.3010692614E3*phi+
0.3010692614E3*old;
    t107 = t83*t62;
    t109 = t59*t62;
    t112 = 1/t61/t42;
    t113 = t38*t112;
    t114 = t88*t67;
    t120 = t86-0.1749821205E4*t6+0.2916368674E3*phi+0.2916368674E3*old;
    t122 = t105*t43-t107*t67-t109*t88+2.0*t113*t114-t63*t120;
    t123 = t122*t45;
    t125 = t90*t69;
    t126 = t125*t45;
    t128 = t123*t47;
    t130 = t125*t48;
    t132 = 0.882498173E1*t123+0.3776949278E2*t126+0.289445110522E2*t128+
0.289445110522E2*t130;
    t136 = 0.4540504209+0.193942126E2*t45+0.950205321859E1*t48;
    t138 = -t55;
    t141 = 0.2889626582E2*t70+0.950205321859E1*t72;
    t143 = -t79;
    t146 = 0.2889626582E2*t91+0.950205321859E1*t93;
    t148 = 1.0-t97+t98-t99;
    t153 = 0.2889626582E2*t123+0.3839831904E2*t126+0.950205321859E1*t128+
0.950205321859E1*t130;
    t160 = 0.18E3*phi+0.18E3*old-0.18E3;
    t162 = -t160;
    t168 = 360.0*t6-0.18E3*phi-0.18E3*old+60.0;
    t173 = 0.18E3*t6-0.9E2*phi-0.9E2*old+0.3E2;
    t180 = 120.0*t16-180.0*t6+0.3E2*phi+0.3E2*old;
    t187 = 0.240855409E4*t16-0.3612831138E4*t6+0.6021385225E3*phi+
0.6021385225E3*old;
    t191 = t88*t88;
    t195 = 0.349964241E4*t6;
    t198 = 0.233309494E4*t16-t195+0.583273735E3*phi+0.583273735E3*old;
    t200 = t187*t43-2.0*t107*t88+2.0*t113*t191-t63*t198;
    t201 = t200*t45;
    t203 = t90*t90;
    t204 = t203*t45;
    t206 = t201*t47;
    t208 = t204*t47;
    t210 = 0.882498173E1*t201+0.3776949278E2*t204+0.289445110522E2*t206+
0.289445110522E2*t208;
    t216 = 0.3612831135E4*t6-0.1806415569E4*phi-0.1806415569E4*old+
0.6021385225E3;
    t218 = t187*t62;
    t220 = t105*t62;
    t223 = t83*t112;
    t228 = t59*t112;
    t231 = t61*t61;
    t232 = 1/t231;
    t233 = t38*t232;
    t234 = t191*t67;
    t237 = t88*t120;
    t241 = t198*t67;
    t246 = 0.583273735E3+t195-0.1749821205E4*phi-0.1749821205E4*old;
    t248 = t216*t43-t218*t67-2.0*t220*t88+4.0*t223*t114-2.0*t107*t120+2.0*t228*
t191-6.0*t233*t234+4.0*t113*t237-t109*t198+2.0*t113*t241-t63*t246;
    t249 = t248*t45;
    t251 = t200*t69;
    t252 = t251*t45;
    t254 = t91*t122;
    t256 = t203*t69;
    t257 = t256*t45;
    t259 = t249*t47;
    t261 = t251*t48;
    t263 = t47*t122;
    t264 = t91*t263;
    t266 = t256*t48;
    t274 = 0.722566227E4*t6-0.3612831138E4*phi-0.3612831138E4*old+
0.1204277045E4;
    t282 = t191*t88;
    t285 = t88*t198;
    t289 = 0.349964241E4*phi;
    t290 = 0.349964241E4*old;
    t291 = 0.116654747E4+0.699928482E4*t6-t289-t290;
    t293 = t274*t43-3.0*t218*t88+6.0*t223*t191-3.0*t107*t198-6.0*t233*t282+6.0*
t113*t285-t63*t291;
    t294 = t293*t45;
    t296 = t200*t90;
    t297 = t296*t45;
    t299 = t203*t90;
    t300 = t299*t45;
    t302 = t294*t47;
    t304 = t296*t48;
    t306 = t300*t47;
    t308 = 0.882498173E1*t294+0.1133084783E3*t297+0.6671400383E2*t300+
0.289445110522E2*t302+0.8683353315E2*t304+0.289445110522E2*t306;
    t368 = -3.0*t107*t246-6.0*t59*t232*t282+24.0*t38/t231/t42*t282*t67-18.0*
t233*t191*t120+6.0*t228*t285-18.0*t233*t285*t67+6.0*t113*t120*t198+6.0*t113*t88
*t246-t109*t291+2.0*t113*t291*t67-t63*(-0.349964241E4+t289+t290);
    t370 = ((-0.3612831138E4+0.3612831135E4*phi+0.3612831135E4*old)*t43-t274*
t62*t67-3.0*t216*t62*t88+6.0*t187*t112*t114-3.0*t218*t120+6.0*t105*t112*t191
-18.0*t83*t232*t234+12.0*t223*t237-3.0*t220*t198+6.0*t223*t241+t368)*t45;
    t372 = t293*t69;
    t373 = t372*t45;
    t375 = t248*t90;
    t376 = t375*t45;
    t378 = t200*t122;
    t379 = t378*t45;
    t381 = t296*t70;
    t383 = t204*t122;
    t385 = t299*t69;
    t386 = t385*t45;
    t388 = t370*t47;
    t390 = t372*t48;
    t392 = t375*t48;
    t394 = t378*t48;
    t396 = t296*t72;
    t398 = t204*t263;
    t400 = t385*t48;
    t402 = 0.882498173E1*t370+0.3776949278E2*t373+0.1133084783E3*t376+
0.1133084783E3*t379+0.2001420114E3*t381+0.2001420115E3*t383+0.9565851488E2*t386
+0.289445110522E2*t388+0.289445110522E2*t390+0.8683353315E2*t392+0.8683353315E2
*t394+0.8683353315E2*t396+0.8683353315E2*t398+0.289445110522E2*t400;
    t405 = -t168;
    t410 = -t180;
    t417 = 0.2889626582E2*t201+0.3839831904E2*t204+0.950205321859E1*t206+
0.950205321859E1*t208;
    t437 = 0.2889626582E2*t294+0.1151949571E3*t297+0.4790037226E2*t300+
0.950205321859E1*t302+0.2850615966E2*t304+0.950205321859E1*t306;
    t453 = 0.2889626582E2*t370+0.3839831904E2*t373+0.1151949571E3*t376+
0.1151949571E3*t379+0.1437011168E3*t381+0.1437011168E3*t383+0.5740242548E2*t386
+0.950205321859E1*t388+0.950205321859E1*t390+0.2850615966E2*t392+0.2850615966E2
*t394+0.2850615966E2*t396+0.2850615966E2*t398+0.950205321859E1*t400;
    t455 = t160*t50+t168*t74+3.0*t173*t95+3.0*t180*t132+3.0*t21*t210+3.0*t55*(
0.882498173E1*t249+0.3776949278E2*t252+0.7553898556E2*t254+0.6671400383E2*t257+
0.289445110522E2*t259+0.289445110522E2*t261+0.578890221E2*t264+0.289445110522E2
*t266)+t79*t308+t100*t402+t162*t136+t405*t141-3.0*t173*t146+3.0*t410*t153+3.0*
t28*t417+3.0*t138*(0.2889626582E2*t249+0.3839831904E2*t252+0.7679663808E2*t254+
0.4790037226E2*t257+0.950205321859E1*t259+0.950205321859E1*t261+0.1900410644E2*
t264+0.950205321859E1*t266)+t143*t437+t148*t453;
    return(2.0*t2*(0.6E1*t6+0.3E1*(2.0*alpha-2.0)*t5-0.3E1*alpha+0.1E1)+t21*t26
+t28*t31-beta*(t21*t50+t55*t74+t79*t95+t100*t132+t28*t136+t138*t141+t143*t146+
t148*t153)+(0.24E2*t2+t160*t26+t162*t31-beta*t455)*(phi-old)/24.0+t2*(0.12E2*
phi+0.12E2*old+12.0*alpha-12.0)/12.0+t168*t26/24.0+t405*t31/24.0-beta*(t168*t50
+3.0*t180*t95+3.0*t55*t210+t100*t308+t405*t136+3.0*t410*t146+3.0*t138*t417+t148
*t437)/24.0);
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t11;
  double t12;
  double t16;
  double t2;
  double t20;
  double t21;
  double t4;
  double t5;
  double t7;
  double t8;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t8 = t4-t5+t7;
    t11 = 0.5*rho+0.5*old;
    t12 = log(t11);
    t16 = 1.0-t4+t5-t7;
    t20 = t11*t11;
    t21 = 1/t20;
    return(t8*(0.882498173E1+0.289445110522E2*t12)+t16*(0.2889626582E2+
0.950205321859E1*t12)+(-0.2894451105E2*t8*t21-0.9502053219E1*t16*t21)*(rho-old)
/24.0);
  }
}

double drhochemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t11;
  double t12;
  double t15;
  double t18;
  double t2;
  double t20;
  double t29;
  double t4;
  double t5;
  double t7;
  double t8;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t8 = t4-t5+t7;
    t11 = 0.5*rho+0.5*old;
    t12 = 1/t11;
    t15 = 1.0-t4+t5-t7;
    t18 = t11*t11;
    t20 = 1/t18/t11;
    t29 = 1/t18;
    return(0.1447225552E2*t8*t12+0.475102661E1*t15*t12+(0.2894451105E2*t8*t20+
0.9502053219E1*t15*t20)*(rho-old)/24.0-0.1206021294E1*t8*t29-0.3959188841*t15*
t29);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t10;
  double t11;
  double t15;
  double t19;
  double t2;
  double t20;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t7 = 30.0*t2-60.0*t1*phi+30.0*t1;
    t10 = 0.5*rho+0.5*old;
    t11 = log(t10);
    t15 = -t7;
    t19 = t10*t10;
    t20 = 1/t19;
    return(t7*(0.882498173E1+0.289445110522E2*t11)+t15*(0.2889626582E2+
0.950205321859E1*t11)+(-0.2894451105E2*t7*t20-0.9502053219E1*t15*t20)*(rho-old)
/24.0);
  }
}

#include <math.h>
double pressure(double rho,double phi)
{
  double t1;
  double t10;
  double t11;
  double t2;
  double t4;
  double t5;
  double t6;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t6 = t1*phi;
    t7 = 10.0*t6;
    t10 = log(rho);
    t11 = rho*t10;
    return((t4-t5+t7)*(-0.2094881058E2+0.4291934897E2*rho-0.2894451105E2*t11+
rho*(-0.1397483792E2+0.289445110522E2*t10))+(1.0-t4+t5-t7)*(-0.6490446091E-1+
0.3405607049E1*rho-0.9502053219E1*t11+rho*(0.609644617E1+0.950205321859E1*t10))
-2.0/delta*A*(t2+(2.0*alpha-2.0)*t6+(-3.0*alpha+1.0)*t1+alpha));
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t1;
  double t2;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t9 = sqrt(0.116654747E11*t2*phi-0.2916368675E11*t2+0.1944245783E11*t1*phi+
950205322.0);
    return(0.1E-3*t9);
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t15;
  double t16;
  double t17;
  double t18;
  double t19;
  double t21;
  double t22;
  double t26;
  double t3;
  double t4;
  double t41;
  double t43;
  double t44;
  double t7;
  {
    t3 = phi*phi;
    t4 = t3*t3;
    t7 = t3*phi;
    t15 = t4*phi;
    t16 = 6.0*t15;
    t17 = 15.0*t4;
    t18 = 10.0*t7;
    t19 = t16-t17+t18;
    t21 = log(rho);
    t22 = rho*t21;
    t26 = 1.0-t16+t17-t18;
    t41 = exp((-0.2889626582E2+0.1204277045E3*t15-0.3010692614E3*t4+
0.2007128409E3*t7)/(0.9502053219E1+0.116654747E3*t15-0.2916368674E3*t4+
0.1944245783E3*t7));
    t43 = log(t41);
    t44 = t41*t43;
    return(2.0/delta*A*(t4+(2.0*alpha-2.0)*t7+(-3.0*alpha+1.0)*t3+alpha)+t19*(
-0.2011952932E2*rho+0.289445110522E2*t22+0.2133795654E2)+t26*(0.193942126E2*rho
+0.950205321859E1*t22+0.4540504209)-beta*(t19*(0.2133795654E2-0.2011952932E2*
t41+0.289445110522E2*t44)+t26*(0.4540504209+0.193942126E2*t41+0.950205321859E1*
t44)));
  }
}

#include <math.h>
double reactionSource(double rho,double phi,double old)
{
  double t100;
  double t107;
  double t114;
  double t116;
  double t120;
  double t121;
  double t122;
  double t129;
  double t131;
  double t132;
  double t134;
  double t135;
  double t137;
  double t139;
  double t157;
  double t172;
  double t174;
  double t175;
  double t178;
  double t180;
  double t182;
  double t184;
  double t2;
  double t20;
  double t24;
  double t26;
  double t27;
  double t29;
  double t31;
  double t34;
  double t36;
  double t40;
  double t44;
  double t45;
  double t47;
  double t49;
  double t5;
  double t50;
  double t52;
  double t54;
  double t55;
  double t56;
  double t57;
  double t6;
  double t61;
  double t63;
  double t64;
  double t65;
  double t69;
  double t7;
  double t71;
  double t72;
  double t74;
  double t76;
  double t80;
  double t82;
  double t85;
  double t98;
  {
    t2 = 1/delta*A;
    t5 = 0.5*phi+0.5*old;
    t6 = t5*t5;
    t7 = t6*t5;
    t20 = t6*t6;
    t24 = 30.0*t20-60.0*t7+30.0*t6;
    t26 = log(rho);
    t27 = rho*t26;
    t29 = -0.2011952932E2*rho+0.289445110522E2*t27+0.2133795654E2;
    t31 = -t24;
    t34 = 0.193942126E2*rho+0.950205321859E1*t27+0.4540504209;
    t36 = t20*t5;
    t40 = -0.2889626582E2+0.1204277045E3*t36-0.3010692614E3*t20+0.2007128409E3*
t7;
    t44 = 0.9502053219E1+0.116654747E3*t36-0.2916368674E3*t20+0.1944245783E3*t7
;
    t45 = 1/t44;
    t47 = exp(t40*t45);
    t49 = log(t47);
    t50 = t47*t49;
    t52 = 0.2133795654E2-0.2011952932E2*t47+0.289445110522E2*t50;
    t54 = 6.0*t36;
    t55 = 15.0*t20;
    t56 = 10.0*t7;
    t57 = t54-t55+t56;
    t61 = 0.6021385225E3*t20-0.1204277046E4*t7+0.6021385227E3*t6;
    t63 = t44*t44;
    t64 = 1/t63;
    t65 = t40*t64;
    t69 = 0.583273735E3*t20-0.116654747E4*t7+0.5832737349E3*t6;
    t71 = t61*t45-t65*t69;
    t72 = t71*t47;
    t74 = t72*t49;
    t76 = 0.882498173E1*t72+0.289445110522E2*t74;
    t80 = 0.4540504209+0.193942126E2*t47+0.950205321859E1*t50;
    t82 = 1.0-t54+t55-t56;
    t85 = 0.2889626582E2*t72+0.950205321859E1*t74;
    t98 = 360.0*t6-0.18E3*phi-0.18E3*old+60.0;
    t100 = -t98;
    t107 = 120.0*t7-180.0*t6+0.3E2*phi+0.3E2*old;
    t114 = 0.240855409E4*t7-0.3612831138E4*t6+0.6021385225E3*phi+0.6021385225E3
*old;
    t116 = t61*t64;
    t120 = 1/t63/t44;
    t121 = t40*t120;
    t122 = t69*t69;
    t129 = 0.233309494E4*t7-0.349964241E4*t6+0.583273735E3*phi+0.583273735E3*
old;
    t131 = t114*t45-2.0*t116*t69+2.0*t121*t122-t65*t129;
    t132 = t131*t47;
    t134 = t71*t71;
    t135 = t134*t47;
    t137 = t132*t49;
    t139 = t135*t49;
    t157 = t63*t63;
    t172 = ((0.722566227E4*t6-0.3612831138E4*phi-0.3612831138E4*old+
0.1204277045E4)*t45-3.0*t114*t64*t69+6.0*t61*t120*t122-3.0*t116*t129-6.0*t40/
t157*t122*t69+6.0*t121*t69*t129-t65*(0.699928482E4*t6-0.349964241E4*phi
-0.349964241E4*old+0.116654747E4))*t47;
    t174 = t131*t71;
    t175 = t174*t47;
    t178 = t134*t71*t47;
    t180 = t172*t49;
    t182 = t174*t50;
    t184 = t178*t49;
    return(2.0*t2*(4.0*t7+3.0*(2.0*alpha-2.0)*t6+2.0*(-3.0*alpha+1.0)*t5)+t24*
t29+t31*t34-beta*(t24*t52+t57*t76+t31*t80+t82*t85)+(2.0*t2*(0.12E2*phi+0.12E2*
old+12.0*alpha-12.0)+t98*t29+t100*t34-beta*(t98*t52+3.0*t107*t76+3.0*t24*(
0.882498173E1*t132+0.3776949278E2*t135+0.289445110522E2*t137+0.289445110522E2*
t139)+t57*(0.882498173E1*t172+0.1133084783E3*t175+0.6671400383E2*t178+
0.289445110522E2*t180+0.8683353315E2*t182+0.289445110522E2*t184)+t100*t80-3.0*
t107*t85+3.0*t31*(0.2889626582E2*t132+0.3839831904E2*t135+0.950205321859E1*t137
+0.950205321859E1*t139)+t82*(0.2889626582E2*t172+0.1151949571E3*t175+
0.4790037226E2*t178+0.950205321859E1*t180+0.2850615966E2*t182+0.950205321859E1*
t184)))*(phi-old)/24.0);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi,double old,double mid)
{
  double t100;
  double t105;
  double t107;
  double t109;
  double t112;
  double t113;
  double t114;
  double t120;
  double t122;
  double t123;
  double t125;
  double t126;
  double t128;
  double t130;
  double t132;
  double t136;
  double t138;
  double t141;
  double t143;
  double t146;
  double t148;
  double t153;
  double t16;
  double t160;
  double t162;
  double t168;
  double t173;
  double t180;
  double t187;
  double t191;
  double t195;
  double t198;
  double t2;
  double t200;
  double t201;
  double t203;
  double t204;
  double t206;
  double t208;
  double t21;
  double t210;
  double t216;
  double t218;
  double t220;
  double t223;
  double t228;
  double t23;
  double t231;
  double t232;
  double t233;
  double t234;
  double t237;
  double t24;
  double t241;
  double t246;
  double t248;
  double t249;
  double t251;
  double t252;
  double t254;
  double t256;
  double t257;
  double t259;
  double t26;
  double t261;
  double t263;
  double t264;
  double t266;
  double t274;
  double t28;
  double t282;
  double t285;
  double t289;
  double t290;
  double t291;
  double t293;
  double t294;
  double t296;
  double t297;
  double t299;
  double t300;
  double t302;
  double t304;
  double t306;
  double t308;
  double t31;
  double t33;
  double t34;
  double t368;
  double t370;
  double t372;
  double t373;
  double t375;
  double t376;
  double t378;
  double t379;
  double t38;
  double t381;
  double t383;
  double t385;
  double t386;
  double t388;
  double t390;
  double t392;
  double t394;
  double t396;
  double t398;
  double t400;
  double t402;
  double t405;
  double t410;
  double t417;
  double t42;
  double t43;
  double t437;
  double t45;
  double t453;
  double t455;
  double t47;
  double t48;
  double t5;
  double t50;
  double t55;
  double t59;
  double t6;
  double t61;
  double t62;
  double t63;
  double t67;
  double t69;
  double t70;
  double t72;
  double t74;
  double t79;
  double t83;
  double t86;
  double t88;
  double t90;
  double t91;
  double t93;
  double t95;
  double t97;
  double t98;
  double t99;
  {
    t2 = 1/delta*A;
    t5 = 0.5*phi+0.5*old;
    t6 = t5*t5;
    t16 = t6*t5;
    t21 = 0.6E2*t16-0.9E2*t6+0.15E2*phi+0.15E2*old;
    t23 = log(rho);
    t24 = rho*t23;
    t26 = -0.2011952932E2*rho+0.289445110522E2*t24+0.2133795654E2;
    t28 = -t21;
    t31 = 0.193942126E2*rho+0.950205321859E1*t24+0.4540504209;
    t33 = t6*t6;
    t34 = t33*t5;
    t38 = -0.2889626582E2+0.1204277045E3*t34-0.3010692614E3*t33+0.2007128409E3*
t16;
    t42 = 0.9502053219E1+0.116654747E3*t34-0.2916368674E3*t33+0.1944245783E3*
t16;
    t43 = 1/t42;
    t45 = exp(t38*t43);
    t47 = log(t45);
    t48 = t45*t47;
    t50 = 0.2133795654E2-0.2011952932E2*t45+0.289445110522E2*t48;
    t55 = 30.0*t33-60.0*t16+30.0*t6;
    t59 = 0.3010692612E3*t33-0.6021385228E3*t16+0.3010692614E3*t6;
    t61 = t42*t42;
    t62 = 1/t61;
    t63 = t38*t62;
    t67 = 0.2916368675E3*t33-0.5832737348E3*t16+0.2916368674E3*t6;
    t69 = t59*t43-t63*t67;
    t70 = t69*t45;
    t72 = t70*t47;
    t74 = 0.882498173E1*t70+0.289445110522E2*t72;
    t79 = 0.15E2*t33-0.3E2*t16+0.15E2*t6;
    t83 = 0.6021385225E3*t33-0.1204277046E4*t16+0.6021385227E3*t6;
    t86 = 0.116654747E4*t16;
    t88 = 0.583273735E3*t33-t86+0.5832737349E3*t6;
    t90 = t83*t43-t63*t88;
    t91 = t90*t45;
    t93 = t91*t47;
    t95 = 0.882498173E1*t91+0.289445110522E2*t93;
    t97 = 6.0*t34;
    t98 = 15.0*t33;
    t99 = 10.0*t16;
    t100 = t97-t98+t99;
    t105 = 0.1204277045E4*t16-0.1806415569E4*t6+0.3010692614E3*phi+
0.3010692614E3*old;
    t107 = t83*t62;
    t109 = t59*t62;
    t112 = 1/t61/t42;
    t113 = t38*t112;
    t114 = t88*t67;
    t120 = t86-0.1749821205E4*t6+0.2916368674E3*phi+0.2916368674E3*old;
    t122 = t105*t43-t107*t67-t109*t88+2.0*t113*t114-t63*t120;
    t123 = t122*t45;
    t125 = t90*t69;
    t126 = t125*t45;
    t128 = t123*t47;
    t130 = t125*t48;
    t132 = 0.882498173E1*t123+0.3776949278E2*t126+0.289445110522E2*t128+
0.289445110522E2*t130;
    t136 = 0.4540504209+0.193942126E2*t45+0.950205321859E1*t48;
    t138 = -t55;
    t141 = 0.2889626582E2*t70+0.950205321859E1*t72;
    t143 = -t79;
    t146 = 0.2889626582E2*t91+0.950205321859E1*t93;
    t148 = 1.0-t97+t98-t99;
    t153 = 0.2889626582E2*t123+0.3839831904E2*t126+0.950205321859E1*t128+
0.950205321859E1*t130;
    t160 = 0.18E3*phi+0.18E3*old-0.18E3;
    t162 = -t160;
    t168 = 360.0*t6-0.18E3*phi-0.18E3*old+60.0;
    t173 = 0.18E3*t6-0.9E2*phi-0.9E2*old+0.3E2;
    t180 = 120.0*t16-180.0*t6+0.3E2*phi+0.3E2*old;
    t187 = 0.240855409E4*t16-0.3612831138E4*t6+0.6021385225E3*phi+
0.6021385225E3*old;
    t191 = t88*t88;
    t195 = 0.349964241E4*t6;
    t198 = 0.233309494E4*t16-t195+0.583273735E3*phi+0.583273735E3*old;
    t200 = t187*t43-2.0*t107*t88+2.0*t113*t191-t63*t198;
    t201 = t200*t45;
    t203 = t90*t90;
    t204 = t203*t45;
    t206 = t201*t47;
    t208 = t204*t47;
    t210 = 0.882498173E1*t201+0.3776949278E2*t204+0.289445110522E2*t206+
0.289445110522E2*t208;
    t216 = 0.3612831135E4*t6-0.1806415569E4*phi-0.1806415569E4*old+
0.6021385225E3;
    t218 = t187*t62;
    t220 = t105*t62;
    t223 = t83*t112;
    t228 = t59*t112;
    t231 = t61*t61;
    t232 = 1/t231;
    t233 = t38*t232;
    t234 = t191*t67;
    t237 = t88*t120;
    t241 = t198*t67;
    t246 = t195-0.1749821205E4*phi-0.1749821205E4*old+0.583273735E3;
    t248 = t216*t43-t218*t67-2.0*t220*t88+4.0*t223*t114-2.0*t107*t120+2.0*t228*
t191-6.0*t233*t234+4.0*t113*t237-t109*t198+2.0*t113*t241-t63*t246;
    t249 = t248*t45;
    t251 = t200*t69;
    t252 = t251*t45;
    t254 = t91*t122;
    t256 = t203*t69;
    t257 = t256*t45;
    t259 = t249*t47;
    t261 = t251*t48;
    t263 = t47*t122;
    t264 = t91*t263;
    t266 = t256*t48;
    t274 = 0.722566227E4*t6-0.3612831138E4*phi-0.3612831138E4*old+
0.1204277045E4;
    t282 = t191*t88;
    t285 = t88*t198;
    t289 = 0.349964241E4*phi;
    t290 = 0.349964241E4*old;
    t291 = 0.699928482E4*t6-t289-t290+0.116654747E4;
    t293 = t274*t43-3.0*t218*t88+6.0*t223*t191-3.0*t107*t198-6.0*t233*t282+6.0*
t113*t285-t63*t291;
    t294 = t293*t45;
    t296 = t200*t90;
    t297 = t296*t45;
    t299 = t203*t90;
    t300 = t299*t45;
    t302 = t294*t47;
    t304 = t296*t48;
    t306 = t300*t47;
    t308 = 0.882498173E1*t294+0.1133084783E3*t297+0.6671400383E2*t300+
0.289445110522E2*t302+0.8683353315E2*t304+0.289445110522E2*t306;
    t368 = -3.0*t107*t246-6.0*t59*t232*t282+24.0*t38/t231/t42*t282*t67-18.0*
t233*t191*t120+6.0*t228*t285-18.0*t233*t285*t67+6.0*t113*t120*t198+6.0*t113*t88
*t246-t109*t291+2.0*t113*t291*t67-t63*(-0.349964241E4+t289+t290);
    t370 = ((0.3612831135E4*phi+0.3612831135E4*old-0.3612831138E4)*t43-t274*t62
*t67-3.0*t216*t62*t88+6.0*t187*t112*t114-3.0*t218*t120+6.0*t105*t112*t191-18.0*
t83*t232*t234+12.0*t223*t237-3.0*t220*t198+6.0*t223*t241+t368)*t45;
    t372 = t293*t69;
    t373 = t372*t45;
    t375 = t248*t90;
    t376 = t375*t45;
    t378 = t200*t122;
    t379 = t378*t45;
    t381 = t296*t70;
    t383 = t204*t122;
    t385 = t299*t69;
    t386 = t385*t45;
    t388 = t370*t47;
    t390 = t372*t48;
    t392 = t375*t48;
    t394 = t378*t48;
    t396 = t296*t72;
    t398 = t204*t263;
    t400 = t385*t48;
    t402 = 0.882498173E1*t370+0.3776949278E2*t373+0.1133084783E3*t376+
0.1133084783E3*t379+0.2001420114E3*t381+0.2001420115E3*t383+0.9565851488E2*t386
+0.289445110522E2*t388+0.289445110522E2*t390+0.8683353315E2*t392+0.8683353315E2
*t394+0.8683353315E2*t396+0.8683353315E2*t398+0.289445110522E2*t400;
    t405 = -t168;
    t410 = -t180;
    t417 = 0.2889626582E2*t201+0.3839831904E2*t204+0.950205321859E1*t206+
0.950205321859E1*t208;
    t437 = 0.2889626582E2*t294+0.1151949571E3*t297+0.4790037226E2*t300+
0.950205321859E1*t302+0.2850615966E2*t304+0.950205321859E1*t306;
    t453 = 0.2889626582E2*t370+0.3839831904E2*t373+0.1151949571E3*t376+
0.1151949571E3*t379+0.1437011168E3*t381+0.1437011168E3*t383+0.5740242548E2*t386
+0.950205321859E1*t388+0.950205321859E1*t390+0.2850615966E2*t392+0.2850615966E2
*t394+0.2850615966E2*t396+0.2850615966E2*t398+0.950205321859E1*t400;
    t455 = t160*t50+t168*t74+3.0*t173*t95+3.0*t180*t132+3.0*t21*t210+3.0*t55*(
0.882498173E1*t249+0.3776949278E2*t252+0.7553898556E2*t254+0.6671400383E2*t257+
0.289445110522E2*t259+0.289445110522E2*t261+0.578890221E2*t264+0.289445110522E2
*t266)+t79*t308+t100*t402+t162*t136+t405*t141-3.0*t173*t146+3.0*t410*t153+3.0*
t28*t417+3.0*t138*(0.2889626582E2*t249+0.3839831904E2*t252+0.7679663808E2*t254+
0.4790037226E2*t257+0.950205321859E1*t259+0.950205321859E1*t261+0.1900410644E2*
t264+0.950205321859E1*t266)+t143*t437+t148*t453;
    return(2.0*t2*(0.6E1*t6+0.3E1*(2.0*alpha-2.0)*t5-0.3E1*alpha+0.1E1)+t21*t26
+t28*t31-beta*(t21*t50+t55*t74+t79*t95+t100*t132+t28*t136+t138*t141+t143*t146+
t148*t153)+(0.24E2*t2+t160*t26+t162*t31-beta*t455)*(phi-old)/24.0+t2*(0.12E2*
phi+0.12E2*old+12.0*alpha-12.0)/12.0+t168*t26/24.0+t405*t31/24.0-beta*(t168*t50
+3.0*t180*t95+3.0*t55*t210+t100*t308+t405*t136+3.0*t410*t146+3.0*t138*t417+t148
*t437)/24.0);
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t11;
  double t12;
  double t16;
  double t2;
  double t20;
  double t21;
  double t4;
  double t5;
  double t7;
  double t8;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t8 = t4-t5+t7;
    t11 = 0.5*rho+0.5*old;
    t12 = log(t11);
    t16 = 1.0-t4+t5-t7;
    t20 = t11*t11;
    t21 = 1/t20;
    return(t8*(0.882498173E1+0.289445110522E2*t12)+t16*(0.2889626582E2+
0.950205321859E1*t12)+(-0.2894451105E2*t8*t21-0.9502053219E1*t16*t21)*(rho-old)
/24.0);
  }
}

double drhochemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t11;
  double t12;
  double t15;
  double t18;
  double t2;
  double t20;
  double t29;
  double t4;
  double t5;
  double t7;
  double t8;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t8 = t4-t5+t7;
    t11 = 0.5*rho+0.5*old;
    t12 = 1/t11;
    t15 = 1.0-t4+t5-t7;
    t18 = t11*t11;
    t20 = 1/t18/t11;
    t29 = 1/t18;
    return(0.1447225552E2*t8*t12+0.475102661E1*t15*t12+(0.2894451105E2*t8*t20+
0.9502053219E1*t15*t20)*(rho-old)/24.0-0.1206021294E1*t8*t29-0.3959188841*t15*
t29);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t10;
  double t11;
  double t15;
  double t19;
  double t2;
  double t20;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t7 = 30.0*t2-60.0*t1*phi+30.0*t1;
    t10 = 0.5*rho+0.5*old;
    t11 = log(t10);
    t15 = -t7;
    t19 = t10*t10;
    t20 = 1/t19;
    return(t7*(0.882498173E1+0.289445110522E2*t11)+t15*(0.2889626582E2+
0.950205321859E1*t11)+(-0.2894451105E2*t7*t20-0.9502053219E1*t15*t20)*(rho-old)
/24.0);
  }
}

#include <math.h>
double pressure(double rho,double phi)
{
  double t1;
  double t10;
  double t11;
  double t2;
  double t4;
  double t5;
  double t6;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t6 = t1*phi;
    t7 = 10.0*t6;
    t10 = log(rho);
    t11 = rho*t10;
    return((t4-t5+t7)*(-0.2094881058E2+0.4291934897E2*rho-0.2894451105E2*t11+
rho*(-0.1397483792E2+0.289445110522E2*t10))+(1.0-t4+t5-t7)*(-0.6490446091E-1+
0.3405607049E1*rho-0.9502053219E1*t11+rho*(0.609644617E1+0.950205321859E1*t10))
-2.0/delta*A*(t2+(2.0*alpha-2.0)*t6+(-3.0*alpha+1.0)*t1+alpha));
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t1;
  double t2;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t9 = sqrt(0.116654747E11*t2*phi-0.2916368675E11*t2+0.1944245783E11*t1*phi+
950205322.0);
    return(0.1E-3*t9);
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t15;
  double t16;
  double t17;
  double t18;
  double t19;
  double t21;
  double t22;
  double t26;
  double t3;
  double t4;
  double t41;
  double t43;
  double t44;
  double t7;
  {
    t3 = phi*phi;
    t4 = t3*t3;
    t7 = t3*phi;
    t15 = t4*phi;
    t16 = 6.0*t15;
    t17 = 15.0*t4;
    t18 = 10.0*t7;
    t19 = t16-t17+t18;
    t21 = log(rho);
    t22 = rho*t21;
    t26 = 1.0-t16+t17-t18;
    t41 = exp((-0.2889626582E2+0.1204277045E3*t15-0.3010692614E3*t4+
0.2007128409E3*t7)/(0.9502053219E1+0.116654747E3*t15-0.2916368674E3*t4+
0.1944245783E3*t7));
    t43 = log(t41);
    t44 = t41*t43;
    return(2.0/delta*A*(t4+(2.0*alpha-2.0)*t7+(-3.0*alpha+1.0)*t3+alpha)+t19*(
-0.2011952932E2*rho+0.289445110522E2*t22+0.2133795654E2)+t26*(0.193942126E2*rho
+0.950205321859E1*t22+0.4540504209)-beta*(t19*(0.2133795654E2-0.2011952932E2*
t41+0.289445110522E2*t44)+t26*(0.4540504209+0.193942126E2*t41+0.950205321859E1*
t44)));
  }
}

#include <math.h>
double reactionSource(double rho,double phi,double old)
{
  double t100;
  double t107;
  double t114;
  double t116;
  double t120;
  double t121;
  double t122;
  double t129;
  double t131;
  double t132;
  double t134;
  double t135;
  double t137;
  double t139;
  double t157;
  double t172;
  double t174;
  double t175;
  double t178;
  double t180;
  double t182;
  double t184;
  double t2;
  double t20;
  double t24;
  double t26;
  double t27;
  double t29;
  double t31;
  double t34;
  double t36;
  double t40;
  double t44;
  double t45;
  double t47;
  double t49;
  double t5;
  double t50;
  double t52;
  double t54;
  double t55;
  double t56;
  double t57;
  double t6;
  double t61;
  double t63;
  double t64;
  double t65;
  double t69;
  double t7;
  double t71;
  double t72;
  double t74;
  double t76;
  double t80;
  double t82;
  double t85;
  double t98;
  {
    t2 = 1/delta*A;
    t5 = 0.5*phi+0.5*old;
    t6 = t5*t5;
    t7 = t6*t5;
    t20 = t6*t6;
    t24 = 30.0*t20-60.0*t7+30.0*t6;
    t26 = log(rho);
    t27 = rho*t26;
    t29 = -0.2011952932E2*rho+0.289445110522E2*t27+0.2133795654E2;
    t31 = -t24;
    t34 = 0.193942126E2*rho+0.950205321859E1*t27+0.4540504209;
    t36 = t20*t5;
    t40 = -0.2889626582E2+0.1204277045E3*t36-0.3010692614E3*t20+0.2007128409E3*
t7;
    t44 = 0.9502053219E1+0.116654747E3*t36-0.2916368674E3*t20+0.1944245783E3*t7
;
    t45 = 1/t44;
    t47 = exp(t40*t45);
    t49 = log(t47);
    t50 = t47*t49;
    t52 = 0.2133795654E2-0.2011952932E2*t47+0.289445110522E2*t50;
    t54 = 6.0*t36;
    t55 = 15.0*t20;
    t56 = 10.0*t7;
    t57 = t54-t55+t56;
    t61 = 0.6021385225E3*t20-0.1204277046E4*t7+0.6021385227E3*t6;
    t63 = t44*t44;
    t64 = 1/t63;
    t65 = t40*t64;
    t69 = 0.583273735E3*t20-0.116654747E4*t7+0.5832737349E3*t6;
    t71 = t61*t45-t65*t69;
    t72 = t71*t47;
    t74 = t72*t49;
    t76 = 0.882498173E1*t72+0.289445110522E2*t74;
    t80 = 0.4540504209+0.193942126E2*t47+0.950205321859E1*t50;
    t82 = 1.0-t54+t55-t56;
    t85 = 0.2889626582E2*t72+0.950205321859E1*t74;
    t98 = 360.0*t6-0.18E3*phi-0.18E3*old+60.0;
    t100 = -t98;
    t107 = 120.0*t7-180.0*t6+0.3E2*phi+0.3E2*old;
    t114 = 0.240855409E4*t7-0.3612831138E4*t6+0.6021385225E3*phi+0.6021385225E3
*old;
    t116 = t61*t64;
    t120 = 1/t63/t44;
    t121 = t40*t120;
    t122 = t69*t69;
    t129 = 0.233309494E4*t7-0.349964241E4*t6+0.583273735E3*phi+0.583273735E3*
old;
    t131 = t114*t45-2.0*t116*t69+2.0*t121*t122-t65*t129;
    t132 = t131*t47;
    t134 = t71*t71;
    t135 = t134*t47;
    t137 = t132*t49;
    t139 = t135*t49;
    t157 = t63*t63;
    t172 = ((0.722566227E4*t6-0.3612831138E4*phi-0.3612831138E4*old+
0.1204277045E4)*t45-3.0*t114*t64*t69+6.0*t61*t120*t122-3.0*t116*t129-6.0*t40/
t157*t122*t69+6.0*t121*t69*t129-t65*(0.699928482E4*t6-0.349964241E4*phi
-0.349964241E4*old+0.116654747E4))*t47;
    t174 = t131*t71;
    t175 = t174*t47;
    t178 = t134*t71*t47;
    t180 = t172*t49;
    t182 = t174*t50;
    t184 = t178*t49;
    return(2.0*t2*(4.0*t7+3.0*(2.0*alpha-2.0)*t6+2.0*(-3.0*alpha+1.0)*t5)+t24*
t29+t31*t34-beta*(t24*t52+t57*t76+t31*t80+t82*t85)+(2.0*t2*(0.12E2*phi+0.12E2*
old+12.0*alpha-12.0)+t98*t29+t100*t34-beta*(t98*t52+3.0*t107*t76+3.0*t24*(
0.882498173E1*t132+0.3776949278E2*t135+0.289445110522E2*t137+0.289445110522E2*
t139)+t57*(0.882498173E1*t172+0.1133084783E3*t175+0.6671400383E2*t178+
0.289445110522E2*t180+0.8683353315E2*t182+0.289445110522E2*t184)+t100*t80-3.0*
t107*t85+3.0*t31*(0.2889626582E2*t132+0.3839831904E2*t135+0.950205321859E1*t137
+0.950205321859E1*t139)+t82*(0.2889626582E2*t172+0.1151949571E3*t175+
0.4790037226E2*t178+0.950205321859E1*t180+0.2850615966E2*t182+0.950205321859E1*
t184)))*(phi-old)/24.0);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi,double old,double mid)
{
  double t100;
  double t105;
  double t107;
  double t109;
  double t112;
  double t113;
  double t114;
  double t120;
  double t122;
  double t123;
  double t125;
  double t126;
  double t128;
  double t130;
  double t132;
  double t136;
  double t138;
  double t141;
  double t143;
  double t146;
  double t148;
  double t153;
  double t16;
  double t160;
  double t162;
  double t168;
  double t173;
  double t180;
  double t187;
  double t191;
  double t195;
  double t198;
  double t2;
  double t200;
  double t201;
  double t203;
  double t204;
  double t206;
  double t208;
  double t21;
  double t210;
  double t216;
  double t218;
  double t220;
  double t223;
  double t228;
  double t23;
  double t231;
  double t232;
  double t233;
  double t234;
  double t237;
  double t24;
  double t241;
  double t246;
  double t248;
  double t249;
  double t251;
  double t252;
  double t254;
  double t256;
  double t257;
  double t259;
  double t26;
  double t261;
  double t263;
  double t264;
  double t266;
  double t274;
  double t28;
  double t282;
  double t285;
  double t289;
  double t290;
  double t291;
  double t293;
  double t294;
  double t296;
  double t297;
  double t299;
  double t300;
  double t302;
  double t304;
  double t306;
  double t308;
  double t31;
  double t33;
  double t34;
  double t368;
  double t370;
  double t372;
  double t373;
  double t375;
  double t376;
  double t378;
  double t379;
  double t38;
  double t381;
  double t383;
  double t385;
  double t386;
  double t388;
  double t390;
  double t392;
  double t394;
  double t396;
  double t398;
  double t400;
  double t402;
  double t405;
  double t410;
  double t417;
  double t42;
  double t43;
  double t437;
  double t45;
  double t453;
  double t455;
  double t47;
  double t48;
  double t5;
  double t50;
  double t55;
  double t59;
  double t6;
  double t61;
  double t62;
  double t63;
  double t67;
  double t69;
  double t70;
  double t72;
  double t74;
  double t79;
  double t83;
  double t86;
  double t88;
  double t90;
  double t91;
  double t93;
  double t95;
  double t97;
  double t98;
  double t99;
  {
    t2 = 1/delta*A;
    t5 = 0.5*phi+0.5*old;
    t6 = t5*t5;
    t16 = t6*t5;
    t21 = 0.6E2*t16-0.9E2*t6+0.15E2*phi+0.15E2*old;
    t23 = log(rho);
    t24 = rho*t23;
    t26 = -0.2011952932E2*rho+0.289445110522E2*t24+0.2133795654E2;
    t28 = -t21;
    t31 = 0.193942126E2*rho+0.950205321859E1*t24+0.4540504209;
    t33 = t6*t6;
    t34 = t33*t5;
    t38 = -0.2889626582E2+0.1204277045E3*t34-0.3010692614E3*t33+0.2007128409E3*
t16;
    t42 = 0.9502053219E1+0.116654747E3*t34-0.2916368674E3*t33+0.1944245783E3*
t16;
    t43 = 1/t42;
    t45 = exp(t38*t43);
    t47 = log(t45);
    t48 = t45*t47;
    t50 = 0.2133795654E2-0.2011952932E2*t45+0.289445110522E2*t48;
    t55 = 30.0*t33-60.0*t16+30.0*t6;
    t59 = 0.3010692612E3*t33-0.6021385228E3*t16+0.3010692614E3*t6;
    t61 = t42*t42;
    t62 = 1/t61;
    t63 = t38*t62;
    t67 = 0.2916368675E3*t33-0.5832737348E3*t16+0.2916368674E3*t6;
    t69 = t59*t43-t63*t67;
    t70 = t69*t45;
    t72 = t70*t47;
    t74 = 0.882498173E1*t70+0.289445110522E2*t72;
    t79 = 0.15E2*t33-0.3E2*t16+0.15E2*t6;
    t83 = 0.6021385225E3*t33-0.1204277046E4*t16+0.6021385227E3*t6;
    t86 = 0.116654747E4*t16;
    t88 = 0.583273735E3*t33-t86+0.5832737349E3*t6;
    t90 = t83*t43-t63*t88;
    t91 = t90*t45;
    t93 = t91*t47;
    t95 = 0.882498173E1*t91+0.289445110522E2*t93;
    t97 = 6.0*t34;
    t98 = 15.0*t33;
    t99 = 10.0*t16;
    t100 = t97-t98+t99;
    t105 = 0.1204277045E4*t16-0.1806415569E4*t6+0.3010692614E3*phi+
0.3010692614E3*old;
    t107 = t83*t62;
    t109 = t59*t62;
    t112 = 1/t61/t42;
    t113 = t38*t112;
    t114 = t88*t67;
    t120 = t86-0.1749821205E4*t6+0.2916368674E3*phi+0.2916368674E3*old;
    t122 = t105*t43-t107*t67-t109*t88+2.0*t113*t114-t63*t120;
    t123 = t122*t45;
    t125 = t90*t69;
    t126 = t125*t45;
    t128 = t123*t47;
    t130 = t125*t48;
    t132 = 0.882498173E1*t123+0.3776949278E2*t126+0.289445110522E2*t128+
0.289445110522E2*t130;
    t136 = 0.4540504209+0.193942126E2*t45+0.950205321859E1*t48;
    t138 = -t55;
    t141 = 0.2889626582E2*t70+0.950205321859E1*t72;
    t143 = -t79;
    t146 = 0.2889626582E2*t91+0.950205321859E1*t93;
    t148 = 1.0-t97+t98-t99;
    t153 = 0.2889626582E2*t123+0.3839831904E2*t126+0.950205321859E1*t128+
0.950205321859E1*t130;
    t160 = 0.18E3*phi+0.18E3*old-0.18E3;
    t162 = -t160;
    t168 = 360.0*t6-0.18E3*phi-0.18E3*old+60.0;
    t173 = 0.18E3*t6-0.9E2*phi-0.9E2*old+0.3E2;
    t180 = 120.0*t16-180.0*t6+0.3E2*phi+0.3E2*old;
    t187 = 0.240855409E4*t16-0.3612831138E4*t6+0.6021385225E3*phi+
0.6021385225E3*old;
    t191 = t88*t88;
    t195 = 0.349964241E4*t6;
    t198 = 0.233309494E4*t16-t195+0.583273735E3*phi+0.583273735E3*old;
    t200 = t187*t43-2.0*t107*t88+2.0*t113*t191-t63*t198;
    t201 = t200*t45;
    t203 = t90*t90;
    t204 = t203*t45;
    t206 = t201*t47;
    t208 = t204*t47;
    t210 = 0.882498173E1*t201+0.3776949278E2*t204+0.289445110522E2*t206+
0.289445110522E2*t208;
    t216 = 0.3612831135E4*t6-0.1806415569E4*phi-0.1806415569E4*old+
0.6021385225E3;
    t218 = t187*t62;
    t220 = t105*t62;
    t223 = t83*t112;
    t228 = t59*t112;
    t231 = t61*t61;
    t232 = 1/t231;
    t233 = t38*t232;
    t234 = t191*t67;
    t237 = t88*t120;
    t241 = t198*t67;
    t246 = t195-0.1749821205E4*phi-0.1749821205E4*old+0.583273735E3;
    t248 = t216*t43-t218*t67-2.0*t220*t88+4.0*t223*t114-2.0*t107*t120+2.0*t228*
t191-6.0*t233*t234+4.0*t113*t237-t109*t198+2.0*t113*t241-t63*t246;
    t249 = t248*t45;
    t251 = t200*t69;
    t252 = t251*t45;
    t254 = t91*t122;
    t256 = t203*t69;
    t257 = t256*t45;
    t259 = t249*t47;
    t261 = t251*t48;
    t263 = t47*t122;
    t264 = t91*t263;
    t266 = t256*t48;
    t274 = 0.722566227E4*t6-0.3612831138E4*phi-0.3612831138E4*old+
0.1204277045E4;
    t282 = t191*t88;
    t285 = t88*t198;
    t289 = 0.349964241E4*phi;
    t290 = 0.349964241E4*old;
    t291 = 0.116654747E4+0.699928482E4*t6-t289-t290;
    t293 = t274*t43-3.0*t218*t88+6.0*t223*t191-3.0*t107*t198-6.0*t233*t282+6.0*
t113*t285-t63*t291;
    t294 = t293*t45;
    t296 = t200*t90;
    t297 = t296*t45;
    t299 = t203*t90;
    t300 = t299*t45;
    t302 = t294*t47;
    t304 = t296*t48;
    t306 = t300*t47;
    t308 = 0.882498173E1*t294+0.1133084783E3*t297+0.6671400383E2*t300+
0.289445110522E2*t302+0.8683353315E2*t304+0.289445110522E2*t306;
    t368 = -3.0*t107*t246-6.0*t59*t232*t282+24.0*t38/t231/t42*t282*t67-18.0*
t233*t191*t120+6.0*t228*t285-18.0*t233*t285*t67+6.0*t113*t120*t198+6.0*t113*t88
*t246-t109*t291+2.0*t113*t291*t67-t63*(-0.349964241E4+t289+t290);
    t370 = ((-0.3612831138E4+0.3612831135E4*phi+0.3612831135E4*old)*t43-t274*
t62*t67-3.0*t216*t62*t88+6.0*t187*t112*t114-3.0*t218*t120+6.0*t105*t112*t191
-18.0*t83*t232*t234+12.0*t223*t237-3.0*t220*t198+6.0*t223*t241+t368)*t45;
    t372 = t293*t69;
    t373 = t372*t45;
    t375 = t248*t90;
    t376 = t375*t45;
    t378 = t200*t122;
    t379 = t378*t45;
    t381 = t296*t70;
    t383 = t204*t122;
    t385 = t299*t69;
    t386 = t385*t45;
    t388 = t370*t47;
    t390 = t372*t48;
    t392 = t375*t48;
    t394 = t378*t48;
    t396 = t296*t72;
    t398 = t204*t263;
    t400 = t385*t48;
    t402 = 0.882498173E1*t370+0.3776949278E2*t373+0.1133084783E3*t376+
0.1133084783E3*t379+0.2001420114E3*t381+0.2001420115E3*t383+0.9565851488E2*t386
+0.289445110522E2*t388+0.289445110522E2*t390+0.8683353315E2*t392+0.8683353315E2
*t394+0.8683353315E2*t396+0.8683353315E2*t398+0.289445110522E2*t400;
    t405 = -t168;
    t410 = -t180;
    t417 = 0.2889626582E2*t201+0.3839831904E2*t204+0.950205321859E1*t206+
0.950205321859E1*t208;
    t437 = 0.2889626582E2*t294+0.1151949571E3*t297+0.4790037226E2*t300+
0.950205321859E1*t302+0.2850615966E2*t304+0.950205321859E1*t306;
    t453 = 0.2889626582E2*t370+0.3839831904E2*t373+0.1151949571E3*t376+
0.1151949571E3*t379+0.1437011168E3*t381+0.1437011168E3*t383+0.5740242548E2*t386
+0.950205321859E1*t388+0.950205321859E1*t390+0.2850615966E2*t392+0.2850615966E2
*t394+0.2850615966E2*t396+0.2850615966E2*t398+0.950205321859E1*t400;
    t455 = t160*t50+t168*t74+3.0*t173*t95+3.0*t180*t132+3.0*t21*t210+3.0*t55*(
0.882498173E1*t249+0.3776949278E2*t252+0.7553898556E2*t254+0.6671400383E2*t257+
0.289445110522E2*t259+0.289445110522E2*t261+0.578890221E2*t264+0.289445110522E2
*t266)+t79*t308+t100*t402+t162*t136+t405*t141-3.0*t173*t146+3.0*t410*t153+3.0*
t28*t417+3.0*t138*(0.2889626582E2*t249+0.3839831904E2*t252+0.7679663808E2*t254+
0.4790037226E2*t257+0.950205321859E1*t259+0.950205321859E1*t261+0.1900410644E2*
t264+0.950205321859E1*t266)+t143*t437+t148*t453;
    return(2.0*t2*(0.6E1*t6+0.3E1*(2.0*alpha-2.0)*t5-0.3E1*alpha+0.1E1)+t21*t26
+t28*t31-beta*(t21*t50+t55*t74+t79*t95+t100*t132+t28*t136+t138*t141+t143*t146+
t148*t153)+(0.24E2*t2+t160*t26+t162*t31-beta*t455)*(phi-old)/24.0+t2*(0.12E2*
phi+0.12E2*old+12.0*alpha-12.0)/12.0+t168*t26/24.0+t405*t31/24.0-beta*(t168*t50
+3.0*t180*t95+3.0*t55*t210+t100*t308+t405*t136+3.0*t410*t146+3.0*t138*t417+t148
*t437)/24.0);
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t11;
  double t12;
  double t16;
  double t2;
  double t20;
  double t21;
  double t4;
  double t5;
  double t7;
  double t8;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t8 = t4-t5+t7;
    t11 = 0.5*rho+0.5*old;
    t12 = log(t11);
    t16 = 1.0-t4+t5-t7;
    t20 = t11*t11;
    t21 = 1/t20;
    return(t8*(0.882498173E1+0.289445110522E2*t12)+t16*(0.2889626582E2+
0.950205321859E1*t12)+(-0.2894451105E2*t8*t21-0.9502053219E1*t16*t21)*(rho-old)
/24.0);
  }
}

double drhochemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t11;
  double t12;
  double t15;
  double t18;
  double t2;
  double t20;
  double t29;
  double t4;
  double t5;
  double t7;
  double t8;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t8 = t4-t5+t7;
    t11 = 0.5*rho+0.5*old;
    t12 = 1/t11;
    t15 = 1.0-t4+t5-t7;
    t18 = t11*t11;
    t20 = 1/t18/t11;
    t29 = 1/t18;
    return(0.1447225552E2*t8*t12+0.475102661E1*t15*t12+(0.2894451105E2*t8*t20+
0.9502053219E1*t15*t20)*(rho-old)/24.0-0.1206021294E1*t8*t29-0.3959188841*t15*
t29);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t10;
  double t11;
  double t15;
  double t19;
  double t2;
  double t20;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t7 = 30.0*t2-60.0*t1*phi+30.0*t1;
    t10 = 0.5*rho+0.5*old;
    t11 = log(t10);
    t15 = -t7;
    t19 = t10*t10;
    t20 = 1/t19;
    return(t7*(0.882498173E1+0.289445110522E2*t11)+t15*(0.2889626582E2+
0.950205321859E1*t11)+(-0.2894451105E2*t7*t20-0.9502053219E1*t15*t20)*(rho-old)
/24.0);
  }
}

#include <math.h>
double pressure(double rho,double phi)
{
  double t1;
  double t10;
  double t11;
  double t2;
  double t4;
  double t5;
  double t6;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t6 = t1*phi;
    t7 = 10.0*t6;
    t10 = log(rho);
    t11 = rho*t10;
    return((t4-t5+t7)*(-0.2094881058E2+0.4291934897E2*rho-0.2894451105E2*t11+
rho*(-0.1397483792E2+0.289445110522E2*t10))+(1.0-t4+t5-t7)*(-0.6490446091E-1+
0.3405607049E1*rho-0.9502053219E1*t11+rho*(0.609644617E1+0.950205321859E1*t10))
-2.0/delta*A*(t2+(2.0*alpha-2.0)*t6+(-3.0*alpha+1.0)*t1+alpha));
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t1;
  double t2;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t9 = sqrt(0.116654747E11*t2*phi-0.2916368675E11*t2+0.1944245783E11*t1*phi+
950205322.0);
    return(0.1E-3*t9);
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t15;
  double t16;
  double t17;
  double t18;
  double t19;
  double t21;
  double t22;
  double t26;
  double t3;
  double t4;
  double t41;
  double t43;
  double t44;
  double t7;
  {
    t3 = phi*phi;
    t4 = t3*t3;
    t7 = t3*phi;
    t15 = t4*phi;
    t16 = 6.0*t15;
    t17 = 15.0*t4;
    t18 = 10.0*t7;
    t19 = t16-t17+t18;
    t21 = log(rho);
    t22 = rho*t21;
    t26 = 1.0-t16+t17-t18;
    t41 = exp((-0.2889626582E2+0.1204277045E3*t15-0.3010692614E3*t4+
0.2007128409E3*t7)/(0.9502053219E1+0.116654747E3*t15-0.2916368674E3*t4+
0.1944245783E3*t7));
    t43 = log(t41);
    t44 = t41*t43;
    return(2.0/delta*A*(t4+(2.0*alpha-2.0)*t7+(-3.0*alpha+1.0)*t3+alpha)+t19*(
-0.2011952932E2*rho+0.289445110522E2*t22+0.2133795654E2)+t26*(0.193942126E2*rho
+0.950205321859E1*t22+0.4540504209)-beta*(t19*(0.2133795654E2-0.2011952932E2*
t41+0.289445110522E2*t44)+t26*(0.4540504209+0.193942126E2*t41+0.950205321859E1*
t44)));
  }
}

#include <math.h>
double reactionSource(double rho,double phi,double old)
{
  double t100;
  double t107;
  double t114;
  double t116;
  double t120;
  double t121;
  double t122;
  double t129;
  double t131;
  double t132;
  double t134;
  double t135;
  double t137;
  double t139;
  double t157;
  double t172;
  double t174;
  double t175;
  double t178;
  double t180;
  double t182;
  double t184;
  double t2;
  double t20;
  double t24;
  double t26;
  double t27;
  double t29;
  double t31;
  double t34;
  double t36;
  double t40;
  double t44;
  double t45;
  double t47;
  double t49;
  double t5;
  double t50;
  double t52;
  double t54;
  double t55;
  double t56;
  double t57;
  double t6;
  double t61;
  double t63;
  double t64;
  double t65;
  double t69;
  double t7;
  double t71;
  double t72;
  double t74;
  double t76;
  double t80;
  double t82;
  double t85;
  double t98;
  {
    t2 = 1/delta*A;
    t5 = 0.5*phi+0.5*old;
    t6 = t5*t5;
    t7 = t6*t5;
    t20 = t6*t6;
    t24 = 30.0*t20-60.0*t7+30.0*t6;
    t26 = log(rho);
    t27 = rho*t26;
    t29 = -0.2011952932E2*rho+0.289445110522E2*t27+0.2133795654E2;
    t31 = -t24;
    t34 = 0.193942126E2*rho+0.950205321859E1*t27+0.4540504209;
    t36 = t20*t5;
    t40 = -0.2889626582E2+0.1204277045E3*t36-0.3010692614E3*t20+0.2007128409E3*
t7;
    t44 = 0.9502053219E1+0.116654747E3*t36-0.2916368674E3*t20+0.1944245783E3*t7
;
    t45 = 1/t44;
    t47 = exp(t40*t45);
    t49 = log(t47);
    t50 = t47*t49;
    t52 = 0.2133795654E2-0.2011952932E2*t47+0.289445110522E2*t50;
    t54 = 6.0*t36;
    t55 = 15.0*t20;
    t56 = 10.0*t7;
    t57 = t54-t55+t56;
    t61 = 0.6021385225E3*t20-0.1204277046E4*t7+0.6021385227E3*t6;
    t63 = t44*t44;
    t64 = 1/t63;
    t65 = t40*t64;
    t69 = 0.583273735E3*t20-0.116654747E4*t7+0.5832737349E3*t6;
    t71 = t61*t45-t65*t69;
    t72 = t71*t47;
    t74 = t72*t49;
    t76 = 0.882498173E1*t72+0.289445110522E2*t74;
    t80 = 0.4540504209+0.193942126E2*t47+0.950205321859E1*t50;
    t82 = 1.0-t54+t55-t56;
    t85 = 0.2889626582E2*t72+0.950205321859E1*t74;
    t98 = 360.0*t6-0.18E3*phi-0.18E3*old+60.0;
    t100 = -t98;
    t107 = 120.0*t7-180.0*t6+0.3E2*phi+0.3E2*old;
    t114 = 0.240855409E4*t7-0.3612831138E4*t6+0.6021385225E3*phi+0.6021385225E3
*old;
    t116 = t61*t64;
    t120 = 1/t63/t44;
    t121 = t40*t120;
    t122 = t69*t69;
    t129 = 0.233309494E4*t7-0.349964241E4*t6+0.583273735E3*phi+0.583273735E3*
old;
    t131 = t114*t45-2.0*t116*t69+2.0*t121*t122-t65*t129;
    t132 = t131*t47;
    t134 = t71*t71;
    t135 = t134*t47;
    t137 = t132*t49;
    t139 = t135*t49;
    t157 = t63*t63;
    t172 = ((0.722566227E4*t6-0.3612831138E4*phi-0.3612831138E4*old+
0.1204277045E4)*t45-3.0*t114*t64*t69+6.0*t61*t120*t122-3.0*t116*t129-6.0*t40/
t157*t122*t69+6.0*t121*t69*t129-t65*(0.699928482E4*t6-0.349964241E4*phi
-0.349964241E4*old+0.116654747E4))*t47;
    t174 = t131*t71;
    t175 = t174*t47;
    t178 = t134*t71*t47;
    t180 = t172*t49;
    t182 = t174*t50;
    t184 = t178*t49;
    return(2.0*t2*(4.0*t7+3.0*(2.0*alpha-2.0)*t6+2.0*(-3.0*alpha+1.0)*t5)+t24*
t29+t31*t34-beta*(t24*t52+t57*t76+t31*t80+t82*t85)+(2.0*t2*(0.12E2*phi+0.12E2*
old+12.0*alpha-12.0)+t98*t29+t100*t34-beta*(t98*t52+3.0*t107*t76+3.0*t24*(
0.882498173E1*t132+0.3776949278E2*t135+0.289445110522E2*t137+0.289445110522E2*
t139)+t57*(0.882498173E1*t172+0.1133084783E3*t175+0.6671400383E2*t178+
0.289445110522E2*t180+0.8683353315E2*t182+0.289445110522E2*t184)+t100*t80-3.0*
t107*t85+3.0*t31*(0.2889626582E2*t132+0.3839831904E2*t135+0.950205321859E1*t137
+0.950205321859E1*t139)+t82*(0.2889626582E2*t172+0.1151949571E3*t175+
0.4790037226E2*t178+0.950205321859E1*t180+0.2850615966E2*t182+0.950205321859E1*
t184)))*(phi-old)/24.0);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi,double old,double mid)
{
  double t100;
  double t105;
  double t107;
  double t109;
  double t112;
  double t113;
  double t114;
  double t120;
  double t122;
  double t123;
  double t125;
  double t126;
  double t128;
  double t130;
  double t132;
  double t136;
  double t138;
  double t141;
  double t143;
  double t146;
  double t148;
  double t153;
  double t16;
  double t160;
  double t162;
  double t168;
  double t173;
  double t180;
  double t187;
  double t191;
  double t195;
  double t198;
  double t2;
  double t200;
  double t201;
  double t203;
  double t204;
  double t206;
  double t208;
  double t21;
  double t210;
  double t216;
  double t218;
  double t220;
  double t223;
  double t228;
  double t23;
  double t231;
  double t232;
  double t233;
  double t234;
  double t237;
  double t24;
  double t241;
  double t246;
  double t248;
  double t249;
  double t251;
  double t252;
  double t254;
  double t256;
  double t257;
  double t259;
  double t26;
  double t261;
  double t263;
  double t264;
  double t266;
  double t274;
  double t28;
  double t282;
  double t285;
  double t289;
  double t290;
  double t291;
  double t293;
  double t294;
  double t296;
  double t297;
  double t299;
  double t300;
  double t302;
  double t304;
  double t306;
  double t308;
  double t31;
  double t33;
  double t34;
  double t368;
  double t370;
  double t372;
  double t373;
  double t375;
  double t376;
  double t378;
  double t379;
  double t38;
  double t381;
  double t383;
  double t385;
  double t386;
  double t388;
  double t390;
  double t392;
  double t394;
  double t396;
  double t398;
  double t400;
  double t402;
  double t405;
  double t410;
  double t417;
  double t42;
  double t43;
  double t437;
  double t45;
  double t453;
  double t455;
  double t47;
  double t48;
  double t5;
  double t50;
  double t55;
  double t59;
  double t6;
  double t61;
  double t62;
  double t63;
  double t67;
  double t69;
  double t70;
  double t72;
  double t74;
  double t79;
  double t83;
  double t86;
  double t88;
  double t90;
  double t91;
  double t93;
  double t95;
  double t97;
  double t98;
  double t99;
  {
    t2 = 1/delta*A;
    t5 = 0.5*phi+0.5*old;
    t6 = t5*t5;
    t16 = t6*t5;
    t21 = 0.6E2*t16-0.9E2*t6+0.15E2*phi+0.15E2*old;
    t23 = log(rho);
    t24 = rho*t23;
    t26 = -0.2011952932E2*rho+0.289445110522E2*t24+0.2133795654E2;
    t28 = -t21;
    t31 = 0.193942126E2*rho+0.950205321859E1*t24+0.4540504209;
    t33 = t6*t6;
    t34 = t33*t5;
    t38 = -0.2889626582E2+0.1204277045E3*t34-0.3010692614E3*t33+0.2007128409E3*
t16;
    t42 = 0.9502053219E1+0.116654747E3*t34-0.2916368674E3*t33+0.1944245783E3*
t16;
    t43 = 1/t42;
    t45 = exp(t38*t43);
    t47 = log(t45);
    t48 = t45*t47;
    t50 = 0.2133795654E2-0.2011952932E2*t45+0.289445110522E2*t48;
    t55 = 30.0*t33-60.0*t16+30.0*t6;
    t59 = 0.3010692612E3*t33-0.6021385228E3*t16+0.3010692614E3*t6;
    t61 = t42*t42;
    t62 = 1/t61;
    t63 = t38*t62;
    t67 = 0.2916368675E3*t33-0.5832737348E3*t16+0.2916368674E3*t6;
    t69 = t59*t43-t63*t67;
    t70 = t69*t45;
    t72 = t70*t47;
    t74 = 0.882498173E1*t70+0.289445110522E2*t72;
    t79 = 0.15E2*t33-0.3E2*t16+0.15E2*t6;
    t83 = 0.6021385225E3*t33-0.1204277046E4*t16+0.6021385227E3*t6;
    t86 = 0.116654747E4*t16;
    t88 = 0.583273735E3*t33-t86+0.5832737349E3*t6;
    t90 = t83*t43-t63*t88;
    t91 = t90*t45;
    t93 = t91*t47;
    t95 = 0.882498173E1*t91+0.289445110522E2*t93;
    t97 = 6.0*t34;
    t98 = 15.0*t33;
    t99 = 10.0*t16;
    t100 = t97-t98+t99;
    t105 = 0.1204277045E4*t16-0.1806415569E4*t6+0.3010692614E3*phi+
0.3010692614E3*old;
    t107 = t83*t62;
    t109 = t59*t62;
    t112 = 1/t61/t42;
    t113 = t38*t112;
    t114 = t88*t67;
    t120 = t86-0.1749821205E4*t6+0.2916368674E3*phi+0.2916368674E3*old;
    t122 = t105*t43-t107*t67-t109*t88+2.0*t113*t114-t63*t120;
    t123 = t122*t45;
    t125 = t90*t69;
    t126 = t125*t45;
    t128 = t123*t47;
    t130 = t125*t48;
    t132 = 0.882498173E1*t123+0.3776949278E2*t126+0.289445110522E2*t128+
0.289445110522E2*t130;
    t136 = 0.4540504209+0.193942126E2*t45+0.950205321859E1*t48;
    t138 = -t55;
    t141 = 0.2889626582E2*t70+0.950205321859E1*t72;
    t143 = -t79;
    t146 = 0.2889626582E2*t91+0.950205321859E1*t93;
    t148 = 1.0-t97+t98-t99;
    t153 = 0.2889626582E2*t123+0.3839831904E2*t126+0.950205321859E1*t128+
0.950205321859E1*t130;
    t160 = 0.18E3*phi+0.18E3*old-0.18E3;
    t162 = -t160;
    t168 = 360.0*t6-0.18E3*phi-0.18E3*old+60.0;
    t173 = 0.18E3*t6-0.9E2*phi-0.9E2*old+0.3E2;
    t180 = 120.0*t16-180.0*t6+0.3E2*phi+0.3E2*old;
    t187 = 0.240855409E4*t16-0.3612831138E4*t6+0.6021385225E3*phi+
0.6021385225E3*old;
    t191 = t88*t88;
    t195 = 0.349964241E4*t6;
    t198 = 0.233309494E4*t16-t195+0.583273735E3*phi+0.583273735E3*old;
    t200 = t187*t43-2.0*t107*t88+2.0*t113*t191-t63*t198;
    t201 = t200*t45;
    t203 = t90*t90;
    t204 = t203*t45;
    t206 = t201*t47;
    t208 = t204*t47;
    t210 = 0.882498173E1*t201+0.3776949278E2*t204+0.289445110522E2*t206+
0.289445110522E2*t208;
    t216 = 0.3612831135E4*t6-0.1806415569E4*phi-0.1806415569E4*old+
0.6021385225E3;
    t218 = t187*t62;
    t220 = t105*t62;
    t223 = t83*t112;
    t228 = t59*t112;
    t231 = t61*t61;
    t232 = 1/t231;
    t233 = t38*t232;
    t234 = t191*t67;
    t237 = t88*t120;
    t241 = t198*t67;
    t246 = t195-0.1749821205E4*phi-0.1749821205E4*old+0.583273735E3;
    t248 = t216*t43-t218*t67-2.0*t220*t88+4.0*t223*t114-2.0*t107*t120+2.0*t228*
t191-6.0*t233*t234+4.0*t113*t237-t109*t198+2.0*t113*t241-t63*t246;
    t249 = t248*t45;
    t251 = t200*t69;
    t252 = t251*t45;
    t254 = t91*t122;
    t256 = t203*t69;
    t257 = t256*t45;
    t259 = t249*t47;
    t261 = t251*t48;
    t263 = t47*t122;
    t264 = t91*t263;
    t266 = t256*t48;
    t274 = 0.722566227E4*t6-0.3612831138E4*phi-0.3612831138E4*old+
0.1204277045E4;
    t282 = t191*t88;
    t285 = t88*t198;
    t289 = 0.349964241E4*phi;
    t290 = 0.349964241E4*old;
    t291 = 0.116654747E4+0.699928482E4*t6-t289-t290;
    t293 = t274*t43-3.0*t218*t88+6.0*t223*t191-3.0*t107*t198-6.0*t233*t282+6.0*
t113*t285-t63*t291;
    t294 = t293*t45;
    t296 = t200*t90;
    t297 = t296*t45;
    t299 = t203*t90;
    t300 = t299*t45;
    t302 = t294*t47;
    t304 = t296*t48;
    t306 = t300*t47;
    t308 = 0.882498173E1*t294+0.1133084783E3*t297+0.6671400383E2*t300+
0.289445110522E2*t302+0.8683353315E2*t304+0.289445110522E2*t306;
    t368 = -3.0*t107*t246-6.0*t59*t232*t282+24.0*t38/t231/t42*t282*t67-18.0*
t233*t191*t120+6.0*t228*t285-18.0*t233*t285*t67+6.0*t113*t120*t198+6.0*t113*t88
*t246-t109*t291+2.0*t113*t291*t67-t63*(-0.349964241E4+t289+t290);
    t370 = ((-0.3612831138E4+0.3612831135E4*phi+0.3612831135E4*old)*t43-t274*
t62*t67-3.0*t216*t62*t88+6.0*t187*t112*t114-3.0*t218*t120+6.0*t105*t112*t191
-18.0*t83*t232*t234+12.0*t223*t237-3.0*t220*t198+6.0*t223*t241+t368)*t45;
    t372 = t293*t69;
    t373 = t372*t45;
    t375 = t248*t90;
    t376 = t375*t45;
    t378 = t200*t122;
    t379 = t378*t45;
    t381 = t296*t70;
    t383 = t204*t122;
    t385 = t299*t69;
    t386 = t385*t45;
    t388 = t370*t47;
    t390 = t372*t48;
    t392 = t375*t48;
    t394 = t378*t48;
    t396 = t296*t72;
    t398 = t204*t263;
    t400 = t385*t48;
    t402 = 0.882498173E1*t370+0.3776949278E2*t373+0.1133084783E3*t376+
0.1133084783E3*t379+0.2001420114E3*t381+0.2001420115E3*t383+0.9565851488E2*t386
+0.289445110522E2*t388+0.289445110522E2*t390+0.8683353315E2*t392+0.8683353315E2
*t394+0.8683353315E2*t396+0.8683353315E2*t398+0.289445110522E2*t400;
    t405 = -t168;
    t410 = -t180;
    t417 = 0.2889626582E2*t201+0.3839831904E2*t204+0.950205321859E1*t206+
0.950205321859E1*t208;
    t437 = 0.2889626582E2*t294+0.1151949571E3*t297+0.4790037226E2*t300+
0.950205321859E1*t302+0.2850615966E2*t304+0.950205321859E1*t306;
    t453 = 0.2889626582E2*t370+0.3839831904E2*t373+0.1151949571E3*t376+
0.1151949571E3*t379+0.1437011168E3*t381+0.1437011168E3*t383+0.5740242548E2*t386
+0.950205321859E1*t388+0.950205321859E1*t390+0.2850615966E2*t392+0.2850615966E2
*t394+0.2850615966E2*t396+0.2850615966E2*t398+0.950205321859E1*t400;
    t455 = t160*t50+t168*t74+3.0*t173*t95+3.0*t180*t132+3.0*t21*t210+3.0*t55*(
0.882498173E1*t249+0.3776949278E2*t252+0.7553898556E2*t254+0.6671400383E2*t257+
0.289445110522E2*t259+0.289445110522E2*t261+0.578890221E2*t264+0.289445110522E2
*t266)+t79*t308+t100*t402+t162*t136+t405*t141-3.0*t173*t146+3.0*t410*t153+3.0*
t28*t417+3.0*t138*(0.2889626582E2*t249+0.3839831904E2*t252+0.7679663808E2*t254+
0.4790037226E2*t257+0.950205321859E1*t259+0.950205321859E1*t261+0.1900410644E2*
t264+0.950205321859E1*t266)+t143*t437+t148*t453;
    return(2.0*t2*(0.6E1*t6+0.3E1*(2.0*alpha-2.0)*t5-0.3E1*alpha+0.1E1)+t21*t26
+t28*t31-beta*(t21*t50+t55*t74+t79*t95+t100*t132+t28*t136+t138*t141+t143*t146+
t148*t153)+(0.24E2*t2+t160*t26+t162*t31-beta*t455)*(phi-old)/24.0+t2*(0.12E2*
phi+0.12E2*old+12.0*alpha-12.0)/12.0+t168*t26/24.0+t405*t31/24.0-beta*(t168*t50
+3.0*t180*t95+3.0*t55*t210+t100*t308+t405*t136+3.0*t410*t146+3.0*t138*t417+t148
*t437)/24.0);
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t11;
  double t12;
  double t16;
  double t2;
  double t20;
  double t21;
  double t4;
  double t5;
  double t7;
  double t8;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t8 = t4-t5+t7;
    t11 = 0.5*rho+0.5*old;
    t12 = log(t11);
    t16 = 1.0-t4+t5-t7;
    t20 = t11*t11;
    t21 = 1/t20;
    return(t8*(0.882498173E1+0.289445110522E2*t12)+t16*(0.2889626582E2+
0.950205321859E1*t12)+(-0.2894451105E2*t8*t21-0.9502053219E1*t16*t21)*(rho-old)
/24.0);
  }
}

double drhochemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t11;
  double t12;
  double t15;
  double t18;
  double t2;
  double t20;
  double t29;
  double t4;
  double t5;
  double t7;
  double t8;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t8 = t4-t5+t7;
    t11 = 0.5*rho+0.5*old;
    t12 = 1/t11;
    t15 = 1.0-t4+t5-t7;
    t18 = t11*t11;
    t20 = 1/t18/t11;
    t29 = 1/t18;
    return(0.1447225552E2*t8*t12+0.475102661E1*t15*t12+(0.2894451105E2*t8*t20+
0.9502053219E1*t15*t20)*(rho-old)/24.0-0.1206021294E1*t8*t29-0.3959188841*t15*
t29);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t10;
  double t11;
  double t15;
  double t19;
  double t2;
  double t20;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t7 = 30.0*t2-60.0*t1*phi+30.0*t1;
    t10 = 0.5*rho+0.5*old;
    t11 = log(t10);
    t15 = -t7;
    t19 = t10*t10;
    t20 = 1/t19;
    return(t7*(0.882498173E1+0.289445110522E2*t11)+t15*(0.2889626582E2+
0.950205321859E1*t11)+(-0.2894451105E2*t7*t20-0.9502053219E1*t15*t20)*(rho-old)
/24.0);
  }
}

#include <math.h>
double pressure(double rho,double phi)
{
  double t1;
  double t10;
  double t11;
  double t2;
  double t4;
  double t5;
  double t6;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t6 = t1*phi;
    t7 = 10.0*t6;
    t10 = log(rho);
    t11 = rho*t10;
    return((t4-t5+t7)*(-0.2094881058E2+0.4291934897E2*rho-0.2894451105E2*t11+
rho*(-0.1397483792E2+0.289445110522E2*t10))+(1.0-t4+t5-t7)*(-0.6490446091E-1+
0.3405607049E1*rho-0.9502053219E1*t11+rho*(0.609644617E1+0.950205321859E1*t10))
-2.0/delta*A*(t2+(2.0*alpha-2.0)*t6+(-3.0*alpha+1.0)*t1+alpha));
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t1;
  double t2;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t9 = sqrt(0.116654747E11*t2*phi-0.2916368675E11*t2+0.1944245783E11*t1*phi+
950205322.0);
    return(0.1E-3*t9);
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t15;
  double t16;
  double t17;
  double t18;
  double t19;
  double t21;
  double t22;
  double t26;
  double t3;
  double t4;
  double t41;
  double t43;
  double t44;
  double t7;
  {
    t3 = phi*phi;
    t4 = t3*t3;
    t7 = t3*phi;
    t15 = t4*phi;
    t16 = 6.0*t15;
    t17 = 15.0*t4;
    t18 = 10.0*t7;
    t19 = t16-t17+t18;
    t21 = log(rho);
    t22 = rho*t21;
    t26 = 1.0-t16+t17-t18;
    t41 = exp((-0.2889626582E2+0.1204277045E3*t15-0.3010692614E3*t4+
0.2007128409E3*t7)/(0.9502053219E1+0.116654747E3*t15-0.2916368674E3*t4+
0.1944245783E3*t7));
    t43 = log(t41);
    t44 = t41*t43;
    return(2.0/delta*A*(t4+(2.0*alpha-2.0)*t7+(-3.0*alpha+1.0)*t3+alpha)+t19*(
0.2133795654E2-0.2011952932E2*rho+0.289445110522E2*t22)+t26*(0.4540504209+
0.193942126E2*rho+0.950205321859E1*t22)-beta*(t19*(0.2133795654E2
-0.2011952932E2*t41+0.289445110522E2*t44)+t26*(0.4540504209+0.193942126E2*t41+
0.950205321859E1*t44)));
  }
}

#include <math.h>
double reactionSource(double rho,double phi,double old)
{
  double t100;
  double t107;
  double t114;
  double t116;
  double t120;
  double t121;
  double t122;
  double t129;
  double t131;
  double t132;
  double t134;
  double t135;
  double t137;
  double t139;
  double t157;
  double t172;
  double t174;
  double t175;
  double t178;
  double t180;
  double t182;
  double t184;
  double t2;
  double t20;
  double t24;
  double t26;
  double t27;
  double t29;
  double t31;
  double t34;
  double t36;
  double t40;
  double t44;
  double t45;
  double t47;
  double t49;
  double t5;
  double t50;
  double t52;
  double t54;
  double t55;
  double t56;
  double t57;
  double t6;
  double t61;
  double t63;
  double t64;
  double t65;
  double t69;
  double t7;
  double t71;
  double t72;
  double t74;
  double t76;
  double t80;
  double t82;
  double t85;
  double t98;
  {
    t2 = 1/delta*A;
    t5 = 0.5*phi+0.5*old;
    t6 = t5*t5;
    t7 = t6*t5;
    t20 = t6*t6;
    t24 = 30.0*t20-60.0*t7+30.0*t6;
    t26 = log(rho);
    t27 = rho*t26;
    t29 = 0.2133795654E2-0.2011952932E2*rho+0.289445110522E2*t27;
    t31 = -t24;
    t34 = 0.4540504209+0.193942126E2*rho+0.950205321859E1*t27;
    t36 = t20*t5;
    t40 = -0.2889626582E2+0.1204277045E3*t36-0.3010692614E3*t20+0.2007128409E3*
t7;
    t44 = 0.9502053219E1+0.116654747E3*t36-0.2916368674E3*t20+0.1944245783E3*t7
;
    t45 = 1/t44;
    t47 = exp(t40*t45);
    t49 = log(t47);
    t50 = t47*t49;
    t52 = 0.2133795654E2-0.2011952932E2*t47+0.289445110522E2*t50;
    t54 = 6.0*t36;
    t55 = 15.0*t20;
    t56 = 10.0*t7;
    t57 = t54-t55+t56;
    t61 = 0.6021385225E3*t20-0.1204277046E4*t7+0.6021385227E3*t6;
    t63 = t44*t44;
    t64 = 1/t63;
    t65 = t40*t64;
    t69 = 0.583273735E3*t20-0.116654747E4*t7+0.5832737349E3*t6;
    t71 = t61*t45-t65*t69;
    t72 = t71*t47;
    t74 = t72*t49;
    t76 = 0.882498173E1*t72+0.289445110522E2*t74;
    t80 = 0.4540504209+0.193942126E2*t47+0.950205321859E1*t50;
    t82 = 1.0-t54+t55-t56;
    t85 = 0.2889626582E2*t72+0.950205321859E1*t74;
    t98 = 360.0*t6-0.18E3*phi-0.18E3*old+60.0;
    t100 = -t98;
    t107 = 120.0*t7-180.0*t6+0.3E2*phi+0.3E2*old;
    t114 = 0.240855409E4*t7-0.3612831138E4*t6+0.6021385225E3*phi+0.6021385225E3
*old;
    t116 = t61*t64;
    t120 = 1/t63/t44;
    t121 = t40*t120;
    t122 = t69*t69;
    t129 = 0.233309494E4*t7-0.349964241E4*t6+0.583273735E3*phi+0.583273735E3*
old;
    t131 = t114*t45-2.0*t116*t69+2.0*t121*t122-t65*t129;
    t132 = t131*t47;
    t134 = t71*t71;
    t135 = t134*t47;
    t137 = t132*t49;
    t139 = t135*t49;
    t157 = t63*t63;
    t172 = ((0.722566227E4*t6-0.3612831138E4*phi-0.3612831138E4*old+
0.1204277045E4)*t45-3.0*t114*t64*t69+6.0*t61*t120*t122-3.0*t116*t129-6.0*t40/
t157*t122*t69+6.0*t121*t69*t129-t65*(0.699928482E4*t6-0.349964241E4*phi
-0.349964241E4*old+0.116654747E4))*t47;
    t174 = t131*t71;
    t175 = t174*t47;
    t178 = t134*t71*t47;
    t180 = t172*t49;
    t182 = t174*t50;
    t184 = t178*t49;
    return(2.0*t2*(4.0*t7+3.0*(2.0*alpha-2.0)*t6+2.0*(-3.0*alpha+1.0)*t5)+t24*
t29+t31*t34-beta*(t24*t52+t57*t76+t31*t80+t82*t85)+(2.0*t2*(0.12E2*phi+0.12E2*
old+12.0*alpha-12.0)+t98*t29+t100*t34-beta*(t98*t52+3.0*t107*t76+3.0*t24*(
0.882498173E1*t132+0.3776949278E2*t135+0.289445110522E2*t137+0.289445110522E2*
t139)+t57*(0.882498173E1*t172+0.1133084783E3*t175+0.6671400383E2*t178+
0.289445110522E2*t180+0.8683353315E2*t182+0.289445110522E2*t184)+t100*t80-3.0*
t107*t85+3.0*t31*(0.2889626582E2*t132+0.3839831904E2*t135+0.950205321859E1*t137
+0.950205321859E1*t139)+t82*(0.2889626582E2*t172+0.1151949571E3*t175+
0.4790037226E2*t178+0.950205321859E1*t180+0.2850615966E2*t182+0.950205321859E1*
t184)))*(phi-old)/24.0);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi,double old,double mid)
{
  double t100;
  double t105;
  double t107;
  double t109;
  double t112;
  double t113;
  double t114;
  double t120;
  double t122;
  double t123;
  double t125;
  double t126;
  double t128;
  double t130;
  double t132;
  double t136;
  double t138;
  double t141;
  double t143;
  double t146;
  double t148;
  double t153;
  double t16;
  double t160;
  double t162;
  double t168;
  double t173;
  double t180;
  double t187;
  double t191;
  double t195;
  double t198;
  double t2;
  double t200;
  double t201;
  double t203;
  double t204;
  double t206;
  double t208;
  double t21;
  double t210;
  double t216;
  double t218;
  double t220;
  double t223;
  double t228;
  double t23;
  double t231;
  double t232;
  double t233;
  double t234;
  double t237;
  double t24;
  double t241;
  double t246;
  double t248;
  double t249;
  double t251;
  double t252;
  double t254;
  double t256;
  double t257;
  double t259;
  double t26;
  double t261;
  double t263;
  double t264;
  double t266;
  double t274;
  double t28;
  double t282;
  double t285;
  double t289;
  double t290;
  double t291;
  double t293;
  double t294;
  double t296;
  double t297;
  double t299;
  double t300;
  double t302;
  double t304;
  double t306;
  double t308;
  double t31;
  double t33;
  double t34;
  double t368;
  double t370;
  double t372;
  double t373;
  double t375;
  double t376;
  double t378;
  double t379;
  double t38;
  double t381;
  double t383;
  double t385;
  double t386;
  double t388;
  double t390;
  double t392;
  double t394;
  double t396;
  double t398;
  double t400;
  double t402;
  double t405;
  double t410;
  double t417;
  double t42;
  double t43;
  double t437;
  double t45;
  double t453;
  double t455;
  double t47;
  double t48;
  double t5;
  double t50;
  double t55;
  double t59;
  double t6;
  double t61;
  double t62;
  double t63;
  double t67;
  double t69;
  double t70;
  double t72;
  double t74;
  double t79;
  double t83;
  double t86;
  double t88;
  double t90;
  double t91;
  double t93;
  double t95;
  double t97;
  double t98;
  double t99;
  {
    t2 = 1/delta*A;
    t5 = 0.5*phi+0.5*old;
    t6 = t5*t5;
    t16 = t6*t5;
    t21 = 0.6E2*t16-0.9E2*t6+0.15E2*phi+0.15E2*old;
    t23 = log(rho);
    t24 = rho*t23;
    t26 = 0.2133795654E2-0.2011952932E2*rho+0.289445110522E2*t24;
    t28 = -t21;
    t31 = 0.4540504209+0.193942126E2*rho+0.950205321859E1*t24;
    t33 = t6*t6;
    t34 = t33*t5;
    t38 = -0.2889626582E2+0.1204277045E3*t34-0.3010692614E3*t33+0.2007128409E3*
t16;
    t42 = 0.9502053219E1+0.116654747E3*t34-0.2916368674E3*t33+0.1944245783E3*
t16;
    t43 = 1/t42;
    t45 = exp(t38*t43);
    t47 = log(t45);
    t48 = t45*t47;
    t50 = 0.2133795654E2-0.2011952932E2*t45+0.289445110522E2*t48;
    t55 = 30.0*t33-60.0*t16+30.0*t6;
    t59 = 0.3010692612E3*t33-0.6021385228E3*t16+0.3010692614E3*t6;
    t61 = t42*t42;
    t62 = 1/t61;
    t63 = t38*t62;
    t67 = 0.2916368675E3*t33-0.5832737348E3*t16+0.2916368674E3*t6;
    t69 = t59*t43-t63*t67;
    t70 = t69*t45;
    t72 = t70*t47;
    t74 = 0.882498173E1*t70+0.289445110522E2*t72;
    t79 = 0.15E2*t33-0.3E2*t16+0.15E2*t6;
    t83 = 0.6021385225E3*t33-0.1204277046E4*t16+0.6021385227E3*t6;
    t86 = 0.116654747E4*t16;
    t88 = 0.583273735E3*t33-t86+0.5832737349E3*t6;
    t90 = t83*t43-t63*t88;
    t91 = t90*t45;
    t93 = t91*t47;
    t95 = 0.882498173E1*t91+0.289445110522E2*t93;
    t97 = 6.0*t34;
    t98 = 15.0*t33;
    t99 = 10.0*t16;
    t100 = t97-t98+t99;
    t105 = 0.1204277045E4*t16-0.1806415569E4*t6+0.3010692614E3*phi+
0.3010692614E3*old;
    t107 = t83*t62;
    t109 = t59*t62;
    t112 = 1/t61/t42;
    t113 = t38*t112;
    t114 = t88*t67;
    t120 = t86-0.1749821205E4*t6+0.2916368674E3*phi+0.2916368674E3*old;
    t122 = t105*t43-t107*t67-t109*t88+2.0*t113*t114-t63*t120;
    t123 = t122*t45;
    t125 = t90*t69;
    t126 = t125*t45;
    t128 = t123*t47;
    t130 = t125*t48;
    t132 = 0.882498173E1*t123+0.3776949278E2*t126+0.289445110522E2*t128+
0.289445110522E2*t130;
    t136 = 0.4540504209+0.193942126E2*t45+0.950205321859E1*t48;
    t138 = -t55;
    t141 = 0.2889626582E2*t70+0.950205321859E1*t72;
    t143 = -t79;
    t146 = 0.2889626582E2*t91+0.950205321859E1*t93;
    t148 = 1.0-t97+t98-t99;
    t153 = 0.2889626582E2*t123+0.3839831904E2*t126+0.950205321859E1*t128+
0.950205321859E1*t130;
    t160 = 0.18E3*phi+0.18E3*old-0.18E3;
    t162 = -t160;
    t168 = 360.0*t6-0.18E3*phi-0.18E3*old+60.0;
    t173 = 0.18E3*t6-0.9E2*phi-0.9E2*old+0.3E2;
    t180 = 120.0*t16-180.0*t6+0.3E2*phi+0.3E2*old;
    t187 = 0.240855409E4*t16-0.3612831138E4*t6+0.6021385225E3*phi+
0.6021385225E3*old;
    t191 = t88*t88;
    t195 = 0.349964241E4*t6;
    t198 = 0.233309494E4*t16-t195+0.583273735E3*phi+0.583273735E3*old;
    t200 = t187*t43-2.0*t107*t88+2.0*t113*t191-t63*t198;
    t201 = t200*t45;
    t203 = t90*t90;
    t204 = t203*t45;
    t206 = t201*t47;
    t208 = t204*t47;
    t210 = 0.882498173E1*t201+0.3776949278E2*t204+0.289445110522E2*t206+
0.289445110522E2*t208;
    t216 = 0.3612831135E4*t6-0.1806415569E4*phi-0.1806415569E4*old+
0.6021385225E3;
    t218 = t187*t62;
    t220 = t105*t62;
    t223 = t83*t112;
    t228 = t59*t112;
    t231 = t61*t61;
    t232 = 1/t231;
    t233 = t38*t232;
    t234 = t191*t67;
    t237 = t88*t120;
    t241 = t198*t67;
    t246 = t195-0.1749821205E4*phi-0.1749821205E4*old+0.583273735E3;
    t248 = t216*t43-t218*t67-2.0*t220*t88+4.0*t223*t114-2.0*t107*t120+2.0*t228*
t191-6.0*t233*t234+4.0*t113*t237-t109*t198+2.0*t113*t241-t63*t246;
    t249 = t248*t45;
    t251 = t200*t69;
    t252 = t251*t45;
    t254 = t91*t122;
    t256 = t203*t69;
    t257 = t256*t45;
    t259 = t249*t47;
    t261 = t251*t48;
    t263 = t47*t122;
    t264 = t91*t263;
    t266 = t256*t48;
    t274 = 0.722566227E4*t6-0.3612831138E4*phi-0.3612831138E4*old+
0.1204277045E4;
    t282 = t191*t88;
    t285 = t88*t198;
    t289 = 0.349964241E4*phi;
    t290 = 0.349964241E4*old;
    t291 = 0.699928482E4*t6-t289-t290+0.116654747E4;
    t293 = t274*t43-3.0*t218*t88+6.0*t223*t191-3.0*t107*t198-6.0*t233*t282+6.0*
t113*t285-t63*t291;
    t294 = t293*t45;
    t296 = t200*t90;
    t297 = t296*t45;
    t299 = t203*t90;
    t300 = t299*t45;
    t302 = t294*t47;
    t304 = t296*t48;
    t306 = t300*t47;
    t308 = 0.882498173E1*t294+0.1133084783E3*t297+0.6671400383E2*t300+
0.289445110522E2*t302+0.8683353315E2*t304+0.289445110522E2*t306;
    t368 = -3.0*t107*t246-6.0*t59*t232*t282+24.0*t38/t231/t42*t282*t67-18.0*
t233*t191*t120+6.0*t228*t285-18.0*t233*t285*t67+6.0*t113*t120*t198+6.0*t113*t88
*t246-t109*t291+2.0*t113*t291*t67-t63*(-0.349964241E4+t289+t290);
    t370 = ((-0.3612831138E4+0.3612831135E4*phi+0.3612831135E4*old)*t43-t274*
t62*t67-3.0*t216*t62*t88+6.0*t187*t112*t114-3.0*t218*t120+6.0*t105*t112*t191
-18.0*t83*t232*t234+12.0*t223*t237-3.0*t220*t198+6.0*t223*t241+t368)*t45;
    t372 = t293*t69;
    t373 = t372*t45;
    t375 = t248*t90;
    t376 = t375*t45;
    t378 = t200*t122;
    t379 = t378*t45;
    t381 = t296*t70;
    t383 = t204*t122;
    t385 = t299*t69;
    t386 = t385*t45;
    t388 = t370*t47;
    t390 = t372*t48;
    t392 = t375*t48;
    t394 = t378*t48;
    t396 = t296*t72;
    t398 = t204*t263;
    t400 = t385*t48;
    t402 = 0.882498173E1*t370+0.3776949278E2*t373+0.1133084783E3*t376+
0.1133084783E3*t379+0.2001420114E3*t381+0.2001420115E3*t383+0.9565851488E2*t386
+0.289445110522E2*t388+0.289445110522E2*t390+0.8683353315E2*t392+0.8683353315E2
*t394+0.8683353315E2*t396+0.8683353315E2*t398+0.289445110522E2*t400;
    t405 = -t168;
    t410 = -t180;
    t417 = 0.2889626582E2*t201+0.3839831904E2*t204+0.950205321859E1*t206+
0.950205321859E1*t208;
    t437 = 0.2889626582E2*t294+0.1151949571E3*t297+0.4790037226E2*t300+
0.950205321859E1*t302+0.2850615966E2*t304+0.950205321859E1*t306;
    t453 = 0.2889626582E2*t370+0.3839831904E2*t373+0.1151949571E3*t376+
0.1151949571E3*t379+0.1437011168E3*t381+0.1437011168E3*t383+0.5740242548E2*t386
+0.950205321859E1*t388+0.950205321859E1*t390+0.2850615966E2*t392+0.2850615966E2
*t394+0.2850615966E2*t396+0.2850615966E2*t398+0.950205321859E1*t400;
    t455 = t160*t50+t168*t74+3.0*t173*t95+3.0*t180*t132+3.0*t21*t210+3.0*t55*(
0.882498173E1*t249+0.3776949278E2*t252+0.7553898556E2*t254+0.6671400383E2*t257+
0.289445110522E2*t259+0.289445110522E2*t261+0.578890221E2*t264+0.289445110522E2
*t266)+t79*t308+t100*t402+t162*t136+t405*t141-3.0*t173*t146+3.0*t410*t153+3.0*
t28*t417+3.0*t138*(0.2889626582E2*t249+0.3839831904E2*t252+0.7679663808E2*t254+
0.4790037226E2*t257+0.950205321859E1*t259+0.950205321859E1*t261+0.1900410644E2*
t264+0.950205321859E1*t266)+t143*t437+t148*t453;
    return(2.0*t2*(0.6E1*t6+0.3E1*(2.0*alpha-2.0)*t5-0.3E1*alpha+0.1E1)+t21*t26
+t28*t31-beta*(t21*t50+t55*t74+t79*t95+t100*t132+t28*t136+t138*t141+t143*t146+
t148*t153)+(0.24E2*t2+t160*t26+t162*t31-beta*t455)*(phi-old)/24.0+t2*(0.12E2*
phi+0.12E2*old+12.0*alpha-12.0)/12.0+t168*t26/24.0+t405*t31/24.0-beta*(t168*t50
+3.0*t180*t95+3.0*t55*t210+t100*t308+t405*t136+3.0*t410*t146+3.0*t138*t417+t148
*t437)/24.0);
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t11;
  double t12;
  double t16;
  double t2;
  double t20;
  double t21;
  double t4;
  double t5;
  double t7;
  double t8;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t8 = t4-t5+t7;
    t11 = 0.5*rho+0.5*old;
    t12 = log(t11);
    t16 = 1.0-t4+t5-t7;
    t20 = t11*t11;
    t21 = 1/t20;
    return(t8*(0.882498173E1+0.289445110522E2*t12)+t16*(0.2889626582E2+
0.950205321859E1*t12)+(-0.2894451105E2*t8*t21-0.9502053219E1*t16*t21)*(rho-old)
/24.0);
  }
}

double drhochemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t11;
  double t12;
  double t15;
  double t18;
  double t2;
  double t20;
  double t29;
  double t4;
  double t5;
  double t7;
  double t8;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t8 = t4-t5+t7;
    t11 = 0.5*rho+0.5*old;
    t12 = 1/t11;
    t15 = 1.0-t4+t5-t7;
    t18 = t11*t11;
    t20 = 1/t18/t11;
    t29 = 1/t18;
    return(0.1447225552E2*t8*t12+0.475102661E1*t15*t12+(0.2894451105E2*t8*t20+
0.9502053219E1*t15*t20)*(rho-old)/24.0-0.1206021294E1*t8*t29-0.3959188841*t15*
t29);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t10;
  double t11;
  double t15;
  double t19;
  double t2;
  double t20;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t7 = 30.0*t2-60.0*t1*phi+30.0*t1;
    t10 = 0.5*rho+0.5*old;
    t11 = log(t10);
    t15 = -t7;
    t19 = t10*t10;
    t20 = 1/t19;
    return(t7*(0.882498173E1+0.289445110522E2*t11)+t15*(0.2889626582E2+
0.950205321859E1*t11)+(-0.2894451105E2*t7*t20-0.9502053219E1*t15*t20)*(rho-old)
/24.0);
  }
}

#include <math.h>
double pressure(double rho,double phi)
{
  double t1;
  double t10;
  double t11;
  double t2;
  double t4;
  double t5;
  double t6;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t6 = t1*phi;
    t7 = 10.0*t6;
    t10 = log(rho);
    t11 = rho*t10;
    return((t4-t5+t7)*(-0.2094881058E2+0.4291934897E2*rho-0.2894451105E2*t11+
rho*(-0.1397483792E2+0.289445110522E2*t10))+(1.0-t4+t5-t7)*(-0.6490446091E-1+
0.3405607049E1*rho-0.9502053219E1*t11+rho*(0.609644617E1+0.950205321859E1*t10))
-2.0/delta*A*(t2+(2.0*alpha-2.0)*t6+(-3.0*alpha+1.0)*t1+alpha));
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t1;
  double t2;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t9 = sqrt(0.116654747E11*t2*phi-0.2916368675E11*t2+0.1944245783E11*t1*phi+
950205322.0);
    return(0.1E-3*t9);
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t15;
  double t16;
  double t17;
  double t18;
  double t19;
  double t21;
  double t22;
  double t26;
  double t3;
  double t4;
  double t41;
  double t43;
  double t44;
  double t7;
  {
    t3 = phi*phi;
    t4 = t3*t3;
    t7 = t3*phi;
    t15 = t4*phi;
    t16 = 6.0*t15;
    t17 = 15.0*t4;
    t18 = 10.0*t7;
    t19 = t16-t17+t18;
    t21 = log(rho);
    t22 = rho*t21;
    t26 = 1.0-t16+t17-t18;
    t41 = exp((-0.2889626582E2+0.1204277045E3*t15-0.3010692614E3*t4+
0.2007128409E3*t7)/(0.9502053219E1+0.116654747E3*t15-0.2916368674E3*t4+
0.1944245783E3*t7));
    t43 = log(t41);
    t44 = t41*t43;
    return(2.0/delta*A*(t4+(2.0*alpha-2.0)*t7+(-3.0*alpha+1.0)*t3+alpha)+t19*(
-0.2011952932E2*rho+0.289445110522E2*t22+0.2133795654E2)+t26*(0.193942126E2*rho
+0.950205321859E1*t22+0.4540504209)-beta*(t19*(0.2133795654E2-0.2011952932E2*
t41+0.289445110522E2*t44)+t26*(0.4540504209+0.193942126E2*t41+0.950205321859E1*
t44)));
  }
}

#include <math.h>
double reactionSource(double rho,double phi,double old)
{
  double t100;
  double t107;
  double t114;
  double t116;
  double t120;
  double t121;
  double t122;
  double t129;
  double t131;
  double t132;
  double t134;
  double t135;
  double t137;
  double t139;
  double t157;
  double t172;
  double t174;
  double t175;
  double t178;
  double t180;
  double t182;
  double t184;
  double t2;
  double t20;
  double t24;
  double t26;
  double t27;
  double t29;
  double t31;
  double t34;
  double t36;
  double t40;
  double t44;
  double t45;
  double t47;
  double t49;
  double t5;
  double t50;
  double t52;
  double t54;
  double t55;
  double t56;
  double t57;
  double t6;
  double t61;
  double t63;
  double t64;
  double t65;
  double t69;
  double t7;
  double t71;
  double t72;
  double t74;
  double t76;
  double t80;
  double t82;
  double t85;
  double t98;
  {
    t2 = 1/delta*A;
    t5 = 0.5*phi+0.5*old;
    t6 = t5*t5;
    t7 = t6*t5;
    t20 = t6*t6;
    t24 = 30.0*t20-60.0*t7+30.0*t6;
    t26 = log(rho);
    t27 = rho*t26;
    t29 = -0.2011952932E2*rho+0.289445110522E2*t27+0.2133795654E2;
    t31 = -t24;
    t34 = 0.193942126E2*rho+0.950205321859E1*t27+0.4540504209;
    t36 = t20*t5;
    t40 = -0.2889626582E2+0.1204277045E3*t36-0.3010692614E3*t20+0.2007128409E3*
t7;
    t44 = 0.9502053219E1+0.116654747E3*t36-0.2916368674E3*t20+0.1944245783E3*t7
;
    t45 = 1/t44;
    t47 = exp(t40*t45);
    t49 = log(t47);
    t50 = t47*t49;
    t52 = 0.2133795654E2-0.2011952932E2*t47+0.289445110522E2*t50;
    t54 = 6.0*t36;
    t55 = 15.0*t20;
    t56 = 10.0*t7;
    t57 = t54-t55+t56;
    t61 = 0.6021385225E3*t20-0.1204277046E4*t7+0.6021385227E3*t6;
    t63 = t44*t44;
    t64 = 1/t63;
    t65 = t40*t64;
    t69 = 0.583273735E3*t20-0.116654747E4*t7+0.5832737349E3*t6;
    t71 = t61*t45-t65*t69;
    t72 = t71*t47;
    t74 = t72*t49;
    t76 = 0.882498173E1*t72+0.289445110522E2*t74;
    t80 = 0.4540504209+0.193942126E2*t47+0.950205321859E1*t50;
    t82 = 1.0-t54+t55-t56;
    t85 = 0.2889626582E2*t72+0.950205321859E1*t74;
    t98 = 360.0*t6-0.18E3*phi-0.18E3*old+60.0;
    t100 = -t98;
    t107 = 120.0*t7-180.0*t6+0.3E2*phi+0.3E2*old;
    t114 = 0.240855409E4*t7-0.3612831138E4*t6+0.6021385225E3*phi+0.6021385225E3
*old;
    t116 = t61*t64;
    t120 = 1/t63/t44;
    t121 = t40*t120;
    t122 = t69*t69;
    t129 = 0.233309494E4*t7-0.349964241E4*t6+0.583273735E3*phi+0.583273735E3*
old;
    t131 = t114*t45-2.0*t116*t69+2.0*t121*t122-t65*t129;
    t132 = t131*t47;
    t134 = t71*t71;
    t135 = t134*t47;
    t137 = t132*t49;
    t139 = t135*t49;
    t157 = t63*t63;
    t172 = ((0.722566227E4*t6-0.3612831138E4*phi-0.3612831138E4*old+
0.1204277045E4)*t45-3.0*t114*t64*t69+6.0*t61*t120*t122-3.0*t116*t129-6.0*t40/
t157*t122*t69+6.0*t121*t69*t129-t65*(0.699928482E4*t6-0.349964241E4*phi
-0.349964241E4*old+0.116654747E4))*t47;
    t174 = t131*t71;
    t175 = t174*t47;
    t178 = t134*t71*t47;
    t180 = t172*t49;
    t182 = t174*t50;
    t184 = t178*t49;
    return(2.0*t2*(4.0*t7+3.0*(2.0*alpha-2.0)*t6+2.0*(-3.0*alpha+1.0)*t5)+t24*
t29+t31*t34-beta*(t24*t52+t57*t76+t31*t80+t82*t85)+(2.0*t2*(0.12E2*phi+0.12E2*
old+12.0*alpha-12.0)+t98*t29+t100*t34-beta*(t98*t52+3.0*t107*t76+3.0*t24*(
0.882498173E1*t132+0.3776949278E2*t135+0.289445110522E2*t137+0.289445110522E2*
t139)+t57*(0.882498173E1*t172+0.1133084783E3*t175+0.6671400383E2*t178+
0.289445110522E2*t180+0.8683353315E2*t182+0.289445110522E2*t184)+t100*t80-3.0*
t107*t85+3.0*t31*(0.2889626582E2*t132+0.3839831904E2*t135+0.950205321859E1*t137
+0.950205321859E1*t139)+t82*(0.2889626582E2*t172+0.1151949571E3*t175+
0.4790037226E2*t178+0.950205321859E1*t180+0.2850615966E2*t182+0.950205321859E1*
t184)))*(phi-old)/24.0);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi,double old,double mid)
{
  double t100;
  double t105;
  double t107;
  double t109;
  double t112;
  double t113;
  double t114;
  double t120;
  double t122;
  double t123;
  double t125;
  double t126;
  double t128;
  double t130;
  double t132;
  double t136;
  double t138;
  double t141;
  double t143;
  double t146;
  double t148;
  double t153;
  double t16;
  double t160;
  double t162;
  double t168;
  double t173;
  double t180;
  double t187;
  double t191;
  double t195;
  double t198;
  double t2;
  double t200;
  double t201;
  double t203;
  double t204;
  double t206;
  double t208;
  double t21;
  double t210;
  double t216;
  double t218;
  double t220;
  double t223;
  double t228;
  double t23;
  double t231;
  double t232;
  double t233;
  double t234;
  double t237;
  double t24;
  double t241;
  double t246;
  double t248;
  double t249;
  double t251;
  double t252;
  double t254;
  double t256;
  double t257;
  double t259;
  double t26;
  double t261;
  double t263;
  double t264;
  double t266;
  double t274;
  double t28;
  double t282;
  double t285;
  double t289;
  double t290;
  double t291;
  double t293;
  double t294;
  double t296;
  double t297;
  double t299;
  double t300;
  double t302;
  double t304;
  double t306;
  double t308;
  double t31;
  double t33;
  double t34;
  double t368;
  double t370;
  double t372;
  double t373;
  double t375;
  double t376;
  double t378;
  double t379;
  double t38;
  double t381;
  double t383;
  double t385;
  double t386;
  double t388;
  double t390;
  double t392;
  double t394;
  double t396;
  double t398;
  double t400;
  double t402;
  double t405;
  double t410;
  double t417;
  double t42;
  double t43;
  double t437;
  double t45;
  double t453;
  double t455;
  double t47;
  double t48;
  double t5;
  double t50;
  double t55;
  double t59;
  double t6;
  double t61;
  double t62;
  double t63;
  double t67;
  double t69;
  double t70;
  double t72;
  double t74;
  double t79;
  double t83;
  double t86;
  double t88;
  double t90;
  double t91;
  double t93;
  double t95;
  double t97;
  double t98;
  double t99;
  {
    t2 = 1/delta*A;
    t5 = 0.5*phi+0.5*old;
    t6 = t5*t5;
    t16 = t6*t5;
    t21 = 0.6E2*t16-0.9E2*t6+0.15E2*phi+0.15E2*old;
    t23 = log(rho);
    t24 = rho*t23;
    t26 = -0.2011952932E2*rho+0.289445110522E2*t24+0.2133795654E2;
    t28 = -t21;
    t31 = 0.193942126E2*rho+0.950205321859E1*t24+0.4540504209;
    t33 = t6*t6;
    t34 = t33*t5;
    t38 = -0.2889626582E2+0.1204277045E3*t34-0.3010692614E3*t33+0.2007128409E3*
t16;
    t42 = 0.9502053219E1+0.116654747E3*t34-0.2916368674E3*t33+0.1944245783E3*
t16;
    t43 = 1/t42;
    t45 = exp(t38*t43);
    t47 = log(t45);
    t48 = t45*t47;
    t50 = 0.2133795654E2-0.2011952932E2*t45+0.289445110522E2*t48;
    t55 = 30.0*t33-60.0*t16+30.0*t6;
    t59 = 0.3010692612E3*t33-0.6021385228E3*t16+0.3010692614E3*t6;
    t61 = t42*t42;
    t62 = 1/t61;
    t63 = t38*t62;
    t67 = 0.2916368675E3*t33-0.5832737348E3*t16+0.2916368674E3*t6;
    t69 = t59*t43-t63*t67;
    t70 = t69*t45;
    t72 = t70*t47;
    t74 = 0.882498173E1*t70+0.289445110522E2*t72;
    t79 = 0.15E2*t33-0.3E2*t16+0.15E2*t6;
    t83 = 0.6021385225E3*t33-0.1204277046E4*t16+0.6021385227E3*t6;
    t86 = 0.116654747E4*t16;
    t88 = 0.583273735E3*t33-t86+0.5832737349E3*t6;
    t90 = t83*t43-t63*t88;
    t91 = t90*t45;
    t93 = t91*t47;
    t95 = 0.882498173E1*t91+0.289445110522E2*t93;
    t97 = 6.0*t34;
    t98 = 15.0*t33;
    t99 = 10.0*t16;
    t100 = t97-t98+t99;
    t105 = 0.1204277045E4*t16-0.1806415569E4*t6+0.3010692614E3*phi+
0.3010692614E3*old;
    t107 = t83*t62;
    t109 = t59*t62;
    t112 = 1/t61/t42;
    t113 = t38*t112;
    t114 = t88*t67;
    t120 = t86-0.1749821205E4*t6+0.2916368674E3*phi+0.2916368674E3*old;
    t122 = t105*t43-t107*t67-t109*t88+2.0*t113*t114-t63*t120;
    t123 = t122*t45;
    t125 = t90*t69;
    t126 = t125*t45;
    t128 = t123*t47;
    t130 = t125*t48;
    t132 = 0.882498173E1*t123+0.3776949278E2*t126+0.289445110522E2*t128+
0.289445110522E2*t130;
    t136 = 0.4540504209+0.193942126E2*t45+0.950205321859E1*t48;
    t138 = -t55;
    t141 = 0.2889626582E2*t70+0.950205321859E1*t72;
    t143 = -t79;
    t146 = 0.2889626582E2*t91+0.950205321859E1*t93;
    t148 = 1.0-t97+t98-t99;
    t153 = 0.2889626582E2*t123+0.3839831904E2*t126+0.950205321859E1*t128+
0.950205321859E1*t130;
    t160 = -0.18E3+0.18E3*phi+0.18E3*old;
    t162 = -t160;
    t168 = 360.0*t6-0.18E3*phi-0.18E3*old+60.0;
    t173 = 0.3E2+0.18E3*t6-0.9E2*phi-0.9E2*old;
    t180 = 120.0*t16-180.0*t6+0.3E2*phi+0.3E2*old;
    t187 = 0.240855409E4*t16-0.3612831138E4*t6+0.6021385225E3*phi+
0.6021385225E3*old;
    t191 = t88*t88;
    t195 = 0.349964241E4*t6;
    t198 = 0.233309494E4*t16-t195+0.583273735E3*phi+0.583273735E3*old;
    t200 = t187*t43-2.0*t107*t88+2.0*t113*t191-t63*t198;
    t201 = t200*t45;
    t203 = t90*t90;
    t204 = t203*t45;
    t206 = t201*t47;
    t208 = t204*t47;
    t210 = 0.882498173E1*t201+0.3776949278E2*t204+0.289445110522E2*t206+
0.289445110522E2*t208;
    t216 = 0.6021385225E3+0.3612831135E4*t6-0.1806415569E4*phi-0.1806415569E4*
old;
    t218 = t187*t62;
    t220 = t105*t62;
    t223 = t83*t112;
    t228 = t59*t112;
    t231 = t61*t61;
    t232 = 1/t231;
    t233 = t38*t232;
    t234 = t191*t67;
    t237 = t88*t120;
    t241 = t198*t67;
    t246 = 0.583273735E3+t195-0.1749821205E4*phi-0.1749821205E4*old;
    t248 = t216*t43-t218*t67-2.0*t220*t88+4.0*t223*t114-2.0*t107*t120+2.0*t228*
t191-6.0*t233*t234+4.0*t113*t237-t109*t198+2.0*t113*t241-t63*t246;
    t249 = t248*t45;
    t251 = t200*t69;
    t252 = t251*t45;
    t254 = t91*t122;
    t256 = t203*t69;
    t257 = t256*t45;
    t259 = t249*t47;
    t261 = t251*t48;
    t263 = t47*t122;
    t264 = t91*t263;
    t266 = t256*t48;
    t274 = 0.722566227E4*t6-0.3612831138E4*phi-0.3612831138E4*old+
0.1204277045E4;
    t282 = t191*t88;
    t285 = t88*t198;
    t289 = 0.349964241E4*phi;
    t290 = 0.349964241E4*old;
    t291 = 0.116654747E4+0.699928482E4*t6-t289-t290;
    t293 = t274*t43-3.0*t218*t88+6.0*t223*t191-3.0*t107*t198-6.0*t233*t282+6.0*
t113*t285-t63*t291;
    t294 = t293*t45;
    t296 = t200*t90;
    t297 = t296*t45;
    t299 = t203*t90;
    t300 = t299*t45;
    t302 = t294*t47;
    t304 = t296*t48;
    t306 = t300*t47;
    t308 = 0.882498173E1*t294+0.1133084783E3*t297+0.6671400383E2*t300+
0.289445110522E2*t302+0.8683353315E2*t304+0.289445110522E2*t306;
    t368 = -3.0*t107*t246-6.0*t59*t232*t282+24.0*t38/t231/t42*t282*t67-18.0*
t233*t191*t120+6.0*t228*t285-18.0*t233*t285*t67+6.0*t113*t120*t198+6.0*t113*t88
*t246-t109*t291+2.0*t113*t291*t67-t63*(-0.349964241E4+t289+t290);
    t370 = ((-0.3612831138E4+0.3612831135E4*phi+0.3612831135E4*old)*t43-t274*
t62*t67-3.0*t216*t62*t88+6.0*t187*t112*t114-3.0*t218*t120+6.0*t105*t112*t191
-18.0*t83*t232*t234+12.0*t223*t237-3.0*t220*t198+6.0*t223*t241+t368)*t45;
    t372 = t293*t69;
    t373 = t372*t45;
    t375 = t248*t90;
    t376 = t375*t45;
    t378 = t200*t122;
    t379 = t378*t45;
    t381 = t296*t70;
    t383 = t204*t122;
    t385 = t299*t69;
    t386 = t385*t45;
    t388 = t370*t47;
    t390 = t372*t48;
    t392 = t375*t48;
    t394 = t378*t48;
    t396 = t296*t72;
    t398 = t204*t263;
    t400 = t385*t48;
    t402 = 0.882498173E1*t370+0.3776949278E2*t373+0.1133084783E3*t376+
0.1133084783E3*t379+0.2001420114E3*t381+0.2001420115E3*t383+0.9565851488E2*t386
+0.289445110522E2*t388+0.289445110522E2*t390+0.8683353315E2*t392+0.8683353315E2
*t394+0.8683353315E2*t396+0.8683353315E2*t398+0.289445110522E2*t400;
    t405 = -t168;
    t410 = -t180;
    t417 = 0.2889626582E2*t201+0.3839831904E2*t204+0.950205321859E1*t206+
0.950205321859E1*t208;
    t437 = 0.2889626582E2*t294+0.1151949571E3*t297+0.4790037226E2*t300+
0.950205321859E1*t302+0.2850615966E2*t304+0.950205321859E1*t306;
    t453 = 0.2889626582E2*t370+0.3839831904E2*t373+0.1151949571E3*t376+
0.1151949571E3*t379+0.1437011168E3*t381+0.1437011168E3*t383+0.5740242548E2*t386
+0.950205321859E1*t388+0.950205321859E1*t390+0.2850615966E2*t392+0.2850615966E2
*t394+0.2850615966E2*t396+0.2850615966E2*t398+0.950205321859E1*t400;
    t455 = t160*t50+t168*t74+3.0*t173*t95+3.0*t180*t132+3.0*t21*t210+3.0*t55*(
0.882498173E1*t249+0.3776949278E2*t252+0.7553898556E2*t254+0.6671400383E2*t257+
0.289445110522E2*t259+0.289445110522E2*t261+0.578890221E2*t264+0.289445110522E2
*t266)+t79*t308+t100*t402+t162*t136+t405*t141-3.0*t173*t146+3.0*t410*t153+3.0*
t28*t417+3.0*t138*(0.2889626582E2*t249+0.3839831904E2*t252+0.7679663808E2*t254+
0.4790037226E2*t257+0.950205321859E1*t259+0.950205321859E1*t261+0.1900410644E2*
t264+0.950205321859E1*t266)+t143*t437+t148*t453;
    return(2.0*t2*(0.6E1*t6+0.3E1*(2.0*alpha-2.0)*t5-0.3E1*alpha+0.1E1)+t21*t26
+t28*t31-beta*(t21*t50+t55*t74+t79*t95+t100*t132+t28*t136+t138*t141+t143*t146+
t148*t153)+(0.24E2*t2+t160*t26+t162*t31-beta*t455)*(phi-old)/24.0+t2*(0.12E2*
phi+0.12E2*old+12.0*alpha-12.0)/12.0+t168*t26/24.0+t405*t31/24.0-beta*(t168*t50
+3.0*t180*t95+3.0*t55*t210+t100*t308+t405*t136+3.0*t410*t146+3.0*t138*t417+t148
*t437)/24.0);
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t11;
  double t12;
  double t16;
  double t2;
  double t20;
  double t21;
  double t4;
  double t5;
  double t7;
  double t8;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t8 = t4-t5+t7;
    t11 = 0.5*rho+0.5*old;
    t12 = log(t11);
    t16 = 1.0-t4+t5-t7;
    t20 = t11*t11;
    t21 = 1/t20;
    return(t8*(0.882498173E1+0.289445110522E2*t12)+t16*(0.2889626582E2+
0.950205321859E1*t12)+(-0.2894451105E2*t8*t21-0.9502053219E1*t16*t21)*(rho-old)
/24.0);
  }
}

double drhochemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t11;
  double t12;
  double t15;
  double t18;
  double t2;
  double t20;
  double t29;
  double t4;
  double t5;
  double t7;
  double t8;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t8 = t4-t5+t7;
    t11 = 0.5*rho+0.5*old;
    t12 = 1/t11;
    t15 = 1.0-t4+t5-t7;
    t18 = t11*t11;
    t20 = 1/t18/t11;
    t29 = 1/t18;
    return(0.1447225552E2*t8*t12+0.475102661E1*t15*t12+(0.2894451105E2*t8*t20+
0.9502053219E1*t15*t20)*(rho-old)/24.0-0.1206021294E1*t8*t29-0.3959188841*t15*
t29);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t10;
  double t11;
  double t15;
  double t19;
  double t2;
  double t20;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t7 = 30.0*t2-60.0*t1*phi+30.0*t1;
    t10 = 0.5*rho+0.5*old;
    t11 = log(t10);
    t15 = -t7;
    t19 = t10*t10;
    t20 = 1/t19;
    return(t7*(0.882498173E1+0.289445110522E2*t11)+t15*(0.2889626582E2+
0.950205321859E1*t11)+(-0.2894451105E2*t7*t20-0.9502053219E1*t15*t20)*(rho-old)
/24.0);
  }
}

#include <math.h>
double pressure(double rho,double phi)
{
  double t1;
  double t10;
  double t11;
  double t2;
  double t4;
  double t5;
  double t6;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t6 = t1*phi;
    t7 = 10.0*t6;
    t10 = log(rho);
    t11 = rho*t10;
    return((t4-t5+t7)*(-0.2094881058E2+0.4291934897E2*rho-0.2894451105E2*t11+
rho*(-0.1397483792E2+0.289445110522E2*t10))+(1.0-t4+t5-t7)*(-0.6490446091E-1+
0.3405607049E1*rho-0.9502053219E1*t11+rho*(0.609644617E1+0.950205321859E1*t10))
-2.0/delta*A*(t2+(2.0*alpha-2.0)*t6+(-3.0*alpha+1.0)*t1+alpha));
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t1;
  double t2;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t9 = sqrt(0.116654747E11*t2*phi-0.2916368675E11*t2+0.1944245783E11*t1*phi+
950205322.0);
    return(0.1E-3*t9);
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t15;
  double t16;
  double t17;
  double t18;
  double t19;
  double t21;
  double t22;
  double t26;
  double t3;
  double t4;
  double t41;
  double t43;
  double t44;
  double t7;
  {
    t3 = phi*phi;
    t4 = t3*t3;
    t7 = t3*phi;
    t15 = t4*phi;
    t16 = 6.0*t15;
    t17 = 15.0*t4;
    t18 = 10.0*t7;
    t19 = t16-t17+t18;
    t21 = log(rho);
    t22 = rho*t21;
    t26 = 1.0-t16+t17-t18;
    t41 = exp((-0.2889626582E2+0.1204277045E3*t15-0.3010692614E3*t4+
0.2007128409E3*t7)/(0.9502053219E1+0.116654747E3*t15-0.2916368674E3*t4+
0.1944245783E3*t7));
    t43 = log(t41);
    t44 = t41*t43;
    return(2.0/delta*A*(t4+(2.0*alpha-2.0)*t7+(-3.0*alpha+1.0)*t3+alpha)+t19*(
-0.2011952932E2*rho+0.289445110522E2*t22+0.2133795654E2)+t26*(0.193942126E2*rho
+0.950205321859E1*t22+0.4540504209)-beta*(t19*(0.2133795654E2-0.2011952932E2*
t41+0.289445110522E2*t44)+t26*(0.4540504209+0.193942126E2*t41+0.950205321859E1*
t44)));
  }
}

#include <math.h>
double reactionSource(double rho,double phi,double old)
{
  double t100;
  double t107;
  double t114;
  double t116;
  double t120;
  double t121;
  double t122;
  double t129;
  double t131;
  double t132;
  double t134;
  double t135;
  double t137;
  double t139;
  double t157;
  double t172;
  double t174;
  double t175;
  double t178;
  double t180;
  double t182;
  double t184;
  double t2;
  double t20;
  double t24;
  double t26;
  double t27;
  double t29;
  double t31;
  double t34;
  double t36;
  double t40;
  double t44;
  double t45;
  double t47;
  double t49;
  double t5;
  double t50;
  double t52;
  double t54;
  double t55;
  double t56;
  double t57;
  double t6;
  double t61;
  double t63;
  double t64;
  double t65;
  double t69;
  double t7;
  double t71;
  double t72;
  double t74;
  double t76;
  double t80;
  double t82;
  double t85;
  double t98;
  {
    t2 = 1/delta*A;
    t5 = 0.5*phi+0.5*old;
    t6 = t5*t5;
    t7 = t6*t5;
    t20 = t6*t6;
    t24 = 30.0*t20-60.0*t7+30.0*t6;
    t26 = log(rho);
    t27 = rho*t26;
    t29 = -0.2011952932E2*rho+0.289445110522E2*t27+0.2133795654E2;
    t31 = -t24;
    t34 = 0.193942126E2*rho+0.950205321859E1*t27+0.4540504209;
    t36 = t20*t5;
    t40 = -0.2889626582E2+0.1204277045E3*t36-0.3010692614E3*t20+0.2007128409E3*
t7;
    t44 = 0.9502053219E1+0.116654747E3*t36-0.2916368674E3*t20+0.1944245783E3*t7
;
    t45 = 1/t44;
    t47 = exp(t40*t45);
    t49 = log(t47);
    t50 = t47*t49;
    t52 = 0.2133795654E2-0.2011952932E2*t47+0.289445110522E2*t50;
    t54 = 6.0*t36;
    t55 = 15.0*t20;
    t56 = 10.0*t7;
    t57 = t54-t55+t56;
    t61 = 0.6021385225E3*t20-0.1204277046E4*t7+0.6021385227E3*t6;
    t63 = t44*t44;
    t64 = 1/t63;
    t65 = t40*t64;
    t69 = 0.583273735E3*t20-0.116654747E4*t7+0.5832737349E3*t6;
    t71 = t61*t45-t65*t69;
    t72 = t71*t47;
    t74 = t72*t49;
    t76 = 0.882498173E1*t72+0.289445110522E2*t74;
    t80 = 0.4540504209+0.193942126E2*t47+0.950205321859E1*t50;
    t82 = 1.0-t54+t55-t56;
    t85 = 0.2889626582E2*t72+0.950205321859E1*t74;
    t98 = 60.0+360.0*t6-0.18E3*phi-0.18E3*old;
    t100 = -t98;
    t107 = 120.0*t7-180.0*t6+0.3E2*phi+0.3E2*old;
    t114 = 0.240855409E4*t7-0.3612831138E4*t6+0.6021385225E3*phi+0.6021385225E3
*old;
    t116 = t61*t64;
    t120 = 1/t63/t44;
    t121 = t40*t120;
    t122 = t69*t69;
    t129 = 0.233309494E4*t7-0.349964241E4*t6+0.583273735E3*phi+0.583273735E3*
old;
    t131 = t114*t45-2.0*t116*t69+2.0*t122*t121-t65*t129;
    t132 = t131*t47;
    t134 = t71*t71;
    t135 = t134*t47;
    t137 = t132*t49;
    t139 = t135*t49;
    t157 = t63*t63;
    t172 = ((0.1204277045E4+0.722566227E4*t6-0.3612831138E4*phi-0.3612831138E4*
old)*t45-3.0*t114*t64*t69+6.0*t61*t120*t122-3.0*t116*t129-6.0*t40/t157*t122*t69
+6.0*t121*t69*t129-t65*(0.116654747E4+0.699928482E4*t6-0.349964241E4*phi
-0.349964241E4*old))*t47;
    t174 = t131*t71;
    t175 = t174*t47;
    t178 = t134*t71*t47;
    t180 = t172*t49;
    t182 = t174*t50;
    t184 = t178*t49;
    return(2.0*t2*(4.0*t7+3.0*(2.0*alpha-2.0)*t6+2.0*(-3.0*alpha+1.0)*t5)+t24*
t29+t31*t34-beta*(t24*t52+t57*t76+t31*t80+t82*t85)+(2.0*t2*(-12.0+0.12E2*phi+
0.12E2*old+12.0*alpha)+t98*t29+t100*t34-beta*(t98*t52+3.0*t107*t76+3.0*t24*(
0.882498173E1*t132+0.3776949278E2*t135+0.289445110522E2*t137+0.289445110522E2*
t139)+t57*(0.882498173E1*t172+0.1133084783E3*t175+0.6671400383E2*t178+
0.289445110522E2*t180+0.8683353315E2*t182+0.289445110522E2*t184)+t100*t80-3.0*
t107*t85+3.0*t31*(0.2889626582E2*t132+0.3839831904E2*t135+0.950205321859E1*t137
+0.950205321859E1*t139)+t82*(0.2889626582E2*t172+0.1151949571E3*t175+
0.4790037226E2*t178+0.950205321859E1*t180+0.2850615966E2*t182+0.950205321859E1*
t184)))*(phi-old)/24.0);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi,double old,double mid)
{
  double t100;
  double t105;
  double t107;
  double t109;
  double t112;
  double t113;
  double t114;
  double t120;
  double t122;
  double t123;
  double t125;
  double t126;
  double t128;
  double t130;
  double t132;
  double t136;
  double t138;
  double t141;
  double t143;
  double t146;
  double t148;
  double t153;
  double t16;
  double t160;
  double t162;
  double t168;
  double t173;
  double t180;
  double t187;
  double t191;
  double t195;
  double t198;
  double t2;
  double t200;
  double t201;
  double t203;
  double t204;
  double t206;
  double t208;
  double t21;
  double t210;
  double t216;
  double t218;
  double t220;
  double t223;
  double t228;
  double t23;
  double t231;
  double t232;
  double t233;
  double t234;
  double t237;
  double t24;
  double t241;
  double t246;
  double t248;
  double t249;
  double t251;
  double t252;
  double t254;
  double t256;
  double t257;
  double t259;
  double t26;
  double t261;
  double t263;
  double t264;
  double t266;
  double t274;
  double t28;
  double t282;
  double t285;
  double t289;
  double t290;
  double t291;
  double t293;
  double t294;
  double t296;
  double t297;
  double t299;
  double t300;
  double t302;
  double t304;
  double t306;
  double t308;
  double t31;
  double t33;
  double t34;
  double t368;
  double t370;
  double t372;
  double t373;
  double t375;
  double t376;
  double t378;
  double t379;
  double t38;
  double t381;
  double t383;
  double t385;
  double t386;
  double t388;
  double t390;
  double t392;
  double t394;
  double t396;
  double t398;
  double t400;
  double t402;
  double t405;
  double t410;
  double t417;
  double t42;
  double t43;
  double t437;
  double t45;
  double t453;
  double t455;
  double t47;
  double t48;
  double t5;
  double t50;
  double t55;
  double t59;
  double t6;
  double t61;
  double t62;
  double t63;
  double t67;
  double t69;
  double t70;
  double t72;
  double t74;
  double t79;
  double t83;
  double t86;
  double t88;
  double t90;
  double t91;
  double t93;
  double t95;
  double t97;
  double t98;
  double t99;
  {
    t2 = 1/delta*A;
    t5 = 0.5*phi+0.5*old;
    t6 = t5*t5;
    t16 = t6*t5;
    t21 = 0.6E2*t16-0.9E2*t6+0.15E2*phi+0.15E2*old;
    t23 = log(rho);
    t24 = rho*t23;
    t26 = -0.2011952932E2*rho+0.289445110522E2*t24+0.2133795654E2;
    t28 = -t21;
    t31 = 0.193942126E2*rho+0.950205321859E1*t24+0.4540504209;
    t33 = t6*t6;
    t34 = t33*t5;
    t38 = -0.2889626582E2+0.1204277045E3*t34-0.3010692614E3*t33+0.2007128409E3*
t16;
    t42 = 0.9502053219E1+0.116654747E3*t34-0.2916368674E3*t33+0.1944245783E3*
t16;
    t43 = 1/t42;
    t45 = exp(t38*t43);
    t47 = log(t45);
    t48 = t45*t47;
    t50 = 0.2133795654E2-0.2011952932E2*t45+0.289445110522E2*t48;
    t55 = 30.0*t33-60.0*t16+30.0*t6;
    t59 = 0.3010692612E3*t33-0.6021385228E3*t16+0.3010692614E3*t6;
    t61 = t42*t42;
    t62 = 1/t61;
    t63 = t38*t62;
    t67 = 0.2916368675E3*t33-0.5832737348E3*t16+0.2916368674E3*t6;
    t69 = t59*t43-t63*t67;
    t70 = t69*t45;
    t72 = t70*t47;
    t74 = 0.882498173E1*t70+0.289445110522E2*t72;
    t79 = 0.15E2*t33-0.3E2*t16+0.15E2*t6;
    t83 = 0.6021385225E3*t33-0.1204277046E4*t16+0.6021385227E3*t6;
    t86 = 0.116654747E4*t16;
    t88 = 0.583273735E3*t33-t86+0.5832737349E3*t6;
    t90 = t83*t43-t63*t88;
    t91 = t90*t45;
    t93 = t91*t47;
    t95 = 0.882498173E1*t91+0.289445110522E2*t93;
    t97 = 6.0*t34;
    t98 = 15.0*t33;
    t99 = 10.0*t16;
    t100 = t97-t98+t99;
    t105 = 0.1204277045E4*t16-0.1806415569E4*t6+0.3010692614E3*phi+
0.3010692614E3*old;
    t107 = t83*t62;
    t109 = t59*t62;
    t112 = 1/t61/t42;
    t113 = t38*t112;
    t114 = t88*t67;
    t120 = t86-0.1749821205E4*t6+0.2916368674E3*phi+0.2916368674E3*old;
    t122 = t105*t43-t107*t67-t109*t88+2.0*t113*t114-t63*t120;
    t123 = t122*t45;
    t125 = t90*t69;
    t126 = t125*t45;
    t128 = t123*t47;
    t130 = t125*t48;
    t132 = 0.882498173E1*t123+0.3776949278E2*t126+0.289445110522E2*t128+
0.289445110522E2*t130;
    t136 = 0.4540504209+0.193942126E2*t45+0.950205321859E1*t48;
    t138 = -t55;
    t141 = 0.2889626582E2*t70+0.950205321859E1*t72;
    t143 = -t79;
    t146 = 0.2889626582E2*t91+0.950205321859E1*t93;
    t148 = 1.0-t97+t98-t99;
    t153 = 0.2889626582E2*t123+0.3839831904E2*t126+0.950205321859E1*t128+
0.950205321859E1*t130;
    t160 = -0.18E3+0.18E3*phi+0.18E3*old;
    t162 = -t160;
    t168 = 60.0+360.0*t6-0.18E3*phi-0.18E3*old;
    t173 = 0.3E2+0.18E3*t6-0.9E2*phi-0.9E2*old;
    t180 = 120.0*t16-180.0*t6+0.3E2*phi+0.3E2*old;
    t187 = 0.240855409E4*t16-0.3612831138E4*t6+0.6021385225E3*phi+
0.6021385225E3*old;
    t191 = t88*t88;
    t195 = 0.349964241E4*t6;
    t198 = 0.233309494E4*t16-t195+0.583273735E3*phi+0.583273735E3*old;
    t200 = t187*t43-2.0*t107*t88+2.0*t113*t191-t63*t198;
    t201 = t200*t45;
    t203 = t90*t90;
    t204 = t203*t45;
    t206 = t201*t47;
    t208 = t204*t47;
    t210 = 0.882498173E1*t201+0.3776949278E2*t204+0.289445110522E2*t206+
0.289445110522E2*t208;
    t216 = 0.6021385225E3+0.3612831135E4*t6-0.1806415569E4*phi-0.1806415569E4*
old;
    t218 = t187*t62;
    t220 = t105*t62;
    t223 = t83*t112;
    t228 = t59*t112;
    t231 = t61*t61;
    t232 = 1/t231;
    t233 = t38*t232;
    t234 = t191*t67;
    t237 = t88*t120;
    t241 = t198*t67;
    t246 = 0.583273735E3+t195-0.1749821205E4*phi-0.1749821205E4*old;
    t248 = t216*t43-t218*t67-2.0*t220*t88+4.0*t223*t114-2.0*t107*t120+2.0*t228*
t191-6.0*t233*t234+4.0*t113*t237-t109*t198+2.0*t113*t241-t63*t246;
    t249 = t248*t45;
    t251 = t200*t69;
    t252 = t251*t45;
    t254 = t91*t122;
    t256 = t203*t69;
    t257 = t256*t45;
    t259 = t249*t47;
    t261 = t251*t48;
    t263 = t47*t122;
    t264 = t91*t263;
    t266 = t256*t48;
    t274 = 0.1204277045E4+0.722566227E4*t6-0.3612831138E4*phi-0.3612831138E4*
old;
    t282 = t191*t88;
    t285 = t88*t198;
    t289 = 0.349964241E4*phi;
    t290 = 0.349964241E4*old;
    t291 = 0.116654747E4+0.699928482E4*t6-t289-t290;
    t293 = t274*t43-3.0*t218*t88+6.0*t223*t191-3.0*t107*t198-6.0*t233*t282+6.0*
t113*t285-t63*t291;
    t294 = t293*t45;
    t296 = t200*t90;
    t297 = t296*t45;
    t299 = t203*t90;
    t300 = t299*t45;
    t302 = t294*t47;
    t304 = t296*t48;
    t306 = t300*t47;
    t308 = 0.882498173E1*t294+0.1133084783E3*t297+0.6671400383E2*t300+
0.289445110522E2*t302+0.8683353315E2*t304+0.289445110522E2*t306;
    t368 = -3.0*t107*t246-6.0*t59*t232*t282+24.0*t38/t231/t42*t282*t67-18.0*
t233*t191*t120+6.0*t228*t285-18.0*t233*t285*t67+6.0*t113*t120*t198+6.0*t113*t88
*t246-t109*t291+2.0*t113*t291*t67-t63*(-0.349964241E4+t289+t290);
    t370 = ((-0.3612831138E4+0.3612831135E4*phi+0.3612831135E4*old)*t43-t274*
t62*t67-3.0*t216*t62*t88+6.0*t187*t112*t114-3.0*t218*t120+6.0*t105*t112*t191
-18.0*t83*t232*t234+12.0*t223*t237-3.0*t220*t198+6.0*t223*t241+t368)*t45;
    t372 = t293*t69;
    t373 = t372*t45;
    t375 = t248*t90;
    t376 = t375*t45;
    t378 = t200*t122;
    t379 = t378*t45;
    t381 = t296*t70;
    t383 = t204*t122;
    t385 = t299*t69;
    t386 = t385*t45;
    t388 = t370*t47;
    t390 = t372*t48;
    t392 = t375*t48;
    t394 = t378*t48;
    t396 = t296*t72;
    t398 = t204*t263;
    t400 = t385*t48;
    t402 = 0.882498173E1*t370+0.3776949278E2*t373+0.1133084783E3*t376+
0.1133084783E3*t379+0.2001420114E3*t381+0.2001420115E3*t383+0.9565851488E2*t386
+0.289445110522E2*t388+0.289445110522E2*t390+0.8683353315E2*t392+0.8683353315E2
*t394+0.8683353315E2*t396+0.8683353315E2*t398+0.289445110522E2*t400;
    t405 = -t168;
    t410 = -t180;
    t417 = 0.2889626582E2*t201+0.3839831904E2*t204+0.950205321859E1*t206+
0.950205321859E1*t208;
    t437 = 0.2889626582E2*t294+0.1151949571E3*t297+0.4790037226E2*t300+
0.950205321859E1*t302+0.2850615966E2*t304+0.950205321859E1*t306;
    t453 = 0.2889626582E2*t370+0.3839831904E2*t373+0.1151949571E3*t376+
0.1151949571E3*t379+0.1437011168E3*t381+0.1437011168E3*t383+0.5740242548E2*t386
+0.950205321859E1*t388+0.950205321859E1*t390+0.2850615966E2*t392+0.2850615966E2
*t394+0.2850615966E2*t396+0.2850615966E2*t398+0.950205321859E1*t400;
    t455 = t160*t50+t168*t74+3.0*t173*t95+3.0*t180*t132+3.0*t21*t210+3.0*t55*(
0.882498173E1*t249+0.3776949278E2*t252+0.7553898556E2*t254+0.6671400383E2*t257+
0.289445110522E2*t259+0.289445110522E2*t261+0.578890221E2*t264+0.289445110522E2
*t266)+t79*t308+t100*t402+t162*t136+t405*t141-3.0*t173*t146+3.0*t410*t153+3.0*
t28*t417+3.0*t138*(0.2889626582E2*t249+0.3839831904E2*t252+0.7679663808E2*t254+
0.4790037226E2*t257+0.950205321859E1*t259+0.950205321859E1*t261+0.1900410644E2*
t264+0.950205321859E1*t266)+t143*t437+t148*t453;
    return(2.0*t2*(0.1E1+0.6E1*t6+0.3E1*(2.0*alpha-2.0)*t5-0.3E1*alpha)+t21*t26
+t28*t31-beta*(t21*t50+t55*t74+t79*t95+t100*t132+t28*t136+t138*t141+t143*t146+
t148*t153)+(0.24E2*t2+t160*t26+t162*t31-beta*t455)*(phi-old)/24.0+t2*(-12.0+
0.12E2*phi+0.12E2*old+12.0*alpha)/12.0+t168*t26/24.0+t405*t31/24.0-beta*(t168*
t50+3.0*t180*t95+3.0*t55*t210+t100*t308+t405*t136+3.0*t410*t146+3.0*t138*t417+
t148*t437)/24.0);
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t11;
  double t12;
  double t16;
  double t2;
  double t20;
  double t21;
  double t4;
  double t5;
  double t7;
  double t8;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t8 = t4-t5+t7;
    t11 = 0.5*rho+0.5*old;
    t12 = log(t11);
    t16 = 1.0-t4+t5-t7;
    t20 = t11*t11;
    t21 = 1/t20;
    return(t8*(0.882498173E1+0.289445110522E2*t12)+t16*(0.2889626582E2+
0.950205321859E1*t12)+(-0.2894451105E2*t8*t21-0.9502053219E1*t16*t21)*(rho-old)
/24.0);
  }
}

double drhochemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t11;
  double t12;
  double t15;
  double t18;
  double t2;
  double t20;
  double t29;
  double t4;
  double t5;
  double t7;
  double t8;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t8 = t4-t5+t7;
    t11 = 0.5*rho+0.5*old;
    t12 = 1/t11;
    t15 = 1.0-t4+t5-t7;
    t18 = t11*t11;
    t20 = 1/t18/t11;
    t29 = 1/t18;
    return(0.1447225552E2*t8*t12+0.475102661E1*t15*t12+(0.2894451105E2*t8*t20+
0.9502053219E1*t15*t20)*(rho-old)/24.0-0.1206021294E1*t8*t29-0.3959188841*t15*
t29);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t10;
  double t11;
  double t15;
  double t19;
  double t2;
  double t20;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t7 = 30.0*t2-60.0*t1*phi+30.0*t1;
    t10 = 0.5*rho+0.5*old;
    t11 = log(t10);
    t15 = -t7;
    t19 = t10*t10;
    t20 = 1/t19;
    return(t7*(0.882498173E1+0.289445110522E2*t11)+t15*(0.2889626582E2+
0.950205321859E1*t11)+(-0.2894451105E2*t7*t20-0.9502053219E1*t15*t20)*(rho-old)
/24.0);
  }
}

#include <math.h>
double pressure(double rho,double phi)
{
  double t1;
  double t10;
  double t11;
  double t2;
  double t4;
  double t5;
  double t6;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t6 = t1*phi;
    t7 = 10.0*t6;
    t10 = log(rho);
    t11 = rho*t10;
    return((t4-t5+t7)*(-0.2094881058E2+0.4291934897E2*rho-0.2894451105E2*t11+
rho*(-0.1397483792E2+0.289445110522E2*t10))+(1.0-t4+t5-t7)*(-0.6490446091E-1+
0.3405607049E1*rho-0.9502053219E1*t11+rho*(0.609644617E1+0.950205321859E1*t10))
-2.0/delta*A*(t2+(2.0*alpha-2.0)*t6+(-3.0*alpha+1.0)*t1+alpha));
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t1;
  double t2;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t9 = sqrt(0.116654747E11*t2*phi-0.2916368675E11*t2+0.1944245783E11*t1*phi+
950205322.0);
    return(0.1E-3*t9);
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t15;
  double t16;
  double t17;
  double t18;
  double t19;
  double t21;
  double t22;
  double t26;
  double t3;
  double t4;
  double t41;
  double t43;
  double t44;
  double t7;
  {
    t3 = phi*phi;
    t4 = t3*t3;
    t7 = t3*phi;
    t15 = t4*phi;
    t16 = 6.0*t15;
    t17 = 15.0*t4;
    t18 = 10.0*t7;
    t19 = t16-t17+t18;
    t21 = log(rho);
    t22 = rho*t21;
    t26 = 1.0-t16+t17-t18;
    t41 = exp((-0.2889626582E2+0.1204277045E3*t15-0.3010692614E3*t4+
0.2007128409E3*t7)/(0.9502053219E1+0.116654747E3*t15-0.2916368674E3*t4+
0.1944245783E3*t7));
    t43 = log(t41);
    t44 = t41*t43;
    return(2.0/delta*A*(t4+(2.0*alpha-2.0)*t7+(-3.0*alpha+1.0)*t3+alpha)+t19*(
-0.2011952932E2*rho+0.289445110522E2*t22+0.2133795654E2)+t26*(0.193942126E2*rho
+0.950205321859E1*t22+0.4540504209)-beta*(t19*(0.2133795654E2-0.2011952932E2*
t41+0.289445110522E2*t44)+t26*(0.4540504209+0.193942126E2*t41+0.950205321859E1*
t44)));
  }
}

#include <math.h>
double reactionSource(double rho,double phi,double old)
{
  double t100;
  double t107;
  double t114;
  double t116;
  double t120;
  double t121;
  double t122;
  double t129;
  double t131;
  double t132;
  double t134;
  double t135;
  double t137;
  double t139;
  double t157;
  double t172;
  double t174;
  double t175;
  double t178;
  double t180;
  double t182;
  double t184;
  double t2;
  double t20;
  double t24;
  double t26;
  double t27;
  double t29;
  double t31;
  double t34;
  double t36;
  double t40;
  double t44;
  double t45;
  double t47;
  double t49;
  double t5;
  double t50;
  double t52;
  double t54;
  double t55;
  double t56;
  double t57;
  double t6;
  double t61;
  double t63;
  double t64;
  double t65;
  double t69;
  double t7;
  double t71;
  double t72;
  double t74;
  double t76;
  double t80;
  double t82;
  double t85;
  double t98;
  {
    t2 = 1/delta*A;
    t5 = 0.5*phi+0.5*old;
    t6 = t5*t5;
    t7 = t6*t5;
    t20 = t6*t6;
    t24 = 30.0*t20-60.0*t7+30.0*t6;
    t26 = log(rho);
    t27 = rho*t26;
    t29 = -0.2011952932E2*rho+0.289445110522E2*t27+0.2133795654E2;
    t31 = -t24;
    t34 = 0.193942126E2*rho+0.950205321859E1*t27+0.4540504209;
    t36 = t20*t5;
    t40 = -0.2889626582E2+0.1204277045E3*t36-0.3010692614E3*t20+0.2007128409E3*
t7;
    t44 = 0.9502053219E1+0.116654747E3*t36-0.2916368674E3*t20+0.1944245783E3*t7
;
    t45 = 1/t44;
    t47 = exp(t40*t45);
    t49 = log(t47);
    t50 = t47*t49;
    t52 = 0.2133795654E2-0.2011952932E2*t47+0.289445110522E2*t50;
    t54 = 6.0*t36;
    t55 = 15.0*t20;
    t56 = 10.0*t7;
    t57 = t54-t55+t56;
    t61 = 0.6021385225E3*t20-0.1204277046E4*t7+0.6021385227E3*t6;
    t63 = t44*t44;
    t64 = 1/t63;
    t65 = t40*t64;
    t69 = 0.583273735E3*t20-0.116654747E4*t7+0.5832737349E3*t6;
    t71 = t61*t45-t65*t69;
    t72 = t71*t47;
    t74 = t72*t49;
    t76 = 0.882498173E1*t72+0.289445110522E2*t74;
    t80 = 0.4540504209+0.193942126E2*t47+0.950205321859E1*t50;
    t82 = 1.0-t54+t55-t56;
    t85 = 0.2889626582E2*t72+0.950205321859E1*t74;
    t98 = 60.0+360.0*t6-0.18E3*phi-0.18E3*old;
    t100 = -t98;
    t107 = 120.0*t7-180.0*t6+0.3E2*phi+0.3E2*old;
    t114 = 0.240855409E4*t7-0.3612831138E4*t6+0.6021385225E3*phi+0.6021385225E3
*old;
    t116 = t61*t64;
    t120 = 1/t63/t44;
    t121 = t40*t120;
    t122 = t69*t69;
    t129 = 0.233309494E4*t7-0.349964241E4*t6+0.583273735E3*phi+0.583273735E3*
old;
    t131 = t114*t45-2.0*t116*t69+2.0*t122*t121-t65*t129;
    t132 = t131*t47;
    t134 = t71*t71;
    t135 = t134*t47;
    t137 = t132*t49;
    t139 = t135*t49;
    t157 = t63*t63;
    t172 = ((0.1204277045E4+0.722566227E4*t6-0.3612831138E4*phi-0.3612831138E4*
old)*t45-3.0*t114*t64*t69+6.0*t61*t120*t122-3.0*t116*t129-6.0*t40/t157*t122*t69
+6.0*t121*t69*t129-t65*(0.116654747E4+0.699928482E4*t6-0.349964241E4*phi
-0.349964241E4*old))*t47;
    t174 = t131*t71;
    t175 = t174*t47;
    t178 = t134*t71*t47;
    t180 = t172*t49;
    t182 = t174*t50;
    t184 = t178*t49;
    return(2.0*t2*(4.0*t7+3.0*(2.0*alpha-2.0)*t6+2.0*(-3.0*alpha+1.0)*t5)+t24*
t29+t31*t34-beta*(t24*t52+t57*t76+t31*t80+t82*t85)+(2.0*t2*(-12.0+0.12E2*phi+
0.12E2*old+12.0*alpha)+t98*t29+t100*t34-beta*(t98*t52+3.0*t107*t76+3.0*t24*(
0.882498173E1*t132+0.3776949278E2*t135+0.289445110522E2*t137+0.289445110522E2*
t139)+t57*(0.882498173E1*t172+0.1133084783E3*t175+0.6671400383E2*t178+
0.289445110522E2*t180+0.8683353315E2*t182+0.289445110522E2*t184)+t100*t80-3.0*
t107*t85+3.0*t31*(0.2889626582E2*t132+0.3839831904E2*t135+0.950205321859E1*t137
+0.950205321859E1*t139)+t82*(0.2889626582E2*t172+0.1151949571E3*t175+
0.4790037226E2*t178+0.950205321859E1*t180+0.2850615966E2*t182+0.950205321859E1*
t184)))*(phi-old)/24.0);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi,double old,double mid)
{
  double t100;
  double t105;
  double t107;
  double t109;
  double t112;
  double t113;
  double t114;
  double t120;
  double t122;
  double t123;
  double t125;
  double t126;
  double t128;
  double t130;
  double t132;
  double t136;
  double t138;
  double t141;
  double t143;
  double t146;
  double t148;
  double t153;
  double t16;
  double t160;
  double t162;
  double t168;
  double t173;
  double t180;
  double t187;
  double t191;
  double t195;
  double t198;
  double t2;
  double t200;
  double t201;
  double t203;
  double t204;
  double t206;
  double t208;
  double t21;
  double t210;
  double t216;
  double t218;
  double t220;
  double t223;
  double t228;
  double t23;
  double t231;
  double t232;
  double t233;
  double t234;
  double t237;
  double t24;
  double t241;
  double t246;
  double t248;
  double t249;
  double t251;
  double t252;
  double t254;
  double t256;
  double t257;
  double t259;
  double t26;
  double t261;
  double t263;
  double t264;
  double t266;
  double t274;
  double t28;
  double t282;
  double t285;
  double t289;
  double t290;
  double t291;
  double t293;
  double t294;
  double t296;
  double t297;
  double t299;
  double t300;
  double t302;
  double t304;
  double t306;
  double t308;
  double t31;
  double t33;
  double t34;
  double t368;
  double t370;
  double t372;
  double t373;
  double t375;
  double t376;
  double t378;
  double t379;
  double t38;
  double t381;
  double t383;
  double t385;
  double t386;
  double t388;
  double t390;
  double t392;
  double t394;
  double t396;
  double t398;
  double t400;
  double t402;
  double t405;
  double t410;
  double t417;
  double t42;
  double t43;
  double t437;
  double t45;
  double t453;
  double t455;
  double t47;
  double t48;
  double t5;
  double t50;
  double t55;
  double t59;
  double t6;
  double t61;
  double t62;
  double t63;
  double t67;
  double t69;
  double t70;
  double t72;
  double t74;
  double t79;
  double t83;
  double t86;
  double t88;
  double t90;
  double t91;
  double t93;
  double t95;
  double t97;
  double t98;
  double t99;
  {
    t2 = 1/delta*A;
    t5 = 0.5*phi+0.5*old;
    t6 = t5*t5;
    t16 = t6*t5;
    t21 = 0.6E2*t16-0.9E2*t6+0.15E2*phi+0.15E2*old;
    t23 = log(rho);
    t24 = rho*t23;
    t26 = -0.2011952932E2*rho+0.289445110522E2*t24+0.2133795654E2;
    t28 = -t21;
    t31 = 0.193942126E2*rho+0.950205321859E1*t24+0.4540504209;
    t33 = t6*t6;
    t34 = t33*t5;
    t38 = -0.2889626582E2+0.1204277045E3*t34-0.3010692614E3*t33+0.2007128409E3*
t16;
    t42 = 0.9502053219E1+0.116654747E3*t34-0.2916368674E3*t33+0.1944245783E3*
t16;
    t43 = 1/t42;
    t45 = exp(t38*t43);
    t47 = log(t45);
    t48 = t45*t47;
    t50 = 0.2133795654E2-0.2011952932E2*t45+0.289445110522E2*t48;
    t55 = 30.0*t33-60.0*t16+30.0*t6;
    t59 = 0.3010692612E3*t33-0.6021385228E3*t16+0.3010692614E3*t6;
    t61 = t42*t42;
    t62 = 1/t61;
    t63 = t38*t62;
    t67 = 0.2916368675E3*t33-0.5832737348E3*t16+0.2916368674E3*t6;
    t69 = t59*t43-t63*t67;
    t70 = t69*t45;
    t72 = t70*t47;
    t74 = 0.882498173E1*t70+0.289445110522E2*t72;
    t79 = 0.15E2*t33-0.3E2*t16+0.15E2*t6;
    t83 = 0.6021385225E3*t33-0.1204277046E4*t16+0.6021385227E3*t6;
    t86 = 0.116654747E4*t16;
    t88 = 0.583273735E3*t33-t86+0.5832737349E3*t6;
    t90 = t83*t43-t63*t88;
    t91 = t90*t45;
    t93 = t91*t47;
    t95 = 0.882498173E1*t91+0.289445110522E2*t93;
    t97 = 6.0*t34;
    t98 = 15.0*t33;
    t99 = 10.0*t16;
    t100 = t97-t98+t99;
    t105 = 0.1204277045E4*t16-0.1806415569E4*t6+0.3010692614E3*phi+
0.3010692614E3*old;
    t107 = t83*t62;
    t109 = t59*t62;
    t112 = 1/t61/t42;
    t113 = t38*t112;
    t114 = t88*t67;
    t120 = t86-0.1749821205E4*t6+0.2916368674E3*phi+0.2916368674E3*old;
    t122 = t105*t43-t107*t67-t109*t88+2.0*t113*t114-t63*t120;
    t123 = t122*t45;
    t125 = t90*t69;
    t126 = t125*t45;
    t128 = t123*t47;
    t130 = t125*t48;
    t132 = 0.882498173E1*t123+0.3776949278E2*t126+0.289445110522E2*t128+
0.289445110522E2*t130;
    t136 = 0.4540504209+0.193942126E2*t45+0.950205321859E1*t48;
    t138 = -t55;
    t141 = 0.2889626582E2*t70+0.950205321859E1*t72;
    t143 = -t79;
    t146 = 0.2889626582E2*t91+0.950205321859E1*t93;
    t148 = 1.0-t97+t98-t99;
    t153 = 0.2889626582E2*t123+0.3839831904E2*t126+0.950205321859E1*t128+
0.950205321859E1*t130;
    t160 = -0.18E3+0.18E3*phi+0.18E3*old;
    t162 = -t160;
    t168 = 60.0+360.0*t6-0.18E3*phi-0.18E3*old;
    t173 = 0.3E2+0.18E3*t6-0.9E2*phi-0.9E2*old;
    t180 = 120.0*t16-180.0*t6+0.3E2*phi+0.3E2*old;
    t187 = 0.240855409E4*t16-0.3612831138E4*t6+0.6021385225E3*phi+
0.6021385225E3*old;
    t191 = t88*t88;
    t195 = 0.349964241E4*t6;
    t198 = 0.233309494E4*t16-t195+0.583273735E3*phi+0.583273735E3*old;
    t200 = t187*t43-2.0*t107*t88+2.0*t113*t191-t63*t198;
    t201 = t200*t45;
    t203 = t90*t90;
    t204 = t203*t45;
    t206 = t201*t47;
    t208 = t204*t47;
    t210 = 0.882498173E1*t201+0.3776949278E2*t204+0.289445110522E2*t206+
0.289445110522E2*t208;
    t216 = 0.6021385225E3+0.3612831135E4*t6-0.1806415569E4*phi-0.1806415569E4*
old;
    t218 = t187*t62;
    t220 = t105*t62;
    t223 = t83*t112;
    t228 = t59*t112;
    t231 = t61*t61;
    t232 = 1/t231;
    t233 = t38*t232;
    t234 = t191*t67;
    t237 = t88*t120;
    t241 = t198*t67;
    t246 = 0.583273735E3+t195-0.1749821205E4*phi-0.1749821205E4*old;
    t248 = t216*t43-t218*t67-2.0*t220*t88+4.0*t223*t114-2.0*t107*t120+2.0*t228*
t191-6.0*t233*t234+4.0*t113*t237-t109*t198+2.0*t113*t241-t63*t246;
    t249 = t248*t45;
    t251 = t200*t69;
    t252 = t251*t45;
    t254 = t91*t122;
    t256 = t203*t69;
    t257 = t256*t45;
    t259 = t249*t47;
    t261 = t251*t48;
    t263 = t47*t122;
    t264 = t91*t263;
    t266 = t256*t48;
    t274 = 0.1204277045E4+0.722566227E4*t6-0.3612831138E4*phi-0.3612831138E4*
old;
    t282 = t191*t88;
    t285 = t88*t198;
    t289 = 0.349964241E4*phi;
    t290 = 0.349964241E4*old;
    t291 = 0.116654747E4+0.699928482E4*t6-t289-t290;
    t293 = t274*t43-3.0*t218*t88+6.0*t223*t191-3.0*t107*t198-6.0*t233*t282+6.0*
t113*t285-t63*t291;
    t294 = t293*t45;
    t296 = t200*t90;
    t297 = t296*t45;
    t299 = t203*t90;
    t300 = t299*t45;
    t302 = t294*t47;
    t304 = t296*t48;
    t306 = t300*t47;
    t308 = 0.882498173E1*t294+0.1133084783E3*t297+0.6671400383E2*t300+
0.289445110522E2*t302+0.8683353315E2*t304+0.289445110522E2*t306;
    t368 = -3.0*t107*t246-6.0*t59*t232*t282+24.0*t38/t231/t42*t282*t67-18.0*
t233*t191*t120+6.0*t228*t285-18.0*t233*t285*t67+6.0*t113*t120*t198+6.0*t113*t88
*t246-t109*t291+2.0*t113*t291*t67-t63*(-0.349964241E4+t289+t290);
    t370 = ((-0.3612831138E4+0.3612831135E4*phi+0.3612831135E4*old)*t43-t274*
t62*t67-3.0*t216*t62*t88+6.0*t187*t112*t114-3.0*t218*t120+6.0*t105*t112*t191
-18.0*t83*t232*t234+12.0*t223*t237-3.0*t220*t198+6.0*t223*t241+t368)*t45;
    t372 = t293*t69;
    t373 = t372*t45;
    t375 = t248*t90;
    t376 = t375*t45;
    t378 = t200*t122;
    t379 = t378*t45;
    t381 = t296*t70;
    t383 = t204*t122;
    t385 = t299*t69;
    t386 = t385*t45;
    t388 = t370*t47;
    t390 = t372*t48;
    t392 = t375*t48;
    t394 = t378*t48;
    t396 = t296*t72;
    t398 = t204*t263;
    t400 = t385*t48;
    t402 = 0.882498173E1*t370+0.3776949278E2*t373+0.1133084783E3*t376+
0.1133084783E3*t379+0.2001420114E3*t381+0.2001420115E3*t383+0.9565851488E2*t386
+0.289445110522E2*t388+0.289445110522E2*t390+0.8683353315E2*t392+0.8683353315E2
*t394+0.8683353315E2*t396+0.8683353315E2*t398+0.289445110522E2*t400;
    t405 = -t168;
    t410 = -t180;
    t417 = 0.2889626582E2*t201+0.3839831904E2*t204+0.950205321859E1*t206+
0.950205321859E1*t208;
    t437 = 0.2889626582E2*t294+0.1151949571E3*t297+0.4790037226E2*t300+
0.950205321859E1*t302+0.2850615966E2*t304+0.950205321859E1*t306;
    t453 = 0.2889626582E2*t370+0.3839831904E2*t373+0.1151949571E3*t376+
0.1151949571E3*t379+0.1437011168E3*t381+0.1437011168E3*t383+0.5740242548E2*t386
+0.950205321859E1*t388+0.950205321859E1*t390+0.2850615966E2*t392+0.2850615966E2
*t394+0.2850615966E2*t396+0.2850615966E2*t398+0.950205321859E1*t400;
    t455 = t160*t50+t168*t74+3.0*t173*t95+3.0*t180*t132+3.0*t21*t210+3.0*t55*(
0.882498173E1*t249+0.3776949278E2*t252+0.7553898556E2*t254+0.6671400383E2*t257+
0.289445110522E2*t259+0.289445110522E2*t261+0.578890221E2*t264+0.289445110522E2
*t266)+t79*t308+t100*t402+t162*t136+t405*t141-3.0*t173*t146+3.0*t410*t153+3.0*
t28*t417+3.0*t138*(0.2889626582E2*t249+0.3839831904E2*t252+0.7679663808E2*t254+
0.4790037226E2*t257+0.950205321859E1*t259+0.950205321859E1*t261+0.1900410644E2*
t264+0.950205321859E1*t266)+t143*t437+t148*t453;
    return(2.0*t2*(0.1E1+0.6E1*t6+0.3E1*(2.0*alpha-2.0)*t5-0.3E1*alpha)+t21*t26
+t28*t31-beta*(t21*t50+t55*t74+t79*t95+t100*t132+t28*t136+t138*t141+t143*t146+
t148*t153)+(0.24E2*t2+t160*t26+t162*t31-beta*t455)*(phi-old)/24.0+t2*(-12.0+
0.12E2*phi+0.12E2*old+12.0*alpha)/12.0+t168*t26/24.0+t405*t31/24.0-beta*(t168*
t50+3.0*t180*t95+3.0*t55*t210+t100*t308+t405*t136+3.0*t410*t146+3.0*t138*t417+
t148*t437)/24.0);
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t11;
  double t12;
  double t16;
  double t2;
  double t20;
  double t21;
  double t4;
  double t5;
  double t7;
  double t8;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t8 = t4-t5+t7;
    t11 = 0.5*rho+0.5*old;
    t12 = log(t11);
    t16 = 1.0-t4+t5-t7;
    t20 = t11*t11;
    t21 = 1/t20;
    return(t8*(0.882498173E1+0.289445110522E2*t12)+t16*(0.2889626582E2+
0.950205321859E1*t12)+(-0.2894451105E2*t8*t21-0.9502053219E1*t16*t21)*(rho-old)
/24.0);
  }
}

double drhochemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t11;
  double t12;
  double t15;
  double t18;
  double t2;
  double t20;
  double t29;
  double t4;
  double t5;
  double t7;
  double t8;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t8 = t4-t5+t7;
    t11 = 0.5*rho+0.5*old;
    t12 = 1/t11;
    t15 = 1.0-t4+t5-t7;
    t18 = t11*t11;
    t20 = 1/t18/t11;
    t29 = 1/t18;
    return(0.1447225552E2*t8*t12+0.475102661E1*t15*t12+(0.2894451105E2*t8*t20+
0.9502053219E1*t15*t20)*(rho-old)/24.0-0.1206021294E1*t8*t29-0.3959188841*t15*
t29);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t10;
  double t11;
  double t15;
  double t19;
  double t2;
  double t20;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t7 = 30.0*t2-60.0*t1*phi+30.0*t1;
    t10 = 0.5*rho+0.5*old;
    t11 = log(t10);
    t15 = -t7;
    t19 = t10*t10;
    t20 = 1/t19;
    return(t7*(0.882498173E1+0.289445110522E2*t11)+t15*(0.2889626582E2+
0.950205321859E1*t11)+(-0.2894451105E2*t7*t20-0.9502053219E1*t15*t20)*(rho-old)
/24.0);
  }
}

#include <math.h>
double pressure(double rho,double phi)
{
  double t1;
  double t10;
  double t11;
  double t2;
  double t4;
  double t5;
  double t6;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t6 = t1*phi;
    t7 = 10.0*t6;
    t10 = log(rho);
    t11 = rho*t10;
    return((t4-t5+t7)*(-0.2094881058E2+0.4291934897E2*rho-0.2894451105E2*t11+
rho*(-0.1397483792E2+0.289445110522E2*t10))+(1.0-t4+t5-t7)*(-0.6490446091E-1+
0.3405607049E1*rho-0.9502053219E1*t11+rho*(0.609644617E1+0.950205321859E1*t10))
-2.0/delta*A*(t2+(2.0*alpha-2.0)*t6+(-3.0*alpha+1.0)*t1+alpha));
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t1;
  double t2;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t9 = sqrt(0.116654747E11*t2*phi-0.2916368675E11*t2+0.1944245783E11*t1*phi+
950205322.0);
    return(0.1E-3*t9);
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t15;
  double t16;
  double t17;
  double t18;
  double t19;
  double t21;
  double t22;
  double t26;
  double t3;
  double t4;
  double t41;
  double t43;
  double t44;
  double t7;
  {
    t3 = phi*phi;
    t4 = t3*t3;
    t7 = t3*phi;
    t15 = t4*phi;
    t16 = 6.0*t15;
    t17 = 15.0*t4;
    t18 = 10.0*t7;
    t19 = t16-t17+t18;
    t21 = log(rho);
    t22 = rho*t21;
    t26 = 1.0-t16+t17-t18;
    t41 = exp((-0.2889626582E2+0.1204277045E3*t15-0.3010692614E3*t4+
0.2007128409E3*t7)/(0.9502053219E1+0.116654747E3*t15-0.2916368674E3*t4+
0.1944245783E3*t7));
    t43 = log(t41);
    t44 = t41*t43;
    return(2.0/delta*A*(t4+(2.0*alpha-2.0)*t7+(-3.0*alpha+1.0)*t3+alpha)+t19*(
-0.2011952932E2*rho+0.289445110522E2*t22+0.2133795654E2)+t26*(0.193942126E2*rho
+0.950205321859E1*t22+0.4540504209)-beta*(t19*(0.2133795654E2-0.2011952932E2*
t41+0.289445110522E2*t44)+t26*(0.4540504209+0.193942126E2*t41+0.950205321859E1*
t44)));
  }
}

#include <math.h>
double reactionSource(double rho,double phi,double old)
{
  double t100;
  double t107;
  double t114;
  double t116;
  double t120;
  double t121;
  double t122;
  double t129;
  double t131;
  double t132;
  double t134;
  double t135;
  double t137;
  double t139;
  double t157;
  double t172;
  double t174;
  double t175;
  double t178;
  double t180;
  double t182;
  double t184;
  double t2;
  double t20;
  double t24;
  double t26;
  double t27;
  double t29;
  double t31;
  double t34;
  double t36;
  double t40;
  double t44;
  double t45;
  double t47;
  double t49;
  double t5;
  double t50;
  double t52;
  double t54;
  double t55;
  double t56;
  double t57;
  double t6;
  double t61;
  double t63;
  double t64;
  double t65;
  double t69;
  double t7;
  double t71;
  double t72;
  double t74;
  double t76;
  double t80;
  double t82;
  double t85;
  double t98;
  {
    t2 = 1/delta*A;
    t5 = 0.5*phi+0.5*old;
    t6 = t5*t5;
    t7 = t6*t5;
    t20 = t6*t6;
    t24 = 30.0*t20-60.0*t7+30.0*t6;
    t26 = log(rho);
    t27 = rho*t26;
    t29 = -0.2011952932E2*rho+0.289445110522E2*t27+0.2133795654E2;
    t31 = -t24;
    t34 = 0.193942126E2*rho+0.950205321859E1*t27+0.4540504209;
    t36 = t20*t5;
    t40 = -0.2889626582E2+0.1204277045E3*t36-0.3010692614E3*t20+0.2007128409E3*
t7;
    t44 = 0.9502053219E1+0.116654747E3*t36-0.2916368674E3*t20+0.1944245783E3*t7
;
    t45 = 1/t44;
    t47 = exp(t40*t45);
    t49 = log(t47);
    t50 = t47*t49;
    t52 = 0.2133795654E2-0.2011952932E2*t47+0.289445110522E2*t50;
    t54 = 6.0*t36;
    t55 = 15.0*t20;
    t56 = 10.0*t7;
    t57 = t54-t55+t56;
    t61 = 0.6021385225E3*t20-0.1204277046E4*t7+0.6021385227E3*t6;
    t63 = t44*t44;
    t64 = 1/t63;
    t65 = t40*t64;
    t69 = 0.583273735E3*t20-0.116654747E4*t7+0.5832737349E3*t6;
    t71 = t61*t45-t65*t69;
    t72 = t71*t47;
    t74 = t72*t49;
    t76 = 0.882498173E1*t72+0.289445110522E2*t74;
    t80 = 0.4540504209+0.193942126E2*t47+0.950205321859E1*t50;
    t82 = 1.0-t54+t55-t56;
    t85 = 0.2889626582E2*t72+0.950205321859E1*t74;
    t98 = 360.0*t6-0.18E3*phi-0.18E3*old+60.0;
    t100 = -t98;
    t107 = 120.0*t7-180.0*t6+0.3E2*phi+0.3E2*old;
    t114 = 0.240855409E4*t7-0.3612831138E4*t6+0.6021385225E3*phi+0.6021385225E3
*old;
    t116 = t61*t64;
    t120 = 1/t63/t44;
    t121 = t40*t120;
    t122 = t69*t69;
    t129 = 0.233309494E4*t7-0.349964241E4*t6+0.583273735E3*phi+0.583273735E3*
old;
    t131 = t114*t45-2.0*t116*t69+2.0*t121*t122-t65*t129;
    t132 = t131*t47;
    t134 = t71*t71;
    t135 = t134*t47;
    t137 = t132*t49;
    t139 = t135*t49;
    t157 = t63*t63;
    t172 = ((0.722566227E4*t6-0.3612831138E4*phi-0.3612831138E4*old+
0.1204277045E4)*t45-3.0*t114*t64*t69+6.0*t61*t120*t122-3.0*t116*t129-6.0*t40/
t157*t122*t69+6.0*t121*t69*t129-t65*(0.699928482E4*t6-0.349964241E4*phi
-0.349964241E4*old+0.116654747E4))*t47;
    t174 = t131*t71;
    t175 = t174*t47;
    t178 = t134*t71*t47;
    t180 = t172*t49;
    t182 = t174*t50;
    t184 = t178*t49;
    return(2.0*t2*(4.0*t7+3.0*(2.0*alpha-2.0)*t6+2.0*(-3.0*alpha+1.0)*t5)+t24*
t29+t31*t34-beta*(t24*t52+t57*t76+t31*t80+t82*t85)+(2.0*t2*(0.12E2*phi+0.12E2*
old+12.0*alpha-12.0)+t98*t29+t100*t34-beta*(t98*t52+3.0*t107*t76+3.0*t24*(
0.882498173E1*t132+0.3776949278E2*t135+0.289445110522E2*t137+0.289445110522E2*
t139)+t57*(0.882498173E1*t172+0.1133084783E3*t175+0.6671400383E2*t178+
0.289445110522E2*t180+0.8683353315E2*t182+0.289445110522E2*t184)+t100*t80-3.0*
t107*t85+3.0*t31*(0.2889626582E2*t132+0.3839831904E2*t135+0.950205321859E1*t137
+0.950205321859E1*t139)+t82*(0.2889626582E2*t172+0.1151949571E3*t175+
0.4790037226E2*t178+0.950205321859E1*t180+0.2850615966E2*t182+0.950205321859E1*
t184)))*(phi-old)/24.0);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi,double old,double mid)
{
  double t100;
  double t105;
  double t107;
  double t109;
  double t112;
  double t113;
  double t114;
  double t120;
  double t122;
  double t123;
  double t125;
  double t126;
  double t128;
  double t130;
  double t132;
  double t136;
  double t138;
  double t141;
  double t143;
  double t146;
  double t148;
  double t153;
  double t16;
  double t160;
  double t162;
  double t168;
  double t173;
  double t180;
  double t187;
  double t191;
  double t195;
  double t198;
  double t2;
  double t200;
  double t201;
  double t203;
  double t204;
  double t206;
  double t208;
  double t21;
  double t210;
  double t216;
  double t218;
  double t220;
  double t223;
  double t228;
  double t23;
  double t231;
  double t232;
  double t233;
  double t234;
  double t237;
  double t24;
  double t241;
  double t246;
  double t248;
  double t249;
  double t251;
  double t252;
  double t254;
  double t256;
  double t257;
  double t259;
  double t26;
  double t261;
  double t263;
  double t264;
  double t266;
  double t274;
  double t28;
  double t282;
  double t285;
  double t289;
  double t290;
  double t291;
  double t293;
  double t294;
  double t296;
  double t297;
  double t299;
  double t300;
  double t302;
  double t304;
  double t306;
  double t308;
  double t31;
  double t33;
  double t34;
  double t368;
  double t370;
  double t372;
  double t373;
  double t375;
  double t376;
  double t378;
  double t379;
  double t38;
  double t381;
  double t383;
  double t385;
  double t386;
  double t388;
  double t390;
  double t392;
  double t394;
  double t396;
  double t398;
  double t400;
  double t402;
  double t405;
  double t410;
  double t417;
  double t42;
  double t43;
  double t437;
  double t45;
  double t453;
  double t455;
  double t47;
  double t48;
  double t5;
  double t50;
  double t55;
  double t59;
  double t6;
  double t61;
  double t62;
  double t63;
  double t67;
  double t69;
  double t70;
  double t72;
  double t74;
  double t79;
  double t83;
  double t86;
  double t88;
  double t90;
  double t91;
  double t93;
  double t95;
  double t97;
  double t98;
  double t99;
  {
    t2 = 1/delta*A;
    t5 = 0.5*phi+0.5*old;
    t6 = t5*t5;
    t16 = t6*t5;
    t21 = 0.6E2*t16-0.9E2*t6+0.15E2*phi+0.15E2*old;
    t23 = log(rho);
    t24 = rho*t23;
    t26 = -0.2011952932E2*rho+0.289445110522E2*t24+0.2133795654E2;
    t28 = -t21;
    t31 = 0.193942126E2*rho+0.950205321859E1*t24+0.4540504209;
    t33 = t6*t6;
    t34 = t33*t5;
    t38 = -0.2889626582E2+0.1204277045E3*t34-0.3010692614E3*t33+0.2007128409E3*
t16;
    t42 = 0.9502053219E1+0.116654747E3*t34-0.2916368674E3*t33+0.1944245783E3*
t16;
    t43 = 1/t42;
    t45 = exp(t38*t43);
    t47 = log(t45);
    t48 = t45*t47;
    t50 = 0.2133795654E2-0.2011952932E2*t45+0.289445110522E2*t48;
    t55 = 30.0*t33-60.0*t16+30.0*t6;
    t59 = 0.3010692612E3*t33-0.6021385228E3*t16+0.3010692614E3*t6;
    t61 = t42*t42;
    t62 = 1/t61;
    t63 = t38*t62;
    t67 = 0.2916368675E3*t33-0.5832737348E3*t16+0.2916368674E3*t6;
    t69 = t59*t43-t63*t67;
    t70 = t69*t45;
    t72 = t70*t47;
    t74 = 0.882498173E1*t70+0.289445110522E2*t72;
    t79 = 0.15E2*t33-0.3E2*t16+0.15E2*t6;
    t83 = 0.6021385225E3*t33-0.1204277046E4*t16+0.6021385227E3*t6;
    t86 = 0.116654747E4*t16;
    t88 = 0.583273735E3*t33-t86+0.5832737349E3*t6;
    t90 = t83*t43-t63*t88;
    t91 = t90*t45;
    t93 = t91*t47;
    t95 = 0.882498173E1*t91+0.289445110522E2*t93;
    t97 = 6.0*t34;
    t98 = 15.0*t33;
    t99 = 10.0*t16;
    t100 = t97-t98+t99;
    t105 = 0.1204277045E4*t16-0.1806415569E4*t6+0.3010692614E3*phi+
0.3010692614E3*old;
    t107 = t83*t62;
    t109 = t59*t62;
    t112 = 1/t61/t42;
    t113 = t38*t112;
    t114 = t88*t67;
    t120 = t86-0.1749821205E4*t6+0.2916368674E3*phi+0.2916368674E3*old;
    t122 = t105*t43-t107*t67-t109*t88+2.0*t113*t114-t63*t120;
    t123 = t122*t45;
    t125 = t90*t69;
    t126 = t125*t45;
    t128 = t123*t47;
    t130 = t125*t48;
    t132 = 0.882498173E1*t123+0.3776949278E2*t126+0.289445110522E2*t128+
0.289445110522E2*t130;
    t136 = 0.4540504209+0.193942126E2*t45+0.950205321859E1*t48;
    t138 = -t55;
    t141 = 0.2889626582E2*t70+0.950205321859E1*t72;
    t143 = -t79;
    t146 = 0.2889626582E2*t91+0.950205321859E1*t93;
    t148 = 1.0-t97+t98-t99;
    t153 = 0.2889626582E2*t123+0.3839831904E2*t126+0.950205321859E1*t128+
0.950205321859E1*t130;
    t160 = 0.18E3*phi+0.18E3*old-0.18E3;
    t162 = -t160;
    t168 = 360.0*t6-0.18E3*phi-0.18E3*old+60.0;
    t173 = 0.18E3*t6-0.9E2*phi-0.9E2*old+0.3E2;
    t180 = 120.0*t16-180.0*t6+0.3E2*phi+0.3E2*old;
    t187 = 0.240855409E4*t16-0.3612831138E4*t6+0.6021385225E3*phi+
0.6021385225E3*old;
    t191 = t88*t88;
    t195 = 0.349964241E4*t6;
    t198 = 0.233309494E4*t16-t195+0.583273735E3*phi+0.583273735E3*old;
    t200 = t187*t43-2.0*t107*t88+2.0*t113*t191-t63*t198;
    t201 = t200*t45;
    t203 = t90*t90;
    t204 = t203*t45;
    t206 = t201*t47;
    t208 = t204*t47;
    t210 = 0.882498173E1*t201+0.3776949278E2*t204+0.289445110522E2*t206+
0.289445110522E2*t208;
    t216 = 0.3612831135E4*t6-0.1806415569E4*phi-0.1806415569E4*old+
0.6021385225E3;
    t218 = t187*t62;
    t220 = t105*t62;
    t223 = t83*t112;
    t228 = t59*t112;
    t231 = t61*t61;
    t232 = 1/t231;
    t233 = t38*t232;
    t234 = t191*t67;
    t237 = t88*t120;
    t241 = t198*t67;
    t246 = 0.583273735E3+t195-0.1749821205E4*phi-0.1749821205E4*old;
    t248 = t216*t43-t218*t67-2.0*t220*t88+4.0*t223*t114-2.0*t107*t120+2.0*t228*
t191-6.0*t233*t234+4.0*t113*t237-t109*t198+2.0*t113*t241-t63*t246;
    t249 = t248*t45;
    t251 = t200*t69;
    t252 = t251*t45;
    t254 = t91*t122;
    t256 = t203*t69;
    t257 = t256*t45;
    t259 = t249*t47;
    t261 = t251*t48;
    t263 = t47*t122;
    t264 = t91*t263;
    t266 = t256*t48;
    t274 = 0.722566227E4*t6-0.3612831138E4*phi-0.3612831138E4*old+
0.1204277045E4;
    t282 = t191*t88;
    t285 = t88*t198;
    t289 = 0.349964241E4*phi;
    t290 = 0.349964241E4*old;
    t291 = 0.116654747E4+0.699928482E4*t6-t289-t290;
    t293 = t274*t43-3.0*t218*t88+6.0*t223*t191-3.0*t107*t198-6.0*t233*t282+6.0*
t113*t285-t63*t291;
    t294 = t293*t45;
    t296 = t200*t90;
    t297 = t296*t45;
    t299 = t203*t90;
    t300 = t299*t45;
    t302 = t294*t47;
    t304 = t296*t48;
    t306 = t300*t47;
    t308 = 0.882498173E1*t294+0.1133084783E3*t297+0.6671400383E2*t300+
0.289445110522E2*t302+0.8683353315E2*t304+0.289445110522E2*t306;
    t368 = -3.0*t107*t246-6.0*t59*t232*t282+24.0*t38/t231/t42*t282*t67-18.0*
t233*t191*t120+6.0*t228*t285-18.0*t233*t285*t67+6.0*t113*t120*t198+6.0*t113*t88
*t246-t109*t291+2.0*t113*t291*t67-t63*(-0.349964241E4+t289+t290);
    t370 = ((-0.3612831138E4+0.3612831135E4*phi+0.3612831135E4*old)*t43-t274*
t62*t67-3.0*t216*t62*t88+6.0*t187*t112*t114-3.0*t218*t120+6.0*t105*t112*t191
-18.0*t83*t232*t234+12.0*t223*t237-3.0*t220*t198+6.0*t223*t241+t368)*t45;
    t372 = t293*t69;
    t373 = t372*t45;
    t375 = t248*t90;
    t376 = t375*t45;
    t378 = t200*t122;
    t379 = t378*t45;
    t381 = t296*t70;
    t383 = t204*t122;
    t385 = t299*t69;
    t386 = t385*t45;
    t388 = t370*t47;
    t390 = t372*t48;
    t392 = t375*t48;
    t394 = t378*t48;
    t396 = t296*t72;
    t398 = t204*t263;
    t400 = t385*t48;
    t402 = 0.882498173E1*t370+0.3776949278E2*t373+0.1133084783E3*t376+
0.1133084783E3*t379+0.2001420114E3*t381+0.2001420115E3*t383+0.9565851488E2*t386
+0.289445110522E2*t388+0.289445110522E2*t390+0.8683353315E2*t392+0.8683353315E2
*t394+0.8683353315E2*t396+0.8683353315E2*t398+0.289445110522E2*t400;
    t405 = -t168;
    t410 = -t180;
    t417 = 0.2889626582E2*t201+0.3839831904E2*t204+0.950205321859E1*t206+
0.950205321859E1*t208;
    t437 = 0.2889626582E2*t294+0.1151949571E3*t297+0.4790037226E2*t300+
0.950205321859E1*t302+0.2850615966E2*t304+0.950205321859E1*t306;
    t453 = 0.2889626582E2*t370+0.3839831904E2*t373+0.1151949571E3*t376+
0.1151949571E3*t379+0.1437011168E3*t381+0.1437011168E3*t383+0.5740242548E2*t386
+0.950205321859E1*t388+0.950205321859E1*t390+0.2850615966E2*t392+0.2850615966E2
*t394+0.2850615966E2*t396+0.2850615966E2*t398+0.950205321859E1*t400;
    t455 = t160*t50+t168*t74+3.0*t173*t95+3.0*t180*t132+3.0*t21*t210+3.0*t55*(
0.882498173E1*t249+0.3776949278E2*t252+0.7553898556E2*t254+0.6671400383E2*t257+
0.289445110522E2*t259+0.289445110522E2*t261+0.578890221E2*t264+0.289445110522E2
*t266)+t79*t308+t100*t402+t162*t136+t405*t141-3.0*t173*t146+3.0*t410*t153+3.0*
t28*t417+3.0*t138*(0.2889626582E2*t249+0.3839831904E2*t252+0.7679663808E2*t254+
0.4790037226E2*t257+0.950205321859E1*t259+0.950205321859E1*t261+0.1900410644E2*
t264+0.950205321859E1*t266)+t143*t437+t148*t453;
    return(2.0*t2*(0.6E1*t6+0.3E1*(2.0*alpha-2.0)*t5-0.3E1*alpha+0.1E1)+t21*t26
+t28*t31-beta*(t21*t50+t55*t74+t79*t95+t100*t132+t28*t136+t138*t141+t143*t146+
t148*t153)+(0.24E2*t2+t160*t26+t162*t31-beta*t455)*(phi-old)/24.0+t2*(0.12E2*
phi+0.12E2*old+12.0*alpha-12.0)/12.0+t168*t26/24.0+t405*t31/24.0-beta*(t168*t50
+3.0*t180*t95+3.0*t55*t210+t100*t308+t405*t136+3.0*t410*t146+3.0*t138*t417+t148
*t437)/24.0);
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t11;
  double t12;
  double t16;
  double t2;
  double t20;
  double t21;
  double t4;
  double t5;
  double t7;
  double t8;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t8 = t4-t5+t7;
    t11 = 0.5*rho+0.5*old;
    t12 = log(t11);
    t16 = 1.0-t4+t5-t7;
    t20 = t11*t11;
    t21 = 1/t20;
    return(t8*(0.882498173E1+0.289445110522E2*t12)+t16*(0.2889626582E2+
0.950205321859E1*t12)+(-0.2894451105E2*t8*t21-0.9502053219E1*t16*t21)*(rho-old)
/24.0);
  }
}

double drhochemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t11;
  double t12;
  double t15;
  double t18;
  double t2;
  double t20;
  double t29;
  double t4;
  double t5;
  double t7;
  double t8;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t8 = t4-t5+t7;
    t11 = 0.5*rho+0.5*old;
    t12 = 1/t11;
    t15 = 1.0-t4+t5-t7;
    t18 = t11*t11;
    t20 = 1/t18/t11;
    t29 = 1/t18;
    return(0.1447225552E2*t8*t12+0.475102661E1*t15*t12+(0.2894451105E2*t8*t20+
0.9502053219E1*t15*t20)*(rho-old)/24.0-0.1206021294E1*t8*t29-0.3959188841*t15*
t29);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t10;
  double t11;
  double t15;
  double t19;
  double t2;
  double t20;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t7 = 30.0*t2-60.0*t1*phi+30.0*t1;
    t10 = 0.5*rho+0.5*old;
    t11 = log(t10);
    t15 = -t7;
    t19 = t10*t10;
    t20 = 1/t19;
    return(t7*(0.882498173E1+0.289445110522E2*t11)+t15*(0.2889626582E2+
0.950205321859E1*t11)+(-0.2894451105E2*t7*t20-0.9502053219E1*t15*t20)*(rho-old)
/24.0);
  }
}

#include <math.h>
double pressure(double rho,double phi)
{
  double t1;
  double t10;
  double t11;
  double t2;
  double t4;
  double t5;
  double t6;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t6 = t1*phi;
    t7 = 10.0*t6;
    t10 = log(rho);
    t11 = rho*t10;
    return((t4-t5+t7)*(-0.2094881058E2+0.4291934897E2*rho-0.2894451105E2*t11+
rho*(-0.1397483792E2+0.289445110522E2*t10))+(1.0-t4+t5-t7)*(-0.6490446091E-1+
0.3405607049E1*rho-0.9502053219E1*t11+rho*(0.609644617E1+0.950205321859E1*t10))
-2.0/delta*A*(t2+(2.0*alpha-2.0)*t6+(-3.0*alpha+1.0)*t1+alpha));
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t1;
  double t2;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t9 = sqrt(0.116654747E11*t2*phi-0.2916368675E11*t2+0.1944245783E11*t1*phi+
950205322.0);
    return(0.1E-3*t9);
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t15;
  double t16;
  double t17;
  double t18;
  double t19;
  double t21;
  double t22;
  double t26;
  double t3;
  double t4;
  double t41;
  double t43;
  double t44;
  double t7;
  {
    t3 = phi*phi;
    t4 = t3*t3;
    t7 = t3*phi;
    t15 = t4*phi;
    t16 = 6.0*t15;
    t17 = 15.0*t4;
    t18 = 10.0*t7;
    t19 = t16-t17+t18;
    t21 = log(rho);
    t22 = rho*t21;
    t26 = 1.0-t16+t17-t18;
    t41 = exp((-0.2889626582E2+0.1204277045E3*t15-0.3010692614E3*t4+
0.2007128409E3*t7)/(0.9502053219E1+0.116654747E3*t15-0.2916368674E3*t4+
0.1944245783E3*t7));
    t43 = log(t41);
    t44 = t41*t43;
    return(2.0/delta*A*(t4+(2.0*alpha-2.0)*t7+(-3.0*alpha+1.0)*t3+alpha)+t19*(
-0.2011952932E2*rho+0.289445110522E2*t22+0.2133795654E2)+t26*(0.193942126E2*rho
+0.950205321859E1*t22+0.4540504209)-beta*(t19*(0.2133795654E2-0.2011952932E2*
t41+0.289445110522E2*t44)+t26*(0.4540504209+0.193942126E2*t41+0.950205321859E1*
t44)));
  }
}

#include <math.h>
double reactionSource(double rho,double phi,double old)
{
  double t100;
  double t107;
  double t114;
  double t116;
  double t120;
  double t121;
  double t122;
  double t129;
  double t131;
  double t132;
  double t134;
  double t135;
  double t137;
  double t139;
  double t157;
  double t172;
  double t174;
  double t175;
  double t178;
  double t180;
  double t182;
  double t184;
  double t2;
  double t20;
  double t24;
  double t26;
  double t27;
  double t29;
  double t31;
  double t34;
  double t36;
  double t40;
  double t44;
  double t45;
  double t47;
  double t49;
  double t5;
  double t50;
  double t52;
  double t54;
  double t55;
  double t56;
  double t57;
  double t6;
  double t61;
  double t63;
  double t64;
  double t65;
  double t69;
  double t7;
  double t71;
  double t72;
  double t74;
  double t76;
  double t80;
  double t82;
  double t85;
  double t98;
  {
    t2 = 1/delta*A;
    t5 = 0.5*phi+0.5*old;
    t6 = t5*t5;
    t7 = t6*t5;
    t20 = t6*t6;
    t24 = 30.0*t20-60.0*t7+30.0*t6;
    t26 = log(rho);
    t27 = rho*t26;
    t29 = -0.2011952932E2*rho+0.289445110522E2*t27+0.2133795654E2;
    t31 = -t24;
    t34 = 0.193942126E2*rho+0.950205321859E1*t27+0.4540504209;
    t36 = t20*t5;
    t40 = -0.2889626582E2+0.1204277045E3*t36-0.3010692614E3*t20+0.2007128409E3*
t7;
    t44 = 0.9502053219E1+0.116654747E3*t36-0.2916368674E3*t20+0.1944245783E3*t7
;
    t45 = 1/t44;
    t47 = exp(t40*t45);
    t49 = log(t47);
    t50 = t47*t49;
    t52 = 0.2133795654E2-0.2011952932E2*t47+0.289445110522E2*t50;
    t54 = 6.0*t36;
    t55 = 15.0*t20;
    t56 = 10.0*t7;
    t57 = t54-t55+t56;
    t61 = 0.6021385225E3*t20-0.1204277046E4*t7+0.6021385227E3*t6;
    t63 = t44*t44;
    t64 = 1/t63;
    t65 = t40*t64;
    t69 = 0.583273735E3*t20-0.116654747E4*t7+0.5832737349E3*t6;
    t71 = t61*t45-t65*t69;
    t72 = t71*t47;
    t74 = t72*t49;
    t76 = 0.882498173E1*t72+0.289445110522E2*t74;
    t80 = 0.4540504209+0.193942126E2*t47+0.950205321859E1*t50;
    t82 = 1.0-t54+t55-t56;
    t85 = 0.2889626582E2*t72+0.950205321859E1*t74;
    t98 = 360.0*t6-0.18E3*phi-0.18E3*old+60.0;
    t100 = -t98;
    t107 = 120.0*t7-180.0*t6+0.3E2*phi+0.3E2*old;
    t114 = 0.240855409E4*t7-0.3612831138E4*t6+0.6021385225E3*phi+0.6021385225E3
*old;
    t116 = t61*t64;
    t120 = 1/t63/t44;
    t121 = t40*t120;
    t122 = t69*t69;
    t129 = 0.233309494E4*t7-0.349964241E4*t6+0.583273735E3*phi+0.583273735E3*
old;
    t131 = t114*t45-2.0*t116*t69+2.0*t121*t122-t65*t129;
    t132 = t131*t47;
    t134 = t71*t71;
    t135 = t134*t47;
    t137 = t132*t49;
    t139 = t135*t49;
    t157 = t63*t63;
    t172 = ((0.722566227E4*t6-0.3612831138E4*phi-0.3612831138E4*old+
0.1204277045E4)*t45-3.0*t114*t64*t69+6.0*t61*t120*t122-3.0*t116*t129-6.0*t40/
t157*t122*t69+6.0*t121*t69*t129-t65*(0.699928482E4*t6-0.349964241E4*phi
-0.349964241E4*old+0.116654747E4))*t47;
    t174 = t131*t71;
    t175 = t174*t47;
    t178 = t134*t71*t47;
    t180 = t172*t49;
    t182 = t174*t50;
    t184 = t178*t49;
    return(2.0*t2*(4.0*t7+3.0*(2.0*alpha-2.0)*t6+2.0*(-3.0*alpha+1.0)*t5)+t24*
t29+t31*t34-beta*(t24*t52+t57*t76+t31*t80+t82*t85)+(2.0*t2*(0.12E2*phi+0.12E2*
old+12.0*alpha-12.0)+t98*t29+t100*t34-beta*(t98*t52+3.0*t107*t76+3.0*t24*(
0.882498173E1*t132+0.3776949278E2*t135+0.289445110522E2*t137+0.289445110522E2*
t139)+t57*(0.882498173E1*t172+0.1133084783E3*t175+0.6671400383E2*t178+
0.289445110522E2*t180+0.8683353315E2*t182+0.289445110522E2*t184)+t100*t80-3.0*
t107*t85+3.0*t31*(0.2889626582E2*t132+0.3839831904E2*t135+0.950205321859E1*t137
+0.950205321859E1*t139)+t82*(0.2889626582E2*t172+0.1151949571E3*t175+
0.4790037226E2*t178+0.950205321859E1*t180+0.2850615966E2*t182+0.950205321859E1*
t184)))*(phi-old)/24.0);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi,double old,double mid)
{
  double t100;
  double t105;
  double t107;
  double t109;
  double t112;
  double t113;
  double t114;
  double t120;
  double t122;
  double t123;
  double t125;
  double t126;
  double t128;
  double t130;
  double t132;
  double t136;
  double t138;
  double t141;
  double t143;
  double t146;
  double t148;
  double t153;
  double t16;
  double t160;
  double t162;
  double t168;
  double t173;
  double t180;
  double t187;
  double t191;
  double t195;
  double t198;
  double t2;
  double t200;
  double t201;
  double t203;
  double t204;
  double t206;
  double t208;
  double t21;
  double t210;
  double t216;
  double t218;
  double t220;
  double t223;
  double t228;
  double t23;
  double t231;
  double t232;
  double t233;
  double t234;
  double t237;
  double t24;
  double t241;
  double t246;
  double t248;
  double t249;
  double t251;
  double t252;
  double t254;
  double t256;
  double t257;
  double t259;
  double t26;
  double t261;
  double t263;
  double t264;
  double t266;
  double t274;
  double t28;
  double t282;
  double t285;
  double t289;
  double t290;
  double t291;
  double t293;
  double t294;
  double t296;
  double t297;
  double t299;
  double t300;
  double t302;
  double t304;
  double t306;
  double t308;
  double t31;
  double t33;
  double t34;
  double t368;
  double t370;
  double t372;
  double t373;
  double t375;
  double t376;
  double t378;
  double t379;
  double t38;
  double t381;
  double t383;
  double t385;
  double t386;
  double t388;
  double t390;
  double t392;
  double t394;
  double t396;
  double t398;
  double t400;
  double t402;
  double t405;
  double t410;
  double t417;
  double t42;
  double t43;
  double t437;
  double t45;
  double t453;
  double t455;
  double t47;
  double t48;
  double t5;
  double t50;
  double t55;
  double t59;
  double t6;
  double t61;
  double t62;
  double t63;
  double t67;
  double t69;
  double t70;
  double t72;
  double t74;
  double t79;
  double t83;
  double t86;
  double t88;
  double t90;
  double t91;
  double t93;
  double t95;
  double t97;
  double t98;
  double t99;
  {
    t2 = 1/delta*A;
    t5 = 0.5*phi+0.5*old;
    t6 = t5*t5;
    t16 = t6*t5;
    t21 = 0.6E2*t16-0.9E2*t6+0.15E2*phi+0.15E2*old;
    t23 = log(rho);
    t24 = rho*t23;
    t26 = -0.2011952932E2*rho+0.289445110522E2*t24+0.2133795654E2;
    t28 = -t21;
    t31 = 0.193942126E2*rho+0.950205321859E1*t24+0.4540504209;
    t33 = t6*t6;
    t34 = t33*t5;
    t38 = -0.2889626582E2+0.1204277045E3*t34-0.3010692614E3*t33+0.2007128409E3*
t16;
    t42 = 0.9502053219E1+0.116654747E3*t34-0.2916368674E3*t33+0.1944245783E3*
t16;
    t43 = 1/t42;
    t45 = exp(t38*t43);
    t47 = log(t45);
    t48 = t45*t47;
    t50 = 0.2133795654E2-0.2011952932E2*t45+0.289445110522E2*t48;
    t55 = 30.0*t33-60.0*t16+30.0*t6;
    t59 = 0.3010692612E3*t33-0.6021385228E3*t16+0.3010692614E3*t6;
    t61 = t42*t42;
    t62 = 1/t61;
    t63 = t38*t62;
    t67 = 0.2916368675E3*t33-0.5832737348E3*t16+0.2916368674E3*t6;
    t69 = t59*t43-t63*t67;
    t70 = t69*t45;
    t72 = t70*t47;
    t74 = 0.882498173E1*t70+0.289445110522E2*t72;
    t79 = 0.15E2*t33-0.3E2*t16+0.15E2*t6;
    t83 = 0.6021385225E3*t33-0.1204277046E4*t16+0.6021385227E3*t6;
    t86 = 0.116654747E4*t16;
    t88 = 0.583273735E3*t33-t86+0.5832737349E3*t6;
    t90 = t83*t43-t63*t88;
    t91 = t90*t45;
    t93 = t91*t47;
    t95 = 0.882498173E1*t91+0.289445110522E2*t93;
    t97 = 6.0*t34;
    t98 = 15.0*t33;
    t99 = 10.0*t16;
    t100 = t97-t98+t99;
    t105 = 0.1204277045E4*t16-0.1806415569E4*t6+0.3010692614E3*phi+
0.3010692614E3*old;
    t107 = t83*t62;
    t109 = t59*t62;
    t112 = 1/t61/t42;
    t113 = t38*t112;
    t114 = t88*t67;
    t120 = t86-0.1749821205E4*t6+0.2916368674E3*phi+0.2916368674E3*old;
    t122 = t105*t43-t107*t67-t109*t88+2.0*t113*t114-t63*t120;
    t123 = t122*t45;
    t125 = t90*t69;
    t126 = t125*t45;
    t128 = t123*t47;
    t130 = t125*t48;
    t132 = 0.882498173E1*t123+0.3776949278E2*t126+0.289445110522E2*t128+
0.289445110522E2*t130;
    t136 = 0.4540504209+0.193942126E2*t45+0.950205321859E1*t48;
    t138 = -t55;
    t141 = 0.2889626582E2*t70+0.950205321859E1*t72;
    t143 = -t79;
    t146 = 0.2889626582E2*t91+0.950205321859E1*t93;
    t148 = 1.0-t97+t98-t99;
    t153 = 0.2889626582E2*t123+0.3839831904E2*t126+0.950205321859E1*t128+
0.950205321859E1*t130;
    t160 = 0.18E3*phi+0.18E3*old-0.18E3;
    t162 = -t160;
    t168 = 360.0*t6-0.18E3*phi-0.18E3*old+60.0;
    t173 = 0.18E3*t6-0.9E2*phi-0.9E2*old+0.3E2;
    t180 = 120.0*t16-180.0*t6+0.3E2*phi+0.3E2*old;
    t187 = 0.240855409E4*t16-0.3612831138E4*t6+0.6021385225E3*phi+
0.6021385225E3*old;
    t191 = t88*t88;
    t195 = 0.349964241E4*t6;
    t198 = 0.233309494E4*t16-t195+0.583273735E3*phi+0.583273735E3*old;
    t200 = t187*t43-2.0*t107*t88+2.0*t113*t191-t63*t198;
    t201 = t200*t45;
    t203 = t90*t90;
    t204 = t203*t45;
    t206 = t201*t47;
    t208 = t204*t47;
    t210 = 0.882498173E1*t201+0.3776949278E2*t204+0.289445110522E2*t206+
0.289445110522E2*t208;
    t216 = 0.3612831135E4*t6-0.1806415569E4*phi-0.1806415569E4*old+
0.6021385225E3;
    t218 = t187*t62;
    t220 = t105*t62;
    t223 = t83*t112;
    t228 = t59*t112;
    t231 = t61*t61;
    t232 = 1/t231;
    t233 = t38*t232;
    t234 = t191*t67;
    t237 = t88*t120;
    t241 = t198*t67;
    t246 = t195-0.1749821205E4*phi-0.1749821205E4*old+0.583273735E3;
    t248 = t216*t43-t218*t67-2.0*t220*t88+4.0*t223*t114-2.0*t107*t120+2.0*t228*
t191-6.0*t233*t234+4.0*t113*t237-t109*t198+2.0*t113*t241-t63*t246;
    t249 = t248*t45;
    t251 = t200*t69;
    t252 = t251*t45;
    t254 = t91*t122;
    t256 = t203*t69;
    t257 = t256*t45;
    t259 = t249*t47;
    t261 = t251*t48;
    t263 = t47*t122;
    t264 = t91*t263;
    t266 = t256*t48;
    t274 = 0.722566227E4*t6-0.3612831138E4*phi-0.3612831138E4*old+
0.1204277045E4;
    t282 = t191*t88;
    t285 = t88*t198;
    t289 = 0.349964241E4*phi;
    t290 = 0.349964241E4*old;
    t291 = 0.116654747E4+0.699928482E4*t6-t289-t290;
    t293 = t274*t43-3.0*t218*t88+6.0*t223*t191-3.0*t107*t198-6.0*t233*t282+6.0*
t113*t285-t63*t291;
    t294 = t293*t45;
    t296 = t200*t90;
    t297 = t296*t45;
    t299 = t203*t90;
    t300 = t299*t45;
    t302 = t294*t47;
    t304 = t296*t48;
    t306 = t300*t47;
    t308 = 0.882498173E1*t294+0.1133084783E3*t297+0.6671400383E2*t300+
0.289445110522E2*t302+0.8683353315E2*t304+0.289445110522E2*t306;
    t368 = -3.0*t107*t246-6.0*t59*t232*t282+24.0*t38/t231/t42*t282*t67-18.0*
t233*t191*t120+6.0*t228*t285-18.0*t233*t285*t67+6.0*t113*t120*t198+6.0*t113*t88
*t246-t109*t291+2.0*t113*t291*t67-t63*(-0.349964241E4+t289+t290);
    t370 = ((-0.3612831138E4+0.3612831135E4*phi+0.3612831135E4*old)*t43-t274*
t62*t67-3.0*t216*t62*t88+6.0*t187*t112*t114-3.0*t218*t120+6.0*t105*t112*t191
-18.0*t83*t232*t234+12.0*t223*t237-3.0*t220*t198+6.0*t223*t241+t368)*t45;
    t372 = t293*t69;
    t373 = t372*t45;
    t375 = t248*t90;
    t376 = t375*t45;
    t378 = t200*t122;
    t379 = t378*t45;
    t381 = t296*t70;
    t383 = t204*t122;
    t385 = t299*t69;
    t386 = t385*t45;
    t388 = t370*t47;
    t390 = t372*t48;
    t392 = t375*t48;
    t394 = t378*t48;
    t396 = t296*t72;
    t398 = t204*t263;
    t400 = t385*t48;
    t402 = 0.882498173E1*t370+0.3776949278E2*t373+0.1133084783E3*t376+
0.1133084783E3*t379+0.2001420114E3*t381+0.2001420115E3*t383+0.9565851488E2*t386
+0.289445110522E2*t388+0.289445110522E2*t390+0.8683353315E2*t392+0.8683353315E2
*t394+0.8683353315E2*t396+0.8683353315E2*t398+0.289445110522E2*t400;
    t405 = -t168;
    t410 = -t180;
    t417 = 0.2889626582E2*t201+0.3839831904E2*t204+0.950205321859E1*t206+
0.950205321859E1*t208;
    t437 = 0.2889626582E2*t294+0.1151949571E3*t297+0.4790037226E2*t300+
0.950205321859E1*t302+0.2850615966E2*t304+0.950205321859E1*t306;
    t453 = 0.2889626582E2*t370+0.3839831904E2*t373+0.1151949571E3*t376+
0.1151949571E3*t379+0.1437011168E3*t381+0.1437011168E3*t383+0.5740242548E2*t386
+0.950205321859E1*t388+0.950205321859E1*t390+0.2850615966E2*t392+0.2850615966E2
*t394+0.2850615966E2*t396+0.2850615966E2*t398+0.950205321859E1*t400;
    t455 = t160*t50+t168*t74+3.0*t173*t95+3.0*t180*t132+3.0*t21*t210+3.0*t55*(
0.882498173E1*t249+0.3776949278E2*t252+0.7553898556E2*t254+0.6671400383E2*t257+
0.289445110522E2*t259+0.289445110522E2*t261+0.578890221E2*t264+0.289445110522E2
*t266)+t79*t308+t100*t402+t162*t136+t405*t141-3.0*t173*t146+3.0*t410*t153+3.0*
t28*t417+3.0*t138*(0.2889626582E2*t249+0.3839831904E2*t252+0.7679663808E2*t254+
0.4790037226E2*t257+0.950205321859E1*t259+0.950205321859E1*t261+0.1900410644E2*
t264+0.950205321859E1*t266)+t143*t437+t148*t453;
    return(2.0*t2*(0.6E1*t6+0.3E1*(2.0*alpha-2.0)*t5-0.3E1*alpha+0.1E1)+t21*t26
+t28*t31-beta*(t21*t50+t55*t74+t79*t95+t100*t132+t28*t136+t138*t141+t143*t146+
t148*t153)+(0.24E2*t2+t160*t26+t162*t31-beta*t455)*(phi-old)/24.0+t2*(0.12E2*
phi+0.12E2*old+12.0*alpha-12.0)/12.0+t168*t26/24.0+t405*t31/24.0-beta*(t168*t50
+3.0*t180*t95+3.0*t55*t210+t100*t308+t405*t136+3.0*t410*t146+3.0*t138*t417+t148
*t437)/24.0);
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t11;
  double t12;
  double t16;
  double t2;
  double t20;
  double t21;
  double t4;
  double t5;
  double t7;
  double t8;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t8 = t4-t5+t7;
    t11 = 0.5*rho+0.5*old;
    t12 = log(t11);
    t16 = 1.0-t4+t5-t7;
    t20 = t11*t11;
    t21 = 1/t20;
    return(t8*(0.882498173E1+0.289445110522E2*t12)+t16*(0.2889626582E2+
0.950205321859E1*t12)+(-0.2894451105E2*t8*t21-0.9502053219E1*t16*t21)*(rho-old)
/24.0);
  }
}

double drhochemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t11;
  double t12;
  double t15;
  double t18;
  double t2;
  double t20;
  double t29;
  double t4;
  double t5;
  double t7;
  double t8;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t8 = t4-t5+t7;
    t11 = 0.5*rho+0.5*old;
    t12 = 1/t11;
    t15 = 1.0-t4+t5-t7;
    t18 = t11*t11;
    t20 = 1/t18/t11;
    t29 = 1/t18;
    return(0.1447225552E2*t8*t12+0.475102661E1*t15*t12+(0.2894451105E2*t8*t20+
0.9502053219E1*t15*t20)*(rho-old)/24.0-0.1206021294E1*t8*t29-0.3959188841*t15*
t29);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t10;
  double t11;
  double t15;
  double t19;
  double t2;
  double t20;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t7 = 30.0*t2-60.0*t1*phi+30.0*t1;
    t10 = 0.5*rho+0.5*old;
    t11 = log(t10);
    t15 = -t7;
    t19 = t10*t10;
    t20 = 1/t19;
    return(t7*(0.882498173E1+0.289445110522E2*t11)+t15*(0.2889626582E2+
0.950205321859E1*t11)+(-0.2894451105E2*t7*t20-0.9502053219E1*t15*t20)*(rho-old)
/24.0);
  }
}

#include <math.h>
double pressure(double rho,double phi)
{
  double t1;
  double t10;
  double t11;
  double t2;
  double t4;
  double t5;
  double t6;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t6 = t1*phi;
    t7 = 10.0*t6;
    t10 = log(rho);
    t11 = rho*t10;
    return((t4-t5+t7)*(-0.2094881058E2+0.4291934897E2*rho-0.2894451105E2*t11+
rho*(-0.1397483792E2+0.289445110522E2*t10))+(1.0-t4+t5-t7)*(-0.6490446091E-1+
0.3405607049E1*rho-0.9502053219E1*t11+rho*(0.609644617E1+0.950205321859E1*t10))
-2.0/delta*A*(t2+(2.0*alpha-2.0)*t6+(-3.0*alpha+1.0)*t1+alpha));
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t1;
  double t2;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t9 = sqrt(0.116654747E11*t2*phi-0.2916368675E11*t2+0.1944245783E11*t1*phi+
950205322.0);
    return(0.1E-3*t9);
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t15;
  double t16;
  double t17;
  double t18;
  double t19;
  double t21;
  double t22;
  double t26;
  double t3;
  double t4;
  double t41;
  double t43;
  double t44;
  double t7;
  {
    t3 = phi*phi;
    t4 = t3*t3;
    t7 = t3*phi;
    t15 = t4*phi;
    t16 = 6.0*t15;
    t17 = 15.0*t4;
    t18 = 10.0*t7;
    t19 = t16-t17+t18;
    t21 = log(rho);
    t22 = rho*t21;
    t26 = 1.0-t16+t17-t18;
    t41 = exp((-0.2889626582E2+0.1204277045E3*t15-0.3010692614E3*t4+
0.2007128409E3*t7)/(0.9502053219E1+0.116654747E3*t15-0.2916368674E3*t4+
0.1944245783E3*t7));
    t43 = log(t41);
    t44 = t41*t43;
    return(2.0/delta*A*(t4+(2.0*alpha-2.0)*t7+(-3.0*alpha+1.0)*t3+alpha)+t19*(
0.2133795654E2-0.2011952932E2*rho+0.289445110522E2*t22)+t26*(0.4540504209+
0.193942126E2*rho+0.950205321859E1*t22)-beta*(t19*(0.2133795654E2
-0.2011952932E2*t41+0.289445110522E2*t44)+t26*(0.4540504209+0.193942126E2*t41+
0.950205321859E1*t44)));
  }
}

#include <math.h>
double reactionSource(double rho,double phi,double old)
{
  double t100;
  double t107;
  double t114;
  double t116;
  double t120;
  double t121;
  double t122;
  double t129;
  double t131;
  double t132;
  double t134;
  double t135;
  double t137;
  double t139;
  double t157;
  double t172;
  double t174;
  double t175;
  double t178;
  double t180;
  double t182;
  double t184;
  double t2;
  double t20;
  double t24;
  double t26;
  double t27;
  double t29;
  double t31;
  double t34;
  double t36;
  double t40;
  double t44;
  double t45;
  double t47;
  double t49;
  double t5;
  double t50;
  double t52;
  double t54;
  double t55;
  double t56;
  double t57;
  double t6;
  double t61;
  double t63;
  double t64;
  double t65;
  double t69;
  double t7;
  double t71;
  double t72;
  double t74;
  double t76;
  double t80;
  double t82;
  double t85;
  double t98;
  {
    t2 = 1/delta*A;
    t5 = 0.5*phi+0.5*old;
    t6 = t5*t5;
    t7 = t6*t5;
    t20 = t6*t6;
    t24 = 30.0*t20-60.0*t7+30.0*t6;
    t26 = log(rho);
    t27 = rho*t26;
    t29 = 0.2133795654E2-0.2011952932E2*rho+0.289445110522E2*t27;
    t31 = -t24;
    t34 = 0.4540504209+0.193942126E2*rho+0.950205321859E1*t27;
    t36 = t20*t5;
    t40 = -0.2889626582E2+0.1204277045E3*t36-0.3010692614E3*t20+0.2007128409E3*
t7;
    t44 = 0.9502053219E1+0.116654747E3*t36-0.2916368674E3*t20+0.1944245783E3*t7
;
    t45 = 1/t44;
    t47 = exp(t40*t45);
    t49 = log(t47);
    t50 = t47*t49;
    t52 = 0.2133795654E2-0.2011952932E2*t47+0.289445110522E2*t50;
    t54 = 6.0*t36;
    t55 = 15.0*t20;
    t56 = 10.0*t7;
    t57 = t54-t55+t56;
    t61 = 0.6021385225E3*t20-0.1204277046E4*t7+0.6021385227E3*t6;
    t63 = t44*t44;
    t64 = 1/t63;
    t65 = t40*t64;
    t69 = 0.583273735E3*t20-0.116654747E4*t7+0.5832737349E3*t6;
    t71 = t61*t45-t65*t69;
    t72 = t71*t47;
    t74 = t72*t49;
    t76 = 0.882498173E1*t72+0.289445110522E2*t74;
    t80 = 0.4540504209+0.193942126E2*t47+0.950205321859E1*t50;
    t82 = 1.0-t54+t55-t56;
    t85 = 0.2889626582E2*t72+0.950205321859E1*t74;
    t98 = 360.0*t6-0.18E3*phi-0.18E3*old+60.0;
    t100 = -t98;
    t107 = 120.0*t7-180.0*t6+0.3E2*phi+0.3E2*old;
    t114 = 0.240855409E4*t7-0.3612831138E4*t6+0.6021385225E3*phi+0.6021385225E3
*old;
    t116 = t61*t64;
    t120 = 1/t63/t44;
    t121 = t40*t120;
    t122 = t69*t69;
    t129 = 0.233309494E4*t7-0.349964241E4*t6+0.583273735E3*phi+0.583273735E3*
old;
    t131 = t114*t45-2.0*t116*t69+2.0*t121*t122-t65*t129;
    t132 = t131*t47;
    t134 = t71*t71;
    t135 = t134*t47;
    t137 = t132*t49;
    t139 = t135*t49;
    t157 = t63*t63;
    t172 = ((0.722566227E4*t6-0.3612831138E4*phi-0.3612831138E4*old+
0.1204277045E4)*t45-3.0*t114*t64*t69+6.0*t61*t120*t122-3.0*t116*t129-6.0*t40/
t157*t122*t69+6.0*t121*t69*t129-t65*(0.699928482E4*t6-0.349964241E4*phi
-0.349964241E4*old+0.116654747E4))*t47;
    t174 = t131*t71;
    t175 = t174*t47;
    t178 = t134*t71*t47;
    t180 = t172*t49;
    t182 = t174*t50;
    t184 = t178*t49;
    return(2.0*t2*(4.0*t7+3.0*(2.0*alpha-2.0)*t6+2.0*(-3.0*alpha+1.0)*t5)+t24*
t29+t31*t34-beta*(t24*t52+t57*t76+t31*t80+t82*t85)+(2.0*t2*(0.12E2*phi+0.12E2*
old+12.0*alpha-12.0)+t98*t29+t100*t34-beta*(t98*t52+3.0*t107*t76+3.0*t24*(
0.882498173E1*t132+0.3776949278E2*t135+0.289445110522E2*t137+0.289445110522E2*
t139)+t57*(0.882498173E1*t172+0.1133084783E3*t175+0.6671400383E2*t178+
0.289445110522E2*t180+0.8683353315E2*t182+0.289445110522E2*t184)+t100*t80-3.0*
t107*t85+3.0*t31*(0.2889626582E2*t132+0.3839831904E2*t135+0.950205321859E1*t137
+0.950205321859E1*t139)+t82*(0.2889626582E2*t172+0.1151949571E3*t175+
0.4790037226E2*t178+0.950205321859E1*t180+0.2850615966E2*t182+0.950205321859E1*
t184)))*(phi-old)/24.0);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi,double old,double mid)
{
  double t100;
  double t105;
  double t107;
  double t109;
  double t112;
  double t113;
  double t114;
  double t120;
  double t122;
  double t123;
  double t125;
  double t126;
  double t128;
  double t130;
  double t132;
  double t136;
  double t138;
  double t141;
  double t143;
  double t146;
  double t148;
  double t153;
  double t16;
  double t160;
  double t162;
  double t168;
  double t173;
  double t180;
  double t187;
  double t191;
  double t195;
  double t198;
  double t2;
  double t200;
  double t201;
  double t203;
  double t204;
  double t206;
  double t208;
  double t21;
  double t210;
  double t216;
  double t218;
  double t220;
  double t223;
  double t228;
  double t23;
  double t231;
  double t232;
  double t233;
  double t234;
  double t237;
  double t24;
  double t241;
  double t246;
  double t248;
  double t249;
  double t251;
  double t252;
  double t254;
  double t256;
  double t257;
  double t259;
  double t26;
  double t261;
  double t263;
  double t264;
  double t266;
  double t274;
  double t28;
  double t282;
  double t285;
  double t289;
  double t290;
  double t291;
  double t293;
  double t294;
  double t296;
  double t297;
  double t299;
  double t300;
  double t302;
  double t304;
  double t306;
  double t308;
  double t31;
  double t33;
  double t34;
  double t368;
  double t370;
  double t372;
  double t373;
  double t375;
  double t376;
  double t378;
  double t379;
  double t38;
  double t381;
  double t383;
  double t385;
  double t386;
  double t388;
  double t390;
  double t392;
  double t394;
  double t396;
  double t398;
  double t400;
  double t402;
  double t405;
  double t410;
  double t417;
  double t42;
  double t43;
  double t437;
  double t45;
  double t453;
  double t455;
  double t47;
  double t48;
  double t5;
  double t50;
  double t55;
  double t59;
  double t6;
  double t61;
  double t62;
  double t63;
  double t67;
  double t69;
  double t70;
  double t72;
  double t74;
  double t79;
  double t83;
  double t86;
  double t88;
  double t90;
  double t91;
  double t93;
  double t95;
  double t97;
  double t98;
  double t99;
  {
    t2 = 1/delta*A;
    t5 = 0.5*phi+0.5*old;
    t6 = t5*t5;
    t16 = t6*t5;
    t21 = 0.6E2*t16-0.9E2*t6+0.15E2*phi+0.15E2*old;
    t23 = log(rho);
    t24 = rho*t23;
    t26 = 0.2133795654E2-0.2011952932E2*rho+0.289445110522E2*t24;
    t28 = -t21;
    t31 = 0.4540504209+0.193942126E2*rho+0.950205321859E1*t24;
    t33 = t6*t6;
    t34 = t33*t5;
    t38 = -0.2889626582E2+0.1204277045E3*t34-0.3010692614E3*t33+0.2007128409E3*
t16;
    t42 = 0.9502053219E1+0.116654747E3*t34-0.2916368674E3*t33+0.1944245783E3*
t16;
    t43 = 1/t42;
    t45 = exp(t38*t43);
    t47 = log(t45);
    t48 = t45*t47;
    t50 = 0.2133795654E2-0.2011952932E2*t45+0.289445110522E2*t48;
    t55 = 30.0*t33-60.0*t16+30.0*t6;
    t59 = 0.3010692612E3*t33-0.6021385228E3*t16+0.3010692614E3*t6;
    t61 = t42*t42;
    t62 = 1/t61;
    t63 = t38*t62;
    t67 = 0.2916368675E3*t33-0.5832737348E3*t16+0.2916368674E3*t6;
    t69 = t59*t43-t63*t67;
    t70 = t69*t45;
    t72 = t70*t47;
    t74 = 0.882498173E1*t70+0.289445110522E2*t72;
    t79 = 0.15E2*t33-0.3E2*t16+0.15E2*t6;
    t83 = 0.6021385225E3*t33-0.1204277046E4*t16+0.6021385227E3*t6;
    t86 = 0.116654747E4*t16;
    t88 = 0.583273735E3*t33-t86+0.5832737349E3*t6;
    t90 = t83*t43-t63*t88;
    t91 = t90*t45;
    t93 = t91*t47;
    t95 = 0.882498173E1*t91+0.289445110522E2*t93;
    t97 = 6.0*t34;
    t98 = 15.0*t33;
    t99 = 10.0*t16;
    t100 = t97-t98+t99;
    t105 = 0.1204277045E4*t16-0.1806415569E4*t6+0.3010692614E3*phi+
0.3010692614E3*old;
    t107 = t83*t62;
    t109 = t59*t62;
    t112 = 1/t61/t42;
    t113 = t38*t112;
    t114 = t88*t67;
    t120 = t86-0.1749821205E4*t6+0.2916368674E3*phi+0.2916368674E3*old;
    t122 = t105*t43-t107*t67-t109*t88+2.0*t113*t114-t63*t120;
    t123 = t122*t45;
    t125 = t90*t69;
    t126 = t125*t45;
    t128 = t123*t47;
    t130 = t125*t48;
    t132 = 0.882498173E1*t123+0.3776949278E2*t126+0.289445110522E2*t128+
0.289445110522E2*t130;
    t136 = 0.4540504209+0.193942126E2*t45+0.950205321859E1*t48;
    t138 = -t55;
    t141 = 0.2889626582E2*t70+0.950205321859E1*t72;
    t143 = -t79;
    t146 = 0.2889626582E2*t91+0.950205321859E1*t93;
    t148 = 1.0-t97+t98-t99;
    t153 = 0.2889626582E2*t123+0.3839831904E2*t126+0.950205321859E1*t128+
0.950205321859E1*t130;
    t160 = 0.18E3*phi+0.18E3*old-0.18E3;
    t162 = -t160;
    t168 = 360.0*t6-0.18E3*phi-0.18E3*old+60.0;
    t173 = 0.18E3*t6-0.9E2*phi-0.9E2*old+0.3E2;
    t180 = 120.0*t16-180.0*t6+0.3E2*phi+0.3E2*old;
    t187 = 0.240855409E4*t16-0.3612831138E4*t6+0.6021385225E3*phi+
0.6021385225E3*old;
    t191 = t88*t88;
    t195 = 0.349964241E4*t6;
    t198 = 0.233309494E4*t16-t195+0.583273735E3*phi+0.583273735E3*old;
    t200 = t187*t43-2.0*t107*t88+2.0*t113*t191-t63*t198;
    t201 = t200*t45;
    t203 = t90*t90;
    t204 = t203*t45;
    t206 = t201*t47;
    t208 = t204*t47;
    t210 = 0.882498173E1*t201+0.3776949278E2*t204+0.289445110522E2*t206+
0.289445110522E2*t208;
    t216 = 0.3612831135E4*t6-0.1806415569E4*phi-0.1806415569E4*old+
0.6021385225E3;
    t218 = t187*t62;
    t220 = t105*t62;
    t223 = t83*t112;
    t228 = t59*t112;
    t231 = t61*t61;
    t232 = 1/t231;
    t233 = t38*t232;
    t234 = t191*t67;
    t237 = t88*t120;
    t241 = t198*t67;
    t246 = t195-0.1749821205E4*phi-0.1749821205E4*old+0.583273735E3;
    t248 = t216*t43-t218*t67-2.0*t220*t88+4.0*t223*t114-2.0*t107*t120+2.0*t228*
t191-6.0*t233*t234+4.0*t113*t237-t109*t198+2.0*t113*t241-t63*t246;
    t249 = t248*t45;
    t251 = t200*t69;
    t252 = t251*t45;
    t254 = t91*t122;
    t256 = t203*t69;
    t257 = t256*t45;
    t259 = t249*t47;
    t261 = t251*t48;
    t263 = t47*t122;
    t264 = t91*t263;
    t266 = t256*t48;
    t274 = 0.722566227E4*t6-0.3612831138E4*phi-0.3612831138E4*old+
0.1204277045E4;
    t282 = t191*t88;
    t285 = t88*t198;
    t289 = 0.349964241E4*phi;
    t290 = 0.349964241E4*old;
    t291 = 0.699928482E4*t6-t289-t290+0.116654747E4;
    t293 = t274*t43-3.0*t218*t88+6.0*t223*t191-3.0*t107*t198-6.0*t233*t282+6.0*
t113*t285-t63*t291;
    t294 = t293*t45;
    t296 = t200*t90;
    t297 = t296*t45;
    t299 = t203*t90;
    t300 = t299*t45;
    t302 = t294*t47;
    t304 = t296*t48;
    t306 = t300*t47;
    t308 = 0.882498173E1*t294+0.1133084783E3*t297+0.6671400383E2*t300+
0.289445110522E2*t302+0.8683353315E2*t304+0.289445110522E2*t306;
    t368 = -3.0*t107*t246-6.0*t59*t232*t282+24.0*t38/t231/t42*t282*t67-18.0*
t233*t191*t120+6.0*t228*t285-18.0*t233*t285*t67+6.0*t113*t120*t198+6.0*t113*t88
*t246-t109*t291+2.0*t113*t291*t67-t63*(-0.349964241E4+t289+t290);
    t370 = ((-0.3612831138E4+0.3612831135E4*phi+0.3612831135E4*old)*t43-t274*
t62*t67-3.0*t216*t62*t88+6.0*t187*t112*t114-3.0*t218*t120+6.0*t105*t112*t191
-18.0*t83*t232*t234+12.0*t223*t237-3.0*t220*t198+6.0*t223*t241+t368)*t45;
    t372 = t293*t69;
    t373 = t372*t45;
    t375 = t248*t90;
    t376 = t375*t45;
    t378 = t200*t122;
    t379 = t378*t45;
    t381 = t296*t70;
    t383 = t204*t122;
    t385 = t299*t69;
    t386 = t385*t45;
    t388 = t370*t47;
    t390 = t372*t48;
    t392 = t375*t48;
    t394 = t378*t48;
    t396 = t296*t72;
    t398 = t204*t263;
    t400 = t385*t48;
    t402 = 0.882498173E1*t370+0.3776949278E2*t373+0.1133084783E3*t376+
0.1133084783E3*t379+0.2001420114E3*t381+0.2001420115E3*t383+0.9565851488E2*t386
+0.289445110522E2*t388+0.289445110522E2*t390+0.8683353315E2*t392+0.8683353315E2
*t394+0.8683353315E2*t396+0.8683353315E2*t398+0.289445110522E2*t400;
    t405 = -t168;
    t410 = -t180;
    t417 = 0.2889626582E2*t201+0.3839831904E2*t204+0.950205321859E1*t206+
0.950205321859E1*t208;
    t437 = 0.2889626582E2*t294+0.1151949571E3*t297+0.4790037226E2*t300+
0.950205321859E1*t302+0.2850615966E2*t304+0.950205321859E1*t306;
    t453 = 0.2889626582E2*t370+0.3839831904E2*t373+0.1151949571E3*t376+
0.1151949571E3*t379+0.1437011168E3*t381+0.1437011168E3*t383+0.5740242548E2*t386
+0.950205321859E1*t388+0.950205321859E1*t390+0.2850615966E2*t392+0.2850615966E2
*t394+0.2850615966E2*t396+0.2850615966E2*t398+0.950205321859E1*t400;
    t455 = t160*t50+t168*t74+3.0*t173*t95+3.0*t180*t132+3.0*t21*t210+3.0*t55*(
0.882498173E1*t249+0.3776949278E2*t252+0.7553898556E2*t254+0.6671400383E2*t257+
0.289445110522E2*t259+0.289445110522E2*t261+0.578890221E2*t264+0.289445110522E2
*t266)+t79*t308+t100*t402+t162*t136+t405*t141-3.0*t173*t146+3.0*t410*t153+3.0*
t28*t417+3.0*t138*(0.2889626582E2*t249+0.3839831904E2*t252+0.7679663808E2*t254+
0.4790037226E2*t257+0.950205321859E1*t259+0.950205321859E1*t261+0.1900410644E2*
t264+0.950205321859E1*t266)+t143*t437+t148*t453;
    return(2.0*t2*(0.6E1*t6+0.3E1*(2.0*alpha-2.0)*t5-0.3E1*alpha+0.1E1)+t21*t26
+t28*t31-beta*(t21*t50+t55*t74+t79*t95+t100*t132+t28*t136+t138*t141+t143*t146+
t148*t153)+(0.24E2*t2+t160*t26+t162*t31-beta*t455)*(phi-old)/24.0+t2*(0.12E2*
phi+0.12E2*old+12.0*alpha-12.0)/12.0+t168*t26/24.0+t405*t31/24.0-beta*(t168*t50
+3.0*t180*t95+3.0*t55*t210+t100*t308+t405*t136+3.0*t410*t146+3.0*t138*t417+t148
*t437)/24.0);
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t11;
  double t12;
  double t16;
  double t2;
  double t20;
  double t21;
  double t4;
  double t5;
  double t7;
  double t8;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t8 = t4-t5+t7;
    t11 = 0.5*rho+0.5*old;
    t12 = log(t11);
    t16 = 1.0-t4+t5-t7;
    t20 = t11*t11;
    t21 = 1/t20;
    return(t8*(0.882498173E1+0.289445110522E2*t12)+t16*(0.2889626582E2+
0.950205321859E1*t12)+(-0.2894451105E2*t8*t21-0.9502053219E1*t16*t21)*(rho-old)
/24.0);
  }
}

double drhochemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t11;
  double t12;
  double t15;
  double t18;
  double t2;
  double t20;
  double t29;
  double t4;
  double t5;
  double t7;
  double t8;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t8 = t4-t5+t7;
    t11 = 0.5*rho+0.5*old;
    t12 = 1/t11;
    t15 = 1.0-t4+t5-t7;
    t18 = t11*t11;
    t20 = 1/t18/t11;
    t29 = 1/t18;
    return(0.1447225552E2*t8*t12+0.475102661E1*t15*t12+(0.2894451105E2*t8*t20+
0.9502053219E1*t15*t20)*(rho-old)/24.0-0.1206021294E1*t8*t29-0.3959188841*t15*
t29);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t10;
  double t11;
  double t15;
  double t19;
  double t2;
  double t20;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t7 = 30.0*t2-60.0*t1*phi+30.0*t1;
    t10 = 0.5*rho+0.5*old;
    t11 = log(t10);
    t15 = -t7;
    t19 = t10*t10;
    t20 = 1/t19;
    return(t7*(0.882498173E1+0.289445110522E2*t11)+t15*(0.2889626582E2+
0.950205321859E1*t11)+(-0.2894451105E2*t7*t20-0.9502053219E1*t15*t20)*(rho-old)
/24.0);
  }
}

#include <math.h>
double pressure(double rho,double phi)
{
  double t1;
  double t10;
  double t11;
  double t2;
  double t4;
  double t5;
  double t6;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t6 = t1*phi;
    t7 = 10.0*t6;
    t10 = log(rho);
    t11 = rho*t10;
    return((t4-t5+t7)*(-0.2094881058E2+0.4291934897E2*rho-0.2894451105E2*t11+
rho*(-0.1397483792E2+0.289445110522E2*t10))+(1.0-t4+t5-t7)*(-0.6490446091E-1+
0.3405607049E1*rho-0.9502053219E1*t11+rho*(0.609644617E1+0.950205321859E1*t10))
-2.0/delta*A*(t2+(2.0*alpha-2.0)*t6+(-3.0*alpha+1.0)*t1+alpha));
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t1;
  double t2;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t9 = sqrt(0.116654747E11*t2*phi-0.2916368675E11*t2+0.1944245783E11*t1*phi+
950205322.0);
    return(0.1E-3*t9);
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t15;
  double t16;
  double t17;
  double t18;
  double t19;
  double t21;
  double t22;
  double t26;
  double t3;
  double t4;
  double t41;
  double t43;
  double t44;
  double t7;
  {
    t3 = phi*phi;
    t4 = t3*t3;
    t7 = t3*phi;
    t15 = t4*phi;
    t16 = 6.0*t15;
    t17 = 15.0*t4;
    t18 = 10.0*t7;
    t19 = t16-t17+t18;
    t21 = log(rho);
    t22 = rho*t21;
    t26 = 1.0-t16+t17-t18;
    t41 = exp((-0.2889626582E2+0.1204277045E3*t15-0.3010692614E3*t4+
0.2007128409E3*t7)/(0.9502053219E1+0.116654747E3*t15-0.2916368674E3*t4+
0.1944245783E3*t7));
    t43 = log(t41);
    t44 = t41*t43;
    return(2.0/delta*A*(t4+(2.0*alpha-2.0)*t7+(-3.0*alpha+1.0)*t3+alpha)+t19*(
0.2133795654E2-0.2011952932E2*rho+0.289445110522E2*t22)+t26*(0.4540504209+
0.193942126E2*rho+0.950205321859E1*t22)-beta*(t19*(0.2133795654E2
-0.2011952932E2*t41+0.289445110522E2*t44)+t26*(0.4540504209+0.193942126E2*t41+
0.950205321859E1*t44)));
  }
}

#include <math.h>
double reactionSource(double rho,double phi,double old)
{
  double t100;
  double t107;
  double t114;
  double t116;
  double t120;
  double t121;
  double t122;
  double t129;
  double t131;
  double t132;
  double t134;
  double t135;
  double t137;
  double t139;
  double t157;
  double t172;
  double t174;
  double t175;
  double t178;
  double t180;
  double t182;
  double t184;
  double t2;
  double t20;
  double t24;
  double t26;
  double t27;
  double t29;
  double t31;
  double t34;
  double t36;
  double t40;
  double t44;
  double t45;
  double t47;
  double t49;
  double t5;
  double t50;
  double t52;
  double t54;
  double t55;
  double t56;
  double t57;
  double t6;
  double t61;
  double t63;
  double t64;
  double t65;
  double t69;
  double t7;
  double t71;
  double t72;
  double t74;
  double t76;
  double t80;
  double t82;
  double t85;
  double t98;
  {
    t2 = 1/delta*A;
    t5 = 0.5*phi+0.5*old;
    t6 = t5*t5;
    t7 = t6*t5;
    t20 = t6*t6;
    t24 = 30.0*t20-60.0*t7+30.0*t6;
    t26 = log(rho);
    t27 = rho*t26;
    t29 = -0.2011952932E2*rho+0.289445110522E2*t27+0.2133795654E2;
    t31 = -t24;
    t34 = 0.193942126E2*rho+0.950205321859E1*t27+0.4540504209;
    t36 = t20*t5;
    t40 = -0.2889626582E2+0.1204277045E3*t36-0.3010692614E3*t20+0.2007128409E3*
t7;
    t44 = 0.9502053219E1+0.116654747E3*t36-0.2916368674E3*t20+0.1944245783E3*t7
;
    t45 = 1/t44;
    t47 = exp(t40*t45);
    t49 = log(t47);
    t50 = t47*t49;
    t52 = 0.2133795654E2-0.2011952932E2*t47+0.289445110522E2*t50;
    t54 = 6.0*t36;
    t55 = 15.0*t20;
    t56 = 10.0*t7;
    t57 = t54-t55+t56;
    t61 = 0.6021385225E3*t20-0.1204277046E4*t7+0.6021385227E3*t6;
    t63 = t44*t44;
    t64 = 1/t63;
    t65 = t40*t64;
    t69 = 0.583273735E3*t20-0.116654747E4*t7+0.5832737349E3*t6;
    t71 = t61*t45-t65*t69;
    t72 = t71*t47;
    t74 = t72*t49;
    t76 = 0.882498173E1*t72+0.289445110522E2*t74;
    t80 = 0.4540504209+0.193942126E2*t47+0.950205321859E1*t50;
    t82 = 1.0-t54+t55-t56;
    t85 = 0.2889626582E2*t72+0.950205321859E1*t74;
    t98 = 360.0*t6-0.18E3*phi-0.18E3*old+60.0;
    t100 = -t98;
    t107 = 120.0*t7-180.0*t6+0.3E2*phi+0.3E2*old;
    t114 = 0.240855409E4*t7-0.3612831138E4*t6+0.6021385225E3*phi+0.6021385225E3
*old;
    t116 = t61*t64;
    t120 = 1/t63/t44;
    t121 = t40*t120;
    t122 = t69*t69;
    t129 = 0.233309494E4*t7-0.349964241E4*t6+0.583273735E3*phi+0.583273735E3*
old;
    t131 = t114*t45-2.0*t116*t69+2.0*t121*t122-t65*t129;
    t132 = t131*t47;
    t134 = t71*t71;
    t135 = t134*t47;
    t137 = t132*t49;
    t139 = t135*t49;
    t157 = t63*t63;
    t172 = ((0.722566227E4*t6-0.3612831138E4*phi-0.3612831138E4*old+
0.1204277045E4)*t45-3.0*t114*t64*t69+6.0*t61*t120*t122-3.0*t116*t129-6.0*t40/
t157*t122*t69+6.0*t121*t69*t129-t65*(0.699928482E4*t6-0.349964241E4*phi
-0.349964241E4*old+0.116654747E4))*t47;
    t174 = t131*t71;
    t175 = t174*t47;
    t178 = t134*t71*t47;
    t180 = t172*t49;
    t182 = t174*t50;
    t184 = t178*t49;
    return(2.0*t2*(4.0*t7+3.0*(2.0*alpha-2.0)*t6+2.0*(-3.0*alpha+1.0)*t5)+t24*
t29+t31*t34-beta*(t24*t52+t57*t76+t31*t80+t82*t85)+(2.0*t2*(0.12E2*phi+0.12E2*
old+12.0*alpha-12.0)+t98*t29+t100*t34-beta*(t98*t52+3.0*t107*t76+3.0*t24*(
0.882498173E1*t132+0.3776949278E2*t135+0.289445110522E2*t137+0.289445110522E2*
t139)+t57*(0.882498173E1*t172+0.1133084783E3*t175+0.6671400383E2*t178+
0.289445110522E2*t180+0.8683353315E2*t182+0.289445110522E2*t184)+t100*t80-3.0*
t107*t85+3.0*t31*(0.2889626582E2*t132+0.3839831904E2*t135+0.950205321859E1*t137
+0.950205321859E1*t139)+t82*(0.2889626582E2*t172+0.1151949571E3*t175+
0.4790037226E2*t178+0.950205321859E1*t180+0.2850615966E2*t182+0.950205321859E1*
t184)))*(phi-old)/24.0);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi,double old,double mid)
{
  double t100;
  double t105;
  double t107;
  double t109;
  double t112;
  double t113;
  double t114;
  double t120;
  double t122;
  double t123;
  double t125;
  double t126;
  double t128;
  double t130;
  double t132;
  double t136;
  double t138;
  double t141;
  double t143;
  double t146;
  double t148;
  double t153;
  double t16;
  double t160;
  double t162;
  double t168;
  double t173;
  double t180;
  double t187;
  double t191;
  double t195;
  double t198;
  double t2;
  double t200;
  double t201;
  double t203;
  double t204;
  double t206;
  double t208;
  double t21;
  double t210;
  double t216;
  double t218;
  double t220;
  double t223;
  double t228;
  double t23;
  double t231;
  double t232;
  double t233;
  double t234;
  double t237;
  double t24;
  double t241;
  double t246;
  double t248;
  double t249;
  double t251;
  double t252;
  double t254;
  double t256;
  double t257;
  double t259;
  double t26;
  double t261;
  double t263;
  double t264;
  double t266;
  double t274;
  double t28;
  double t282;
  double t285;
  double t289;
  double t290;
  double t291;
  double t293;
  double t294;
  double t296;
  double t297;
  double t299;
  double t300;
  double t302;
  double t304;
  double t306;
  double t308;
  double t31;
  double t33;
  double t34;
  double t368;
  double t370;
  double t372;
  double t373;
  double t375;
  double t376;
  double t378;
  double t379;
  double t38;
  double t381;
  double t383;
  double t385;
  double t386;
  double t388;
  double t390;
  double t392;
  double t394;
  double t396;
  double t398;
  double t400;
  double t402;
  double t405;
  double t410;
  double t417;
  double t42;
  double t43;
  double t437;
  double t45;
  double t453;
  double t455;
  double t47;
  double t48;
  double t5;
  double t50;
  double t55;
  double t59;
  double t6;
  double t61;
  double t62;
  double t63;
  double t67;
  double t69;
  double t70;
  double t72;
  double t74;
  double t79;
  double t83;
  double t86;
  double t88;
  double t90;
  double t91;
  double t93;
  double t95;
  double t97;
  double t98;
  double t99;
  {
    t2 = 1/delta*A;
    t5 = 0.5*phi+0.5*old;
    t6 = t5*t5;
    t16 = t6*t5;
    t21 = 0.6E2*t16-0.9E2*t6+0.15E2*phi+0.15E2*old;
    t23 = log(rho);
    t24 = rho*t23;
    t26 = -0.2011952932E2*rho+0.289445110522E2*t24+0.2133795654E2;
    t28 = -t21;
    t31 = 0.193942126E2*rho+0.950205321859E1*t24+0.4540504209;
    t33 = t6*t6;
    t34 = t33*t5;
    t38 = -0.2889626582E2+0.1204277045E3*t34-0.3010692614E3*t33+0.2007128409E3*
t16;
    t42 = 0.9502053219E1+0.116654747E3*t34-0.2916368674E3*t33+0.1944245783E3*
t16;
    t43 = 1/t42;
    t45 = exp(t38*t43);
    t47 = log(t45);
    t48 = t45*t47;
    t50 = 0.2133795654E2-0.2011952932E2*t45+0.289445110522E2*t48;
    t55 = 30.0*t33-60.0*t16+30.0*t6;
    t59 = 0.3010692612E3*t33-0.6021385228E3*t16+0.3010692614E3*t6;
    t61 = t42*t42;
    t62 = 1/t61;
    t63 = t38*t62;
    t67 = 0.2916368675E3*t33-0.5832737348E3*t16+0.2916368674E3*t6;
    t69 = t59*t43-t63*t67;
    t70 = t69*t45;
    t72 = t70*t47;
    t74 = 0.882498173E1*t70+0.289445110522E2*t72;
    t79 = 0.15E2*t33-0.3E2*t16+0.15E2*t6;
    t83 = 0.6021385225E3*t33-0.1204277046E4*t16+0.6021385227E3*t6;
    t86 = 0.116654747E4*t16;
    t88 = 0.583273735E3*t33-t86+0.5832737349E3*t6;
    t90 = t83*t43-t63*t88;
    t91 = t90*t45;
    t93 = t91*t47;
    t95 = 0.882498173E1*t91+0.289445110522E2*t93;
    t97 = 6.0*t34;
    t98 = 15.0*t33;
    t99 = 10.0*t16;
    t100 = t97-t98+t99;
    t105 = 0.1204277045E4*t16-0.1806415569E4*t6+0.3010692614E3*phi+
0.3010692614E3*old;
    t107 = t83*t62;
    t109 = t59*t62;
    t112 = 1/t61/t42;
    t113 = t38*t112;
    t114 = t88*t67;
    t120 = t86-0.1749821205E4*t6+0.2916368674E3*phi+0.2916368674E3*old;
    t122 = t105*t43-t107*t67-t109*t88+2.0*t113*t114-t63*t120;
    t123 = t122*t45;
    t125 = t90*t69;
    t126 = t125*t45;
    t128 = t123*t47;
    t130 = t125*t48;
    t132 = 0.882498173E1*t123+0.3776949278E2*t126+0.289445110522E2*t128+
0.289445110522E2*t130;
    t136 = 0.4540504209+0.193942126E2*t45+0.950205321859E1*t48;
    t138 = -t55;
    t141 = 0.2889626582E2*t70+0.950205321859E1*t72;
    t143 = -t79;
    t146 = 0.2889626582E2*t91+0.950205321859E1*t93;
    t148 = 1.0-t97+t98-t99;
    t153 = 0.2889626582E2*t123+0.3839831904E2*t126+0.950205321859E1*t128+
0.950205321859E1*t130;
    t160 = 0.18E3*phi+0.18E3*old-0.18E3;
    t162 = -t160;
    t168 = 360.0*t6-0.18E3*phi-0.18E3*old+60.0;
    t173 = 0.18E3*t6-0.9E2*phi-0.9E2*old+0.3E2;
    t180 = 120.0*t16-180.0*t6+0.3E2*phi+0.3E2*old;
    t187 = 0.240855409E4*t16-0.3612831138E4*t6+0.6021385225E3*phi+
0.6021385225E3*old;
    t191 = t88*t88;
    t195 = 0.349964241E4*t6;
    t198 = 0.233309494E4*t16-t195+0.583273735E3*phi+0.583273735E3*old;
    t200 = t187*t43-2.0*t107*t88+2.0*t113*t191-t63*t198;
    t201 = t200*t45;
    t203 = t90*t90;
    t204 = t203*t45;
    t206 = t201*t47;
    t208 = t204*t47;
    t210 = 0.882498173E1*t201+0.3776949278E2*t204+0.289445110522E2*t206+
0.289445110522E2*t208;
    t216 = 0.3612831135E4*t6-0.1806415569E4*phi-0.1806415569E4*old+
0.6021385225E3;
    t218 = t187*t62;
    t220 = t105*t62;
    t223 = t83*t112;
    t228 = t59*t112;
    t231 = t61*t61;
    t232 = 1/t231;
    t233 = t38*t232;
    t234 = t191*t67;
    t237 = t88*t120;
    t241 = t198*t67;
    t246 = t195-0.1749821205E4*phi-0.1749821205E4*old+0.583273735E3;
    t248 = t216*t43-t218*t67-2.0*t220*t88+4.0*t223*t114-2.0*t107*t120+2.0*t228*
t191-6.0*t233*t234+4.0*t113*t237-t109*t198+2.0*t113*t241-t63*t246;
    t249 = t248*t45;
    t251 = t200*t69;
    t252 = t251*t45;
    t254 = t91*t122;
    t256 = t203*t69;
    t257 = t256*t45;
    t259 = t249*t47;
    t261 = t251*t48;
    t263 = t47*t122;
    t264 = t91*t263;
    t266 = t256*t48;
    t274 = 0.722566227E4*t6-0.3612831138E4*phi-0.3612831138E4*old+
0.1204277045E4;
    t282 = t191*t88;
    t285 = t88*t198;
    t289 = 0.349964241E4*phi;
    t290 = 0.349964241E4*old;
    t291 = 0.699928482E4*t6-t289-t290+0.116654747E4;
    t293 = t274*t43-3.0*t218*t88+6.0*t223*t191-3.0*t107*t198-6.0*t233*t282+6.0*
t113*t285-t63*t291;
    t294 = t293*t45;
    t296 = t200*t90;
    t297 = t296*t45;
    t299 = t203*t90;
    t300 = t299*t45;
    t302 = t294*t47;
    t304 = t296*t48;
    t306 = t300*t47;
    t308 = 0.882498173E1*t294+0.1133084783E3*t297+0.6671400383E2*t300+
0.289445110522E2*t302+0.8683353315E2*t304+0.289445110522E2*t306;
    t368 = -3.0*t107*t246-6.0*t59*t232*t282+24.0*t38/t231/t42*t282*t67-18.0*
t233*t191*t120+6.0*t228*t285-18.0*t233*t285*t67+6.0*t113*t120*t198+6.0*t113*t88
*t246-t109*t291+2.0*t113*t291*t67-t63*(-0.349964241E4+t289+t290);
    t370 = ((0.3612831135E4*phi+0.3612831135E4*old-0.3612831138E4)*t43-t274*t62
*t67-3.0*t216*t62*t88+6.0*t187*t112*t114-3.0*t218*t120+6.0*t105*t112*t191-18.0*
t83*t232*t234+12.0*t223*t237-3.0*t220*t198+6.0*t223*t241+t368)*t45;
    t372 = t293*t69;
    t373 = t372*t45;
    t375 = t248*t90;
    t376 = t375*t45;
    t378 = t200*t122;
    t379 = t378*t45;
    t381 = t296*t70;
    t383 = t204*t122;
    t385 = t299*t69;
    t386 = t385*t45;
    t388 = t370*t47;
    t390 = t372*t48;
    t392 = t375*t48;
    t394 = t378*t48;
    t396 = t296*t72;
    t398 = t204*t263;
    t400 = t385*t48;
    t402 = 0.882498173E1*t370+0.3776949278E2*t373+0.1133084783E3*t376+
0.1133084783E3*t379+0.2001420114E3*t381+0.2001420115E3*t383+0.9565851488E2*t386
+0.289445110522E2*t388+0.289445110522E2*t390+0.8683353315E2*t392+0.8683353315E2
*t394+0.8683353315E2*t396+0.8683353315E2*t398+0.289445110522E2*t400;
    t405 = -t168;
    t410 = -t180;
    t417 = 0.2889626582E2*t201+0.3839831904E2*t204+0.950205321859E1*t206+
0.950205321859E1*t208;
    t437 = 0.2889626582E2*t294+0.1151949571E3*t297+0.4790037226E2*t300+
0.950205321859E1*t302+0.2850615966E2*t304+0.950205321859E1*t306;
    t453 = 0.2889626582E2*t370+0.3839831904E2*t373+0.1151949571E3*t376+
0.1151949571E3*t379+0.1437011168E3*t381+0.1437011168E3*t383+0.5740242548E2*t386
+0.950205321859E1*t388+0.950205321859E1*t390+0.2850615966E2*t392+0.2850615966E2
*t394+0.2850615966E2*t396+0.2850615966E2*t398+0.950205321859E1*t400;
    t455 = t160*t50+t168*t74+3.0*t173*t95+3.0*t180*t132+3.0*t21*t210+3.0*t55*(
0.882498173E1*t249+0.3776949278E2*t252+0.7553898556E2*t254+0.6671400383E2*t257+
0.289445110522E2*t259+0.289445110522E2*t261+0.578890221E2*t264+0.289445110522E2
*t266)+t79*t308+t100*t402+t162*t136+t405*t141-3.0*t173*t146+3.0*t410*t153+3.0*
t28*t417+3.0*t138*(0.2889626582E2*t249+0.3839831904E2*t252+0.7679663808E2*t254+
0.4790037226E2*t257+0.950205321859E1*t259+0.950205321859E1*t261+0.1900410644E2*
t264+0.950205321859E1*t266)+t143*t437+t148*t453;
    return(2.0*t2*(0.6E1*t6+0.3E1*(2.0*alpha-2.0)*t5-0.3E1*alpha+0.1E1)+t21*t26
+t28*t31-beta*(t21*t50+t55*t74+t79*t95+t100*t132+t28*t136+t138*t141+t143*t146+
t148*t153)+(0.24E2*t2+t160*t26+t162*t31-beta*t455)*(phi-old)/24.0+t2*(0.12E2*
phi+0.12E2*old+12.0*alpha-12.0)/12.0+t168*t26/24.0+t405*t31/24.0-beta*(t168*t50
+3.0*t180*t95+3.0*t55*t210+t100*t308+t405*t136+3.0*t410*t146+3.0*t138*t417+t148
*t437)/24.0);
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t11;
  double t12;
  double t16;
  double t2;
  double t20;
  double t21;
  double t4;
  double t5;
  double t7;
  double t8;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t8 = t4-t5+t7;
    t11 = 0.5*rho+0.5*old;
    t12 = log(t11);
    t16 = 1.0-t4+t5-t7;
    t20 = t11*t11;
    t21 = 1/t20;
    return(t8*(0.882498173E1+0.289445110522E2*t12)+t16*(0.2889626582E2+
0.950205321859E1*t12)+(-0.2894451105E2*t8*t21-0.9502053219E1*t16*t21)*(rho-old)
/24.0);
  }
}

double drhochemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t11;
  double t12;
  double t15;
  double t18;
  double t2;
  double t20;
  double t29;
  double t4;
  double t5;
  double t7;
  double t8;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t8 = t4-t5+t7;
    t11 = 0.5*rho+0.5*old;
    t12 = 1/t11;
    t15 = 1.0-t4+t5-t7;
    t18 = t11*t11;
    t20 = 1/t18/t11;
    t29 = 1/t18;
    return(0.1447225552E2*t8*t12+0.475102661E1*t15*t12+(0.2894451105E2*t8*t20+
0.9502053219E1*t15*t20)*(rho-old)/24.0-0.1206021294E1*t8*t29-0.3959188841*t15*
t29);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t10;
  double t11;
  double t15;
  double t19;
  double t2;
  double t20;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t7 = 30.0*t2-60.0*t1*phi+30.0*t1;
    t10 = 0.5*rho+0.5*old;
    t11 = log(t10);
    t15 = -t7;
    t19 = t10*t10;
    t20 = 1/t19;
    return(t7*(0.882498173E1+0.289445110522E2*t11)+t15*(0.2889626582E2+
0.950205321859E1*t11)+(-0.2894451105E2*t7*t20-0.9502053219E1*t15*t20)*(rho-old)
/24.0);
  }
}

#include <math.h>
double pressure(double rho,double phi)
{
  double t1;
  double t10;
  double t11;
  double t2;
  double t4;
  double t5;
  double t6;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t6 = t1*phi;
    t7 = 10.0*t6;
    t10 = log(rho);
    t11 = rho*t10;
    return((t4-t5+t7)*(-0.2094881058E2+0.4291934897E2*rho-0.2894451105E2*t11+
rho*(-0.1397483792E2+0.289445110522E2*t10))+(1.0-t4+t5-t7)*(-0.6490446091E-1+
0.3405607049E1*rho-0.9502053219E1*t11+rho*(0.609644617E1+0.950205321859E1*t10))
-2.0/delta*A*(t2+(2.0*alpha-2.0)*t6+(-3.0*alpha+1.0)*t1+alpha));
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t1;
  double t2;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t9 = sqrt(0.116654747E11*t2*phi-0.2916368675E11*t2+0.1944245783E11*t1*phi+
950205322.0);
    return(0.1E-3*t9);
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t15;
  double t16;
  double t17;
  double t18;
  double t19;
  double t21;
  double t22;
  double t26;
  double t3;
  double t4;
  double t41;
  double t43;
  double t44;
  double t7;
  {
    t3 = phi*phi;
    t4 = t3*t3;
    t7 = t3*phi;
    t15 = t4*phi;
    t16 = 6.0*t15;
    t17 = 15.0*t4;
    t18 = 10.0*t7;
    t19 = t16-t17+t18;
    t21 = log(rho);
    t22 = rho*t21;
    t26 = 1.0-t16+t17-t18;
    t41 = exp((-0.2889626582E2+0.1204277045E3*t15-0.3010692614E3*t4+
0.2007128409E3*t7)/(0.9502053219E1+0.116654747E3*t15-0.2916368674E3*t4+
0.1944245783E3*t7));
    t43 = log(t41);
    t44 = t41*t43;
    return(2.0/delta*A*(t4+(2.0*alpha-2.0)*t7+(-3.0*alpha+1.0)*t3+alpha)+t19*(
-0.2011952932E2*rho+0.289445110522E2*t22+0.2133795654E2)+t26*(0.193942126E2*rho
+0.950205321859E1*t22+0.4540504209)-beta*(t19*(0.2133795654E2-0.2011952932E2*
t41+0.289445110522E2*t44)+t26*(0.4540504209+0.193942126E2*t41+0.950205321859E1*
t44)));
  }
}

#include <math.h>
double reactionSource(double rho,double phi,double old)
{
  double t100;
  double t107;
  double t114;
  double t116;
  double t120;
  double t121;
  double t122;
  double t129;
  double t131;
  double t132;
  double t134;
  double t135;
  double t137;
  double t139;
  double t157;
  double t172;
  double t174;
  double t175;
  double t178;
  double t180;
  double t182;
  double t184;
  double t2;
  double t20;
  double t24;
  double t26;
  double t27;
  double t29;
  double t31;
  double t34;
  double t36;
  double t40;
  double t44;
  double t45;
  double t47;
  double t49;
  double t5;
  double t50;
  double t52;
  double t54;
  double t55;
  double t56;
  double t57;
  double t6;
  double t61;
  double t63;
  double t64;
  double t65;
  double t69;
  double t7;
  double t71;
  double t72;
  double t74;
  double t76;
  double t80;
  double t82;
  double t85;
  double t98;
  {
    t2 = 1/delta*A;
    t5 = 0.5*phi+0.5*old;
    t6 = t5*t5;
    t7 = t6*t5;
    t20 = t6*t6;
    t24 = 30.0*t20-60.0*t7+30.0*t6;
    t26 = log(rho);
    t27 = rho*t26;
    t29 = -0.2011952932E2*rho+0.289445110522E2*t27+0.2133795654E2;
    t31 = -t24;
    t34 = 0.193942126E2*rho+0.950205321859E1*t27+0.4540504209;
    t36 = t20*t5;
    t40 = -0.2889626582E2+0.1204277045E3*t36-0.3010692614E3*t20+0.2007128409E3*
t7;
    t44 = 0.9502053219E1+0.116654747E3*t36-0.2916368674E3*t20+0.1944245783E3*t7
;
    t45 = 1/t44;
    t47 = exp(t40*t45);
    t49 = log(t47);
    t50 = t47*t49;
    t52 = 0.2133795654E2-0.2011952932E2*t47+0.289445110522E2*t50;
    t54 = 6.0*t36;
    t55 = 15.0*t20;
    t56 = 10.0*t7;
    t57 = t54-t55+t56;
    t61 = 0.6021385225E3*t20-0.1204277046E4*t7+0.6021385227E3*t6;
    t63 = t44*t44;
    t64 = 1/t63;
    t65 = t40*t64;
    t69 = 0.583273735E3*t20-0.116654747E4*t7+0.5832737349E3*t6;
    t71 = t61*t45-t65*t69;
    t72 = t71*t47;
    t74 = t72*t49;
    t76 = 0.882498173E1*t72+0.289445110522E2*t74;
    t80 = 0.4540504209+0.193942126E2*t47+0.950205321859E1*t50;
    t82 = 1.0-t54+t55-t56;
    t85 = 0.2889626582E2*t72+0.950205321859E1*t74;
    t98 = 360.0*t6-0.18E3*phi-0.18E3*old+60.0;
    t100 = -t98;
    t107 = 120.0*t7-180.0*t6+0.3E2*phi+0.3E2*old;
    t114 = 0.240855409E4*t7-0.3612831138E4*t6+0.6021385225E3*phi+0.6021385225E3
*old;
    t116 = t61*t64;
    t120 = 1/t63/t44;
    t121 = t40*t120;
    t122 = t69*t69;
    t129 = 0.233309494E4*t7-0.349964241E4*t6+0.583273735E3*phi+0.583273735E3*
old;
    t131 = t114*t45-2.0*t116*t69+2.0*t121*t122-t65*t129;
    t132 = t131*t47;
    t134 = t71*t71;
    t135 = t134*t47;
    t137 = t132*t49;
    t139 = t135*t49;
    t157 = t63*t63;
    t172 = ((0.722566227E4*t6-0.3612831138E4*phi-0.3612831138E4*old+
0.1204277045E4)*t45-3.0*t114*t64*t69+6.0*t61*t120*t122-3.0*t116*t129-6.0*t40/
t157*t122*t69+6.0*t121*t69*t129-t65*(0.699928482E4*t6-0.349964241E4*phi
-0.349964241E4*old+0.116654747E4))*t47;
    t174 = t131*t71;
    t175 = t174*t47;
    t178 = t134*t71*t47;
    t180 = t172*t49;
    t182 = t174*t50;
    t184 = t178*t49;
    return(2.0*t2*(4.0*t7+3.0*(2.0*alpha-2.0)*t6+2.0*(-3.0*alpha+1.0)*t5)+t24*
t29+t31*t34-beta*(t24*t52+t57*t76+t31*t80+t82*t85)+(2.0*t2*(0.12E2*phi+0.12E2*
old+12.0*alpha-12.0)+t98*t29+t100*t34-beta*(t98*t52+3.0*t107*t76+3.0*t24*(
0.882498173E1*t132+0.3776949278E2*t135+0.289445110522E2*t137+0.289445110522E2*
t139)+t57*(0.882498173E1*t172+0.1133084783E3*t175+0.6671400383E2*t178+
0.289445110522E2*t180+0.8683353315E2*t182+0.289445110522E2*t184)+t100*t80-3.0*
t107*t85+3.0*t31*(0.2889626582E2*t132+0.3839831904E2*t135+0.950205321859E1*t137
+0.950205321859E1*t139)+t82*(0.2889626582E2*t172+0.1151949571E3*t175+
0.4790037226E2*t178+0.950205321859E1*t180+0.2850615966E2*t182+0.950205321859E1*
t184)))*(phi-old)/24.0);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi,double old,double mid)
{
  double t100;
  double t105;
  double t107;
  double t109;
  double t112;
  double t113;
  double t114;
  double t120;
  double t122;
  double t123;
  double t125;
  double t126;
  double t128;
  double t130;
  double t132;
  double t136;
  double t138;
  double t141;
  double t143;
  double t146;
  double t148;
  double t153;
  double t16;
  double t160;
  double t162;
  double t168;
  double t173;
  double t180;
  double t187;
  double t191;
  double t195;
  double t198;
  double t2;
  double t200;
  double t201;
  double t203;
  double t204;
  double t206;
  double t208;
  double t21;
  double t210;
  double t216;
  double t218;
  double t220;
  double t223;
  double t228;
  double t23;
  double t231;
  double t232;
  double t233;
  double t234;
  double t237;
  double t24;
  double t241;
  double t246;
  double t248;
  double t249;
  double t251;
  double t252;
  double t254;
  double t256;
  double t257;
  double t259;
  double t26;
  double t261;
  double t263;
  double t264;
  double t266;
  double t274;
  double t28;
  double t282;
  double t285;
  double t289;
  double t290;
  double t291;
  double t293;
  double t294;
  double t296;
  double t297;
  double t299;
  double t300;
  double t302;
  double t304;
  double t306;
  double t308;
  double t31;
  double t33;
  double t34;
  double t368;
  double t370;
  double t372;
  double t373;
  double t375;
  double t376;
  double t378;
  double t379;
  double t38;
  double t381;
  double t383;
  double t385;
  double t386;
  double t388;
  double t390;
  double t392;
  double t394;
  double t396;
  double t398;
  double t400;
  double t402;
  double t405;
  double t410;
  double t417;
  double t42;
  double t43;
  double t437;
  double t45;
  double t453;
  double t455;
  double t47;
  double t48;
  double t5;
  double t50;
  double t55;
  double t59;
  double t6;
  double t61;
  double t62;
  double t63;
  double t67;
  double t69;
  double t70;
  double t72;
  double t74;
  double t79;
  double t83;
  double t86;
  double t88;
  double t90;
  double t91;
  double t93;
  double t95;
  double t97;
  double t98;
  double t99;
  {
    t2 = 1/delta*A;
    t5 = 0.5*phi+0.5*old;
    t6 = t5*t5;
    t16 = t6*t5;
    t21 = 0.6E2*t16-0.9E2*t6+0.15E2*phi+0.15E2*old;
    t23 = log(rho);
    t24 = rho*t23;
    t26 = -0.2011952932E2*rho+0.289445110522E2*t24+0.2133795654E2;
    t28 = -t21;
    t31 = 0.193942126E2*rho+0.950205321859E1*t24+0.4540504209;
    t33 = t6*t6;
    t34 = t33*t5;
    t38 = -0.2889626582E2+0.1204277045E3*t34-0.3010692614E3*t33+0.2007128409E3*
t16;
    t42 = 0.9502053219E1+0.116654747E3*t34-0.2916368674E3*t33+0.1944245783E3*
t16;
    t43 = 1/t42;
    t45 = exp(t38*t43);
    t47 = log(t45);
    t48 = t45*t47;
    t50 = 0.2133795654E2-0.2011952932E2*t45+0.289445110522E2*t48;
    t55 = 30.0*t33-60.0*t16+30.0*t6;
    t59 = 0.3010692612E3*t33-0.6021385228E3*t16+0.3010692614E3*t6;
    t61 = t42*t42;
    t62 = 1/t61;
    t63 = t38*t62;
    t67 = 0.2916368675E3*t33-0.5832737348E3*t16+0.2916368674E3*t6;
    t69 = t59*t43-t63*t67;
    t70 = t69*t45;
    t72 = t70*t47;
    t74 = 0.882498173E1*t70+0.289445110522E2*t72;
    t79 = 0.15E2*t33-0.3E2*t16+0.15E2*t6;
    t83 = 0.6021385225E3*t33-0.1204277046E4*t16+0.6021385227E3*t6;
    t86 = 0.116654747E4*t16;
    t88 = 0.583273735E3*t33-t86+0.5832737349E3*t6;
    t90 = t83*t43-t63*t88;
    t91 = t90*t45;
    t93 = t91*t47;
    t95 = 0.882498173E1*t91+0.289445110522E2*t93;
    t97 = 6.0*t34;
    t98 = 15.0*t33;
    t99 = 10.0*t16;
    t100 = t97-t98+t99;
    t105 = 0.1204277045E4*t16-0.1806415569E4*t6+0.3010692614E3*phi+
0.3010692614E3*old;
    t107 = t83*t62;
    t109 = t59*t62;
    t112 = 1/t61/t42;
    t113 = t38*t112;
    t114 = t88*t67;
    t120 = t86-0.1749821205E4*t6+0.2916368674E3*phi+0.2916368674E3*old;
    t122 = t105*t43-t107*t67-t109*t88+2.0*t113*t114-t63*t120;
    t123 = t122*t45;
    t125 = t90*t69;
    t126 = t125*t45;
    t128 = t123*t47;
    t130 = t125*t48;
    t132 = 0.882498173E1*t123+0.3776949278E2*t126+0.289445110522E2*t128+
0.289445110522E2*t130;
    t136 = 0.4540504209+0.193942126E2*t45+0.950205321859E1*t48;
    t138 = -t55;
    t141 = 0.2889626582E2*t70+0.950205321859E1*t72;
    t143 = -t79;
    t146 = 0.2889626582E2*t91+0.950205321859E1*t93;
    t148 = 1.0-t97+t98-t99;
    t153 = 0.2889626582E2*t123+0.3839831904E2*t126+0.950205321859E1*t128+
0.950205321859E1*t130;
    t160 = -0.18E3+0.18E3*phi+0.18E3*old;
    t162 = -t160;
    t168 = 360.0*t6-0.18E3*phi-0.18E3*old+60.0;
    t173 = 0.3E2+0.18E3*t6-0.9E2*phi-0.9E2*old;
    t180 = 120.0*t16-180.0*t6+0.3E2*phi+0.3E2*old;
    t187 = 0.240855409E4*t16-0.3612831138E4*t6+0.6021385225E3*phi+
0.6021385225E3*old;
    t191 = t88*t88;
    t195 = 0.349964241E4*t6;
    t198 = 0.233309494E4*t16-t195+0.583273735E3*phi+0.583273735E3*old;
    t200 = t187*t43-2.0*t107*t88+2.0*t113*t191-t63*t198;
    t201 = t200*t45;
    t203 = t90*t90;
    t204 = t203*t45;
    t206 = t201*t47;
    t208 = t204*t47;
    t210 = 0.882498173E1*t201+0.3776949278E2*t204+0.289445110522E2*t206+
0.289445110522E2*t208;
    t216 = 0.6021385225E3+0.3612831135E4*t6-0.1806415569E4*phi-0.1806415569E4*
old;
    t218 = t187*t62;
    t220 = t105*t62;
    t223 = t83*t112;
    t228 = t59*t112;
    t231 = t61*t61;
    t232 = 1/t231;
    t233 = t38*t232;
    t234 = t191*t67;
    t237 = t88*t120;
    t241 = t198*t67;
    t246 = 0.583273735E3+t195-0.1749821205E4*phi-0.1749821205E4*old;
    t248 = t216*t43-t218*t67-2.0*t220*t88+4.0*t223*t114-2.0*t107*t120+2.0*t228*
t191-6.0*t233*t234+4.0*t113*t237-t109*t198+2.0*t113*t241-t63*t246;
    t249 = t248*t45;
    t251 = t200*t69;
    t252 = t251*t45;
    t254 = t91*t122;
    t256 = t203*t69;
    t257 = t256*t45;
    t259 = t249*t47;
    t261 = t251*t48;
    t263 = t47*t122;
    t264 = t91*t263;
    t266 = t256*t48;
    t274 = 0.722566227E4*t6-0.3612831138E4*phi-0.3612831138E4*old+
0.1204277045E4;
    t282 = t191*t88;
    t285 = t88*t198;
    t289 = 0.349964241E4*phi;
    t290 = 0.349964241E4*old;
    t291 = 0.116654747E4+0.699928482E4*t6-t289-t290;
    t293 = t274*t43-3.0*t218*t88+6.0*t223*t191-3.0*t107*t198-6.0*t233*t282+6.0*
t113*t285-t63*t291;
    t294 = t293*t45;
    t296 = t200*t90;
    t297 = t296*t45;
    t299 = t203*t90;
    t300 = t299*t45;
    t302 = t294*t47;
    t304 = t296*t48;
    t306 = t300*t47;
    t308 = 0.882498173E1*t294+0.1133084783E3*t297+0.6671400383E2*t300+
0.289445110522E2*t302+0.8683353315E2*t304+0.289445110522E2*t306;
    t368 = -3.0*t107*t246-6.0*t59*t232*t282+24.0*t38/t231/t42*t282*t67-18.0*
t233*t191*t120+6.0*t228*t285-18.0*t233*t285*t67+6.0*t113*t120*t198+6.0*t113*t88
*t246-t109*t291+2.0*t113*t291*t67-t63*(-0.349964241E4+t289+t290);
    t370 = ((-0.3612831138E4+0.3612831135E4*phi+0.3612831135E4*old)*t43-t274*
t62*t67-3.0*t216*t62*t88+6.0*t187*t112*t114-3.0*t218*t120+6.0*t105*t112*t191
-18.0*t83*t232*t234+12.0*t223*t237-3.0*t220*t198+6.0*t223*t241+t368)*t45;
    t372 = t293*t69;
    t373 = t372*t45;
    t375 = t248*t90;
    t376 = t375*t45;
    t378 = t200*t122;
    t379 = t378*t45;
    t381 = t296*t70;
    t383 = t204*t122;
    t385 = t299*t69;
    t386 = t385*t45;
    t388 = t370*t47;
    t390 = t372*t48;
    t392 = t375*t48;
    t394 = t378*t48;
    t396 = t296*t72;
    t398 = t204*t263;
    t400 = t385*t48;
    t402 = 0.882498173E1*t370+0.3776949278E2*t373+0.1133084783E3*t376+
0.1133084783E3*t379+0.2001420114E3*t381+0.2001420115E3*t383+0.9565851488E2*t386
+0.289445110522E2*t388+0.289445110522E2*t390+0.8683353315E2*t392+0.8683353315E2
*t394+0.8683353315E2*t396+0.8683353315E2*t398+0.289445110522E2*t400;
    t405 = -t168;
    t410 = -t180;
    t417 = 0.2889626582E2*t201+0.3839831904E2*t204+0.950205321859E1*t206+
0.950205321859E1*t208;
    t437 = 0.2889626582E2*t294+0.1151949571E3*t297+0.4790037226E2*t300+
0.950205321859E1*t302+0.2850615966E2*t304+0.950205321859E1*t306;
    t453 = 0.2889626582E2*t370+0.3839831904E2*t373+0.1151949571E3*t376+
0.1151949571E3*t379+0.1437011168E3*t381+0.1437011168E3*t383+0.5740242548E2*t386
+0.950205321859E1*t388+0.950205321859E1*t390+0.2850615966E2*t392+0.2850615966E2
*t394+0.2850615966E2*t396+0.2850615966E2*t398+0.950205321859E1*t400;
    t455 = t160*t50+t168*t74+3.0*t173*t95+3.0*t180*t132+3.0*t21*t210+3.0*t55*(
0.882498173E1*t249+0.3776949278E2*t252+0.7553898556E2*t254+0.6671400383E2*t257+
0.289445110522E2*t259+0.289445110522E2*t261+0.578890221E2*t264+0.289445110522E2
*t266)+t79*t308+t100*t402+t162*t136+t405*t141-3.0*t173*t146+3.0*t410*t153+3.0*
t28*t417+3.0*t138*(0.2889626582E2*t249+0.3839831904E2*t252+0.7679663808E2*t254+
0.4790037226E2*t257+0.950205321859E1*t259+0.950205321859E1*t261+0.1900410644E2*
t264+0.950205321859E1*t266)+t143*t437+t148*t453;
    return(2.0*t2*(0.6E1*t6+0.3E1*(2.0*alpha-2.0)*t5-0.3E1*alpha+0.1E1)+t21*t26
+t28*t31-beta*(t21*t50+t55*t74+t79*t95+t100*t132+t28*t136+t138*t141+t143*t146+
t148*t153)+(0.24E2*t2+t160*t26+t162*t31-beta*t455)*(phi-old)/24.0+t2*(0.12E2*
phi+0.12E2*old+12.0*alpha-12.0)/12.0+t168*t26/24.0+t405*t31/24.0-beta*(t168*t50
+3.0*t180*t95+3.0*t55*t210+t100*t308+t405*t136+3.0*t410*t146+3.0*t138*t417+t148
*t437)/24.0);
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t11;
  double t12;
  double t16;
  double t2;
  double t20;
  double t21;
  double t4;
  double t5;
  double t7;
  double t8;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t8 = t4-t5+t7;
    t11 = 0.5*rho+0.5*old;
    t12 = log(t11);
    t16 = 1.0-t4+t5-t7;
    t20 = t11*t11;
    t21 = 1/t20;
    return(t8*(0.882498173E1+0.289445110522E2*t12)+t16*(0.2889626582E2+
0.950205321859E1*t12)+(-0.2894451105E2*t8*t21-0.9502053219E1*t16*t21)*(rho-old)
/24.0);
  }
}

double drhochemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t11;
  double t12;
  double t15;
  double t18;
  double t2;
  double t20;
  double t29;
  double t4;
  double t5;
  double t7;
  double t8;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t8 = t4-t5+t7;
    t11 = 0.5*rho+0.5*old;
    t12 = 1/t11;
    t15 = 1.0-t4+t5-t7;
    t18 = t11*t11;
    t20 = 1/t18/t11;
    t29 = 1/t18;
    return(0.1447225552E2*t8*t12+0.475102661E1*t15*t12+(0.2894451105E2*t8*t20+
0.9502053219E1*t15*t20)*(rho-old)/24.0-0.1206021294E1*t8*t29-0.3959188841*t15*
t29);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t10;
  double t11;
  double t15;
  double t19;
  double t2;
  double t20;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t7 = 30.0*t2-60.0*t1*phi+30.0*t1;
    t10 = 0.5*rho+0.5*old;
    t11 = log(t10);
    t15 = -t7;
    t19 = t10*t10;
    t20 = 1/t19;
    return(t7*(0.882498173E1+0.289445110522E2*t11)+t15*(0.2889626582E2+
0.950205321859E1*t11)+(-0.2894451105E2*t7*t20-0.9502053219E1*t15*t20)*(rho-old)
/24.0);
  }
}

#include <math.h>
double pressure(double rho,double phi)
{
  double t1;
  double t10;
  double t11;
  double t2;
  double t4;
  double t5;
  double t6;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t6 = t1*phi;
    t7 = 10.0*t6;
    t10 = log(rho);
    t11 = rho*t10;
    return((t4-t5+t7)*(-0.2094881058E2+0.4291934897E2*rho-0.2894451105E2*t11+
rho*(-0.1397483792E2+0.289445110522E2*t10))+(1.0-t4+t5-t7)*(-0.6490446091E-1+
0.3405607049E1*rho-0.9502053219E1*t11+rho*(0.609644617E1+0.950205321859E1*t10))
-2.0/delta*A*(t2+(2.0*alpha-2.0)*t6+(-3.0*alpha+1.0)*t1+alpha));
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t1;
  double t2;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t9 = sqrt(0.116654747E11*t2*phi-0.2916368675E11*t2+0.1944245783E11*t1*phi+
950205322.0);
    return(0.1E-3*t9);
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t15;
  double t16;
  double t17;
  double t18;
  double t19;
  double t21;
  double t22;
  double t26;
  double t3;
  double t4;
  double t41;
  double t43;
  double t44;
  double t7;
  {
    t3 = phi*phi;
    t4 = t3*t3;
    t7 = t3*phi;
    t15 = t4*phi;
    t16 = 6.0*t15;
    t17 = 15.0*t4;
    t18 = 10.0*t7;
    t19 = t16-t17+t18;
    t21 = log(rho);
    t22 = rho*t21;
    t26 = 1.0-t16+t17-t18;
    t41 = exp((-0.2889626582E2+0.1204277045E3*t15-0.3010692614E3*t4+
0.2007128409E3*t7)/(0.9502053219E1+0.116654747E3*t15-0.2916368674E3*t4+
0.1944245783E3*t7));
    t43 = log(t41);
    t44 = t41*t43;
    return(2.0/delta*A*(t4+(2.0*alpha-2.0)*t7+(-3.0*alpha+1.0)*t3+alpha)+t19*(
-0.2011952932E2*rho+0.289445110522E2*t22+0.2133795654E2)+t26*(0.193942126E2*rho
+0.950205321859E1*t22+0.4540504209)-beta*(t19*(0.2133795654E2-0.2011952932E2*
t41+0.289445110522E2*t44)+t26*(0.4540504209+0.193942126E2*t41+0.950205321859E1*
t44)));
  }
}

#include <math.h>
double reactionSource(double rho,double phi,double old)
{
  double t100;
  double t107;
  double t114;
  double t116;
  double t120;
  double t121;
  double t122;
  double t129;
  double t131;
  double t132;
  double t134;
  double t135;
  double t137;
  double t139;
  double t157;
  double t172;
  double t174;
  double t175;
  double t178;
  double t180;
  double t182;
  double t184;
  double t2;
  double t20;
  double t24;
  double t26;
  double t27;
  double t29;
  double t31;
  double t34;
  double t36;
  double t40;
  double t44;
  double t45;
  double t47;
  double t49;
  double t5;
  double t50;
  double t52;
  double t54;
  double t55;
  double t56;
  double t57;
  double t6;
  double t61;
  double t63;
  double t64;
  double t65;
  double t69;
  double t7;
  double t71;
  double t72;
  double t74;
  double t76;
  double t80;
  double t82;
  double t85;
  double t98;
  {
    t2 = 1/delta*A;
    t5 = 0.5*phi+0.5*old;
    t6 = t5*t5;
    t7 = t6*t5;
    t20 = t6*t6;
    t24 = 30.0*t20-60.0*t7+30.0*t6;
    t26 = log(rho);
    t27 = rho*t26;
    t29 = -0.2011952932E2*rho+0.289445110522E2*t27+0.2133795654E2;
    t31 = -t24;
    t34 = 0.193942126E2*rho+0.950205321859E1*t27+0.4540504209;
    t36 = t20*t5;
    t40 = -0.2889626582E2+0.1204277045E3*t36-0.3010692614E3*t20+0.2007128409E3*
t7;
    t44 = 0.9502053219E1+0.116654747E3*t36-0.2916368674E3*t20+0.1944245783E3*t7
;
    t45 = 1/t44;
    t47 = exp(t40*t45);
    t49 = log(t47);
    t50 = t47*t49;
    t52 = 0.2133795654E2-0.2011952932E2*t47+0.289445110522E2*t50;
    t54 = 6.0*t36;
    t55 = 15.0*t20;
    t56 = 10.0*t7;
    t57 = t54-t55+t56;
    t61 = 0.6021385225E3*t20-0.1204277046E4*t7+0.6021385227E3*t6;
    t63 = t44*t44;
    t64 = 1/t63;
    t65 = t40*t64;
    t69 = 0.583273735E3*t20-0.116654747E4*t7+0.5832737349E3*t6;
    t71 = t61*t45-t65*t69;
    t72 = t71*t47;
    t74 = t72*t49;
    t76 = 0.882498173E1*t72+0.289445110522E2*t74;
    t80 = 0.4540504209+0.193942126E2*t47+0.950205321859E1*t50;
    t82 = 1.0-t54+t55-t56;
    t85 = 0.2889626582E2*t72+0.950205321859E1*t74;
    t98 = 360.0*t6-0.18E3*phi-0.18E3*old+60.0;
    t100 = -t98;
    t107 = 120.0*t7-180.0*t6+0.3E2*phi+0.3E2*old;
    t114 = 0.240855409E4*t7-0.3612831138E4*t6+0.6021385225E3*phi+0.6021385225E3
*old;
    t116 = t61*t64;
    t120 = 1/t63/t44;
    t121 = t40*t120;
    t122 = t69*t69;
    t129 = 0.233309494E4*t7-0.349964241E4*t6+0.583273735E3*phi+0.583273735E3*
old;
    t131 = t114*t45-2.0*t116*t69+2.0*t121*t122-t65*t129;
    t132 = t131*t47;
    t134 = t71*t71;
    t135 = t134*t47;
    t137 = t132*t49;
    t139 = t135*t49;
    t157 = t63*t63;
    t172 = ((0.722566227E4*t6-0.3612831138E4*phi-0.3612831138E4*old+
0.1204277045E4)*t45-3.0*t114*t64*t69+6.0*t61*t120*t122-3.0*t116*t129-6.0*t40/
t157*t122*t69+6.0*t121*t69*t129-t65*(0.699928482E4*t6-0.349964241E4*phi
-0.349964241E4*old+0.116654747E4))*t47;
    t174 = t131*t71;
    t175 = t174*t47;
    t178 = t134*t71*t47;
    t180 = t172*t49;
    t182 = t174*t50;
    t184 = t178*t49;
    return(2.0*t2*(4.0*t7+3.0*(2.0*alpha-2.0)*t6+2.0*(-3.0*alpha+1.0)*t5)+t24*
t29+t31*t34-beta*(t24*t52+t57*t76+t31*t80+t82*t85)+(2.0*t2*(0.12E2*phi+0.12E2*
old+12.0*alpha-12.0)+t98*t29+t100*t34-beta*(t98*t52+3.0*t107*t76+3.0*t24*(
0.882498173E1*t132+0.3776949278E2*t135+0.289445110522E2*t137+0.289445110522E2*
t139)+t57*(0.882498173E1*t172+0.1133084783E3*t175+0.6671400383E2*t178+
0.289445110522E2*t180+0.8683353315E2*t182+0.289445110522E2*t184)+t100*t80-3.0*
t107*t85+3.0*t31*(0.2889626582E2*t132+0.3839831904E2*t135+0.950205321859E1*t137
+0.950205321859E1*t139)+t82*(0.2889626582E2*t172+0.1151949571E3*t175+
0.4790037226E2*t178+0.950205321859E1*t180+0.2850615966E2*t182+0.950205321859E1*
t184)))*(phi-old)/24.0);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi,double old,double mid)
{
  double t100;
  double t105;
  double t107;
  double t109;
  double t112;
  double t113;
  double t114;
  double t120;
  double t122;
  double t123;
  double t125;
  double t126;
  double t128;
  double t130;
  double t132;
  double t136;
  double t138;
  double t141;
  double t143;
  double t146;
  double t148;
  double t153;
  double t16;
  double t160;
  double t162;
  double t168;
  double t173;
  double t180;
  double t187;
  double t191;
  double t195;
  double t198;
  double t2;
  double t200;
  double t201;
  double t203;
  double t204;
  double t206;
  double t208;
  double t21;
  double t210;
  double t216;
  double t218;
  double t220;
  double t223;
  double t228;
  double t23;
  double t231;
  double t232;
  double t233;
  double t234;
  double t237;
  double t24;
  double t241;
  double t246;
  double t248;
  double t249;
  double t251;
  double t252;
  double t254;
  double t256;
  double t257;
  double t259;
  double t26;
  double t261;
  double t263;
  double t264;
  double t266;
  double t274;
  double t28;
  double t282;
  double t285;
  double t289;
  double t290;
  double t291;
  double t293;
  double t294;
  double t296;
  double t297;
  double t299;
  double t300;
  double t302;
  double t304;
  double t306;
  double t308;
  double t31;
  double t33;
  double t34;
  double t368;
  double t370;
  double t372;
  double t373;
  double t375;
  double t376;
  double t378;
  double t379;
  double t38;
  double t381;
  double t383;
  double t385;
  double t386;
  double t388;
  double t390;
  double t392;
  double t394;
  double t396;
  double t398;
  double t400;
  double t402;
  double t405;
  double t410;
  double t417;
  double t42;
  double t43;
  double t437;
  double t45;
  double t453;
  double t455;
  double t47;
  double t48;
  double t5;
  double t50;
  double t55;
  double t59;
  double t6;
  double t61;
  double t62;
  double t63;
  double t67;
  double t69;
  double t70;
  double t72;
  double t74;
  double t79;
  double t83;
  double t86;
  double t88;
  double t90;
  double t91;
  double t93;
  double t95;
  double t97;
  double t98;
  double t99;
  {
    t2 = 1/delta*A;
    t5 = 0.5*phi+0.5*old;
    t6 = t5*t5;
    t16 = t6*t5;
    t21 = 0.6E2*t16-0.9E2*t6+0.15E2*phi+0.15E2*old;
    t23 = log(rho);
    t24 = rho*t23;
    t26 = -0.2011952932E2*rho+0.289445110522E2*t24+0.2133795654E2;
    t28 = -t21;
    t31 = 0.193942126E2*rho+0.950205321859E1*t24+0.4540504209;
    t33 = t6*t6;
    t34 = t33*t5;
    t38 = -0.2889626582E2+0.1204277045E3*t34-0.3010692614E3*t33+0.2007128409E3*
t16;
    t42 = 0.9502053219E1+0.116654747E3*t34-0.2916368674E3*t33+0.1944245783E3*
t16;
    t43 = 1/t42;
    t45 = exp(t38*t43);
    t47 = log(t45);
    t48 = t45*t47;
    t50 = 0.2133795654E2-0.2011952932E2*t45+0.289445110522E2*t48;
    t55 = 30.0*t33-60.0*t16+30.0*t6;
    t59 = 0.3010692612E3*t33-0.6021385228E3*t16+0.3010692614E3*t6;
    t61 = t42*t42;
    t62 = 1/t61;
    t63 = t38*t62;
    t67 = 0.2916368675E3*t33-0.5832737348E3*t16+0.2916368674E3*t6;
    t69 = t59*t43-t63*t67;
    t70 = t69*t45;
    t72 = t70*t47;
    t74 = 0.882498173E1*t70+0.289445110522E2*t72;
    t79 = 0.15E2*t33-0.3E2*t16+0.15E2*t6;
    t83 = 0.6021385225E3*t33-0.1204277046E4*t16+0.6021385227E3*t6;
    t86 = 0.116654747E4*t16;
    t88 = 0.583273735E3*t33-t86+0.5832737349E3*t6;
    t90 = t83*t43-t63*t88;
    t91 = t90*t45;
    t93 = t91*t47;
    t95 = 0.882498173E1*t91+0.289445110522E2*t93;
    t97 = 6.0*t34;
    t98 = 15.0*t33;
    t99 = 10.0*t16;
    t100 = t97-t98+t99;
    t105 = 0.1204277045E4*t16-0.1806415569E4*t6+0.3010692614E3*phi+
0.3010692614E3*old;
    t107 = t83*t62;
    t109 = t59*t62;
    t112 = 1/t61/t42;
    t113 = t38*t112;
    t114 = t88*t67;
    t120 = t86-0.1749821205E4*t6+0.2916368674E3*phi+0.2916368674E3*old;
    t122 = t105*t43-t107*t67-t109*t88+2.0*t113*t114-t63*t120;
    t123 = t122*t45;
    t125 = t90*t69;
    t126 = t125*t45;
    t128 = t123*t47;
    t130 = t125*t48;
    t132 = 0.882498173E1*t123+0.3776949278E2*t126+0.289445110522E2*t128+
0.289445110522E2*t130;
    t136 = 0.4540504209+0.193942126E2*t45+0.950205321859E1*t48;
    t138 = -t55;
    t141 = 0.2889626582E2*t70+0.950205321859E1*t72;
    t143 = -t79;
    t146 = 0.2889626582E2*t91+0.950205321859E1*t93;
    t148 = 1.0-t97+t98-t99;
    t153 = 0.2889626582E2*t123+0.3839831904E2*t126+0.950205321859E1*t128+
0.950205321859E1*t130;
    t160 = 0.18E3*phi+0.18E3*old-0.18E3;
    t162 = -t160;
    t168 = 360.0*t6-0.18E3*phi-0.18E3*old+60.0;
    t173 = 0.18E3*t6-0.9E2*phi-0.9E2*old+0.3E2;
    t180 = 120.0*t16-180.0*t6+0.3E2*phi+0.3E2*old;
    t187 = 0.240855409E4*t16-0.3612831138E4*t6+0.6021385225E3*phi+
0.6021385225E3*old;
    t191 = t88*t88;
    t195 = 0.349964241E4*t6;
    t198 = 0.233309494E4*t16-t195+0.583273735E3*phi+0.583273735E3*old;
    t200 = t187*t43-2.0*t107*t88+2.0*t113*t191-t63*t198;
    t201 = t200*t45;
    t203 = t90*t90;
    t204 = t203*t45;
    t206 = t201*t47;
    t208 = t204*t47;
    t210 = 0.882498173E1*t201+0.3776949278E2*t204+0.289445110522E2*t206+
0.289445110522E2*t208;
    t216 = 0.6021385225E3+0.3612831135E4*t6-0.1806415569E4*phi-0.1806415569E4*
old;
    t218 = t187*t62;
    t220 = t105*t62;
    t223 = t83*t112;
    t228 = t59*t112;
    t231 = t61*t61;
    t232 = 1/t231;
    t233 = t38*t232;
    t234 = t191*t67;
    t237 = t88*t120;
    t241 = t198*t67;
    t246 = 0.583273735E3+t195-0.1749821205E4*phi-0.1749821205E4*old;
    t248 = t216*t43-t218*t67-2.0*t220*t88+4.0*t223*t114-2.0*t107*t120+2.0*t228*
t191-6.0*t233*t234+4.0*t113*t237-t109*t198+2.0*t113*t241-t63*t246;
    t249 = t248*t45;
    t251 = t200*t69;
    t252 = t251*t45;
    t254 = t91*t122;
    t256 = t203*t69;
    t257 = t256*t45;
    t259 = t249*t47;
    t261 = t251*t48;
    t263 = t47*t122;
    t264 = t91*t263;
    t266 = t256*t48;
    t274 = 0.722566227E4*t6-0.3612831138E4*phi-0.3612831138E4*old+
0.1204277045E4;
    t282 = t191*t88;
    t285 = t88*t198;
    t289 = 0.349964241E4*phi;
    t290 = 0.349964241E4*old;
    t291 = 0.116654747E4+0.699928482E4*t6-t289-t290;
    t293 = t274*t43-3.0*t218*t88+6.0*t223*t191-3.0*t107*t198-6.0*t233*t282+6.0*
t113*t285-t63*t291;
    t294 = t293*t45;
    t296 = t200*t90;
    t297 = t296*t45;
    t299 = t203*t90;
    t300 = t299*t45;
    t302 = t294*t47;
    t304 = t296*t48;
    t306 = t300*t47;
    t308 = 0.882498173E1*t294+0.1133084783E3*t297+0.6671400383E2*t300+
0.289445110522E2*t302+0.8683353315E2*t304+0.289445110522E2*t306;
    t368 = -3.0*t107*t246-6.0*t59*t232*t282+24.0*t38/t231/t42*t282*t67-18.0*
t233*t191*t120+6.0*t228*t285-18.0*t233*t285*t67+6.0*t113*t120*t198+6.0*t113*t88
*t246-t109*t291+2.0*t113*t291*t67-t63*(-0.349964241E4+t289+t290);
    t370 = ((-0.3612831138E4+0.3612831135E4*phi+0.3612831135E4*old)*t43-t274*
t62*t67-3.0*t216*t62*t88+6.0*t187*t112*t114-3.0*t218*t120+6.0*t105*t112*t191
-18.0*t83*t232*t234+12.0*t223*t237-3.0*t220*t198+6.0*t223*t241+t368)*t45;
    t372 = t293*t69;
    t373 = t372*t45;
    t375 = t248*t90;
    t376 = t375*t45;
    t378 = t200*t122;
    t379 = t378*t45;
    t381 = t296*t70;
    t383 = t204*t122;
    t385 = t299*t69;
    t386 = t385*t45;
    t388 = t370*t47;
    t390 = t372*t48;
    t392 = t375*t48;
    t394 = t378*t48;
    t396 = t296*t72;
    t398 = t204*t263;
    t400 = t385*t48;
    t402 = 0.882498173E1*t370+0.3776949278E2*t373+0.1133084783E3*t376+
0.1133084783E3*t379+0.2001420114E3*t381+0.2001420115E3*t383+0.9565851488E2*t386
+0.289445110522E2*t388+0.289445110522E2*t390+0.8683353315E2*t392+0.8683353315E2
*t394+0.8683353315E2*t396+0.8683353315E2*t398+0.289445110522E2*t400;
    t405 = -t168;
    t410 = -t180;
    t417 = 0.2889626582E2*t201+0.3839831904E2*t204+0.950205321859E1*t206+
0.950205321859E1*t208;
    t437 = 0.2889626582E2*t294+0.1151949571E3*t297+0.4790037226E2*t300+
0.950205321859E1*t302+0.2850615966E2*t304+0.950205321859E1*t306;
    t453 = 0.2889626582E2*t370+0.3839831904E2*t373+0.1151949571E3*t376+
0.1151949571E3*t379+0.1437011168E3*t381+0.1437011168E3*t383+0.5740242548E2*t386
+0.950205321859E1*t388+0.950205321859E1*t390+0.2850615966E2*t392+0.2850615966E2
*t394+0.2850615966E2*t396+0.2850615966E2*t398+0.950205321859E1*t400;
    t455 = t160*t50+t168*t74+3.0*t173*t95+3.0*t180*t132+3.0*t21*t210+3.0*t55*(
0.882498173E1*t249+0.3776949278E2*t252+0.7553898556E2*t254+0.6671400383E2*t257+
0.289445110522E2*t259+0.289445110522E2*t261+0.578890221E2*t264+0.289445110522E2
*t266)+t79*t308+t100*t402+t162*t136+t405*t141-3.0*t173*t146+3.0*t410*t153+3.0*
t28*t417+3.0*t138*(0.2889626582E2*t249+0.3839831904E2*t252+0.7679663808E2*t254+
0.4790037226E2*t257+0.950205321859E1*t259+0.950205321859E1*t261+0.1900410644E2*
t264+0.950205321859E1*t266)+t143*t437+t148*t453;
    return(2.0*t2*(0.6E1*t6+0.3E1*(2.0*alpha-2.0)*t5-0.3E1*alpha+0.1E1)+t21*t26
+t28*t31-beta*(t21*t50+t55*t74+t79*t95+t100*t132+t28*t136+t138*t141+t143*t146+
t148*t153)+(0.24E2*t2+t160*t26+t162*t31-beta*t455)*(phi-old)/24.0+t2*(0.12E2*
phi+0.12E2*old+12.0*alpha-12.0)/12.0+t168*t26/24.0+t405*t31/24.0-beta*(t168*t50
+3.0*t180*t95+3.0*t55*t210+t100*t308+t405*t136+3.0*t410*t146+3.0*t138*t417+t148
*t437)/24.0);
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t11;
  double t12;
  double t16;
  double t2;
  double t20;
  double t21;
  double t4;
  double t5;
  double t7;
  double t8;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t8 = t4-t5+t7;
    t11 = 0.5*rho+0.5*old;
    t12 = log(t11);
    t16 = 1.0-t4+t5-t7;
    t20 = t11*t11;
    t21 = 1/t20;
    return(t8*(0.882498173E1+0.289445110522E2*t12)+t16*(0.2889626582E2+
0.950205321859E1*t12)+(-0.2894451105E2*t8*t21-0.9502053219E1*t16*t21)*(rho-old)
/24.0);
  }
}

double drhochemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t11;
  double t12;
  double t15;
  double t18;
  double t2;
  double t20;
  double t29;
  double t4;
  double t5;
  double t7;
  double t8;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t8 = t4-t5+t7;
    t11 = 0.5*rho+0.5*old;
    t12 = 1/t11;
    t15 = 1.0-t4+t5-t7;
    t18 = t11*t11;
    t20 = 1/t18/t11;
    t29 = 1/t18;
    return(0.1447225552E2*t8*t12+0.475102661E1*t15*t12+(0.2894451105E2*t8*t20+
0.9502053219E1*t15*t20)*(rho-old)/24.0-0.1206021294E1*t8*t29-0.3959188841*t15*
t29);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi,double old)
{
  double t1;
  double t10;
  double t11;
  double t15;
  double t19;
  double t2;
  double t20;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t7 = 30.0*t2-60.0*t1*phi+30.0*t1;
    t10 = 0.5*rho+0.5*old;
    t11 = log(t10);
    t15 = -t7;
    t19 = t10*t10;
    t20 = 1/t19;
    return(t7*(0.882498173E1+0.289445110522E2*t11)+t15*(0.2889626582E2+
0.950205321859E1*t11)+(-0.2894451105E2*t7*t20-0.9502053219E1*t15*t20)*(rho-old)
/24.0);
  }
}

#include <math.h>
double pressure(double rho,double phi)
{
  double t1;
  double t10;
  double t11;
  double t2;
  double t4;
  double t5;
  double t6;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t6 = t1*phi;
    t7 = 10.0*t6;
    t10 = log(rho);
    t11 = rho*t10;
    return((t4-t5+t7)*(-0.2094881058E2+0.4291934897E2*rho-0.2894451105E2*t11+
rho*(-0.1397483792E2+0.289445110522E2*t10))+(1.0-t4+t5-t7)*(-0.6490446091E-1+
0.3405607049E1*rho-0.9502053219E1*t11+rho*(0.609644617E1+0.950205321859E1*t10))
-2.0/delta*A*(t2+(2.0*alpha-2.0)*t6+(-3.0*alpha+1.0)*t1+alpha));
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t1;
  double t2;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t9 = sqrt(0.116654747E11*t2*phi-0.2916368675E11*t2+0.1944245783E11*t1*phi+
950205322.0);
    return(0.1E-3*t9);
  }
}

