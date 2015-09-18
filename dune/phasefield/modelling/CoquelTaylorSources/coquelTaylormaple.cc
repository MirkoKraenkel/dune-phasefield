
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
  double t108;
  double t110;
  double t114;
  double t115;
  double t116;
  double t12;
  double t123;
  double t125;
  double t129;
  double t130;
  double t146;
  double t15;
  double t160;
  double t163;
  double t166;
  double t167;
  double t175;
  double t177;
  double t186;
  double t187;
  double t188;
  double t189;
  double t19;
  double t191;
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
  double t97;
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
    t97 = -120.0*t5+180.0*t4-0.3E2*phi-0.3E2*old;
    t108 = 0.8317766168E2*t5-0.1247664925E3*t4+0.2079441542E2*phi+
0.2079441542E2*old;
    t110 = t52*t55;
    t114 = 1/t54/t37;
    t115 = t33*t114;
    t116 = t60*t60;
    t123 = -0.6E2*t5+0.9E2*t4-0.15E2*phi-0.15E2*old;
    t125 = t108*t38-2.0*t110*t60+2.0*t115*t116-t56*t123;
    t129 = t62*t62;
    t130 = t19*t129;
    t146 = t54*t54;
    t160 = (0.249532985E3*t4-0.1247664925E3*phi-0.1247664925E3*old+
0.4158883084E2)*t38-3.0*t108*t55*t60+6.0*t52*t114*t116-3.0*t110*t123-6.0*t33/
t146*t116*t60+6.0*t115*t60*t123-t56*(-0.18E3*t4+0.9E2*phi+0.9E2*old-0.3E2);
    t163 = t48*t125;
    t166 = t129*t62;
    t167 = t48*t166;
    t175 = t125*t40;
    t177 = t129*t40;
    t186 = t160*t40;
    t187 = 0.6931471806*t186;
    t188 = t125*t62;
    t189 = t188*t40;
    t191 = t166*t40;
    t203 = t87*t43+t89*t67+3.0*t97*t74-3.0*t97*t62*t42+3.0*t19*t125*t42+3.0*
t130*t40+t48*t160*t42+3.0*t163*t70+2.0*t167*t40+3.0*t130*t42+3.0*t163*t72+t167*
t42+3.0*t24*(0.6931471806*t175+0.2193147181E1*t177+0.15E1*t175*t41+0.15E1*t177*
t41)+t69*(t187+0.6579441543E1*t189+0.3693147181E1*t191+0.15E1*t186*t41+0.45E1*
t188*t42+0.15E1*t191*t41)-t187-0.2079441542E1*t189-0.6931471806*t191;
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
  double t149;
  double t150;
  double t151;
  double t153;
  double t154;
  double t155;
  double t157;
  double t161;
  double t162;
  double t163;
  double t166;
  double t169;
  double t170;
  double t173;
  double t174;
  double t179;
  double t18;
  double t180;
  double t181;
  double t182;
  double t185;
  double t19;
  double t191;
  double t193;
  double t194;
  double t195;
  double t20;
  double t203;
  double t204;
  double t205;
  double t21;
  double t211;
  double t213;
  double t220;
  double t224;
  double t23;
  double t231;
  double t235;
  double t238;
  double t240;
  double t249;
  double t26;
  double t262;
  double t265;
  double t268;
  double t271;
  double t278;
  double t28;
  double t29;
  double t293;
  double t3;
  double t31;
  double t313;
  double t314;
  double t315;
  double t316;
  double t317;
  double t318;
  double t321;
  double t324;
  double t33;
  double t350;
  double t351;
  double t352;
  double t354;
  double t359;
  double t362;
  double t367;
  double t37;
  double t370;
  double t376;
  double t38;
  double t390;
  double t393;
  double t397;
  double t399;
  double t4;
  double t40;
  double t404;
  double t41;
  double t411;
  double t42;
  double t425;
  double t43;
  double t450;
  double t454;
  double t459;
  double t460;
  double t461;
  double t463;
  double t471;
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
    t149 = t75*t75;
    t150 = t149*t40;
    t151 = t150*t103;
    t153 = t149*t75;
    t154 = t153*t61;
    t155 = t154*t40;
    t157 = t67*t153;
    t161 = 0.1247664925E3*phi;
    t162 = 0.1247664925E3*old;
    t163 = 0.4158883084E2+0.249532985E3*t4-t161-t162;
    t166 = 0.1247664925E3*t4;
    t169 = 0.8317766168E2*t13-t166+0.2079441542E2*phi+0.2079441542E2*old;
    t170 = t169*t54;
    t173 = t71*t92;
    t174 = t73*t73;
    t179 = t53*t53;
    t180 = 1/t179;
    t181 = t33*t180;
    t182 = t174*t73;
    t185 = t73*t23;
    t191 = -0.18E3*t4+0.9E2*phi+0.9E2*old-0.3E2;
    t193 = t163*t38-3.0*t170*t73+6.0*t173*t174-3.0*t87*t23-6.0*t181*t182+6.0*
t93*t185-t55*t191;
    t194 = t193*t61;
    t195 = t194*t40;
    t203 = t169*t38-2.0*t87*t73+2.0*t93*t174-t55*t23;
    t204 = t203*t103;
    t205 = t204*t40;
    t211 = t81*t193;
    t213 = t81*t149;
    t220 = 120.0*t13-180.0*t4+0.3E2*phi+0.3E2*old;
    t224 = t203*t40;
    t231 = 0.6931471806*t224+0.2193147181E1*t150+0.15E1*t224*t41+0.15E1*t150*
t41;
    t235 = t220*t75;
    t238 = -0.2079441542E1*t151-0.6931471806*t155+2.0*t157*t40-0.6931471806*
t195-0.2079441542E1*t205-3.0*t191*t75*t42+t211*t108+6.0*t213*t127+3.0*t220*t103
*t42+3.0*t23*t231+t145*t43+3.0*t235*t107;
    t240 = t18*t149;
    t249 = 0.2079441542E2+t166-0.6238324625E2*phi-0.6238324625E2*old;
    t262 = t174*t59;
    t265 = t73*t101;
    t268 = t85*t54;
    t271 = t23*t59;
    t278 = -0.15E2-0.9E2*t4+0.45E2*phi+0.45E2*old;
    t293 = t51*t92;
    t313 = -3.0*t87*t278-6.0*t51*t180*t182+24.0*t33/t179/t37*t182*t59-18.0*t181
*t174*t101+6.0*t293*t185-18.0*t181*t185*t59+6.0*t93*t101*t23+6.0*t93*t73*t278-
t89*t191+2.0*t93*t191*t59-t55*(0.9E2-0.9E2*phi-0.9E2*old);
    t314 = (-0.1247664925E3+t161+t162)*t38-t163*t54*t59-3.0*t249*t54*t73+6.0*
t169*t92*t94-3.0*t170*t101+6.0*t85*t92*t174-18.0*t71*t180*t262+12.0*t173*t265
-3.0*t268*t23+6.0*t173*t271+t313;
    t315 = t314*t40;
    t316 = 0.6931471806*t315;
    t317 = t203*t75;
    t318 = t317*t107;
    t321 = t67*t203;
    t324 = t81*t203;
    t350 = t249*t38-t170*t59-2.0*t268*t73+4.0*t173*t94-2.0*t87*t101+2.0*t293*
t174-6.0*t181*t262+4.0*t93*t265-t89*t23+2.0*t93*t271-t55*t278;
    t351 = t350*t75;
    t352 = t351*t40;
    t354 = t81*t153;
    t359 = t211*t107+3.0*t240*t42-t316-0.2079441542E1*t318+t147*t113+3.0*t321*
t122+3.0*t324*t75*t108+3.0*t191*t124+3.0*t321*t120-0.2079441542E1*t352+3.0*t354
*t107+3.0*t324*t132;
    t362 = t48*t149;
    t367 = t48*t203;
    t370 = t81*t350;
    t376 = 360.0*t4-0.18E3*phi-0.18E3*old+60.0;
    t390 = t354*t108+3.0*t362*t108+3.0*t240*t40+3.0*t367*t107+3.0*t370*t120+
t376*t61*t42+t157*t42+3.0*t324*t127+3.0*t235*t108+t67*t193*t42+t81*t314*t42+6.0
*t362*t107;
    t393 = -t376;
    t397 = t350*t40;
    t399 = t203*t61;
    t404 = t149*t61;
    t411 = t41*t103;
    t425 = t48*t75;
    t450 = t316+0.2193147181E1*t195+0.6579441543E1*t352+0.6579441543E1*t205+
0.1107944154E2*t318+0.1107944154E2*t151+0.5193147181E1*t155+0.15E1*t315*t41+
0.15E1*t194*t42+0.45E1*t351*t42+0.45E1*t204*t42+0.45E1*t317*t108+0.45E1*t150*
t411+0.15E1*t154*t42;
    t454 = -t220;
    t459 = t193*t40;
    t460 = 0.6931471806*t459;
    t461 = t317*t40;
    t463 = t153*t40;
    t471 = t460+0.6579441543E1*t461+0.3693147181E1*t463+0.15E1*t459*t41+0.45E1*
t317*t42+0.15E1*t463*t41;
    t473 = 3.0*t367*t108+t393*t118+6.0*t324*t130+3.0*t115*(0.6931471806*t397+
0.2193147181E1*t399*t40+0.4386294362E1*t120*t103+0.3693147181E1*t404*t40+0.15E1
*t397*t41+0.15E1*t399*t42+0.3E1*t120*t411+0.15E1*t404*t42)+3.0*t18*t203*t42+3.0
*t48*t350*t42+6.0*t425*t127+3.0*t370*t122+t126*t450+3.0*t213*t132+3.0*t454*t136
+6.0*t425*t132+t73*t471;
    t516 = t376*t43+t393*t113+3.0*t454*t124+3.0*t235*t42+3.0*t367*t42+3.0*t362*
t40+t211*t42+3.0*t324*t120+2.0*t354*t40+3.0*t362*t42+3.0*t324*t122+t354*t42+3.0
*t115*t231+t126*t471-t460-0.2079441542E1*t461-0.6931471806*t463;
    return(2.0*A_*(0.6E1*t4-0.3E1*phi-0.3E1*old+0.1E1)*t10+t18*t21+t23*t26-beta_*
t139+(0.24E2*A_*t10+t145*t21+t147*t26-beta_*(t238+t359+t390+t473))*(phi-old)/24.0
+A_*(0.12E2*phi+0.12E2*old-12.0)*t10/12.0+t376*t21/24.0+t393*t26/24.0-beta_*t516/
24.0);
  }
}


inline double drhoreactionSource ( double rho ,double phi ,double old ) const
{
  double t10;
  double t11;
  double t15;
  double t20;
  double t3;
  double t4;
  double t5;
  {
    t3 = 0.5*phi+0.5*old;
    t4 = t3*t3;
    t5 = t4*t4;
    t10 = 30.0*t5-60.0*t4*t3+30.0*t4;
    t11 = log(rho);
    t15 = 0.6931471806+0.15E1*t11;
    t20 = 360.0*t4-0.18E3*phi-0.18E3*old+60.0;
    return(t10*t11-t10*t15+(t20*t11-t20*t15)*(phi-old)/24.0);
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
    return(0.5*t8*t12+0.75*t15*t12+(0.1E1*t8*t20+0.15E1*t20*t15)*(rho-old)/24.0
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

