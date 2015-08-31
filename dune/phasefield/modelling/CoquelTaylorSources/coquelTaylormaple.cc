
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
  double t103;
  double t105;
  double t109;
  double t110;
  double t111;
  double t118;
  double t12;
  double t120;
  double t129;
  double t131;
  double t132;
  double t15;
  double t154;
  double t168;
  double t169;
  double t170;
  double t171;
  double t172;
  double t174;
  double t175;
  double t185;
  double t19;
  double t190;
  double t193;
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
  double t95;
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
    t95 = 120.0*t5-180.0*t4+0.3E2*phi+0.3E2*old;
    t103 = 0.8317766168E2*t5-0.1247664925E3*t4+0.2079441542E2*phi+
0.2079441542E2*old;
    t105 = t52*t55;
    t109 = 1/t54/t37;
    t110 = t33*t109;
    t111 = t60*t60;
    t118 = -0.6E2*t5+0.9E2*t4-0.15E2*phi-0.15E2*old;
    t120 = t103*t38-2.0*t105*t60+2.0*t110*t111-t56*t118;
    t129 = t120*t40;
    t131 = t62*t62;
    t132 = t131*t40;
    t154 = t54*t54;
    t168 = (0.249532985E3*t4-0.1247664925E3*phi-0.1247664925E3*old+
0.4158883084E2)*t38-3.0*t103*t55*t60+6.0*t52*t109*t111-3.0*t105*t118-6.0*t33/
t154*t111*t60+6.0*t110*t60*t118-t56*(-0.18E3*t4+0.9E2*phi+0.9E2*old-0.3E2);
    t169 = t168*t40;
    t170 = 0.6931471806*t169;
    t171 = t120*t62;
    t172 = t171*t40;
    t174 = t131*t62;
    t175 = t174*t40;
    t185 = t19*t131;
    t190 = t48*t120;
    t193 = t48*t174;
    t203 = 3.0*t95*t62*t42+3.0*t19*t120*t42+t87*t43+t89*t67-3.0*t95*t74+3.0*t24
*(0.6931471806*t129+0.2193147181E1*t132+0.15E1*t129*t41+0.15E1*t132*t41)+t69*(
t170+0.6579441543E1*t172+0.3693147181E1*t175+0.15E1*t169*t41+0.45E1*t171*t42+
0.15E1*t175*t41)+3.0*t185*t40+t48*t168*t42+3.0*t190*t70+2.0*t193*t40+3.0*t185*
t42+3.0*t190*t72-0.2079441542E1*t172+t193*t42-t170-0.6931471806*t175;
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
  double t152;
  double t157;
  double t158;
  double t160;
  double t163;
  double t168;
  double t169;
  double t172;
  double t173;
  double t174;
  double t175;
  double t178;
  double t18;
  double t182;
  double t188;
  double t19;
  double t190;
  double t194;
  double t195;
  double t20;
  double t200;
  double t202;
  double t21;
  double t211;
  double t212;
  double t217;
  double t218;
  double t220;
  double t224;
  double t227;
  double t228;
  double t23;
  double t231;
  double t232;
  double t235;
  double t236;
  double t237;
  double t245;
  double t248;
  double t254;
  double t256;
  double t259;
  double t26;
  double t260;
  double t28;
  double t29;
  double t3;
  double t31;
  double t321;
  double t322;
  double t325;
  double t327;
  double t33;
  double t333;
  double t336;
  double t342;
  double t343;
  double t347;
  double t348;
  double t350;
  double t351;
  double t356;
  double t363;
  double t366;
  double t37;
  double t370;
  double t373;
  double t38;
  double t380;
  double t381;
  double t387;
  double t391;
  double t392;
  double t4;
  double t40;
  double t409;
  double t41;
  double t414;
  double t416;
  double t419;
  double t42;
  double t421;
  double t426;
  double t43;
  double t451;
  double t452;
  double t453;
  double t455;
  double t463;
  double t465;
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
    t149 = 0.1247664925E3*t4;
    t152 = 0.2079441542E2+t149-0.6238324625E2*phi-0.6238324625E2*old;
    t157 = 0.8317766168E2*t13-t149+0.2079441542E2*phi+0.2079441542E2*old;
    t158 = t157*t54;
    t160 = t85*t54;
    t163 = t71*t92;
    t168 = t51*t92;
    t169 = t73*t73;
    t172 = t53*t53;
    t173 = 1/t172;
    t174 = t33*t173;
    t175 = t169*t59;
    t178 = t73*t101;
    t182 = t23*t59;
    t188 = -0.15E2-0.9E2*t4+0.45E2*phi+0.45E2*old;
    t190 = t152*t38-t158*t59-2.0*t160*t73+4.0*t163*t94-2.0*t87*t101+2.0*t168*
t169-6.0*t174*t175+4.0*t93*t178-t89*t23+2.0*t93*t182-t55*t188;
    t194 = t190*t75;
    t195 = t194*t40;
    t200 = -60.0-360.0*t4+0.18E3*phi+0.18E3*old;
    t202 = t48*t75;
    t211 = t157*t38-2.0*t87*t73+2.0*t93*t169-t55*t23;
    t212 = t48*t211;
    t217 = t211*t75;
    t218 = t217*t107;
    t220 = t81*t211;
    t224 = t81*t190;
    t227 = t75*t75;
    t228 = t18*t227;
    t231 = t227*t75;
    t232 = t67*t231;
    t235 = 0.1247664925E3*phi;
    t236 = 0.1247664925E3*old;
    t237 = 0.4158883084E2+0.249532985E3*t4-t235-t236;
    t245 = t169*t73;
    t248 = t73*t23;
    t254 = -0.18E3*t4+0.9E2*phi+0.9E2*old-0.3E2;
    t256 = t237*t38-3.0*t158*t73+6.0*t163*t169-3.0*t87*t23-6.0*t174*t245+6.0*
t93*t248-t55*t254;
    t259 = 3.0*t48*t190*t42-0.2079441542E1*t195+t200*t118+6.0*t202*t127+3.0*
t212*t107+3.0*t212*t108-0.2079441542E1*t218+3.0*t220*t75*t108+3.0*t224*t122+3.0
*t228*t42+t232*t42+t67*t256*t42;
    t260 = t48*t227;
    t321 = -3.0*t87*t188-6.0*t51*t173*t245+24.0*t33/t172/t37*t245*t59-18.0*t174
*t169*t101+6.0*t168*t248-18.0*t174*t248*t59+6.0*t93*t101*t23+6.0*t93*t73*t188-
t89*t254+2.0*t93*t254*t59-t55*(0.9E2-0.9E2*phi-0.9E2*old);
    t322 = (-0.1247664925E3+t235+t236)*t38-t237*t54*t59-3.0*t152*t54*t73+6.0*
t157*t92*t94-3.0*t158*t101+6.0*t85*t92*t169-18.0*t71*t173*t175+12.0*t163*t178
-3.0*t160*t23+6.0*t163*t182+t321;
    t325 = t211*t40;
    t327 = t227*t40;
    t333 = 0.6931471806*t325+0.2193147181E1*t327+0.15E1*t325*t41+0.15E1*t327*
t41;
    t336 = t67*t211;
    t342 = t211*t103;
    t343 = t342*t40;
    t347 = t231*t61;
    t348 = t347*t40;
    t350 = t322*t40;
    t351 = 0.6931471806*t350;
    t356 = 6.0*t260*t107+t81*t322*t42+3.0*t23*t333+3.0*t336*t122+t145*t43+6.0*
t202*t132-0.2079441542E1*t343+3.0*t336*t120-0.6931471806*t348-t351+6.0*t220*
t130+3.0*t220*t132;
    t363 = -120.0*t13+180.0*t4-0.3E2*phi-0.3E2*old;
    t366 = t327*t103;
    t370 = t81*t227;
    t373 = t81*t231;
    t380 = -t363;
    t381 = t380*t75;
    t387 = t81*t256;
    t391 = t256*t61;
    t392 = t391*t40;
    t409 = t41*t103;
    t414 = t351+0.2193147181E1*t392+0.6579441543E1*t195+0.6579441543E1*t343+
0.1107944154E2*t218+0.1107944154E2*t366+0.5193147181E1*t348+0.15E1*t350*t41+
0.15E1*t391*t42+0.45E1*t194*t42+0.45E1*t342*t42+0.45E1*t217*t108+0.45E1*t327*
t409+0.15E1*t347*t42;
    t416 = t147*t113+3.0*t363*t136-0.2079441542E1*t366+3.0*t220*t127+6.0*t370*
t127+3.0*t373*t107-3.0*t254*t75*t42+3.0*t381*t107+3.0*t380*t103*t42+t387*t108+
3.0*t254*t124+t126*t414;
    t419 = t190*t40;
    t421 = t211*t61;
    t426 = t227*t61;
    t451 = t256*t40;
    t452 = 0.6931471806*t451;
    t453 = t217*t40;
    t455 = t231*t40;
    t463 = t452+0.6579441543E1*t453+0.3693147181E1*t455+0.15E1*t451*t41+0.45E1*
t217*t42+0.15E1*t455*t41;
    t465 = -t200;
    t473 = 3.0*t260*t108+3.0*t115*(0.6931471806*t419+0.2193147181E1*t421*t40+
0.4386294362E1*t120*t103+0.3693147181E1*t426*t40+0.15E1*t419*t41+0.15E1*t421*
t42+0.3E1*t120*t409+0.15E1*t426*t42)+3.0*t370*t132+t373*t108+3.0*t381*t108+3.0*
t18*t211*t42+3.0*t224*t120+t387*t107+t73*t463+t465*t61*t42+2.0*t232*t40
-0.6931471806*t392+3.0*t228*t40;
    t516 = 3.0*t381*t42+3.0*t212*t42+t465*t43+t200*t113+3.0*t363*t124+3.0*t115*
t333+t126*t463+3.0*t260*t40+t387*t42+3.0*t220*t120+2.0*t373*t40+3.0*t260*t42+
3.0*t220*t122-0.2079441542E1*t453+t373*t42-t452-0.6931471806*t455;
    return(2.0*A_*(0.6E1*t4-0.3E1*phi-0.3E1*old+0.1E1)*t10+t18*t21+t23*t26-beta_*
t139+(0.24E2*A_*t10+t145*t21+t147*t26-beta_*(t259+t356+t416+t473))*(phi-old)/24.0
+A_*(0.12E2*phi+0.12E2*old-12.0)*t10/12.0+t465*t21/24.0+t200*t26/24.0-beta_*t516/
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

