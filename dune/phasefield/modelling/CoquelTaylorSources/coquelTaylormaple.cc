
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
  double t104;
  double t106;
  double t110;
  double t111;
  double t112;
  double t119;
  double t12;
  double t121;
  double t125;
  double t126;
  double t142;
  double t15;
  double t156;
  double t159;
  double t162;
  double t163;
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
  double t96;
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
    t96 = 120.0*t5-180.0*t4+0.3E2*phi+0.3E2*old;
    t104 = 0.8317766168E2*t5-0.1247664925E3*t4+0.2079441542E2*phi+
0.2079441542E2*old;
    t106 = t52*t55;
    t110 = 1/t54/t37;
    t111 = t33*t110;
    t112 = t60*t60;
    t119 = -0.6E2*t5+0.9E2*t4-0.15E2*phi-0.15E2*old;
    t121 = t104*t38-2.0*t106*t60+2.0*t111*t112-t56*t119;
    t125 = t62*t62;
    t126 = t19*t125;
    t142 = t54*t54;
    t156 = (0.249532985E3*t4-0.1247664925E3*phi-0.1247664925E3*old+
0.4158883084E2)*t38-3.0*t104*t55*t60+6.0*t52*t110*t112-3.0*t106*t119-6.0*t33/
t142*t112*t60+6.0*t111*t60*t119-t56*(-0.18E3*t4+0.9E2*phi+0.9E2*old-0.3E2);
    t159 = t48*t121;
    t162 = t125*t62;
    t163 = t48*t162;
    t175 = t121*t40;
    t177 = t125*t40;
    t186 = t156*t40;
    t187 = 0.6931471806*t186;
    t188 = t121*t62;
    t189 = t188*t40;
    t191 = t162*t40;
    t203 = t87*t43+3.0*t96*t62*t42+3.0*t19*t121*t42+3.0*t126*t40+t48*t156*t42+
3.0*t159*t70+2.0*t163*t40+3.0*t126*t42+3.0*t159*t72+t163*t42+t89*t67-3.0*t96*
t74+3.0*t24*(0.6931471806*t175+0.2193147181E1*t177+0.15E1*t175*t41+0.15E1*t177*
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
  double t154;
  double t157;
  double t161;
  double t165;
  double t166;
  double t169;
  double t170;
  double t172;
  double t175;
  double t18;
  double t182;
  double t187;
  double t189;
  double t19;
  double t191;
  double t194;
  double t199;
  double t20;
  double t202;
  double t203;
  double t204;
  double t205;
  double t208;
  double t21;
  double t212;
  double t218;
  double t220;
  double t221;
  double t222;
  double t224;
  double t225;
  double t228;
  double t23;
  double t233;
  double t236;
  double t237;
  double t242;
  double t245;
  double t256;
  double t257;
  double t26;
  double t260;
  double t264;
  double t265;
  double t266;
  double t274;
  double t277;
  double t28;
  double t281;
  double t284;
  double t286;
  double t29;
  double t291;
  double t298;
  double t3;
  double t309;
  double t31;
  double t311;
  double t321;
  double t33;
  double t37;
  double t38;
  double t387;
  double t388;
  double t389;
  double t390;
  double t391;
  double t392;
  double t396;
  double t397;
  double t4;
  double t40;
  double t400;
  double t401;
  double t41;
  double t417;
  double t419;
  double t42;
  double t429;
  double t43;
  double t436;
  double t438;
  double t446;
  double t447;
  double t448;
  double t450;
  double t458;
  double t463;
  double t470;
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
    t150 = t18*t149;
    t154 = 0.1247664925E3*t4;
    t157 = 0.8317766168E2*t13-t154+0.2079441542E2*phi+0.2079441542E2*old;
    t161 = t73*t73;
    t165 = t157*t38-2.0*t87*t73+2.0*t93*t161-t55*t23;
    t166 = t81*t165;
    t169 = t149*t40;
    t170 = t169*t103;
    t172 = t48*t165;
    t175 = t165*t40;
    t182 = 0.6931471806*t175+0.2193147181E1*t169+0.15E1*t175*t41+0.15E1*t169*
t41;
    t187 = 0.2079441542E2+t154-0.6238324625E2*phi-0.6238324625E2*old;
    t189 = t157*t54;
    t191 = t85*t54;
    t194 = t71*t92;
    t199 = t51*t92;
    t202 = t53*t53;
    t203 = 1/t202;
    t204 = t33*t203;
    t205 = t161*t59;
    t208 = t73*t101;
    t212 = t23*t59;
    t218 = -0.15E2-0.9E2*t4+0.45E2*phi+0.45E2*old;
    t220 = t187*t38-t189*t59-2.0*t191*t73+4.0*t194*t94-2.0*t87*t101+2.0*t199*
t161-6.0*t204*t205+4.0*t93*t208-t89*t23+2.0*t93*t212-t55*t218;
    t221 = t220*t75;
    t222 = t221*t40;
    t224 = t149*t75;
    t225 = t67*t224;
    t228 = t81*t220;
    t233 = t81*t149;
    t236 = t165*t103;
    t237 = t236*t40;
    t242 = -0.18E3*t4+0.9E2*phi+0.9E2*old-0.3E2;
    t245 = 3.0*t150*t40+6.0*t166*t130-0.2079441542E1*t170+3.0*t172*t108+3.0*t23
*t182-0.2079441542E1*t222+2.0*t225*t40+3.0*t228*t120+3.0*t166*t132+6.0*t233*
t127-0.2079441542E1*t237+3.0*t242*t124;
    t256 = 120.0*t13-180.0*t4+0.3E2*phi+0.3E2*old;
    t257 = t256*t75;
    t260 = t48*t149;
    t264 = 0.1247664925E3*phi;
    t265 = 0.1247664925E3*old;
    t266 = 0.4158883084E2+0.249532985E3*t4-t264-t265;
    t274 = t161*t73;
    t277 = t73*t23;
    t281 = t266*t38-3.0*t189*t73+6.0*t194*t161-3.0*t87*t23-6.0*t204*t274+6.0*
t93*t277-t55*t242;
    t284 = t220*t40;
    t286 = t165*t61;
    t291 = t149*t61;
    t298 = t41*t103;
    t309 = -60.0-360.0*t4+0.18E3*phi+0.18E3*old;
    t311 = t81*t281;
    t321 = 3.0*t48*t220*t42+3.0*t228*t122+t147*t113+3.0*t257*t108+6.0*t260*t107
+t67*t281*t42+3.0*t115*(0.6931471806*t284+0.2193147181E1*t286*t40+
0.4386294362E1*t120*t103+0.3693147181E1*t291*t40+0.15E1*t284*t41+0.15E1*t286*
t42+0.3E1*t120*t298+0.15E1*t291*t42)+t309*t118+t311*t107+3.0*t256*t103*t42+3.0*
t257*t107+3.0*t18*t165*t42;
    t387 = -3.0*t87*t218-6.0*t51*t203*t274+24.0*t33/t202/t37*t274*t59-18.0*t204
*t161*t101+6.0*t199*t277-18.0*t204*t277*t59+6.0*t93*t101*t23+6.0*t93*t73*t218-
t89*t242+2.0*t93*t242*t59-t55*(0.9E2-0.9E2*phi-0.9E2*old);
    t388 = (-0.1247664925E3+t264+t265)*t38-t266*t54*t59-3.0*t187*t54*t73+6.0*
t157*t92*t94-3.0*t189*t101+6.0*t85*t92*t161-18.0*t71*t203*t205+12.0*t194*t208
-3.0*t191*t23+6.0*t194*t212+t387;
    t389 = t388*t40;
    t390 = 0.6931471806*t389;
    t391 = t281*t61;
    t392 = t391*t40;
    t396 = t165*t75;
    t397 = t396*t107;
    t400 = t224*t61;
    t401 = t400*t40;
    t417 = t390+0.2193147181E1*t392+0.6579441543E1*t222+0.6579441543E1*t237+
0.1107944154E2*t397+0.1107944154E2*t170+0.5193147181E1*t401+0.15E1*t389*t41+
0.15E1*t391*t42+0.45E1*t221*t42+0.45E1*t236*t42+0.45E1*t396*t108+0.45E1*t169*
t298+0.15E1*t400*t42;
    t419 = t48*t75;
    t429 = t67*t165;
    t436 = 3.0*t166*t127+3.0*t233*t132+3.0*t260*t108+t126*t417+6.0*t419*t132+
6.0*t419*t127+3.0*t150*t42+t225*t42+t81*t388*t42+3.0*t429*t122+3.0*t172*t107+
3.0*t429*t120;
    t438 = -t309;
    t446 = t281*t40;
    t447 = 0.6931471806*t446;
    t448 = t396*t40;
    t450 = t224*t40;
    t458 = t447+0.6579441543E1*t448+0.3693147181E1*t450+0.15E1*t446*t41+0.45E1*
t396*t42+0.15E1*t450*t41;
    t463 = t81*t224;
    t470 = -t256;
    t473 = t311*t108-t390+t438*t61*t42-3.0*t242*t75*t42-0.2079441542E1*t397+t73
*t458+3.0*t166*t75*t108+t463*t108-0.6931471806*t401-0.6931471806*t392+3.0*t463*
t107+t145*t43+3.0*t470*t136;
    t516 = t438*t43+3.0*t257*t42+3.0*t172*t42+3.0*t260*t40+t311*t42+3.0*t166*
t120+2.0*t463*t40+3.0*t260*t42+3.0*t166*t122+t463*t42+t309*t113+3.0*t470*t124+
3.0*t115*t182+t126*t458-t447-0.2079441542E1*t448-0.6931471806*t450;
    return(2.0*A_*(0.6E1*t4-0.3E1*phi-0.3E1*old+0.1E1)*t10+t18*t21+t23*t26-beta_*
t139+(0.24E2*A_*t10+t145*t21+t147*t26-beta_*(t245+t321+t436+t473))*(phi-old)/24.0
+A_*(0.12E2*phi+0.12E2*old-12.0)*t10/12.0+t438*t21/24.0+t309*t26/24.0-beta_*t516/
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

