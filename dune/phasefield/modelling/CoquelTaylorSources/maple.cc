
inline double helmholtz ( double rho ,double phi ) const
{
  double t1;
  double t15;
  double t16;
  double t17;
  double t18;
  double t19;
  double t2;
  double t20;
  double t21;
  double t23;
  double t24;
  double t26;
  double t27;
  double t30;
  double t43;
  double t45;
  double t46;
  double t48;
  double t5;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t5 = t1*phi;
    t15 = t2*phi;
    t16 = 6.0*t15;
    t17 = 15.0*t2;
    t18 = 10.0*t5;
    t19 = t16-t17+t18;
    t20 = log(2.0);
    t21 = t20-0.15E1;
    t23 = log(rho);
    t24 = rho*t23;
    t26 = t20-0.15E-9;
    t27 = t26*rho;
    t30 = 1.0-t16+t17-t18;
    t43 = exp((t20-0.15E-9-t19*t20)/(0.3E1*t15-0.75E1*t2+0.5E1*t5+1.0));
    t45 = log(t43);
    t46 = t43*t45;
    t48 = t26*t43;
    return(2.0*A_*(t2+(2.0*alpha_-2.0)*t5+(-3.0*alpha_+1.0)*t1+alpha_)/delta_+beta_*(
t19*(t21*rho+0.15E1*t24-t27+0.15E1)+t30*(-rho+t24+0.2E1-t27))-t19*(t21*t43+
0.15E1*t46-t48+0.15E1)-t30*(-t43+t46+0.2E1-t48));
  }
}


inline double reactionSource ( double rho ,double phi ,double old ) const
{
  double t108;
  double t111;
  double t113;
  double t117;
  double t118;
  double t119;
  double t126;
  double t128;
  double t129;
  double t131;
  double t134;
  double t135;
  double t137;
  double t138;
  double t142;
  double t143;
  double t145;
  double t159;
  double t17;
  double t173;
  double t178;
  double t181;
  double t182;
  double t185;
  double t187;
  double t189;
  double t190;
  double t195;
  double t197;
  double t199;
  double t20;
  double t200;
  double t24;
  double t25;
  double t26;
  double t28;
  double t29;
  double t3;
  double t31;
  double t32;
  double t33;
  double t35;
  double t36;
  double t4;
  double t40;
  double t41;
  double t42;
  double t43;
  double t44;
  double t46;
  double t5;
  double t50;
  double t51;
  double t53;
  double t55;
  double t56;
  double t58;
  double t59;
  double t61;
  double t63;
  double t64;
  double t65;
  double t69;
  double t71;
  double t74;
  double t75;
  double t79;
  double t80;
  double t82;
  double t84;
  double t85;
  double t97;
  double t99;
  {
    t3 = 0.5*phi+0.5*old;
    t4 = t3*t3;
    t5 = t4*t3;
    t17 = 1/delta_;
    t20 = t4*t4;
    t24 = 30.0*t20-60.0*t5+30.0*t4;
    t25 = log(2.0);
    t26 = t25-0.15E1;
    t28 = log(rho);
    t29 = rho*t28;
    t31 = t25-0.15E-9;
    t32 = t31*rho;
    t33 = t26*rho+0.15E1*t29-t32+0.15E1;
    t35 = -t24;
    t36 = -rho+t29+0.2E1-t32;
    t40 = t20*t3;
    t41 = 6.0*t40;
    t42 = 15.0*t20;
    t43 = 10.0*t5;
    t44 = t41-t42+t43;
    t46 = t25-0.15E-9-t44*t25;
    t50 = 0.3E1*t40-0.75E1*t20+0.5E1*t5+1.0;
    t51 = 1/t50;
    t53 = exp(t46*t51);
    t55 = log(t53);
    t56 = t53*t55;
    t58 = t31*t53;
    t59 = t26*t53+0.15E1*t56-t58+0.15E1;
    t61 = t24*t25;
    t63 = t50*t50;
    t64 = 1/t63;
    t65 = t46*t64;
    t69 = 0.15E2*t20-0.3E2*t5+0.15E2*t4;
    t71 = -t61*t51-t65*t69;
    t74 = t71*t53;
    t75 = t74*t55;
    t79 = t31*t71*t53;
    t80 = t26*t71*t53+0.15E1*t75+0.15E1*t74-t79;
    t82 = -t53+t56+0.2E1-t58;
    t84 = 1.0-t41+t42-t43;
    t85 = t75-t79;
    t97 = 360.0*t4-0.18E3*phi-0.18E3*old+60.0;
    t99 = -t97;
    t108 = 120.0*t5-180.0*t4+0.3E2*phi+0.3E2*old;
    t111 = t108*t25;
    t113 = t64*t69;
    t117 = 1/t63/t50;
    t118 = t46*t117;
    t119 = t69*t69;
    t126 = 0.6E2*t5-0.9E2*t4+0.15E2*phi+0.15E2*old;
    t128 = -t111*t51+2.0*t61*t113+2.0*t118*t119-t65*t126;
    t129 = t26*t128;
    t131 = t71*t71;
    t134 = t128*t53;
    t135 = t134*t55;
    t137 = t131*t53;
    t138 = t137*t55;
    t142 = t31*t128;
    t143 = t142*t53;
    t145 = t31*t131*t53;
    t159 = t63*t63;
    t173 = -t97*t25*t51+3.0*t111*t113-6.0*t61*t117*t119+3.0*t61*t64*t126-6.0*
t46/t159*t119*t69+6.0*t118*t69*t126-t65*(0.18E3*t4-0.9E2*phi-0.9E2*old+0.3E2);
    t178 = t131*t71;
    t181 = t173*t53;
    t182 = t181*t55;
    t185 = t128*t71*t56;
    t187 = t134*t71;
    t189 = t178*t53;
    t190 = t189*t55;
    t195 = t31*t173*t53;
    t197 = 3.0*t142*t74;
    t199 = t31*t178*t53;
    t200 = t26*t173*t53+3.0*t129*t74+t26*t178*t53+0.15E1*t182+0.45E1*t185+0.9E1
*t187+0.15E1*t190+0.45E1*t189+0.15E1*t181-t195-t197-t199;
    return(2.0*A_*(4.0*t5+3.0*(2.0*alpha_-2.0)*t4+2.0*(-3.0*alpha_+1.0)*t3)*t17+
beta_*(t24*t33+t35*t36)-t24*t59-t44*t80-t35*t82-t84*t85+(2.0*A_*(0.12E2*phi+
0.12E2*old+12.0*alpha_-12.0)*t17+beta_*(t97*t33+t99*t36)-t97*t59-3.0*t108*t80-3.0
*t24*(t129*t53+t26*t131*t53+0.15E1*t135+0.15E1*t138+0.3E1*t137+0.15E1*t134-t143
-t145)-t44*t200-t99*t82+3.0*t108*t85-3.0*t35*(t135+t138+t137-t143-t145)-t84*(
t182+3.0*t185+3.0*t187+t190+2.0*t189-t195-t197-t199))*(phi-old)/24.0);
  }
}


inline double dphireactionSource ( double rho ,double phi ,double old ) const
{
  double t100;
  double t102;
  double t104;
  double t106;
  double t110;
  double t111;
  double t119;
  double t121;
  double t125;
  double t126;
  double t129;
  double t13;
  double t131;
  double t135;
  double t136;
  double t137;
  double t139;
  double t141;
  double t142;
  double t144;
  double t145;
  double t147;
  double t148;
  double t155;
  double t157;
  double t16;
  double t165;
  double t170;
  double t177;
  double t180;
  double t182;
  double t188;
  double t189;
  double t191;
  double t192;
  double t194;
  double t195;
  double t197;
  double t198;
  double t202;
  double t203;
  double t204;
  double t205;
  double t206;
  double t209;
  double t21;
  double t214;
  double t215;
  double t218;
  double t22;
  double t221;
  double t225;
  double t226;
  double t227;
  double t23;
  double t240;
  double t242;
  double t243;
  double t249;
  double t25;
  double t250;
  double t253;
  double t255;
  double t257;
  double t258;
  double t26;
  double t261;
  double t263;
  double t265;
  double t268;
  double t269;
  double t270;
  double t272;
  double t273;
  double t274;
  double t277;
  double t28;
  double t281;
  double t289;
  double t29;
  double t293;
  double t294;
  double t298;
  double t299;
  double t3;
  double t30;
  double t301;
  double t302;
  double t304;
  double t305;
  double t307;
  double t309;
  double t310;
  double t314;
  double t315;
  double t317;
  double t318;
  double t319;
  double t32;
  double t320;
  double t322;
  double t324;
  double t33;
  double t344;
  double t355;
  double t37;
  double t38;
  double t384;
  double t387;
  double t39;
  double t390;
  double t391;
  double t394;
  double t397;
  double t4;
  double t40;
  double t400;
  double t401;
  double t402;
  double t404;
  double t407;
  double t409;
  double t41;
  double t412;
  double t416;
  double t417;
  double t42;
  double t421;
  double t423;
  double t424;
  double t427;
  double t432;
  double t434;
  double t438;
  double t44;
  double t442;
  double t447;
  double t450;
  double t456;
  double t462;
  double t473;
  double t475;
  double t48;
  double t49;
  double t507;
  double t51;
  double t53;
  double t54;
  double t56;
  double t57;
  double t62;
  double t66;
  double t69;
  double t70;
  double t71;
  double t75;
  double t77;
  double t80;
  double t81;
  double t85;
  double t86;
  double t88;
  double t91;
  double t92;
  double t94;
  double t95;
  double t98;
  double t99;
  {
    t3 = 0.5*phi+0.5*old;
    t4 = t3*t3;
    t13 = 1/delta_;
    t16 = t4*t3;
    t21 = 0.6E2*t16-0.9E2*t4+0.15E2*phi+0.15E2*old;
    t22 = log(2.0);
    t23 = t22-0.15E1;
    t25 = log(rho);
    t26 = rho*t25;
    t28 = t22-0.15E-9;
    t29 = rho*t28;
    t30 = t23*rho+0.15E1*t26-t29+0.15E1;
    t32 = -t21;
    t33 = -rho+t26+0.2E1-t29;
    t37 = t4*t4;
    t38 = t37*t3;
    t39 = 6.0*t38;
    t40 = 15.0*t37;
    t41 = 10.0*t16;
    t42 = t39-t40+t41;
    t44 = t22-0.15E-9-t42*t22;
    t48 = 0.3E1*t38-0.75E1*t37+0.5E1*t16+1.0;
    t49 = 1/t48;
    t51 = exp(t44*t49);
    t53 = log(t51);
    t54 = t51*t53;
    t56 = t28*t51;
    t57 = t23*t51+0.15E1*t54-t56+0.15E1;
    t62 = 30.0*t37-60.0*t16+30.0*t4;
    t66 = 0.15E2*t37-0.3E2*t16+0.15E2*t4;
    t69 = t48*t48;
    t70 = 1/t69;
    t71 = t44*t70;
    t75 = 0.75E1*t37-0.15E2*t16+0.75E1*t4;
    t77 = -t66*t22*t49-t71*t75;
    t80 = t77*t51;
    t81 = t80*t53;
    t85 = t28*t77*t51;
    t86 = t23*t77*t51+0.15E1*t81+0.15E1*t80-t85;
    t88 = t62*t22;
    t91 = -t88*t49-t71*t66;
    t92 = t23*t91;
    t94 = t91*t51;
    t95 = t94*t53;
    t98 = t28*t91;
    t99 = t98*t51;
    t100 = t92*t51+0.15E1*t95+0.15E1*t94-t99;
    t102 = t21*t22;
    t104 = t70*t75;
    t106 = t66*t66;
    t110 = 1/t69/t48;
    t111 = t44*t110;
    t119 = 0.3E2*t16-0.45E2*t4+0.75E1*phi+0.75E1*old;
    t121 = -t102*t49+t88*t104+t106*t22*t70+2.0*t111*t66*t75-t71*t119;
    t125 = t121*t51;
    t126 = t125*t53;
    t129 = t91*t77*t54;
    t131 = t94*t77;
    t135 = t28*t121*t51;
    t136 = t98*t80;
    t137 = t23*t121*t51+t92*t80+0.15E1*t126+0.15E1*t129+0.3E1*t131+0.15E1*t125-
t135-t136;
    t139 = 0.2E1-t51+t54-t56;
    t141 = -t62;
    t142 = t81-t85;
    t144 = -t66;
    t145 = t95-t99;
    t147 = 1.0-t39+t40-t41;
    t148 = t126+t129+t131-t135-t136;
    t155 = -0.18E3+0.18E3*phi+0.18E3*old;
    t157 = -t155;
    t165 = 360.0*t4-0.18E3*phi-0.18E3*old+60.0;
    t170 = 0.18E3*t4-0.9E2*phi-0.9E2*old+0.3E2;
    t177 = 120.0*t16-180.0*t4+0.3E2*phi+0.3E2*old;
    t180 = t177*t22;
    t182 = t70*t66;
    t188 = -t180*t49+2.0*t88*t182+2.0*t111*t106-t71*t21;
    t189 = t23*t188;
    t191 = t91*t91;
    t192 = t23*t191;
    t194 = t188*t51;
    t195 = t194*t53;
    t197 = t191*t51;
    t198 = t197*t53;
    t202 = t28*t188;
    t203 = t202*t51;
    t204 = t28*t191;
    t205 = t204*t51;
    t206 = t189*t51+t192*t51+0.15E1*t195+0.15E1*t198+0.3E1*t197+0.15E1*t194-
t203-t205;
    t209 = t170*t22;
    t214 = t110*t66;
    t215 = t214*t75;
    t218 = t70*t119;
    t221 = t106*t66;
    t225 = t69*t69;
    t226 = 1/t225;
    t227 = t44*t226;
    t240 = 0.15E2+0.9E2*t4-0.45E2*phi-0.45E2*old;
    t242 = -t209*t49+t180*t104+3.0*t102*t182-4.0*t88*t215+2.0*t88*t218-2.0*t221
*t22*t110-6.0*t227*t106*t75+4.0*t111*t66*t119+2.0*t111*t21*t75-t71*t240;
    t243 = t23*t242;
    t249 = t242*t51;
    t250 = t249*t53;
    t253 = t188*t77*t54;
    t255 = t194*t77;
    t257 = t53*t121;
    t258 = t94*t257;
    t261 = t191*t77*t54;
    t263 = t197*t77;
    t265 = t94*t121;
    t268 = t28*t242;
    t269 = t268*t51;
    t270 = t202*t80;
    t272 = 2.0*t98*t125;
    t273 = t204*t80;
    t274 = t243*t51+t189*t80+2.0*t92*t125+t192*t80+0.15E1*t250+0.15E1*t253+
0.3E1*t255+0.3E1*t258+0.15E1*t261+0.45E1*t263+0.6E1*t265+0.15E1*t249-t269-t270-
t272-t273;
    t277 = t165*t22;
    t281 = t110*t106;
    t289 = t66*t21;
    t293 = -t277*t49+3.0*t180*t182-6.0*t88*t281+3.0*t88*t70*t21-6.0*t227*t221+
6.0*t111*t289-t71*t170;
    t294 = t23*t293;
    t298 = t191*t91;
    t299 = t23*t298;
    t301 = t293*t51;
    t302 = t301*t53;
    t304 = t188*t91;
    t305 = t304*t54;
    t307 = t194*t91;
    t309 = t298*t51;
    t310 = t309*t53;
    t314 = t28*t293;
    t315 = t314*t51;
    t317 = 3.0*t202*t94;
    t318 = t28*t298;
    t319 = t318*t51;
    t320 = t294*t51+3.0*t189*t94+t299*t51+0.15E1*t302+0.45E1*t305+0.9E1*t307+
0.15E1*t310+0.45E1*t309+0.15E1*t301-t315-t317-t319;
    t322 = t197*t121;
    t324 = t309*t77;
    t344 = t21*t21;
    t355 = t106*t106;
    t384 = -t155*t22*t49+t277*t104+4.0*t209*t182-6.0*t180*t215+3.0*t180*t218
-12.0*t102*t281+18.0*t88*t226*t106*t75-12.0*t88*t214*t119+3.0*t344*t22*t70-6.0*
t88*t110*t21*t75+3.0*t88*t70*t240+6.0*t355*t22*t226+24.0*t44/t225/t48*t221*t75
-18.0*t227*t106*t119-18.0*t227*t289*t75+6.0*t111*t119*t21+6.0*t111*t66*t240+2.0
*t111*t170*t75-t71*(-0.9E2+0.9E2*phi+0.9E2*old);
    t387 = t194*t121;
    t390 = 3.0*t268*t94;
    t391 = t249*t91;
    t394 = t293*t77*t54;
    t397 = t242*t91*t54;
    t400 = 3.0*t202*t131;
    t401 = t384*t51;
    t402 = t401*t53;
    t404 = t197*t257;
    t407 = t298*t77*t54;
    t409 = t318*t80;
    t412 = 0.135E2*t322+0.6E1*t324+t23*t384*t51+0.9E1*t387-t390+0.9E1*t391+
0.15E1*t394+0.45E1*t397-t400+0.15E1*t402+0.45E1*t404+0.15E1*t407-t409+3.0*t189*
t125;
    t416 = t28*t384*t51;
    t417 = t301*t77;
    t421 = t314*t80;
    t423 = 3.0*t202*t125;
    t424 = t304*t80;
    t427 = t304*t81;
    t432 = 3.0*t204*t125;
    t434 = t188*t121*t54;
    t438 = 3.0*t243*t94-t416+0.3E1*t417+3.0*t189*t131-t421-t423+0.135E2*t424+
0.15E1*t401+0.45E1*t427+3.0*t192*t125-t432+0.45E1*t434+t294*t80+t299*t80;
    t442 = -t165;
    t447 = -t177;
    t450 = t195+t198+t197-t203-t205;
    t456 = t250+t253+t255+2.0*t258+t261+2.0*t263+2.0*t265-t269-t270-t272-t273;
    t462 = t302+3.0*t305+3.0*t307+t310+2.0*t309-t315-t317-t319;
    t473 = t402+t394+t417+3.0*t397+3.0*t434+3.0*t427+6.0*t424+3.0*t391+3.0*t387
+3.0*t404+t407+3.0*t324+6.0*t322-t416-t421-t390-t423-t400-t432-t409;
    t475 = 0.24E2*A_*t13+beta_*(t155*t30+t157*t33)-t155*t57-t165*t86-3.0*t170*
t100-3.0*t177*t137-3.0*t21*t206-3.0*t62*t274-t66*t320-t42*(t412+t438)-t157*t139
-t442*t142+3.0*t170*t145-3.0*t447*t148-3.0*t32*t450-3.0*t141*t456-t144*t462-
t147*t473;
    t507 = t475*(phi-old)/24.0+A_*(0.12E2*phi+0.12E2*old+12.0*alpha_-12.0)*t13/
12.0+beta_*(t165*t30+t442*t33)/24.0-t165*t57/24.0-t177*t100/8.0-t62*t206/8.0-t42
*t320/24.0-t442*t139/24.0-t447*t145/8.0-t141*t450/8.0-t147*t462/24.0;
    return(2.0*A_*(0.6E1*t4+0.3E1*(2.0*alpha_-2.0)*t3-0.3E1*alpha_+0.1E1)*t13+beta_
*(t21*t30+t32*t33)-t21*t57-t62*t86-t66*t100-t42*t137-t32*t139-t141*t142-t144*
t145-t147*t148+t507);
  }
}


inline double chemicalPotential ( double rho ,double phi ,double old ) const
{
  double t1;
  double t11;
  double t12;
  double t16;
  double t17;
  double t2;
  double t22;
  double t23;
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
    t17 = log(2.0);
    t22 = t11*t11;
    t23 = 1/t22;
    return(beta_*(t8*(0.15E-9+0.15E1*t12)+t16*(0.15E-9+t12-t17))+beta_*(-0.15E1*
t8*t23-t16*t23)*(rho-old)/24.0);
  }
}

inline double drhochemicalPotential ( double rho ,double phi ,double old ) const
{
  double t1;
  double t11;
  double t12;
  double t15;
  double t2;
  double t20;
  double t22;
  double t32;
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
    t20 = t11*t11;
    t22 = 1/t20/t11;
    t32 = 1/t20;
    return(beta_*(0.75*t8*t12+0.5*t15*t12)+beta_*(0.15E1*t8*t22+0.1E1*t15*t22)*(
rho-old)/24.0+beta_*(-0.15E1*t8*t32-t15*t32)/24.0);
  }
}


inline double dphichemicalPotential ( double rho ,double phi ,double old ) const
{
  double t1;
  double t10;
  double t11;
  double t15;
  double t16;
  double t2;
  double t21;
  double t22;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t7 = 30.0*t2-60.0*t1*phi+30.0*t1;
    t10 = 0.5*rho+0.5*old;
    t11 = log(t10);
    t15 = -t7;
    t16 = log(2.0);
    t21 = t10*t10;
    t22 = 1/t21;
    return(beta_*(t7*(0.15E-9+0.15E1*t11)+t15*(0.15E-9+t11-t16))+beta_*(-0.15E1*
t7*t22-t15*t22)*(rho-old)/24.0);
  }
}


inline double pressure ( double rho ,double phi ) const
{
  double t1;
  double t12;
  double t2;
  double t4;
  double t5;
  double t6;
  double t7;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t6 = t1*phi;
    t7 = 10.0*t6;
    t9 = log(2.0);
    t12 = log(rho);
    return(beta_*((t4-t5+t7)*(-(t9-0.15E1)*rho-0.15E1*rho*t12+rho*(t9+0.15E1*t12
))+(1.0-t4+t5-t7)*(rho-0.5))-2.0*A_*(t2+(2.0*alpha_-2.0)*t6+(-3.0*alpha_+1.0)*t1+
alpha_)/delta_);
  }
}


inline double a ( double rho ,double phi ) const
{
  double t1;
  double t10;
  double t2;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t10 = sqrt(beta_*(6.0*t2*phi-15.0*t2+10.0*t1*phi+2.0));
    return(0.7071067812*t10);
  }
}

