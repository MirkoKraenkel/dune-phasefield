
inline double helmholtz ( double rho ,double phi ) const
{
  double t1;
  double t15;
  double t16;
  double t17;
  double t18;
  double t19;
  double t2;
  double t21;
  double t22;
  double t26;
  double t41;
  double t43;
  double t44;
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
    t21 = log(rho);
    t22 = rho*t21;
    t26 = 1.0-t16+t17-t18;
    t41 = exp((-0.2889626582E2+0.1204277045E3*t15-0.3010692614E3*t2+
0.2007128409E3*t5)/(0.9502053219E1+0.116654747E3*t15-0.2916368674E3*t2+
0.1944245783E3*t5));
    t43 = log(t41);
    t44 = t41*t43;
    return(2.0*A_*(t2+(2.0*alpha_-2.0)*t5+(-3.0*alpha_+1.0)*t1+alpha_)/delta_+t19*(
-0.2011952932E2*rho+0.289445110522E2*t22+0.2133795654E2)+t26*(0.193942126E2*rho
+0.950205321859E1*t22+0.4540504209)-beta_*(t19*(0.2133795654E2-0.2011952932E2*
t41+0.289445110522E2*t44)+t26*(0.4540504209+0.193942126E2*t41+0.950205321859E1*
t44)));
  }
}


inline double reactionSource ( double rho ,double phi ,double old ) const
{
  double t101;
  double t108;
  double t115;
  double t117;
  double t121;
  double t122;
  double t123;
  double t130;
  double t132;
  double t133;
  double t135;
  double t136;
  double t138;
  double t140;
  double t158;
  double t17;
  double t173;
  double t175;
  double t176;
  double t179;
  double t181;
  double t183;
  double t185;
  double t20;
  double t24;
  double t26;
  double t27;
  double t29;
  double t3;
  double t31;
  double t34;
  double t36;
  double t4;
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
  double t61;
  double t63;
  double t64;
  double t65;
  double t69;
  double t71;
  double t72;
  double t74;
  double t76;
  double t80;
  double t82;
  double t85;
  double t99;
  {
    t3 = 0.5*phi+0.5*old;
    t4 = t3*t3;
    t5 = t4*t3;
    t17 = 1/delta_;
    t20 = t4*t4;
    t24 = 30.0*t20-60.0*t5+30.0*t4;
    t26 = log(rho);
    t27 = rho*t26;
    t29 = -0.2011952932E2*rho+0.289445110522E2*t27+0.2133795654E2;
    t31 = -t24;
    t34 = 0.193942126E2*rho+0.950205321859E1*t27+0.4540504209;
    t36 = t20*t3;
    t40 = -0.2889626582E2+0.1204277045E3*t36-0.3010692614E3*t20+0.2007128409E3*
t5;
    t44 = 0.9502053219E1+0.116654747E3*t36-0.2916368674E3*t20+0.1944245783E3*t5
;
    t45 = 1/t44;
    t47 = exp(t40*t45);
    t49 = log(t47);
    t50 = t47*t49;
    t52 = 0.2133795654E2-0.2011952932E2*t47+0.289445110522E2*t50;
    t54 = 6.0*t36;
    t55 = 15.0*t20;
    t56 = 10.0*t5;
    t57 = t54-t55+t56;
    t61 = 0.6021385225E3*t20-0.1204277046E4*t5+0.6021385227E3*t4;
    t63 = t44*t44;
    t64 = 1/t63;
    t65 = t40*t64;
    t69 = 0.583273735E3*t20-0.116654747E4*t5+0.5832737349E3*t4;
    t71 = t61*t45-t65*t69;
    t72 = t71*t47;
    t74 = t72*t49;
    t76 = 0.882498173E1*t72+0.289445110522E2*t74;
    t80 = 0.4540504209+0.193942126E2*t47+0.950205321859E1*t50;
    t82 = 1.0-t54+t55-t56;
    t85 = 0.2889626582E2*t72+0.950205321859E1*t74;
    t99 = 360.0*t4-0.18E3*phi-0.18E3*old+60.0;
    t101 = -t99;
    t108 = 120.0*t5-180.0*t4+0.3E2*phi+0.3E2*old;
    t115 = 0.240855409E4*t5-0.3612831138E4*t4+0.6021385225E3*phi+0.6021385225E3
*old;
    t117 = t61*t64;
    t121 = 1/t63/t44;
    t122 = t40*t121;
    t123 = t69*t69;
    t130 = 0.233309494E4*t5-0.349964241E4*t4+0.583273735E3*phi+0.583273735E3*
old;
    t132 = t115*t45-2.0*t117*t69+2.0*t122*t123-t65*t130;
    t133 = t132*t47;
    t135 = t71*t71;
    t136 = t135*t47;
    t138 = t133*t49;
    t140 = t136*t49;
    t158 = t63*t63;
    t173 = ((0.722566227E4*t4-0.3612831138E4*phi-0.3612831138E4*old+
0.1204277045E4)*t45-3.0*t115*t64*t69+6.0*t61*t121*t123-3.0*t117*t130-6.0*t40/
t158*t123*t69+6.0*t122*t69*t130-t65*(0.699928482E4*t4-0.349964241E4*phi
-0.349964241E4*old+0.116654747E4))*t47;
    t175 = t132*t71;
    t176 = t175*t47;
    t179 = t135*t71*t47;
    t181 = t173*t49;
    t183 = t175*t50;
    t185 = t179*t49;
    return(2.0*A_*(4.0*t5+3.0*(2.0*alpha_-2.0)*t4+2.0*(-3.0*alpha_+1.0)*t3)*t17+
t24*t29+t31*t34-beta_*(t24*t52+t57*t76+t31*t80+t82*t85)+(2.0*A_*(0.12E2*phi+
0.12E2*old+12.0*alpha_-12.0)*t17+t99*t29+t101*t34-beta_*(t99*t52+3.0*t108*t76+3.0
*t24*(0.882498173E1*t133+0.3776949278E2*t136+0.289445110522E2*t138+
0.289445110522E2*t140)+t57*(0.882498173E1*t173+0.1133084783E3*t176+
0.6671400383E2*t179+0.289445110522E2*t181+0.8683353315E2*t183+0.289445110522E2*
t185)+t101*t80-3.0*t108*t85+3.0*t31*(0.2889626582E2*t133+0.3839831904E2*t136+
0.950205321859E1*t138+0.950205321859E1*t140)+t82*(0.2889626582E2*t173+
0.1151949571E3*t176+0.4790037226E2*t179+0.950205321859E1*t181+0.2850615966E2*
t183+0.950205321859E1*t185)))*(phi-old)/24.0);
  }
}


inline double dphireactionSource ( double rho ,double phi ,double old ) const
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
  double t13;
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
  double t161;
  double t163;
  double t169;
  double t174;
  double t181;
  double t188;
  double t192;
  double t196;
  double t199;
  double t201;
  double t202;
  double t204;
  double t205;
  double t207;
  double t209;
  double t21;
  double t211;
  double t217;
  double t219;
  double t221;
  double t224;
  double t229;
  double t23;
  double t232;
  double t233;
  double t234;
  double t235;
  double t238;
  double t24;
  double t242;
  double t247;
  double t249;
  double t250;
  double t252;
  double t253;
  double t255;
  double t257;
  double t258;
  double t26;
  double t260;
  double t262;
  double t264;
  double t265;
  double t267;
  double t275;
  double t28;
  double t283;
  double t286;
  double t290;
  double t291;
  double t292;
  double t294;
  double t295;
  double t297;
  double t298;
  double t3;
  double t300;
  double t301;
  double t303;
  double t305;
  double t307;
  double t309;
  double t31;
  double t33;
  double t34;
  double t369;
  double t371;
  double t373;
  double t374;
  double t376;
  double t377;
  double t379;
  double t38;
  double t380;
  double t382;
  double t384;
  double t386;
  double t387;
  double t389;
  double t391;
  double t393;
  double t395;
  double t397;
  double t399;
  double t4;
  double t401;
  double t403;
  double t406;
  double t411;
  double t418;
  double t42;
  double t43;
  double t438;
  double t45;
  double t454;
  double t456;
  double t47;
  double t48;
  double t50;
  double t55;
  double t59;
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
    t3 = 0.5*phi+0.5*old;
    t4 = t3*t3;
    t13 = 1/delta_;
    t16 = t4*t3;
    t21 = 0.6E2*t16-0.9E2*t4+0.15E2*phi+0.15E2*old;
    t23 = log(rho);
    t24 = rho*t23;
    t26 = -0.2011952932E2*rho+0.289445110522E2*t24+0.2133795654E2;
    t28 = -t21;
    t31 = 0.193942126E2*rho+0.950205321859E1*t24+0.4540504209;
    t33 = t4*t4;
    t34 = t33*t3;
    t38 = -0.2889626582E2+0.1204277045E3*t34-0.3010692614E3*t33+0.2007128409E3*
t16;
    t42 = 0.9502053219E1+0.116654747E3*t34-0.2916368674E3*t33+0.1944245783E3*
t16;
    t43 = 1/t42;
    t45 = exp(t38*t43);
    t47 = log(t45);
    t48 = t45*t47;
    t50 = 0.2133795654E2-0.2011952932E2*t45+0.289445110522E2*t48;
    t55 = 30.0*t33-60.0*t16+30.0*t4;
    t59 = 0.3010692612E3*t33-0.6021385228E3*t16+0.3010692614E3*t4;
    t61 = t42*t42;
    t62 = 1/t61;
    t63 = t38*t62;
    t67 = 0.2916368675E3*t33-0.5832737348E3*t16+0.2916368674E3*t4;
    t69 = t59*t43-t63*t67;
    t70 = t69*t45;
    t72 = t70*t47;
    t74 = 0.882498173E1*t70+0.289445110522E2*t72;
    t79 = 0.15E2*t33-0.3E2*t16+0.15E2*t4;
    t83 = 0.6021385225E3*t33-0.1204277046E4*t16+0.6021385227E3*t4;
    t86 = 0.116654747E4*t16;
    t88 = 0.583273735E3*t33-t86+0.5832737349E3*t4;
    t90 = t83*t43-t63*t88;
    t91 = t90*t45;
    t93 = t91*t47;
    t95 = 0.882498173E1*t91+0.289445110522E2*t93;
    t97 = 6.0*t34;
    t98 = 15.0*t33;
    t99 = 10.0*t16;
    t100 = t97-t98+t99;
    t105 = 0.1204277045E4*t16-0.1806415569E4*t4+0.3010692614E3*phi+
0.3010692614E3*old;
    t107 = t83*t62;
    t109 = t59*t62;
    t112 = 1/t61/t42;
    t113 = t38*t112;
    t114 = t88*t67;
    t120 = t86-0.1749821205E4*t4+0.2916368674E3*phi+0.2916368674E3*old;
    t122 = t105*t43-t107*t67-t109*t88+2.0*t114*t113-t63*t120;
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
    t161 = -0.18E3+0.18E3*phi+0.18E3*old;
    t163 = -t161;
    t169 = 360.0*t4-0.18E3*phi-0.18E3*old+60.0;
    t174 = 0.3E2+0.18E3*t4-0.9E2*phi-0.9E2*old;
    t181 = 120.0*t16-180.0*t4+0.3E2*phi+0.3E2*old;
    t188 = 0.240855409E4*t16-0.3612831138E4*t4+0.6021385225E3*phi+
0.6021385225E3*old;
    t192 = t88*t88;
    t196 = 0.349964241E4*t4;
    t199 = 0.233309494E4*t16-t196+0.583273735E3*phi+0.583273735E3*old;
    t201 = t188*t43-2.0*t107*t88+2.0*t113*t192-t63*t199;
    t202 = t201*t45;
    t204 = t90*t90;
    t205 = t204*t45;
    t207 = t202*t47;
    t209 = t205*t47;
    t211 = 0.882498173E1*t202+0.3776949278E2*t205+0.289445110522E2*t207+
0.289445110522E2*t209;
    t217 = 0.6021385225E3+0.3612831135E4*t4-0.1806415569E4*phi-0.1806415569E4*
old;
    t219 = t188*t62;
    t221 = t105*t62;
    t224 = t83*t112;
    t229 = t59*t112;
    t232 = t61*t61;
    t233 = 1/t232;
    t234 = t38*t233;
    t235 = t192*t67;
    t238 = t88*t120;
    t242 = t199*t67;
    t247 = 0.583273735E3+t196-0.1749821205E4*phi-0.1749821205E4*old;
    t249 = t217*t43-t219*t67-2.0*t221*t88+4.0*t224*t114-2.0*t107*t120+2.0*t229*
t192-6.0*t234*t235+4.0*t113*t238-t109*t199+2.0*t113*t242-t63*t247;
    t250 = t249*t45;
    t252 = t201*t69;
    t253 = t252*t45;
    t255 = t91*t122;
    t257 = t204*t69;
    t258 = t257*t45;
    t260 = t250*t47;
    t262 = t252*t48;
    t264 = t47*t122;
    t265 = t91*t264;
    t267 = t257*t48;
    t275 = 0.722566227E4*t4-0.3612831138E4*phi-0.3612831138E4*old+
0.1204277045E4;
    t283 = t192*t88;
    t286 = t88*t199;
    t290 = 0.349964241E4*phi;
    t291 = 0.349964241E4*old;
    t292 = 0.116654747E4+0.699928482E4*t4-t290-t291;
    t294 = t275*t43-3.0*t219*t88+6.0*t224*t192-3.0*t107*t199-6.0*t234*t283+6.0*
t113*t286-t63*t292;
    t295 = t294*t45;
    t297 = t201*t90;
    t298 = t297*t45;
    t300 = t204*t90;
    t301 = t300*t45;
    t303 = t295*t47;
    t305 = t297*t48;
    t307 = t301*t47;
    t309 = 0.882498173E1*t295+0.1133084783E3*t298+0.6671400383E2*t301+
0.289445110522E2*t303+0.8683353315E2*t305+0.289445110522E2*t307;
    t369 = -3.0*t107*t247-6.0*t59*t233*t283+24.0*t38/t232/t42*t283*t67-18.0*
t234*t192*t120+6.0*t229*t286-18.0*t234*t286*t67+6.0*t113*t120*t199+6.0*t113*t88
*t247-t109*t292+2.0*t113*t292*t67-t63*(-0.349964241E4+t290+t291);
    t371 = ((-0.3612831138E4+0.3612831135E4*phi+0.3612831135E4*old)*t43-t275*
t62*t67-3.0*t217*t62*t88+6.0*t188*t112*t114-3.0*t219*t120+6.0*t105*t112*t192
-18.0*t83*t233*t235+12.0*t224*t238-3.0*t221*t199+6.0*t224*t242+t369)*t45;
    t373 = t294*t69;
    t374 = t373*t45;
    t376 = t249*t90;
    t377 = t376*t45;
    t379 = t201*t122;
    t380 = t379*t45;
    t382 = t297*t70;
    t384 = t122*t205;
    t386 = t300*t69;
    t387 = t386*t45;
    t389 = t371*t47;
    t391 = t373*t48;
    t393 = t376*t48;
    t395 = t379*t48;
    t397 = t297*t72;
    t399 = t205*t264;
    t401 = t386*t48;
    t403 = 0.882498173E1*t371+0.3776949278E2*t374+0.1133084783E3*t377+
0.1133084783E3*t380+0.2001420114E3*t382+0.2001420115E3*t384+0.9565851488E2*t387
+0.289445110522E2*t389+0.289445110522E2*t391+0.8683353315E2*t393+0.8683353315E2
*t395+0.8683353315E2*t397+0.8683353315E2*t399+0.289445110522E2*t401;
    t406 = -t169;
    t411 = -t181;
    t418 = 0.2889626582E2*t202+0.3839831904E2*t205+0.950205321859E1*t207+
0.950205321859E1*t209;
    t438 = 0.2889626582E2*t295+0.1151949571E3*t298+0.4790037226E2*t301+
0.950205321859E1*t303+0.2850615966E2*t305+0.950205321859E1*t307;
    t454 = 0.2889626582E2*t371+0.3839831904E2*t374+0.1151949571E3*t377+
0.1151949571E3*t380+0.1437011168E3*t382+0.1437011168E3*t384+0.5740242548E2*t387
+0.950205321859E1*t389+0.950205321859E1*t391+0.2850615966E2*t393+0.2850615966E2
*t395+0.2850615966E2*t397+0.2850615966E2*t399+0.950205321859E1*t401;
    t456 = t161*t50+t169*t74+3.0*t174*t95+3.0*t181*t132+3.0*t21*t211+3.0*t55*(
0.882498173E1*t250+0.3776949278E2*t253+0.7553898556E2*t255+0.6671400383E2*t258+
0.289445110522E2*t260+0.289445110522E2*t262+0.578890221E2*t265+0.289445110522E2
*t267)+t79*t309+t100*t403+t163*t136+t406*t141-3.0*t174*t146+3.0*t411*t153+3.0*
t28*t418+3.0*t138*(0.2889626582E2*t250+0.3839831904E2*t253+0.7679663808E2*t255+
0.4790037226E2*t258+0.950205321859E1*t260+0.950205321859E1*t262+0.1900410644E2*
t265+0.950205321859E1*t267)+t143*t438+t148*t454;
    return(2.0*A_*(0.1E1+0.6E1*t4+0.3E1*(2.0*alpha_-2.0)*t3-0.3E1*alpha_)*t13+t21*
t26+t28*t31-beta_*(t21*t50+t55*t74+t79*t95+t100*t132+t28*t136+t138*t141+t143*
t146+t148*t153)+(0.24E2*A_*t13+t161*t26+t163*t31-beta_*t456)*(phi-old)/24.0+A_*(
0.12E2*phi+0.12E2*old+12.0*alpha_-12.0)*t13/12.0+t169*t26/24.0+t406*t31/24.0-
beta_*(t169*t50+3.0*t181*t95+3.0*t55*t211+t100*t309+t406*t136+3.0*t411*t146+3.0*
t138*t418+t148*t438)/24.0);
  }
}


inline double drhoreactionSource ( double rho ,double phi ,double old ) const
{
  double t10;
  double t11;
  double t13;
  double t17;
  double t22;
  double t3;
  double t4;
  double t5;
  {
    t3 = 0.5*phi+0.5*old;
    t4 = t3*t3;
    t5 = t4*t4;
    t10 = 30.0*t5-60.0*t4*t3+30.0*t4;
    t11 = log(rho);
    t13 = 0.882498173E1+0.289445110522E2*t11;
    t17 = 0.2889626582E2+0.950205321859E1*t11;
    t22 = 360.0*t4-0.18E3*phi-0.18E3*old+60.0;
    return(t10*t13-t10*t17+(t22*t13-t22*t17)*(phi-old)/24.0);
  }
}


inline double chemicalPotential ( double rho ,double phi ,double old ) const
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
    return(0.1447225552E2*t8*t12+0.475102661E1*t15*t12+(0.2894451105E2*t8*t20+
0.9502053219E1*t15*t20)*(rho-old)/24.0-0.1206021294E1*t8*t29-0.3959188841*t15*
t29);
  }
}


inline double dphichemicalPotential ( double rho ,double phi ,double old ) const
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


inline double pressure ( double rho ,double phi ) const
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
-2.0*A_*(t2+(2.0*alpha_-2.0)*t6+(-3.0*alpha_+1.0)*t1+alpha_)/delta_);
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
    t9 = sqrt(0.116654747E11*t2*phi-0.2916368675E11*t2+0.1944245783E11*t1*phi+
950205322.0);
    return(0.1E-3*t9);
  }
}

