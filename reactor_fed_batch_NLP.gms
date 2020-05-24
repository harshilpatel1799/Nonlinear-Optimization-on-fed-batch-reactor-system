$Ontext
Survey of Optimal Production in Fed-Batch Yeast Reactors using Applications of Non-Linear Programming
Created by Harshil Patel on March 5, 2020

objective function is measured in grams
Volume is in liter (L)
mass in grams (g)
concentration is in grams per liter (g/L)

For this instance use tf(final cycle time) as 10 hours- this a data point as a parameter

Use iterative solution to solve when reduced gradient less than internal software tolerance.

$Offtext

*Set the number of interations and save solutions for final evaluation (solution settings)
$if     set n  $set nh %n%
$if not set nh $set nh 75



set nh Number of subintervals / 0*%nh% /;
alias (nh,k);


*System Setting
Scalar    tf   final time           / 10  /
          x1_0 initial value for x1 / 1.0 /
          x2_0 initial value for x2 / 5.0 /
          x3_0 initial value for x3 / 0.0 /
          x4_0 initial value for x4 / 0.0 /
          x5_0 initial value for x5 / 1.0 /
          h ;
          h=tf/%nh%;

Variables x1(nh) biomass concentration
          x2(nh)  glucose concentration
          x3(nh)  total initial protein mass
          x4(nh)  secreted protein concentration
          x5(nh)  culture volume
          u(nh)   control variable
          a1(nh)  feedRate1
          a2(nh)  feedRate2
          a3(nh)  feedRate3
          obj     criterion ;

Equations eobj        criterion definition
          state1(nh)  state equation 1
          state2(nh)  state equation 2
          state3(nh)  state equation 3
          state4(nh)  state equation 4
          state5(nh)  state equation 5
          boundX1
          boundX2
          boundX3
          boundX4
          boundX5
          ea1
          ea2
          ea3;

*Objective function
eobj..obj =e= x4['%nh%']*x5['%nh%'];

*State Equations
state1(nh(k+1))..x1[k+1] =e= x1(k)+(h/2)*( a1(k)*x1(k) - u(k)*x1(k)/x5(k) +a1(k+1)*x1(k+1) - u(k+1)*x1(k+1)/x5(k+1) ) ;

state2(nh(k+1))..x2[k+1] =e= x2(k)+(h/2)*( -7.3*a1(k)*x1(k) - u(k)*(x2(k)-20)/x5(k)-7.3*a1(k+1)*x1(k+1) - u(k+1)*(x2(k+1)-20)/x5(k+1) );

state3(nh(k+1))..x3[k+1] =e= x3(k)+(h/2)*( a2(k)*x1(k) - u(k)*x3(k)/x5(k) +a2(k+1)*x1(k+1) - u(k+1)*x3(k+1)/x5(k+1) );

state4(nh(k+1))..x4[k+1] =e= x4(k)+(h/2)*( a3(k)*(x3(k)-x4(k)) - u(k)*x4(k)/x5(k) +a3(k+1)*(x3(k+1)-x4(k+1)) - u(k+1)*x4(k+1)/x5(k+1) );

state5(nh(k+1))..x5[k+1] =e= x5(k) + (h/2)*( u(k) + u(k+1) );


*Descion Varible Bounds Equations
boundX1(nh(k))..x1(k) =g=0;
boundX2(nh(k))..x2(k) =g=0;
boundX3(nh(k))..x3(k) =g=0;
boundX4(nh(k))..x4(k) =g=0;
boundX5(nh(k))..x5(k) =g=0;


*Co-State Equations
ea1(nh(k)).. a1(k) =e= 21.87*x2(k)/((x2(k)+0.4)*(x2(k)+62.5));
ea2(nh(k)).. a2(k) =e= (x2(k)*exp(-5*x2(k)))/(0.1+x2(k));
ea3(nh(k)).. a3(k) =e= 4.75*a1(k)/(0.12+a1(k));



*Control Feed Rate Bounds for System (Found out how to does this on GAMS website)
u.lo(nh) = 0.0;
u.up(nh) = 5;


*Initial point
x1.l[nh]=1.0;
x2.l[nh]=5.0;
x3.l[nh]=0.0;
x4.l(nh)=0.0;
x5.l(nh)=1.0;
u.l(nh) =0.0;

x1.fx ['0'] = x1_0;
x2.fx ['0'] = x2_0;
x3.fx ['0'] = x3_0;
x4.fx ['0'] = x4_0;
x5.fx ['0'] = x5_0;

Model reactor /all/;
Solve reactor maximizing obj using nlp;

