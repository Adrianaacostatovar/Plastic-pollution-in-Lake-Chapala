function dydt = BaselineModel(t, y, p, Gamma_t)
% BASELINEMODEL Define the 26-compartment plastic pollution model.
%
% Inputs:
%   t       - Simulation time
%   y       - State vector of the 26 plastic compartments
%   p       - Structure containing the model parameters
%   Gamma_t - Function handle for the time-dependent plastic input
%
% Output:
%   dydt    - Time derivatives of the 26 compartments

% Model parameters
lambda_1  = p.lambda_1;
lambda_2  = p.lambda_2;
lambda_3  = p.lambda_3;
lambdah_1 = p.lambdah_1;
lambdah_2 = p.lambdah_2;
lambdah_3 = p.lambdah_3;

q1 = p.q1;
q2 = p.q2;
q3 = p.q3;
q4 = p.q4;

p1 = p.p1;
p2 = p.p2;
p3 = p.p3;

r12  = p.r12;
rh12 = p.rh12;
r13  = p.r13;
rh13 = p.rh13;
r23  = p.r23;
rh23 = p.rh23;

u12 = p.u12;
u13 = p.u13;
u23 = p.u23;

alpha1  = p.alpha1;
alpha2  = p.alpha2;
alpha3  = p.alpha3;
alpha1h = p.alpha1h;
alpha2h = p.alpha2h;
alpha3h = p.alpha3h;

tau1  = p.tau1;
tau2  = p.tau2;
tau3  = p.tau3;
tau4  = p.tau4;
tauh1 = p.tauh1;
tauh2 = p.tauh2;
tauh3 = p.tauh3;
tauh4 = p.tauh4;

% Time-dependent plastic input
Gamma = Gamma_t(t);

Gamma1  = 2.26/100  * Gamma;
Gamma2  = 0.28/100  * Gamma;
Gamma3  = 0.56/100  * Gamma;
Gamma4  = 0.03/100  * Gamma;

Gammah1 = 54.24/100 * Gamma;
Gammah2 = 13.72/100 * Gamma;
Gammah3 = 27.44/100 * Gamma;
Gammah4 = 1.47/100  * Gamma;

% State variables
% Surface macroplastics M1
M_11 = y(1);
M_12 = y(2);
M_13 = y(3);
M_14 = y(4);

% Sediment-associated macroplastics SM1
Mh_11 = y(5);
Mh_12 = y(6);
Mh_13 = y(7);
Mh_14 = y(8);

% Surface macroplastics M2
M_21 = y(9);
M_22 = y(10);
M_23 = y(11);
M_24 = y(12);

% Sediment-associated macroplastics SM2
Mh_21 = y(13);
Mh_22 = y(14);
Mh_23 = y(15);
Mh_24 = y(16);

% Surface mesoplastics M3
M_31 = y(17);
M_32 = y(18);
M_33 = y(19);
M_34 = y(20);

% Sediment-associated mesoplastics SM3
Mh_31 = y(21);
Mh_32 = y(22);
Mh_33 = y(23);
Mh_34 = y(24);

% Microplastics
M_m  = y(25);
Mh_m = y(26);

% Differential equations
    %M1
    DM11=Gamma1-lambda_1*M_11-tau1*M_11;
  	DM12= (1-p1-q1)*lambda_1*M_11 - lambda_1*M_12-tau1*M_12;
    DM13= (1-p1-q1)*lambda_1*M_12 - lambda_1*M_13-tau1*M_13;
    DM14=  (1-p1-q1)*lambda_1*M_13 - lambda_1*M_14-tau1*M_14;
	%Mh1
    DMh11=Gammah1-lambdah_1*Mh_11-tauh1*Mh_11;
    DMh12=p1*lambda_1*M_11+lambdah_1*Mh_11 - lambdah_1*Mh_12-tauh1*Mh_12;
    DMh13=p1*lambda_1*M_12+lambdah_1*Mh_12 - lambdah_1*Mh_13-tauh1*Mh_13;
	DMh14=p1*lambda_1*M_13+lambdah_1*Mh_13 - lambdah_1*Mh_14-tauh1*Mh_14;
    %M2
    DM21=Gamma2+r12*lambda_1*M_14-lambda_2*M_21-tau2*M_21;
    DM22= (1-p2-q2)*lambda_2*M_21 - lambda_2*M_22-tau2*M_22;
    DM23= (1-p2-q2)*lambda_2*M_22 - lambda_2*M_23-tau2*M_23;
    DM24=  (1-p2-q2)*lambda_2*M_23 - lambda_2*M_24-tau2*M_24;
	%Mh2
    DMh21=Gammah2 + (u12)*lambdah_1*Mh_14+rh12*lambda_1*M_14-lambdah_2*Mh_21-tauh2*Mh_21;
    DMh22=p2*lambda_2*M_21+lambdah_2*Mh_21 - lambdah_2*Mh_22-tauh2*Mh_22;
    DMh23=p2*lambda_2*M_22+lambdah_2*Mh_22 - lambdah_2*Mh_23-tauh2*Mh_23;
	DMh24=p2*lambda_2*M_23+lambdah_2*Mh_23 - lambdah_2*Mh_24-tauh2*Mh_24;
    %M3
    DM31=Gamma3+r23*lambda_2*M_24+r13*lambda_1*M_14-lambda_3*M_31-tau3*M_31;
   	DM32= (1-p3-q3)*lambda_3*M_31 - lambda_3*M_32-tau3*M_32;
    DM33= (1-p3-q3)*lambda_3*M_32 - lambda_3*M_33-tau3*M_33;
    DM34=  (1-p3-q3)*lambda_3*M_33 - lambda_3*M_34-tau3*M_34;
    %Mh3
    DMh31=Gammah3+u23*lambdah_2*Mh_24+(u13)*lambdah_1*Mh_14+(1-q2-r23)*lambda_2*M_24+ (1-q1-r12-rh12-r13)*lambda_1*M_14-lambdah_3*Mh_31-tauh3*Mh_31;
    DMh32=p3*lambda_3*M_31+lambdah_3*Mh_31 - lambdah_3*Mh_32-tauh3*Mh_32;
    DMh33=p3*lambda_3*M_32+lambdah_3*Mh_32 - lambdah_3*Mh_33-tauh3*Mh_33;
	DMh34=p3*lambda_3*M_33+lambdah_3*Mh_33 - lambdah_3*Mh_34-tauh3*Mh_34;
    %micro
    DMm= Gamma4-(tau4+q4)*M_m+alpha1*q1*lambda_1*(M_11 + M_12 + M_13 + M_14)+alpha2*q2*lambda_2*(M_21 + M_22 + M_23 + M_24)+alpha3*q3*lambda_3*(M_31 + M_32 + M_33)+alpha3*lambda_3*M_34 +alpha3h*lambdah_3*Mh_34 +alpha1h*(1-u12-u13)*lambdah_1*Mh_14+alpha2h*(1-u23)*lambdah_2*Mh_24;
    DMhm= Gammah4+q4*M_m-tauh4*Mh_m+(1-alpha1)*q1*lambda_1*(M_11 + M_12 + M_13 + M_14)+(1-alpha2)*q2*lambda_2*(M_21 + M_22 + M_23 + M_24)+(1-alpha3)*q3*lambda_3*(M_31 + M_32 + M_33)+(1-alpha3)*lambda_3*M_34 +(1-alpha3h)*lambdah_3*Mh_34+(1-alpha1h)*(1-u12-u13)*lambdah_1*Mh_14 +(1-alpha2h)*(1-u23)*lambdah_2*Mh_24;
    JSed = ...
    Gammah1 + Gammah2 + Gammah3 + Gammah4 ...
    + p1*lambda_1*(M_11 + M_12 + M_13) ...
    + rh12*lambda_1*M_14 ...
    + rh13*lambda_1*M_14 ...
    + (1-alpha1)*q1*lambda_1*(M_11 + M_12 + M_13 + M_14) ...
    + p2*lambda_2*(M_21 + M_22 + M_23) ...
    + rh23*lambda_2*M_24 ...
    + (1-alpha2)*q2*lambda_2*(M_21 + M_22 + M_23 + M_24) ...
    + p3*lambda_3*(M_31 + M_32 + M_33) ...
    + (1-alpha3)*q3*lambda_3*(M_31 + M_32 + M_33) ...
    + (1-alpha3)*lambda_3*M_34 ...
    + q4*M_m;
    dydt = [DM11;DM12;DM13;DM14;DMh11;DMh12;DMh13;DMh14;DM21;DM22;DM23;DM24;DMh21;DMh22;DMh23;DMh24; DM31;DM32;DM33;DM34;DMh31;DMh32;DMh33;DMh34;DMm;DMhm;];
end