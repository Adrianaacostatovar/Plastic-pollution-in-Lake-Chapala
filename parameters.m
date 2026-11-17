function p = parameters()
% PARAMETERS Define the Lake Chapala plastic pollution model parameters.
%
% Output:
%   p - Structure containing baseline parameters, input data,
%       initial conditions, and reference intervention settings.
%
% Reference: Main manuscript and Supplementary Material, Tables S1-S5.

% Fragmentation residence times (years)
p.mu1 = 230;
p.mu2 = 295;
p.mu3 = 590;
p.muh1 = 2*p.mu1;
p.muh2 = 2*p.mu2;
p.muh3 = 2*p.mu3;
p.n = 4;

% Transfer rates (1/year)
p.lambda_1  = 0.0174;%p.n/p.mu1;
p.lambda_2  = 0.0135;%p.n/p.mu2;
p.lambda_3  = 0.0067;%p.n/p.mu3;
p.lambdah_1 = 0.0086;%p.n/p.muh1;
p.lambdah_2 = 0.00677;%p.n/p.muh2;
p.lambdah_3 = 0.00338;%p.n/p.muh3;

% Sedimentation and fragmentation fractions
p.p1 = 0.01;
p.p2 = 0.01;
p.p3 = 0.001;
p.q1 = 0.01;
p.q2 = 0.001;
p.q3 = 0.001;
p.q4 = 0.001;

% Transitions between plastic classes
p.r12  = 0.8;
p.rh12 = 0.002;
p.r13  = 0.1;
p.rh13 = 1 - (p.q1 + p.r12 + p.rh12 + p.r13);
p.r23  = 0.95;
p.rh23 = 1 - p.q2 - p.r23;
p.u12  = 0.49;
p.u13  = 0.50;
p.u23  = 0.50;

% Microplastic redistribution
p.alpha1  = 0.1;
p.alpha2  = 0.5;
p.alpha3  = 0.1;
p.alpha1h = 0.001;
p.alpha2h = 0.001;
p.alpha3h = 0.001;

% Natural loss rates (1/year)
p.tau1 = 3.8e-3;
p.tau2 = 3.8e-3;
p.tau3 = 3.8e-3;
p.tau4 = 3.8e-3;
p.tauh1 = 3.8e-3;
p.tauh2 = 3.8e-3;
p.tauh3 = 3.8e-3;
p.tauh4 = 3.8e-3;

% Plastic input and hydrological parameters
p.input.RSUjal   = 13703055500;    % g/year
p.input.RSUMich  = 5529885000;     % g/year
p.input.fp       = 0.1299;
p.input.alphaagua = 0.4;
p.input.theta    = 0.2;
p.input.Vol      = 7.6e9;          % m^3
p.input.Qcol     = 2.2769e7;       % m^3/year
p.input.betaNT   = 1;
p.input.rho      = 1e-5;
p.input.Cpcuenca = 3;              % g/m^3
p.input.area     = 1.112e9;        % m^2
p.input.kappa    = 0.15;

% Baseline management conditions
p.repBase = 0.58;
p.FugaBase = 0.07;

% Historical input profile (years since 1950)
p.year0 = 1950;
p.t1 = 40;                    % 1990
p.tSwitch = 65;               % 2015
p.GammaStart = 0.1;           % g/(m^3 year)
p.Tcontrol = 76;              % 2026

% Initial concentrations (g/m^3), 26 state variables
p.initial_conditions = [ ...
    1e-4; 2e-5; 2e-5; 2.5e-5; ... % M1
    1e-4; 2e-5; 2e-5; 2.5e-5; ... % SM1
    5e-5; 2e-5; 2e-5; 5e-5; ...   % M2
    5e-5; 2e-5; 2e-5; 5e-5; ...   % SM2
    5e-5; 2e-5; 2.5e-5; 5e-5; ... % M3
    5e-5; 2e-5; 2.5e-5; 5e-5; ... % SM3
    5e-5; 5e-5];                  % Mm, SMm

% Reference removal parameters (Supplementary Table S5)
p.removal.Cmax_sur = 1.61e-2;  % g/(m^3 year)
p.removal.KM_sur = 0.125;      % g/m^3
p.removal.uM = [0, 10, 20, 30];
p.removal.Cmax_sed = 8.04e-3;  % g/(m^3 year)
p.removal.KM_sed = 3.369;      % g/m^3
p.removal.uSM = [0, 10, 20];

end
