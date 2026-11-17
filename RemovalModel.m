function dydt = RemovalModel(t, y, p, Gamma_t, uM_t, uSM_t)
% REMOVALMODEL Define the plastic pollution model with removal interventions.
%
% Inputs:
%   t       - Simulation time (years since 1950)
%   y       - State vector (26 compartments and 3 cumulative variables)
%   p       - Structure containing the model and removal parameters
%   Gamma_t - Function handle for the time-dependent plastic input
%   uM_t    - Function handle for surface removal intensity
%   uSM_t   - Function handle for sediment removal intensity
%
% Output:
%   dydt    - Time derivatives of the 29 state variables
%
% State variables:
%   1:26 - Plastic concentrations, in the same order as BaselineModel
%   27   - Cumulative plastic input into sediment (g/m^3)
%   28   - Cumulative plastic removed from sediment (g/m^3)
%   29   - Cumulative plastic removed from surface (g/m^3)
%
% Removal acts only on M1, M2, SM1, and SM2, as in the manuscript.

%% Baseline dynamics (original 26 equations)
dydt_base = BaselineModel(t, y(1:26), p, Gamma_t);

%% Targeted plastic concentrations
surface_idx  = [1:4, 9:12];
sediment_idx = [5:8, 13:16];

M  = sum(y(surface_idx));
SM = sum(y(sediment_idx));

%% Saturating seasonal removal rates (Supplementary Eqs. S11-S12)
R_sur = uM_t(t)  * p.removal.Cmax_sur * M  / (p.removal.KM_sur + M);
R_sed = uSM_t(t) * p.removal.Cmax_sed * SM / (p.removal.KM_sed + SM);

% Distribute the removal in proportion to the available concentration.
% These expressions are algebraically equivalent to the original script.
dydt_base(surface_idx) = dydt_base(surface_idx) ...
    - uM_t(t) * p.removal.Cmax_sur .* y(surface_idx) ...
    / (p.removal.KM_sur + M);

dydt_base(sediment_idx) = dydt_base(sediment_idx) ...
    - uSM_t(t) * p.removal.Cmax_sed .* y(sediment_idx) ...
    / (p.removal.KM_sed + SM);

%% Cumulative sediment input (same expression as original model)
M1  = sum(y(1:4));
M2  = sum(y(9:12));
M3a = sum(y(17:19));

Gamma = Gamma_t(t);
Gamma_sed = (54.24 + 13.72 + 27.44 + 1.47)/100 * Gamma;

JSed = Gamma_sed ...
    + p.p1*p.lambda_1*sum(y(1:3)) ...
    + p.rh12*p.lambda_1*y(4) ...
    + p.rh13*p.lambda_1*y(4) ...
    + (1-p.alpha1)*p.q1*p.lambda_1*M1 ...
    + p.p2*p.lambda_2*sum(y(9:11)) ...
    + p.rh23*p.lambda_2*y(12) ...
    + (1-p.alpha2)*p.q2*p.lambda_2*M2 ...
    + p.p3*p.lambda_3*M3a ...
    + (1-p.alpha3)*p.q3*p.lambda_3*M3a ...
    + (1-p.alpha3)*p.lambda_3*y(20) ...
    + p.q4*y(25);

%% Full model: 26 concentrations + 3 cumulative quantities
dydt = [dydt_base; JSed; R_sed; R_sur];

end
