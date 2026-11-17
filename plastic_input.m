function Gamma_t = plastic_input(p, epsilon, F_leak, mode)
% PLASTIC_INPUT Define the time-dependent plastic input to Lake Chapala.
%
% Inputs:
%   p       - Structure containing the model parameters
%   epsilon - Fraction of collected wastewater receiving treatment
%   F_leak  - Fraction of plastic waste leaking into the environment
%   mode    - 'historical' (default) or 'constant'
%
% Output:
%   Gamma_t - Function handle for the plastic input (g/m^3/year)
%
% The historical profile follows Supplementary Eq. (S8).
% The contemporary input follows Supplementary Eq. (S7).

if nargin < 2 || isempty(epsilon)
    epsilon = p.repBase;
end
if nargin < 3 || isempty(F_leak)
    F_leak = p.FugaBase;
end
if nargin < 4 || isempty(mode)
    mode = 'historical';
end

%% Contemporary input (Supplementary Eq. S7)
in = p.input;
MSW = in.RSUjal + in.RSUMich;

Qeff = epsilon * in.Qcol;
QUT  = (1 - epsilon) * in.Qcol;
Qlake = Qeff + in.betaNT * QUT;

if Qlake <= 0
    error('plastic_input:InvalidFlow', 'Effective wastewater flow must be positive.');
end

Imun_dir = in.theta * MSW * in.fp * F_leak * in.alphaagua;
Imun_W   = (1 - in.theta) * MSW * in.fp * F_leak * in.alphaagua;

Cp = in.Cpcuenca + Imun_W / Qlake;
Ct = in.rho * Cp;
IW = Ct * Qeff + Cp * in.betaNT * QUT;

heff = 1; % Effective normalization depth (m), Supplementary Section S5.3
Gamma_base = in.kappa * (IW + Imun_dir) / (in.area * heff);

%% Plastic input profile (Supplementary Eq. S8)
switch lower(mode)
    case 'historical'
        peak = 7 * p.GammaStart;
        shift = peak / (1 + (p.tSwitch - p.t1)^2 / 7);
        Gamma_t = @(t) Gamma_base + ...
            (t < p.tSwitch) .* (peak ./ (1 + (t - p.t1).^2 / 7) - shift);

    case 'constant'
        Gamma_t = @(t) Gamma_base + zeros(size(t));

    otherwise
        error('plastic_input:InvalidMode', ...
            'Mode must be ''historical'' or ''constant''.');
end

end
