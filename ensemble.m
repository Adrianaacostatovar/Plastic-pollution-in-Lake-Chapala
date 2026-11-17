function [Ymean, Yall, tspan] = ensemble(Nsim, t_year, variation, mode)
% ENSEMBLE Evaluate sensitivity of plastic pollution trajectories.
%
% Inputs:
%   Nsim      - Number of ensemble simulations
%   t_year    - Simulation duration in years (starting in 1950)
%   variation - Relative uniform parameter variation (e.g., 0.30)
%   mode      - 'loss' or 'transitions'
%
% Outputs:
%   Ymean     - Mean trajectories of eight plastic classes and total (nT x 9)
%   Yall      - Trajectories of all simulations (nT x 9 x Nsim)
%   tspan     - Simulation times in years since 1950
%
% Produces Supplementary Figures S1 and S2 (EPS vector format).

%% Settings
p = parameters();
Gamma_t = plastic_input(p);
tspan = 0:1:t_year;
nT = numel(tspan);
Yall = zeros(nT, 9, Nsim);
rng(2026);

switch lower(mode)
    case 'loss'
        names = {'tau1','tau2','tau3','tau4', ...
                 'tauh1','tauh2','tauh3','tauh4'};
        filename = 'Figure_S1_Loss_Rates.eps';
    case 'transitions'
        names = {'p1','q1','p2','q2','p3','q3','q4', ...
                 'r12','rh12','r13','r23', ...
                 'u12','u13','u23', ...
                 'alpha1','alpha2','alpha3', ...
                 'alpha1h','alpha2h','alpha3h'};
        filename = 'Figure_S2_Transitions.eps';
    otherwise
        error('ensemble:InvalidMode', ...
              'Mode must be ''loss'' or ''transitions''.');
end

%% Ensemble simulations
for i = 1:Nsim
    valid = false;

    while ~valid
        pi = p;
        for k = 1:numel(names)
            name = names{k};
            pi.(name) = p.(name) * (1 + variation * (2*rand - 1));
        end

        if strcmpi(mode, 'transitions')
            % Dependent fractions must satisfy the original constraints.
            pi.rh13 = 1 - (pi.q1 + pi.r12 + pi.rh12 + pi.r13);
            pi.rh23 = 1 - pi.q2 - pi.r23;

            valid = (pi.p1 + pi.q1 <= 1) && ...
                    (pi.p2 + pi.q2 <= 1) && ...
                    (pi.p3 + pi.q3 <= 1) && ...
                    (pi.rh13 >= 0) && (pi.rh23 >= 0) && ...
                    (pi.u12 + pi.u13 <= 1) && (pi.u23 <= 1) && ...
                    all([pi.alpha1,pi.alpha2,pi.alpha3, ...
                         pi.alpha1h,pi.alpha2h,pi.alpha3h] >= 0) && ...
                    all([pi.alpha1,pi.alpha2,pi.alpha3, ...
                         pi.alpha1h,pi.alpha2h,pi.alpha3h] <= 1);
        else
            valid = true;
        end
    end

    [~, Y] = ode45(@(t,y) BaselineModel(t,y,pi,Gamma_t), ...
                   tspan, p.initial_conditions);

    % Consistent order: M1 M2 M3 SM1 SM2 SM3 Mm SMm Total.
    Yall(:,:,i) = [sum(Y(:,1:4),2), sum(Y(:,9:12),2), ...
                   sum(Y(:,17:20),2), sum(Y(:,5:8),2), ...
                   sum(Y(:,13:16),2), sum(Y(:,21:24),2), ...
                   Y(:,25), Y(:,26), sum(Y(:,1:26),2)];
end

Ymean = mean(Yall, 3);

%% Supplementary figure
output_dir = fullfile(fileparts(mfilename('fullpath')), 'Ensemble');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

years = tspan + p.year0;
labels = {'M_1','M_2','M_3','SM_1','SM_2','SM_3', ...
          'M_m','SM_m','Total'};

fig = figure('Name', ['Ensemble - ' mode], ...
             'Color', 'w', 'WindowState', 'maximized');

tl = tiledlayout(fig, 3, 3, ...
    'TileSpacing', 'compact', 'Padding', 'compact');
for j = 1:9
    nexttile;
    hold on;

    % Plot all trajectories at once, avoiding an additional plotting loop.
    curves = reshape(Yall(:,j,:), nT, Nsim);
    h = plot(years, curves, 'Color', [0.75, 0.50, 0.75], ...
             'LineWidth', 0.5);
    set(h, 'HandleVisibility', 'off');

    hRandom = plot(nan, nan, 'Color', [0.75, 0.50, 0.75], ...
                   'LineWidth', 2);
    hMean = plot(years, Ymean(:,j), 'k-.', 'LineWidth', 2);

    grid on;
    box on;
    set(gca, 'FontSize', 18);
    xlim([p.year0, p.year0 + t_year]);
    legend([hRandom, hMean], ...
           {[labels{j}, ' ET'], ...
            [labels{j}, ' MT']}, ...
           'Location', 'southeast', 'FontSize', 12);
end

drawnow;

exportgraphics(tl, fullfile(output_dir, filename), ...
               'ContentType', 'vector');
end
