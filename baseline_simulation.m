function [t, Y, T] = baseline_simulation(t_year)
% BASELINE_SIMULATION Simulate the reference plastic pollution scenario.
%
% Input:
%   t_year - Simulation duration in years, starting in 1950
%
% Outputs:
%   t      - Simulation time vector (years since 1950)
%   Y      - Solutions of the 26 plastic compartments (g/m^3)
%   T      - Baseline class concentrations and distributions in 2026 and 2050
%
% Generates Figure 1 (plastic input), Figure 2 (plastic classes),
% and Table 3 of the manuscript.

%% Baseline model
p = parameters();
Gamma_t = plastic_input(p);
tspan = 0:1:t_year;

[t, Y] = ode45(@(t,y) BaselineModel(t,y,p,Gamma_t), ...
               tspan, p.initial_conditions);

years = t + p.year0;

%% Aggregate the 26 states into eight classes and total concentration
Yplot = [ ...
    sum(Y(:,1:4),2), ...       % M1
    sum(Y(:,9:12),2), ...      % M2
    sum(Y(:,17:20),2), ...     % M3
    sum(Y(:,5:8),2), ...      % SM1
    sum(Y(:,13:16),2), ...    % SM2
    sum(Y(:,21:24),2), ...    % SM3
    Y(:,25), ...              % Mm
    Y(:,26), ...              % SMm
    sum(Y(:,1:26),2) ...      % Total
];

labels = {'M_1','M_2','M_3','SM_1','SM_2','SM_3', ...
          'M_m','SM_m','Total'};

%% Output directory
output_dir = fullfile(fileparts(mfilename('fullpath')), 'Baseline');
if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

%% Figure 1: Historical plastic input
fig = figure('Name','Historical plastic input', ...
             'Color','w','WindowState','maximized');
plot(years, Gamma_t(t), 'LineWidth',2.5);
grid on;
set(gca,'FontSize',30);
xlabel('Year','FontSize',35);
ylabel('\Gamma(t)','FontSize',35);
xlim([p.year0, p.year0 + t_year]);

xline(p.year0 + p.t1, 'k:', '1990', ...
      'LineWidth',2.3,'FontSize',20,'HandleVisibility','off');
xline(p.year0 + p.tSwitch, 'k:', '2015', ...
      'LineWidth',2.3,'FontSize',20,'HandleVisibility','off');
xline(p.year0 + p.Tcontrol, 'k:', '2026', ...
      'LineWidth',2.3,'FontSize',20,'HandleVisibility','off');

exportgraphics(fig, fullfile(output_dir,'gamma.eps'), ...
               'ContentType','vector');

%% Figure 2: Plastic classes and total concentration
fig = figure('Name','Baseline plastic classes', ...
             'Color','w','WindowState','maximized');
tiledlayout(3,3,'TileSpacing','compact','Padding','compact');

for i = 1:9
    nexttile;
    plot(years, Yplot(:,i), 'LineWidth',2.2);
    grid on;
    box on;
    set(gca,'FontSize',20);
    legend(labels{i},'FontSize',20,'Location','southeast');
    xlim([p.year0, p.year0 + t_year]);

    xline(76+p.year0, 'r:', ...
    'LineWidth', 2.5, ...
    'FontSize', 30, ...
    'LabelVerticalAlignment', 'middle', ...
    'HandleVisibility', 'off');

end

exportgraphics(fig, fullfile(output_dir,'2026.eps'), ...
               'ContentType','vector');

%% Table 3: Baseline concentrations and distributions in 2026 and 2050
% Both years must belong to the simulation interval.
reference_years = [2026,2050];
if any(reference_years > p.year0 + t_year)
    error('baseline_simulation:TimeHorizon', ...
          'Simulation duration must include 2026 and 2050 to generate Table 3.');
end

concentrations = zeros(9,2);
distributions = zeros(9,2);
for j = 1:2
    idx = find(years == reference_years(j),1);
    concentrations(:,j) = Yplot(idx,:)';
    distributions(:,j) = 100*concentrations(:,j)/concentrations(end,j);
end

T = table(string({'M1';'M2';'M3';'SM1';'SM2';'SM3';'Mm';'SMm';'Total'}), ...
    concentrations(:,1), distributions(:,1), ...
    concentrations(:,2), distributions(:,2), ...
    'VariableNames', {'Class','Concentration_2026','Distribution_2026', ...
                      'Concentration_2050','Distribution_2050'});

fprintf('\n========== BASELINE SCENARIO (TABLE 3) ==========\n');
disp(T);

end
