function [favorable, unfavorable, T] = management_scenarios(t_year, Ngrid)
% MANAGEMENT_SCENARIOS Evaluate changes in plastic-waste management.
%
% Inputs:
%   t_year - Simulation duration in years, starting in 1950
%   Ngrid  - Number of points along each heatmap axis
%
% Outputs:
%   favorable   - Favorable scenario parameter grids and concentrations
%   unfavorable - Unfavorable scenario parameter grids and concentrations
%   T           - Summary of total concentrations in 2050 (g/m^3)
%
% Generates Figures 3 and 4 of the manuscript in vector EPS format.
% Management changes begin in 2026 and are evaluated in 2050.

%% Parameters and evaluation years
p = parameters();
t_control = p.Tcontrol;   % 2026, measured from 1950
t_eval = 100;             % 2050, measured from 1950

if t_year < t_eval
    error('management_scenarios:TimeHorizon', ...
          't_year must be at least 100 to evaluate 2050.');
end
if Ngrid < 2 || Ngrid ~= round(Ngrid)
    error('management_scenarios:InvalidGrid', ...
          'Ngrid must be an integer greater than 1.');
end

%% Reference state at the start of management changes (2026)
Gamma_base = plastic_input(p);
[~, Ybase] = ode45(@(t,y) BaselineModel(t,y,p,Gamma_base), ...
                   0:1:t_control, p.initial_conditions);
y_control = Ybase(end,:)';

%% Parameter ranges from the manuscript
favorable.epsilon = linspace(p.repBase, 1, Ngrid);
favorable.F_leak = linspace(p.FugaBase, 0, Ngrid);

unfavorable.epsilon = linspace(0, p.repBase, Ngrid);
unfavorable.F_leak = linspace(1, p.FugaBase, Ngrid);

%% Evaluate total concentration in 2050
scenarios = {'favorable', 'unfavorable'};
for k = 1:numel(scenarios)
    scenario = scenarios{k};
    if k == 1
        eps_values = favorable.epsilon;
        leak_values = favorable.F_leak;
    else
        eps_values = unfavorable.epsilon;
        leak_values = unfavorable.F_leak;
    end

    C = zeros(Ngrid, Ngrid);
    fprintf('\nComputing %s management scenarios...\n', scenario);

    for iF = 1:Ngrid
        for iE = 1:Ngrid
            Gamma_t = plastic_input(p, eps_values(iE), leak_values(iF), 'constant');
            [~, Y] = ode45(@(t,y) BaselineModel(t,y,p,Gamma_t), ...
                           [t_control, t_eval], y_control);
            C(iF,iE) = sum(Y(end,1:26));
        end
    end

    if k == 1
        favorable.total = C;
    else
        unfavorable.total = C;
    end
end

%% Reference and extreme scenario comparison
C_reference = favorable.total(1,1);
C_favorable = favorable.total(end,end);
C_unfavorable = unfavorable.total(1,1);

% Both grids must reproduce the same baseline at their common corner.
if abs(C_reference - unfavorable.total(end,end)) > 1e-5
    warning('management_scenarios:BaselineMismatch', ...
            'The reference concentration differs between heatmap grids.');
end

T = table(string({'Reference';'Favorable extreme';'Unfavorable extreme'}), ...
          [C_reference; C_favorable; C_unfavorable], ...
          100 * ([C_reference; C_favorable; C_unfavorable] / C_reference - 1), ...
          'VariableNames', {'Scenario','Concentration_2050','Change_percent'});

fprintf('\n========== MANAGEMENT SCENARIOS (2050) ==========\n');
disp(T);

%% Output directory and figure colors
output_dir = fullfile(fileparts(mfilename('fullpath')), 'Management');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

ncolors = 256;
colors_favorable = [0.05 0.25 0.75; 0.05 0.70 0.80; 0.20 0.65 0.25];
colors_unfavorable = [0.20 0.65 0.25; 0.95 0.80 0.15; 0.80 0.10 0.10];
cmap_favorable = interp1(linspace(0,1,3), colors_favorable, linspace(0,1,ncolors));
cmap_unfavorable = interp1(linspace(0,1,3), colors_unfavorable, linspace(0,1,ncolors));

%% Figure 3: Favorable management changes
fig = figure('Name','Favorable management scenarios', ...
             'Color','w','WindowState','maximized');
imagesc(favorable.epsilon, favorable.F_leak, favorable.total);
set(gca, 'YDir','normal','FontSize',35);
colormap(gca, cmap_favorable);
clim([min(favorable.total(:)), C_reference]);
colorbar;
xlabel('\epsilon','FontSize',45);
ylabel('F_{leak}','FontSize',40);
hold on;
plot(p.repBase, p.FugaBase, 'p', 'MarkerSize',20, ...
     'MarkerFaceColor','k','MarkerEdgeColor','w','LineWidth',2, ...
     'Clipping','off');
plot(1, 0, 's', 'MarkerSize',20, ...
     'MarkerFaceColor','k','MarkerEdgeColor','w','LineWidth',2, ...
     'Clipping','off');
hold off;
exportgraphics(fig, fullfile(output_dir,'Favorable_Management.eps'), ...
               'ContentType','vector');

%% Figure 4: Deteriorating management conditions
fig = figure('Name','Unfavorable management scenarios', ...
             'Color','w','WindowState','maximized');
imagesc(unfavorable.epsilon, unfavorable.F_leak, unfavorable.total);
set(gca, 'YDir','normal','FontSize',35);
colormap(gca, cmap_unfavorable);
clim([unfavorable.total(end,end), max(unfavorable.total(:))]);
colorbar;
xlabel('\epsilon','FontSize',45);
ylabel('F_{leak}','FontSize',40);
hold on;
plot(p.repBase, p.FugaBase, 'p', 'MarkerSize',20, ...
     'MarkerFaceColor','k','MarkerEdgeColor','w','LineWidth',2, ...
     'Clipping','off');
plot(0, 1, 's', 'MarkerSize',20, ...
     'MarkerFaceColor','k','MarkerEdgeColor','w','LineWidth',2, ...
     'Clipping','off');
hold off;
exportgraphics(fig, fullfile(output_dir,'Unfavorable_Management.eps'), ...
               'ContentType','vector');

end
