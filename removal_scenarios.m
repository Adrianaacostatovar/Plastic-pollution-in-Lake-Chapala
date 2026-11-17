function [t_removal, Y_removal, T_removal] = removal_scenarios(t_year)
% REMOVAL_SCENARIOS Simulate seasonal plastic-removal interventions.
%
% Input:
%   t_year    - Simulation duration in years since 1950 (at least 100)
%
% Outputs:
%   t_removal - Simulation time vector (years since 1950, from 2026 onward)
%   Y_removal - Structure containing the 12 scenario trajectories
%   T_removal - Structure containing manuscript Tables 4-7
%
% Generates Figure 5 and Tables 4-7 of the manuscript.
% Interventions start in 2026 and are active for six months per year.

%% Parameters and evaluation dates
p = parameters();
Gamma_t = plastic_input(p);
t0 = p.Tcontrol;               % 2026 (years since 1950)
teval = 2050 - p.year0;        % Evaluation year
if t_year < teval
    error('removal_scenarios:TimeHorizon', ...
        't_year must include 2050 (t_year >= %g).', teval);
end
tend = min(t_year, 2055 - p.year0); % Figure 5 extends through 2055

uM_levels  = p.removal.uM;
uSM_levels = p.removal.uSM;

%% Common initial state at the beginning of intervention
[~, Ypre] = ode45(@(t,y) BaselineModel(t,y,p,Gamma_t), ...
    [0 t0], p.initial_conditions);
y0 = [Ypre(end,:)'; 0; 0; 0];

%% Solve the 12 scenarios once (piecewise-constant seasonal activation)
nM = numel(uM_levels);
nSM = numel(uSM_levels);
Y_removal = struct();
C2050 = zeros(11,nM,nSM);
removed2050 = zeros(nM,nSM);
t_removal = [];

for j = 1:nSM
    for i = 1:nM
        uM = uM_levels(i);
        uSM = uSM_levels(j);
        fprintf('Removal scenario: uM = %g, uSM = %g\n',uM,uSM);

        % Solve each active/inactive half-year separately so the solver
        % never steps across a discontinuity of the removal schedule.
        y_current = y0;
        tt = [];
        YY = [];
        nPeriods = ceil(2*(tend-t0));
        for k = 1:nPeriods
            ta = t0 + (k-1)/2;
            tb = min(t0 + k/2, tend);
            if ta >= tb
                break;
            end
            active = mod(k,2) == 1;
            uSurface = uM * double(active);
            uSediment = uSM * double(active);
            uM_t  = @(t) uSurface + zeros(size(t));
            uSM_t = @(t) uSediment + zeros(size(t));

            tseg = linspace(ta,tb,3);
            [tk,Yk] = ode45(@(t,y) RemovalModel( ...
                t,y,p,Gamma_t,uM_t,uSM_t), tseg, y_current);
            if k > 1
                tk = tk(2:end);
                Yk = Yk(2:end,:);
            end
            tt = [tt; tk]; %#ok<AGROW>
            YY = [YY; Yk]; %#ok<AGROW>
            y_current = Yk(end,:)';
        end

        if isempty(t_removal)
            t_removal = tt;
        end
        Y_removal(i,j).uM = uM;
        Y_removal(i,j).uSM = uSM;
        Y_removal(i,j).Y = YY;

        [~, idx2050] = min(abs(tt-teval));
        if abs(tt(idx2050)-teval)>1e-9
            error('removal_scenarios:EvaluationTime', ...
                'The numerical solution does not include 2050.');
        end
        C2050(:,i,j) = class_concentrations(YY(idx2050,:));
        removed2050(i,j) = sum(YY(idx2050,28:29));
    end
end

%% No-control reference for percentage reductions
idxM0 = find(uM_levels==0,1);
idxSM0 = find(uSM_levels==0,1);
if isempty(idxM0) || isempty(idxSM0)
    error('removal_scenarios:ReferenceMissing', ...
        'Both intensity arrays must contain 0.');
end
Cref = C2050(:,idxM0,idxSM0);

%% Tables 4-6: concentrations and reductions in 2050
row_names = ["M1";"M2";"M3";"SM1";"SM2";"SM3"; ...
             "Mm";"SMm";"Surface";"Sediment";"Total"];
T_removal = struct();
for j = 1:nSM
    if uSM_levels(j)==0
        selected = find(uM_levels~=0); % Table 4 omits the no-removal column
        tableNumber = 4;
    elseif uSM_levels(j)==10
        selected = 1:nM;
        tableNumber = 5;
    elseif uSM_levels(j)==20
        selected = 1:nM;
        tableNumber = 6;
    else
        continue;
    end

    tab = table(row_names,'VariableNames',{'Class'});
    for i = selected
        conc = C2050(:,i,j);
        reduction = 100*(1-conc./Cref);
        tab.(sprintf('uM_%g_Conc',uM_levels(i))) = round(conc, 3);
        tab.(sprintf('uM_%g_Reduction',uM_levels(i))) = round(reduction, 3);
    end
    field = sprintf('Table%d',tableNumber);
    T_removal.(field) = tab;
    fprintf('\n%s: uSM = %g (concentration g/m^3; reduction %%)\n', ...
        field,uSM_levels(j));
    disp(tab);
end

%% Table 7: Cumulative removal between 2026 and 2050
% Equivalent mass uses a one-metre effective water layer, not lake volume.
tonnes_per_gm3 = p.input.area / 1e6;
tab = table(uM_levels(:),'VariableNames',{'uM'});
for j = 1:nSM
    tag = sprintf('uSM_%g',uSM_levels(j));
    tab.([tag '_Conc']) = round(removed2050(:,j),3);
    tab.([tag '_Tonnes']) = round(removed2050(:,j)*tonnes_per_gm3,3);
end
T_removal.Table7 = tab;
fprintf('\nTable7: cumulative removed concentration and equivalent mass, 2026-2050\n');
disp(tab);

%% Removal figures
% One figure for each sediment-removal intensity
Area_Chapala = p.input.area;
h_eff = 1;  % Effective water layer (m)

gm3_to_tonnes = Area_Chapala * h_eff / 1e6;
output_dir = fullfile(fileparts(mfilename('fullpath')), 'Removal');

if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

% Original manuscript colors
colores_control = [
    0.00 0.45 0.74
    0.85 0.33 0.10
    0.93 0.69 0.13
    0.49 0.18 0.56
    0.47 0.67 0.19
];

years = p.year0 + t_removal;

for j = 1:nSM

    fig = figure('Color','w','WindowState','maximized');

    tl = tiledlayout(fig,2,1, ...
        'TileSpacing','compact', ...
        'Padding','compact');

    %% Remaining plastic mass
    ax1 = nexttile(tl);
    hold(ax1,'on');

    for i = 1:nM

        Y = Y_removal(i,j).Y;

        % Total remaining concentration
        total_conc = sum(Y(:,1:26),2);

        % Convert concentration to equivalent mass
        total_tonnes = total_conc*gm3_to_tonnes;

        plot(ax1,years,total_tonnes, ...
            'Color',colores_control(i,:), ...
            'LineWidth',2.2, ...
            'DisplayName',sprintf('u_M = %g',uM_levels(i)));

    end

    xline(ax1,2026,'k:', ...
        'LineWidth',2.3, ...
        'HandleVisibility','off');

    grid(ax1,'on');
    box(ax1,'on');

    set(ax1,'FontSize',20);

    xlabel(ax1,'Year','FontSize',25);
    ylabel(ax1,'tonnes','FontSize',25);

    title(ax1,sprintf('Sediment removal: u_{SM} = %g', ...
        uSM_levels(j)), ...
        'FontSize',35, ...
        'Interpreter','tex');

    legend(ax1,'Location','best','FontSize',10);

    xlim(ax1,[2026 2055]);

    %% Cumulative removed plastic mass
    ax2 = nexttile(tl);
    hold(ax2,'on');

    for i = 1:nM

        Y = Y_removal(i,j).Y;

        % Cumulative surface and sediment removal
        removed_conc = Y(:,28) + Y(:,29);

        % Convert concentration to equivalent mass
        removed_tonnes = removed_conc*gm3_to_tonnes;

        plot(ax2,years,removed_tonnes, ...
            'Color',colores_control(i,:), ...
            'LineWidth',2.2, ...
            'DisplayName',sprintf('u_M = %g',uM_levels(i)));

    end

    xline(ax2,2050,'k:', ...
        'LineWidth',4, ...
        'HandleVisibility','off');

    grid(ax2,'on');
    box(ax2,'on');

    set(ax2,'FontSize',20);

    xlabel(ax2,'Year','FontSize',25);
    ylabel(ax2,'Removed (tonnes)','FontSize',25);

    legend(ax2,'Location','best','FontSize',10);

    xlim(ax2,[2026 2055]);

    %% Export figure
    filename = sprintf('Removal_uSM_%g.eps',uSM_levels(j));

    drawnow;

    exportgraphics(tl,fullfile(output_dir,filename), ...
        'ContentType','vector');

end

end

function C = class_concentrations(y)
% CLASS_CONCENTRATIONS Aggregate the 26 states using manuscript order.
y = y(:);
C = [sum(y(1:4));sum(y(9:12));sum(y(17:20)); ...
     sum(y(5:8));sum(y(13:16));sum(y(21:24)); ...
     y(25);y(26); ...
     sum(y([1:4,9:12,17:20,25])); ...
     sum(y([5:8,13:16,21:24,26])); ...
     sum(y(1:26))];
end
