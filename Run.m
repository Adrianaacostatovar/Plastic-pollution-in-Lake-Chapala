%% LAKE CHAPALA PLASTIC POLLUTION MODEL
%Plastic pollution in Lake Chapala, Mexico: Long-term accumulation, vulnerability to management changes, and removal scenarios
%Adriana Acosta-Tovar, Roberto A. Saenz
%DOI:

clear; clc; close all;

% Simulation settings
year0 = 1950;
t_year = 150;

% Ensemble settings
Nsim = 10; %10000;
variation_tau = 0.30;
variation_transitions = 0.50;

% Heatmap settings
Ngrid = 31;

% Model parameters
T = parameter_table();

% Baseline simulation
[t_baseline, Y_baseline] = baseline_simulation(t_year);

% Management scenarios
management_scenarios(t_year, Ngrid);

%% Plastic removal scenarios

[t_removal, Y_removal, T_removal] = removal_scenarios(t_year);
%% Ensemble simulations
ensemble(Nsim, t_year, variation_tau, 'loss');
ensemble(Nsim, t_year, variation_transitions, 'transitions');

%% Summary tables
%summary_tables();