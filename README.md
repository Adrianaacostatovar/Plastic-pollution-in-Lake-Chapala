# Plastic Pollution in Lake Chapala, Mexico

**MATLAB implementation of a compartmental model for long-term plastic accumulation, vulnerability to waste-management changes, and seasonal removal scenarios in Lake Chapala, Mexico.**

**Authors:** Adriana Acosta-Tovar and Roberto A. Saenz.

This repository contains the numerical model and scripts used to investigate the evolution of plastic concentrations in Lake Chapala from 1950 to 2100. The simulations examine a reference trajectory, changes in wastewater treatment and plastic-waste leakage, seasonal surface and sediment removal, and sensitivity to uncertain loss and transfer parameters. The model supports exploratory scenario comparisons; it is **not** a reconstruction of a measured historical time series.

## Model overview

The model consists of **26 concentration compartments** (g/m³): three surface plastic classes (`M1`, `M2`, `M3`), their sediment-associated counterparts (`SM1`, `SM2`, `SM3`), and surface and sediment-associated microplastics (`Mm`, `SMm`). Each of the six non-microplastic classes is divided into four internal subclasses. The model includes plastic inputs, progression and fragmentation between classes, redistribution between surface and sediment, and natural loss processes.

The baseline equations are implemented in `BaselineModel.m`. The removal model, `RemovalModel.m`, extends the 26 equations with three cumulative accounting variables, for a total of **29 state variables**. Removal acts on the `M1`, `M2`, `SM1`, and `SM2` classes.

The time-dependent input `Gamma(t)` is constructed in `plastic_input.m` from waste-management and hydrological parameters. The historical profile includes a reference peak around 1990 and switches to the contemporary baseline input in 2015. Management changes and removal interventions begin in **2026**.

## Repository structure

```text
Lake_Chapala_Plastic_Pollution/
├── Run.m                  # Main script: run all analyses
├── parameters.m           # Shared model and intervention parameters
├── parameter_table.m      # Display parameter tables (Supplementary Tables S1-S4)
├── plastic_input.m        # Historical and constant plastic-input functions
├── BaselineModel.m        # Baseline 26-compartment ODE system
├── baseline_simulation.m  # Reference trajectories and 2026/2050 summary
├── management_scenarios.m # Favorable/unfavorable management heatmaps
├── RemovalModel.m         # Removal-extended ODE system (29 states)
├── removal_scenarios.m    # Seasonal removal comparisons and summary tables
└── ensemble.m             # Loss-rate and transition-parameter ensembles
```

The `Baseline/`, `Management/`, `Removal/`, and `Ensemble/` directories are created when the scripts export figures.

## Requirements

- MATLAB with `ode45`, `tiledlayout`, `exportgraphics`, and figure `WindowState` support. **MATLAB R2020b or later is recommended.**
- No third-party MATLAB code is called by the supplied scripts.

The scripts have not been independently executed as part of this README preparation; check availability of the MATLAB functions in your installed release.

## Getting started

1. Download or clone the repository.
2. Open MATLAB and set the current folder to the directory containing the `.m` files.
3. Run the main script:

   ```matlab
   Run
   ```

`Run.m` displays the parameter tables, runs the baseline model, generates both management heatmaps, evaluates seasonal removal scenarios, and produces two parameter ensembles. Results are displayed in MATLAB and figures are exported as vector EPS files.

**Filename note:** MATLAB function files must have filenames matching their function names. If your downloaded files contain suffixes such as `ensemble(1).m` or `parameters(1).m`, rename them to `ensemble.m` and `parameters.m` (and similarly for the other files) before running `Run`.

## Experiment settings

The main script currently uses:

```matlab
% Simulation horizon
 t_year = 150;                  % 1950-2100

% Ensemble analysis
 Nsim = 10;                     % Small test configuration
 variation_tau = 0.30;          % +/- 30% loss-rate perturbations
 variation_transitions = 0.50;  % +/- 50% transition perturbations

% Management heatmaps
 Ngrid = 31;                    % Points per parameter axis
```

`Run.m` contains a comment indicating `10000` as the larger ensemble setting. Change `Nsim` if you need a larger ensemble; doing so increases computational and figure-export costs. **Do not assume that the current defaults reproduce every manuscript figure unchanged.**

Shared numerical parameters and initial conditions are centralized in `parameters.m`. In particular, the contemporary management baseline is `epsilon = 0.58` (collected wastewater receiving treatment) and `F_leak = 0.07` (plastic-waste leakage fraction). Removal intensity arrays are `uM = [0, 10, 20, 30]` and `uSM = [0, 10, 20]`.

## Analysis workflow

### 1. Baseline trajectory

`baseline_simulation.m` integrates the 26-compartment model with `ode45` from 1950 through the configured final year. It plots the historical plastic input and the aggregated plastic classes, and prints a table of modeled concentrations and class percentages for **2026 and 2050**. The horizon must include 2050 for this table to be generated.

**Outputs:** `Baseline/gamma.eps`, `Baseline/2026.eps`, and a MATLAB table.

### 2. Management scenarios

`management_scenarios.m` examines **two families of changes beginning in 2026**, using a grid of wastewater-treatment fraction `epsilon` and waste-leakage fraction `F_leak`:

- **Favorable:** treatment increases from its baseline value toward 1, while leakage decreases toward 0.
- **Unfavorable:** treatment decreases toward 0, while leakage increases toward 1.

The simulations compare **total plastic concentration in 2050 (g/m³)**. The routine prints the reference concentration and the extreme favorable/unfavorable concentrations with their percentage changes.

**Outputs:** `Management/Favorable_Management.eps` and `Management/Unfavorable_Management.eps`.

### 3. Seasonal removal strategies

`removal_scenarios.m` compares **12 scenarios** combining four surface-removal intensities (`uM`) and three sediment-removal intensities (`uSM`). Interventions start in **2026** and are active for the first **six months of each year**. Saturating removal terms target surface and sediment macroplastic classes through `RemovalModel.m`.

The routine reports modeled concentrations and percentage reductions by class in **2050**, along with cumulative removed concentration and equivalent mass. Equivalent mass is computed using the lake surface area and an assumed **1 m effective water layer**, rather than multiplying concentration by the lake's entire water volume.

**Outputs:** `Removal/Removal_uSM_0.eps`, `Removal/Removal_uSM_10.eps`, `Removal/Removal_uSM_20.eps`, and MATLAB tables labeled Tables 4–7. Removal-trajectory figures extend through **2055** when the simulation horizon permits.

### 4. Parametric robustness ensembles

`ensemble.m` explores the effect of uncertainty in two separate parameter families:

| Mode | Parameters varied | Default perturbation |
|---|---|---:|
| `'loss'` | Natural loss rates `tau1`–`tau4`, `tauh1`–`tauh4` | ±30% |
| `'transitions'` | Fragmentation, sedimentation, transfer, and redistribution parameters | ±50% |

Parameters are sampled using independent uniform relative perturbations, subject to additional feasibility conditions for the transition mode. The script uses `rng(2026)` for reproducibility. It returns all simulated trajectories, the pointwise ensemble mean, and the time vector. The ensemble statistics cover the eight aggregated classes plus total plastic concentration.

**Outputs:** `Ensemble/Figure_S1_Loss_Rates.eps` and `Ensemble/Figure_S2_Transitions.eps`.

These runs are **deterministic ODE simulations under sampled parameter sets**, not stochastic differential-equation simulations.

## Generated outputs

| Directory | Contents | Files |
|---|---|---|
| `Baseline/` | Historical input and reference model trajectories | `gamma.eps`, `2026.eps` |
| `Management/` | Favorable and unfavorable 2050 heatmaps | `Favorable_Management.eps`, `Unfavorable_Management.eps` |
| `Removal/` | Seasonal intervention comparisons | `Removal_uSM_0.eps`, `Removal_uSM_10.eps`, `Removal_uSM_20.eps` |
| `Ensemble/` | Loss and transition sensitivity ensembles | `Figure_S1_Loss_Rates.eps`, `Figure_S2_Transitions.eps` |
| MATLAB Command Window | Parameter and scenario summary tables | Text/tables |

## Reproducibility notes and scope

- Calendar year is represented as `1950 + t`; the default `t_year = 150` covers 1950–2100.
- The baseline system has 26 state variables, with concentrations in g/m³.
- Management and removal interventions are introduced in 2026 and evaluated primarily in 2050.
- The 2026 and 2050 tables describe **modeled concentrations**, not direct observations.
- The historical plastic input and intervention capacities include explicit modeling assumptions. Consult the associated manuscript and supplementary material for parameter sources, scenario definitions, and limitations.
- A large `Nsim` may create large vector EPS files because the ensemble plots include individual simulation trajectories.

## Citation

If you use this code, please cite the associated manuscript:

> Acosta-Tovar, A., and Saenz, R. A. *Plastic pollution in Lake Chapala, Mexico: Long-term accumulation, vulnerability to management changes, and removal scenarios*. Publication details and DOI ().

## License

No license file is included here. Unless the authors publish explicit licensing terms, the repository should not be interpreted as granting general permission to reuse, modify, or redistribute the code. Contact the authors regarding permissions.
