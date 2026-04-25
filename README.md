# StrainModel

Reduced-order MATLAB model for strain-dependent pressure transmission in a membrane-confined hydrodynamic sensor interface.

## Overview

This repository contains a reduced-order numerical model used to interpret experimentally observed strain-dependent pressure sensitivity in a hydrodynamic sensor module for Autonomous Underwater Vehicles (AUVs).

The current implementation, `sensor_interface_numerical_model_v20.m`, models the membrane-confined cavity/interface system using an iterative force-balance formulation. Version 20 emphasizes physically defensible interpretation by removing active geometry-amplification terms and reporting the model response as an interface-relative transmission normalized by the zero-strain membrane condition.

The model output is reported as:

`S_model(eps) = T_cav(eps) / T_cav(0)`

where:

- `T_cav(eps)` is the local cavity/interface pressure transmission response
- `T_cav(0)` is the zero-strain membrane reference response
- `S_model(eps)` is the interface-relative normalized sensitivity

No empirical bare-port scale factor is applied in the active v20 formulation. The bare-port reference is retained in figures as an experimental reference line.

## Physical Interpretation

The model treats each sensing side as an air-filled cavity coupled to a strained membrane/plate-like interface. A small external pressure perturbation is transmitted through the interface while the internal cavity pressure evolves according to:

- trapped-air compression
- reduced-order bending resistance
- strain-dependent membrane thickness
- strain-activated effective pre-tension
- bounded cavity-volume change
- relaxed iterative pressure updates for stable convergence

The model is intended as a reduced-order physical interpretation tool. It is not a full nonlinear membrane mechanics model, constitutive material model, or spatial fluid-structure interaction simulation.

## v20 Assumptions

The v20 implementation adopts the following assumptions:

1. Axisymmetric reduced-order interface response
2. Effective circular loaded region
3. Small-deflection-style volume estimate regularized by smooth saturation
4. Fixed effective loaded radius independent of installed membrane diameter
5. Strain-dependent thickness update based on incompressibility
6. Pre-strain represented through a smoothed effective pre-tension term
7. Cavity pressure solved iteratively using ideal-gas-style compression
8. Closure terms approximate higher-order mechanics not explicitly resolved in the reduced-order formulation

## Key Changes in v20

Compared with earlier versions, v20:

- removes active geometry-amplification terms from the model response
- uses a fixed effective loaded radius based on the forced loading radius
- removes empirical bare-port scaling from the active model
- interprets model output as interface-relative transmission
- retains material, thickness, pressure-level, and calibration sweeps as diagnostics
- adds a normalized mechanics diagnostic plot
- overlays a physically tuned thickness case on the main model/data comparison figure

## Current Model Formulation

For a given engineering strain `eps`, the model:

1. Interpolates the installed membrane/interface diameter from measured strain states
2. Computes installed geometric radius `a_geom(eps)`
3. Uses a fixed effective loaded radius `a_load = r_forced`
4. Updates membrane thickness using an incompressibility-style relation:

   `h(eps) = h0 / (1 + eps)^2`

5. Computes plate bending stiffness:

   `D_plate(eps) = E*h(eps)^3 / (12*(1 - nu^2))`

6. Computes a smooth pre-tension activation function:

   `phi_pre(eps) = 1 - exp(-(eps / eps_char)^2)`

7. Computes an effective pre-tension closure:

   `Tpre_eff(eps) = kT0 * phi_pre(eps) * E*h(eps)*eps / (1 - nu)`

8. Solves cavity pressure iteratively using a reduced-order cavity/interface force balance
9. Computes local transmission from symmetric positive and negative pressure perturbations
10. Normalizes the response by the zero-strain membrane condition

## Reduced-Order Deflection Basis

The interface response uses a bending-plus-tension denominator of the form:

`w_max ~ q*a^4 / (64*D_plate + 4*Tpre_eff*a^2)`

where:

- `q` is the net membrane load
- `a` is the effective loaded radius
- `D_plate` is the reduced-order bending stiffness
- `Tpre_eff` is the effective pre-tension term

The reduced-order cavity volume change is estimated as:

`dV_raw = pi*a^2*w_max / 3`

To prevent unphysical divergence, the volume change is smoothly bounded using:

`dV = dV_cap * tanh(dV_raw / dV_cap)`

where:

`dV_cap = dV_cap_fraction * V00`

## Key Model Parameters

Important model parameters include:

- `kT0`  
  Effective pre-tension engagement factor. This represents the fraction of the ideal equibiaxial membrane tension realized in the installed interface.

- `eps_char`  
  Characteristic strain controlling smooth activation of pre-tension.

- `t_plate0`  
  Baseline membrane thickness. The nominal hygienic latex value is 0.020 in, or approximately 0.508 mm.

- `E_plate0`  
  Baseline effective Young's modulus.

- `nu_plate`  
  Poisson's ratio.

- `dV_cap_fraction`  
  Saturation limit for cavity-volume change relative to nominal cavity volume.

- `relax`  
  Relaxation factor used in the iterative cavity-pressure update.

- `tolP`  
  Pressure convergence tolerance.

- `c1`, `c2`  
  Optional strain-dependent modulus coefficients. These are currently disabled in the baseline v20 model.

## Current Interface and Cavity Assumptions

The current script uses:

- Interface material: hygienic-grade natural latex rubber
- Nominal membrane thickness: 0.020 in, approximately 0.508 mm
- Baseline effective modulus: `6.0e5 Pa`
- Poisson ratio: `0.49`
- Working fluid: trapped air
- Cavity model: ideal-gas-style pressure-volume update
- Interface representation: reduced-order circular plate/membrane approximation
- Cavity volume: estimated from cylindrical plus spherical-cap geometry per sensing side
- Loading model: uniform loading over an effective circular loaded region

## Experimental Reference Points

The script includes experimental summary values for:

- engineering strain
- normalized sensitivity
- uncertainty/error bars

The current experimental comparison uses:

`eps_data = [0.100, 0.150, 0.200, 0.250]`

`y_data = [1.167, 1.328, 1.270, 1.224]`

`y_err = [0.093, 0.082, 0.045, 0.042]`

These values are used for visual comparison, residual calculation, RMSE reporting, and WRMSE reporting.

## Current v20 Interpretation

The baseline v20 model shows that increasing strain reduces the interface-relative cavity transmission when only fixed loaded radius, thickness evolution, cavity compression, and effective pre-tension are included.

This is useful because it separates two effects:

1. The reduced-order mechanics model predicts how membrane stiffening and pre-tension affect transmission.
2. Experimental data show sensitivity amplification relative to the bare-port reference.

The v20 results therefore support the interpretation that the experimentally observed sensitivity response cannot be explained by pre-tension stiffening alone under the current fixed-radius assumptions. The thickness sensitivity sweep suggests that effective interface thickness and installation-controlled mechanical boundary conditions strongly influence the predicted response.

The main figure includes:

- baseline model using nominal membrane thickness, `h0 = 0.508 mm`
- physically tuned model case using `h0 = 1.50 mm`
- experimental normalized sensitivity measurements
- bare-port reference line

The tuned thickness case is interpreted as a physically controlled parameter sensitivity case, not as a unique fitted material property.

## Parameter Sensitivity Analysis

The v20 implementation includes diagnostic sweeps for:

- effective pre-tension scale factor `kT0`
- membrane modulus `E`
- Poisson's ratio `nu`
- membrane thickness `h0`
- pressure-level dependence

Optional RMSE and WRMSE plots are included but are disabled by default for manuscript-focused runs.

## Example Outputs

The script can generate:

- main normalized sensitivity comparison against experimental data
- interface-relative cavity transmission plot
- normalized mechanics diagnostic plot
- material sensitivity plots for `E`, `nu`, and `h0`
- material fit overlay plots
- optional pressure-level diagnostic plots
- optional RMSE and WRMSE sweep plots

Figures are saved to `results_v20/` in:

- `.png`
- `.fig`
- `.eps`

formats.

## Files

- `sensor_interface_numerical_model_v20.m`  
  Main MATLAB script implementing the iterative reduced-order cavity/interface model, diagnostics, parameter sweeps, and figure export workflow.

## Getting Started

Open MATLAB in this repository folder and run:

```matlab
sensor_interface_numerical_model_v20
```

The script will:

1. Define experimental strain and sensitivity data
2. Evaluate the baseline reduced-order model over a continuous strain range
3. Print model parameters and convergence information
4. Print a comparison table between experimental and modeled values
5. Report RMSE and weighted RMSE
6. Run pressure-level diagnostic calculations
7. Generate enabled plots
8. Save figures to the results_v20/ folder

## Notes

- The model is a reduced-order interpretive tool, not a full nonlinear membrane mechanics solution.
- The pre-tension term is an effective closure term, not a directly measured installed membrane tension.
- The parameter kT0 should be interpreted as an engagement factor that accounts for nonuniform installation strain, clamping compliance, seating effects, and departures from ideal linear-elastic membrane behavior.
- The fixed loaded-radius assumption is intentional in v20 to avoid hidden geometry-amplification fitting.
- The active model response is normalized to the zero-strain membrane condition.
- The bare-port reference is retained for experimental comparison but is not used as an empirical scaling factor in the active v20 formulation.
- The h0 = 1.50 mm case is a physically tuned sensitivity case and should not be interpreted as the actual measured latex thickness.
- Closure terms are included to preserve physical traceability while representing unresolved higher-order effects.

## Author

Tyler J. Inkley

Department of Ocean and Resources Engineering

University of Hawaiʻi at Mānoa