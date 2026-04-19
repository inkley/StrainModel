# StrainModel

Reduced-order MATLAB model for strain-dependent pressure transmission in a membrane-confined hydrodynamic sensor interface.

## Overview

This repository contains a reduced-order numerical model used to interpret experimentally observed strain-dependent pressure sensitivity in a hydrodynamic sensor module for Autonomous Underwater Vehicles (AUVs).

The current implementation, `sensor_interface_numerical_model_v18.m`, is an iterative reduced-order force-balance model for the cavity/interface system. It is developed to assess whether strain-dependent effective pre-tension, cavity compression, and modest geometry evolution are sufficient to explain the observed normalized sensitivity response.

Rather than imposing a standalone empirical coupling curve, the model computes normalized sensitivity directly from the cavity/interface transmission response. The formulation combines an iterative ideal-gas cavity-pressure update with a reduced-order plate/membrane-style interface response, using measured strain-dependent interface diameter states as geometric input.

The model output is reported as:

`S_norm(eps) = S0_interface_bare * T_cav_rel0(eps)`

where:

- `T_cav_rel0(eps) = T_cav(eps) / T_cav(0)` is the cavity/interface transmission response normalized by its zero-strain value
- `S0_interface_bare` is a baseline multiplicative calibration that maps interface-relative transmission to the experimentally referenced bare-port sensitivity scale

In the present formulation, strain affects the model through:

- installed interface diameter
- effective loaded radius
- interface thickness kinematics through pre-stretch
- plate bending stiffness
- a smoothed effective pre-tension term
- cavity compression feedback through the iterative pressure solve

## Physical Interpretation

The model treats each sensing side as an air-filled cavity coupled to a strained membrane/plate-like interface. A small applied surface pressure perturbation is transmitted through the interface, while internal cavity pressure evolves iteratively according to:

- trapped-air compression
- strain-dependent interface geometry
- reduced-order bending resistance
- strain-dependent effective pre-tension
- bounded cavity-volume change
- relaxed iterative pressure updates for stable convergence

The present implementation is intended as a reduced-order physical interpretation tool that preserves physical traceability without attempting to fully resolve nonlinear membrane mechanics, large-deflection membrane tension, or a full spatial fluid-structure interaction field.

## v18 Assumptions

The current v18 implementation adopts the following assumptions:

1. axisymmetric reduced-order response  
2. effective circular loaded region  
3. small-deflection-style volume estimate regularized by saturation  
4. pre-strain represented through a smoothed effective in-plane tension term  
5. cavity response solved iteratively via ideal-gas compression  
6. closure terms approximate higher-order mechanics not explicitly included in the reduced-order formulation  

## Current Model Formulation

For a given engineering strain `eps`, the model:

1. Interpolates the installed interface diameter from measured strain states
2. Computes the installed geometric radius `a_geom(eps)`
3. Computes an effective loaded radius `a_load(eps)` using a weak, capped geometry-based correction relative to the forced loading radius
4. Computes thickness evolution through  
   `t(eps) = t0 / (1 + eps)^2`
5. Computes effective modulus using either a constant value or the optional form  
   `E_eff(eps) = E_plate0 * (1 + c1*eps + c2*eps^2)`
6. Computes plate bending stiffness  
   `D_plate(eps) = E_eff * t^3 / (12 * (1 - nu^2))`
7. Computes a smoothed pre-tension activation term  
   `phi_pre(eps) = 1 - exp(-(eps / eps_char)^2)`
8. Computes an effective pre-tension closure  
   `Tpre_eff(eps) = kT0 * phi_pre * E_eff * t * eps / (1 - nu)`
9. Solves cavity pressure iteratively using a reduced-order cavity/interface force balance
10. Computes normalized transmission from symmetric positive and negative pressure perturbations
11. Scales the normalized cavity/interface transmission to the bare-port experimental reference using `S0_interface_bare`

## Reduced-Order Deflection Basis

The interface response is based on a circular bending-plus-tension denominator of the form:

`w_max ~ q*a^4 / (64*D + 4*Tpre_eff*a^2)`

where:

- `q` is the net membrane load
- `a` is the effective loaded radius
- `D` is the reduced-order bending stiffness
- `Tpre_eff` is the effective pre-tension closure term

The code then estimates reduced-order cavity volume change from:

`dV_raw = pi*a^2*w_max / 3`

To maintain numerical stability and prevent unphysical divergence, the raw volume change is smoothly bounded using:

`dV = dV_cap * tanh(dV_raw / dV_cap)`

where:

`dV_cap = dV_cap_fraction * V00`

This formulation preserves a physically interpretable response structure while regularizing the cavity update for stable iterative convergence.

## Key Model Parameters

Important model parameters include:

- `kT0`  
  Scaling coefficient for the effective pre-tension closure term

- `eps_char`  
  Characteristic strain controlling the smooth activation of pre-tension through `phi_pre`

- `S0_interface_bare`  
  Baseline multiplicative calibration that maps normalized cavity/interface transmission to the experimentally referenced bare-port sensitivity scale

- `alpha_load`  
  Weighting parameter controlling how the installed geometric radius modifies the effective loaded radius

- `a_load_cap_fraction`  
  Cap on the effective loaded-radius correction relative to the forced loading radius

- `dV_cap_fraction`  
  Cap on cavity volume change relative to nominal cavity volume

- `relax`  
  Relaxation factor used to improve solver stability during the iterative cavity-pressure update

- `tolP`  
  Pressure convergence tolerance for the iterative cavity solve

- `c1`, `c2`  
  Optional linear and quadratic coefficients for a strain-dependent effective modulus model

## Current Interpretation of v18 Results

The v18 implementation shifts the reduced-order interpretation away from geometry-dominated gain fitting and toward a mechanics-informed cavity/interface force-balance model with effective pre-tension and bounded geometry evolution.

Current results suggest that:

- the model is numerically stable and convergent across the full strain sweep
- a strain-activated effective pre-tension term provides a plausible first-order explanation for the reduction in normalized transmission with increasing strain
- modest loaded-radius evolution contributes secondary geometry-dependent correction
- cavity compression remains an essential part of the coupled response through the iterative pressure update
- the model should be interpreted as a reduced-order mechanics-informed framework rather than a fully predictive constitutive membrane model
- the closure terms are intended to represent unresolved higher-order effects, not to imply unique parameter identification

In other words, the present v18 model supports the interpretation that increasing pre-strain alters interface compliance and transmitted pressure response through coupled tension, geometry, and cavity-compression effects within a reduced-order force-balance framework.

## Current Interface and Cavity Assumptions

The current script uses:

- **Interface material:** hygienic latex
- **Membrane/interface thickness:** 0.020 in
- **Nominal effective modulus:** 6.0e5 Pa
- **Poisson ratio:** 0.49
- **Working fluid:** trapped air
- **Cavity model:** isothermal ideal-gas-style pressure-volume update
- **Interface representation:** reduced-order circular plate/membrane-style approximation
- **Cavity volume:** estimated from a cylindrical plus spherical-cap geometry per sensing side
- **Load model:** reduced-order uniform loading over an effective circular loaded region

## Experimental Reference Points

The script includes experimental calibration summary values for:

- engineering strain
- normalized sensitivity
- approximate uncertainty / error bars

These are used for direct visual comparison with model predictions and for reporting:

- residuals
- RMSE
- weighted RMSE

## Parameter Sensitivity Analysis

The current v18 implementation includes focused parameter sweeps around the current working region.

Implemented sweeps include:

- `kT0` sweep
- `S0_interface_bare` sweep
- `alpha_load` sweep
- joint contour sweep of `kT0` and `S0_interface_bare`
- RMSE and WRMSE reporting for all sweep cases
- baseline versus minimum-error locations in the joint sweep space

The optional strain-dependent modulus coefficients `c1` and `c2` remain available in the script but are currently disabled in the working v18 formulation.

## Example Outputs

The script generates:

- baseline normalized sensitivity comparison against experimental data
- cavity/interface transmission response normalized to the zero-strain interface state
- mechanics diagnostic plot showing `E_eff`, `D_plate`, `Tpre_eff`, and `phi_pre`
- geometry diagnostic plot showing `a_geom`, `a_load`, and `r_forced`
- joint RMSE contour plot over `kT0` and `S0_interface_bare`
- joint WRMSE contour plot over `kT0` and `S0_interface_bare`

Figures are saved to the results directory in:

- `.png`
- `.fig`
- `.eps`

formats.

## Files

- `sensor_interface_numerical_model_v18.m`  
  Main MATLAB script implementing the iterative reduced-order cavity/interface force-balance model, diagnostics, parameter sweeps, and figure export workflow

## Getting Started

Open MATLAB in this repository folder and run:

```matlab
sensor_interface_numerical_model_v18
```

The script will:

1.	Define experimental strain and sensitivity data
2.	Evaluate the baseline force-balance model over a continuous strain range
3.	Print a comparison table between experimental and modeled values
4.	Report RMSE and weighted RMSE
5.	Generate baseline diagnostic plots
6.	Run focused parameter sweeps
7.	Generate joint kT0 / S0_interface_bare contour plots
8.	Save figures to the results_v18/ folder as .png, .fig, and .eps files

## Notes

- The present model is a reduced-order interpretive tool, not a full nonlinear membrane mechanics solution
- The interface mechanics are represented using a reduced-order bending-plus-effective-tension formulation
- The pre-tension contribution is introduced as a smoothed effective closure term, not as a fully derived constitutive membrane tension law
- The cavity-pressure solution assumes an air-filled cavity with an ideal-gas-style pressure-volume update
- Cavity volume change is smoothly bounded using a hyperbolic tangent saturation to improve numerical stability
- The effective loaded radius is allowed to evolve modestly with strain through a weak, capped geometry-based correction
- The model response is normalized to the zero-strain interface state and then mapped to the bare-port reference scale using S0_interface_bare
- Closure terms are included to approximate unresolved higher-order mechanics while preserving physical traceability in the reduced-order formulation

## Author

Tyler J. Inkley

Department of Ocean and Resources Engineering

University of Hawaiʻi at Mānoa