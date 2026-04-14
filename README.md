# StrainModel

Reduced-order MATLAB model for strain-dependent pressure transmission in a membrane-confined hydrodynamic sensor interface.

## Overview

This repository contains a reduced-order numerical model used to interpret experimentally observed strain-dependent pressure sensitivity in a hydrodynamic sensor module for Autonomous Underwater Vehicles (AUVs).

The current implementation, `sensor_interface_numerical_model_v11.m`, uses an iterative force-balance formulation in which measured strain-dependent interface geometry and a reduced-order plate response modify pressure transmission through an air-filled sensing cavity. Rather than imposing a standalone empirical coupling curve, the model computes normalized sensitivity directly from the cavity/interface response.

The model output is reported as:

`S_norm(eps) = T_cav,norm(eps)`

where `T_cav,norm(eps)` is the normalized cavity/interface transmission response evaluated relative to its zero-strain value.

A geometry-only comparison curve is also included in v11:

`S_geom_only(eps) = G_eff(eps)`

This provides a direct visual comparison between:
- the embedded geometry contribution alone
- the full iterative cavity/interface response

In the present formulation, strain affects the model through:

- installed interface diameter
- effective plate radius
- strain-dependent effective modulus
- plate bending stiffness
- embedded geometry coupling through `beta_geom`

## Physical Interpretation

The model treats each sensing side as an air-filled cavity coupled to a strained membrane/plate-like interface. A small applied surface pressure perturbation is transmitted through the interface, while internal cavity pressure evolves iteratively according to:

- trapped-air compression
- strain-dependent interface geometry
- strain-dependent plate stiffness
- a capped and relaxed net force-balance update

The present implementation is intended as a reduced-order interpretive tool that preserves physical traceability without attempting to fully resolve nonlinear membrane mechanics, large-deflection membrane tension, or a full spatial fluid-structure interaction field.

## Current Model Formulation

For a given engineering strain `eps`, the model:

1. Interpolates the installed interface diameter from measured strain states
2. Computes effective plate radius `a_plate(eps)`
3. Computes effective modulus using  
   `E_eff(eps) = E_plate0 * (1 + c1*eps + c2*eps^2)`
4. Computes plate bending stiffness  
   `D_plate(eps) = E_eff * t^3 / (12 * (1 - nu^2))`
5. Computes a geometry-based gain term from the forced loading area relative to the installed interface area
6. Blends that gain through the geometry coupling parameter `beta_geom`
7. Solves cavity pressure iteratively using a relaxed and capped net-load update
8. Computes normalized transmission from symmetric positive and negative pressure perturbations
9. Reports both the full model response and a geometry-only comparison curve

## Plate-Deflection Basis

The plate response used in the model is based on the classical small-deflection circular plate solution under uniform load. In reduced-order form, the code uses the following closed-form results:

- maximum center deflection  
  `w_max = q*a^4 / (64*D)`

- integrated deflection volume  
  `dV = pi*q*a^6 / (192*D)`

These expressions are the closed-form consequences of the circular plate deflection solution and are used directly to update cavity volume and internal pressure during each iterative solve. The model therefore incorporates the plate-deflection and integrated volume-change logic in reduced-order form, rather than explicitly solving the full spatial deflection field `w(r)` at every step.

## Key Model Parameters

Important model parameters include:

- `beta_geom`  
  Geometry coupling parameter controlling how embedded area-ratio amplification enters the effective forcing

- `c1`, `c2`  
  Linear and quadratic coefficients controlling the strain-dependent effective modulus

- `q_cap`  
  Cap on net plate load during iterative updates

- `dV_cap_fraction`  
  Cap on cavity volume change relative to nominal cavity volume

- `relax`, `relax_q`  
  Relaxation factors used to improve solver stability

- `tolP`  
  Pressure convergence tolerance for the iterative cavity solve

## Current Interpretation of v11 Results

The cleaned v11 implementation preserved stable convergence across the full strain sweep while improving interpretability of the dominant reduced-order mechanisms.

Current results suggest that:

- the model is numerically stable and convergent across the full strain range
- the normalized sensitivity trend is dominated primarily by strain-dependent embedded geometry
- cavity compression and plate deflection act as secondary corrections in the present formulation
- the geometry coupling parameter `beta_geom` remains the most influential fitting parameter over the tested range
- the strain-dependent stiffness coefficients `c1` and `c2` were retained, but earlier sweeps showed that they had minimal effect on fit quality relative to `beta_geom`

In other words, the present model supports a geometry-dominated interpretation of the observed sensitivity enhancement, with cavity mechanics providing a smaller corrective contribution.

## Current Interface and Cavity Assumptions

The current script uses:

- **Interface material:** hygienic latex
- **Membrane/interface thickness:** 0.020 in
- **Nominal effective modulus:** 6.0e5 Pa
- **Working fluid:** trapped air
- **Cavity model:** isothermal ideal-gas-style pressure-volume update
- **Interface representation:** circular plate-style reduced-order approximation
- **Cavity volume:** estimated from a cylindrical plus spherical-cap geometry per sensing side
- **Load model:** uniform reduced-order net load on the interface

## Experimental Reference Points

The script includes experimental calibration summary values for:

- engineering strain
- normalized sensitivity
- approximate uncertainty / error bars

These are used for direct visual comparison with model predictions and for reporting:

- residuals
- RMSE
- weighted RMSE

## Sensitivity Analysis

The current v11 implementation includes sensitivity studies focused on the dominant reduced-order fitting parameter.

Implemented sweeps include:

- `beta_geom` sweep
- 2D contour sweep of `beta_geom` versus engineering strain
- RMSE and WRMSE versus `beta_geom`

Earlier exploratory versions also included `c1` and `c2` sweeps, but these were removed from the current presentation-oriented v11 workflow because they produced minimal variation in model response over the tested strain range.

## Example Outputs

The script generates:

- baseline normalized sensitivity comparison against experimental data
- geometry-only versus full-model comparison
- transmission and embedded geometry diagnostics
- iterative cavity volume / plate deflection diagnostics
- plate mechanics diagnostics
- interface diameter interpolation plot
- solver iteration count and convergence flag plots
- `beta_geom` sensitivity sweep
- 2D contour plot of normalized sensitivity over `beta_geom` and strain
- error metric plot showing RMSE and WRMSE versus `beta_geom`

## Files

- `sensor_interface_numerical_model_v11.m`  
  Main MATLAB script implementing the iterative force-balance model, diagnostics, geometry-only comparison, `beta_geom` sweeps, and 2D contour analysis

## Getting Started

Open MATLAB in this repository folder and run:

```matlab
sensor_interface_numerical_model_v11
```

The script will:

	1.	Define experimental strain and sensitivity data
	2.	Evaluate the baseline force-balance model over a continuous strain range
	3.	Print a comparison table between experimental and modeled values
	4.	Report RMSE and weighted RMSE
	5.	Generate baseline diagnostic plots
	6.	Generate a geometry-only versus full-model comparison
	7.	Run beta_geom sensitivity sweeps
	8.	Generate a 2D beta_geom versus strain contour plot
	9.	Save figures to the results_v11/ folder as both .png and .fig files

## Notes

- The present model is a reduced-order interpretive tool, not a full nonlinear membrane mechanics solution
- The interface mechanics are represented using a plate-style stiffness approximation with strain-dependent effective modulus
- The plate-deflection and integrated deflection-volume relations are incorporated in closed form through the reduced-order update used in loadedCavityState_fromDeltaP_v11
- The cavity-pressure solution assumes an air-filled cavity with an ideal-gas-style pressure-volume update
- Geometry amplification is not imposed as a standalone fit curve; it enters through the embedded area-ratio formulation and the coupling parameter beta_geom
- The v11 geometry-only comparison was added to clarify how much of the predicted sensitivity trend is driven by embedded geometry versus cavity feedback
- Small localized kinks may remain in the response curve due to the capped and relaxed iterative update structure, but the solver is fully convergent across the tested strain range

## Author

Tyler J. Inkley

Department of Ocean and Resources Engineering

University of Hawaiʻi at Mānoa