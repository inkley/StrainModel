# StrainModel

Reduced-order MATLAB model for strain-dependent pressure transmission in a membrane-confined hydrodynamic sensor interface.

## Overview

This repository contains a reduced-order numerical model used to interpret experimentally observed strain-dependent pressure sensitivity in a hydrodynamic sensor module for Autonomous Underwater Vehicles (AUVs).

The current implementation uses an iterative force-balance formulation in which strain-dependent interface geometry and plate stiffness modify the pressure transmission of an air-filled sensing cavity. Rather than imposing a separate empirical coupling curve, the model computes normalized sensitivity directly from the cavity/interface response.

The model output is reported as:

S_norm(eps) = T_cav,norm(eps)

where `T_cav,norm(eps)` is the normalized cavity/interface transmission response, evaluated relative to its zero-strain value.

In the present formulation, strain affects the model through:

- installed interface diameter  
- effective plate radius  
- strain-dependent effective modulus  
- plate bending stiffness  
- embedded geometry coupling through `beta_geom`  

This produces a smooth, monotonic amplification trend that can be compared directly against experimental calibration data.

## Physical Interpretation

The model treats each sensing side as an air-filled cavity coupled to a strained membrane/plate-like interface. A small applied surface pressure differential is transmitted through the interface, while the internal cavity pressure evolves iteratively according to:

- trapped-air compression  
- strain-dependent interface geometry  
- strain-dependent plate stiffness  
- a capped and relaxed net force-balance update  

The model is intended as a reduced-order interpretive tool that preserves physical traceability without attempting to fully resolve nonlinear membrane mechanics.

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

## Key Model Parameters

Important model parameters include:

- `beta_geom`  
  Blending factor for embedded geometry amplification

- `c1`, `c2`  
  Coefficients controlling the strain-dependent effective modulus

- `q_cap`  
  Cap on net plate load during iterative updates

- `dV_cap_fraction`  
  Cap on cavity volume change relative to nominal cavity volume

- `relax`, `relax_q`  
  Relaxation factors used to improve solver stability

- `tolP`  
  Pressure convergence tolerance for the iterative cavity solve

## Current Interface and Cavity Assumptions

The current script uses:

- **Interface material:** hygienic latex  
- **Membrane/interface thickness:** 0.020 in  
- **Nominal effective modulus:** 6.0e5 Pa  
- **Working fluid:** trapped air  
- **Cavity model:** isothermal ideal gas compression  
- **Interface representation:** circular plate-style reduced-order approximation  
- **Cavity volume:** estimated from a cylindrical plus spherical-cap geometry per sensing side  

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

The current version also includes parameter sensitivity studies to assess robustness of the reduced-order formulation.

Implemented sweeps include:

- `beta_geom` sweep  
- `c1` sweep  
- `c2` sweep  
- 2D contour sweep of `beta_geom` versus engineering strain  

These plots help show whether the observed model trend is robust to reasonable parameter variation, rather than being dependent on a single narrowly chosen setting.

## Example Outputs

The script generates:

- baseline normalized sensitivity comparison against experimental data  
- transmission and embedded geometry diagnostics  
- iterative cavity volume / plate deflection diagnostics  
- plate mechanics diagnostics  
- interface diameter interpolation plot  
- solver iteration count and convergence flag plots  
- parameter sensitivity sweep plots for `beta_geom`, `c1`, and `c2`  
- 2D contour plot of normalized sensitivity over `beta_geom` and strain  
- error metric plot showing RMSE and WRMSE versus `beta_geom`  

## Files

- `sensor_interface_numerical_model_v10_sensitivity_contour.m`  
  Main MATLAB script implementing the iterative force-balance model, diagnostics, sensitivity sweeps, and 2D contour analysis

## Getting Started

Open MATLAB in this repository folder and run:

```matlab
sensor_interface_numerical_model_v10_sensitivity_contour
```

The script will:

1. Define experimental strain and sensitivity data  
2. Evaluate the baseline force-balance model over a continuous strain range  
3. Print a comparison table between experimental and modeled values  
4. Report RMSE and weighted RMSE  
5. Generate baseline diagnostic plots  
6. Run 1D parameter sensitivity sweeps  
7. Generate a 2D `beta_geom` versus strain contour plot  
8. Save figures to the `results/` folder as both `.png` and `.fig` files  

## Notes

- The present model is a reduced-order interpretive tool, not a full nonlinear membrane mechanics solution  
- The interface mechanics are represented using a plate-style stiffness approximation with strain-dependent modulus  
- The cavity-pressure solution assumes an isothermal air-filled cavity and ideal-gas compression  
- Geometry amplification is not imposed as a standalone fit curve; it enters through the embedded area-ratio formulation and the blending parameter `beta_geom`  
- Sensitivity sweeps are included to assess robustness and physical defensibility of the overall trend  

## Author

Tyler J. Inkley  
Department of Ocean and Resources Engineering  
University of Hawaiʻi at Mānoa  