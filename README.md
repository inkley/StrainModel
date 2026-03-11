# StrainModel

Reduced-order numerical model for strain-dependent pressure transmission in a membrane-confined hydrodynamic sensor interface.

## Overview

This repository contains a first-pass MATLAB reduced-order model used to interpret experimentally observed strain-dependent pressure sensitivity in a hydrodynamic sensor module for Autonomous Underwater Vehicles (AUVs).

The model represents normalized sensitivity as:

S_norm(eps) = G(eps) * T(eps)

where:
- `T(eps)` captures transmission effects associated with cavity compressibility and interface compliance
- `G(eps)` is an effective strain-dependent gain term introduced to reproduce experimentally observed sensitivity amplification at moderate strain

The script compares model predictions against experimentally derived normalized sensitivity values and generates publication-style and diagnostic plots.

## Files

- `StrainModel.m`  
  Main MATLAB script implementing the reduced-order model, parameter definitions, plotting routines, and local cavity-pressure solver.

- `NSES_Long.fig`  
  MATLAB figure file containing the experimental normalized sensitivity plot used for optional overlay of the reduced-order model.

## Model Features

- Strain-dependent cavity volume
- Strain-dependent membrane/interface compliance
- Effective peaked gain term
- Small-signal pressure sweep for slope-based sensitivity estimation
- Comparison against experimental calibration summary points
- Optional overlay onto saved experimental MATLAB figure

## Experimental Reference Points

The script includes experimental summary points for:
- membrane engineering strain
- normalized sensitivity
- approximate error bars

These values are used for visual comparison between experimental observations and model output.

## Usage

Open MATLAB in this repository folder and run:

```matlab
StrainModel
```

The script will:

1. Define calibration strain points and experimental summary data
2. Evaluate the reduced-order model across a continuous strain range
3. Compute normalized sensitivity based on cavity compressibility and interface compliance
4. Print a comparison between experimental and modeled sensitivity values
5. Generate publication-style and diagnostic plots
6. Optionally overlay the reduced-order model on `NSES_Long.fig` if the figure file is present

## Notes

- The present model is intended as a reduced-order interpretive tool rather than a fully predictive physics model.
- The cavity-pressure solution assumes an **isothermal air-filled cavity**, consistent with the calibration configuration used in the experiments.
- Future implementations may replace the air-filled cavity with mineral oil to reduce compressibility effects at operational depth.

## Author

Tyler J. Inkley  
Department of Ocean and Resources Engineering  
University of Hawaiʻi at Mānoa