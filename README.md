# StrainModel

Reduced-order numerical model for strain-dependent pressure transmission in a membrane-confined hydrodynamic sensor interface.

# StrainModel

Reduced-order numerical model for strain-dependent pressure transmission in a membrane-confined hydrodynamic sensor interface.

## Example Output

![Model comparison](figures/model_result.png)

Example output from the reduced-order model showing normalized sensitivity as a function of membrane engineering strain, compared with experimental calibration data.

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

## Generated Figures

Running the script produces several figures that illustrate how the reduced-order model behaves and how it compares to experimental calibration results.

### 1. Normalized Sensitivity vs Engineering Strain

This is the primary publication-style figure used to compare the reduced-order model with experimental data.

The figure shows:

- **Experimental normalized sensitivity** for each membrane configuration (75%, 80%, 85%, and 90% interface scale) with error bars.
- A **bare-port reference line** at normalized sensitivity = 1, representing the case without a membrane interface.
- The **reduced-order model prediction** of normalized sensitivity as a continuous function of engineering strain.

This figure demonstrates that the model reproduces the experimentally observed trend:
- Sensitivity increases at moderate membrane strain
- Peaks near the **85% interface condition**
- Decreases at higher strain levels due to increasing interface stiffness and compressibility effects.

This figure corresponds directly to the comparison presented in the associated manuscript.

---

### 2. Model Decomposition Plot

This diagnostic figure separates the contributions of the two components of the reduced-order model:

S_norm(ε) = G(ε) · T(ε)

The plot shows three curves:

- **T(ε)** – transmission due to cavity compressibility and membrane/interface compliance
- **G(ε)** – an effective strain-dependent gain term representing additional coupling mechanisms
- **S_norm(ε)** – the resulting modeled normalized sensitivity

The decomposition helps illustrate how:

- **compressibility and interface compliance attenuate pressure transmission**, while  
- **membrane mechanics and geometric coupling can amplify the effective pressure signal**.

The combined effect produces the observed peak in sensitivity at intermediate strain levels.

---

### 3. Diagnostic Parameter Plot

This figure visualizes internal model parameters as functions of engineering strain.

The plot includes:

- **λ(ε) = P_ref · C_m / V_0**  
  A dimensionless parameter representing the relative influence of cavity compressibility and membrane compliance.

- **V₀(ε)**  
  The strain-dependent cavity volume used in the model.

- **C_m(ε)**  
  The effective interface compliance, which decreases with increasing strain due to membrane stiffening.

This diagnostic view helps interpret how strain simultaneously affects:

- cavity volume,
- membrane compliance, and
- the resulting pressure transmission behavior.

---

### 4. Optional Experimental Overlay

If the file `NSES_Long.fig` is present in the repository, the script will open the saved experimental figure and automatically overlay the reduced-order model prediction.

This allows direct comparison between:

- previously generated experimental calibration plots, and  
- the numerical model prediction.

This feature was used to produce the model–experiment comparison figure included in the associated publication.

## Related Publication

This model was developed to support the analysis presented in:

Inkley, T. J., Glass, G., & Krieg, M.  
*Modeling and Experimental Validation of a Hydrodynamic Sensor Module for Enhanced AUV Perception*  
(Under review)

## Author

Tyler J. Inkley  
Department of Ocean and Resources Engineering  
University of Hawaiʻi at Mānoa