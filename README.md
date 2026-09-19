# StrainModel

MATLAB model for strain-dependent pressure transmission through the membrane-covered air cavities of a hydrodynamic sensor module.

## Current Version

The current implementation is `strainModel_v40.m`.

Version 40 is a standalone cleanup of the v39 manuscript model. It embeds the experimental summaries, numerical solvers, figure formatting, and solver checks in one file. It requires no other project files.

The model evaluates how membrane installation, bending, retained tension, deformation-induced tension, and air compression affect differential-pressure transmission. Its purpose is to examine observed trends and the consequences of uncertain installation conditions.

## Requirements and Running

Tested with MATLAB R2025b.

Copy `strainModel_v40.m` into a writable directory, open MATLAB in that directory, and run:

```matlab
strainModel_v40
```

Alternatively, open the script in the MATLAB Editor and click **Run**.

Outputs are written to `results_v40` beside the script. Re-running replaces the generated files in that folder.

There are no run modes or external formatter dependencies. Earlier exploratory sweeps remain in the previous model versions.

## Physical Model

Each sensing pathway consists of a membrane coupled to a sealed air cavity. External pressure deflects the membrane, changing cavity volume and internal pressure until equilibrium is reached.

The two cavities are solved independently. Their internal pressure difference gives the predicted differential sensor response.

Primary parameters are:

| Parameter | Value |
|---|---:|
| Pressure-loaded radius | 5.5 mm |
| Unstretched membrane thickness | 0.508 mm |
| Effective Young's modulus | 0.6 MPa |
| Poisson's ratio | 0.49 |
| Initial cavity volume, per side | Approximately 691 mm³ |
| Initial absolute gas pressure | 101325 Pa |
| Assumed water density | 1000 kg/m³ |
| Calibration depth | 0.3048 m |
| Total applied pressure difference | 800 Pa |
| Pressure offsets about hydrostatic loading | ±400 Pa |

Air compression is quasi-static and isothermal. Additional tubing, sensor, and fitting volumes were not measured and are not included in the cavity volume.

## Model Cases

The script compares three cases:

1. **Full nominal-strain tension:** retains the installed tension calculated from nominal installation strain.
2. **Effective-strain tension, volume-conserving:** calculates center thickness using an illustrative clamping rule, then infers an equivalent strain and retained tension.
3. **Locally relaxed nonlinear response:** removes retained installation tension while preserving bending and deformation-induced stretching resistance.

The nominal and locally relaxed cases retain the Poisson-ratio-based thickness relation used in v39. The effective-strain case uses the volume-conservation mapping below.

### Illustrative Clamping Rule

The effective-strain case assumes incompressible installation stretching and complete inward redistribution of material displaced beneath the washers:

```text
h_nom = h0 / (1 + epsilon_nom)^2
h_center = h_nom + 0.5 * (h0 - h_nom)

A_total * h_nom = A_washer * h_washer + A_center * h_center

epsilon_eff = sqrt(h0 / h_center) - 1
```

The rule restores half the thickness lost during installation. This fraction is prescribed, not measured or fitted to the response data.

Across the four experimental installation states, the implied washer-region thickness reduction is approximately 3.9–13.0% relative to the stretched thickness.

Complete inward redistribution is also an assumption. Outward material displacement would reduce center thickening for a given washer compression.

The effective strain is an inferred model coordinate, not a measured local strain. Converting that coordinate into retained tension is a separate mechanical approximation. All three cases are plotted against the same effective-strain coordinate for comparison.

## Experimental Data

The script embeds the retained experimental means, SEMs, and trial counts at their archived precision.

| Interface scale | Nominal strain | Inferred effective strain | Response-plot trials | Normalized trials | Mean normalized sensitivity ± SEM |
|---|---:|---:|---:|---:|---:|
| 90% | 0.111 | 0.05113 | 11 | 10 | 1.152 ± 0.085 |
| 85% | 0.176 | 0.07736 | 6 | 6 | 1.327 ± 0.082 |
| 80% | 0.250 | 0.10432 | 6 | 3 | 1.270 ± 0.045 |
| 75% | 0.333 | 0.13127 | 10 | 7 | 1.225 ± 0.042 |

Experimental sensitivity is the magnitude of each trial's fitted voltage–pressure slope, divided by the mean sensitivity of retained bare-port trials from the same day.

The normalized analysis contains 26 membrane trials. Thirteen bare-port trials provide their same-day references. The response plots contain 33 membrane trials and 17 bare-port trials.

SEMs describe variation among normalized trial values. Uncertainty in the estimated same-day reference means is not separately propagated.

Version 40 reproduces the model comparison from these embedded summaries. It does not reprocess raw calibration recordings or regenerate manuscript Figures 5–6.

## Bare-Port Normalization

Raw pressure transmission is normalized as:

```text
R_relative = R_cavity / R_bare
```

The displayed cases use separate conditional reference factors:

| Case | R_bare | Selection method |
|---|---:|---|
| Full nominal-strain tension | 0.290 | Selected to place the prediction slightly below the measured 85% mean |
| Effective-strain tension | Approximately 0.364 | Unweighted least-squares fit to response levels at 85%, 80%, and 75% |
| Locally relaxed nonlinear response | Approximately 0.727 | Unweighted least-squares fit to response levels at all four scales |

The nominal value is not the least-squares optimum or a rigorous lower bound. Its least-squares alternative is retained in the fit summary.

These factors are alternatives conditional on each model, not three measured properties of the same bare-port pathway. Fitted agreement is calibration, not independent model validation.

Changing `R_bare` rescales the normalized response and its slope. It does not change membrane equilibrium, raw transmission, or structural energy.

Inverse-SEM-squared weighted fits are also reported for comparison. They do not establish confidence intervals because shared-reference correlations are not included.

## Interpretation

All three cases predict passive raw pressure transmission.

The retained-tension cases capture the direction of the measured postpeak decrease, but predict a larger decrease than observed. The locally relaxed case predicts an increasing response. None reproduces the apparent intermediate-strain maximum.

The selected and fitted normalizations improve agreement in response magnitude without uniquely identifying the installed membrane state or absolute bare-port transmission.

## Potential-Energy Calculation

Structural energy is evaluated at the hydrostatic, high-side, and low-side equilibrium states.

This is a separate nonlinear calculation using the thickness and retained tension of the volume-conserving effective-strain case, together with deformation-induced stretching and a displaced-volume coefficient of `C_V = 1/2`.

It is not a calculation of the structural energy of all three transmission cases. Changing the clamping rule can change these energies; changing `R_bare` does not.

## Outputs

The script creates:

```text
results_v40/
    formatted/
        fig07_potential_energy_vs_tension.fig
        fig07_potential_energy_vs_tension.png
        fig07_potential_energy_vs_tension.pdf
        fig07_potential_energy_vs_tension.eps
        fig08_model_experiment_comparison.fig
        fig08_model_experiment_comparison.png
        fig08_model_experiment_comparison.pdf
        fig08_model_experiment_comparison.eps

    champion_curves.csv
    champion_geometry.csv
    champion_energy.csv
    champion_fit_summary.csv
    champion_experimental_comparison.csv
    embedded_experimental_summary.csv
    champion_workspace.mat
    champion_validation.txt
    solver_comparison.csv
    solver_validation.txt
```

Unformatted source figures and an energy-versus-strain diagnostic are also generated.

The `champion_*` table names are retained for continuity with v39. The formatted outputs include editable MATLAB figures, 600-dpi PNGs, and vector PDF/EPS files.

## Numerical Checks

Each run checks:

- Convergence across the 243-state strain grid.
- Finite, passive raw pressure transmission.
- Volume conservation under the prescribed clamping rule.
- A strictly increasing effective-strain coordinate beginning at zero.
- Inclusion of the exact experimental installation states.
- Agreement between analytical normalization fits and a numerical parameter sweep.
- Agreement between bisection and under-relaxed fixed-point solutions for the two retained-tension cases.

The solver comparison includes zero strain and the four experimental installation states, each evaluated under hydrostatic and ±400 Pa loading: 30 comparisons in total.

For the verified default configuration, both solvers converged in all 30 comparisons. Their maximum cavity-pressure difference was approximately `7.88e-5 Pa`.

The fixed-point solver stops on the change between relaxed iterates, rather than on the equilibrium residual. Its stopping tolerance therefore does not directly bound its difference from bisection.

The bisection solver ignores the supplied starting pressure. Identical direct and hydrostatic-start solutions are consequently expected by construction and do not independently establish physical loading-path independence.

## Version 40 Verification

Version 40 was tested in an isolated folder containing only `strainModel_v40.m`, with MATLAB's path reset before execution.

Verification confirmed:

- Zero MATLAB Code Analyzer messages, without warning suppressions.
- Successful standalone execution.
- Exact agreement with the saved v39 numerical tables.
- Successful completion of all 30 solver comparisons.
- Correctly rendered formatted Figures 7–8.

Release-check records are saved separately as:

```text
results_v40/code_analyzer.txt
results_v40/v39_parity.txt
results_v40/standalone_run.log
```

These records are not runtime inputs and are not regenerated by the model script.

## Scope and Limitations

- The model is quasi-static.
- Dynamic membrane inertia, damping, external fluid inertia, and viscous fluid–structure interactions are excluded.
- Material relaxation, hysteresis, clamp slip, and other history-dependent behavior are not explicitly modeled.
- Washer compression, local strain, and retained tension were not measured.
- The clamping rule and inferred tension require experimental assessment.
- Bare-port normalization does not constitute a validated capillary-pressure correction.
- The default numerical checks do not establish physical validity at untested depths or loading conditions.

Earlier versions retain development history and exploratory diagnostics. They should not be interpreted as the current manuscript formulation.

## Citation

If you use this model, please cite the associated hydrodynamic sensor-module manuscript once its final bibliographic information is available.