# StrainModel

MATLAB tools for studying strain-dependent pressure transmission through the membrane-covered air cavities of a hydrodynamic sensor module.

The repository contains two workflows:

| Workflow | Purpose |
|---|---|
| `strainModel_v40.m` | Evaluate membrane–cavity models and generate manuscript Figures 7–9. |
| `hydrostatic_testing/hydrostatic_testing_v19.m` | Reduce experimental recordings and generate manuscript Figures 5–6. |

The model runs independently using embedded experimental summaries. Experimental reduction requires the raw recordings and trial-register workbook.

## Requirements

Tested with MATLAB R2025b.

The model script includes its numerical solvers and figure-formatting functions. No additional project files are required to run it.

The hydrostatic workflow uses the helper functions in `hydrostatic_testing/`. Keep experimental recordings outside the repository.

## Running the Model

Open MATLAB in the repository directory and run:

```matlab
strainModel_v40
```

Alternatively, open `strainModel_v40.m` in the MATLAB Editor and click **Run**.

Outputs are written to `results_v40` beside the script. Re-running overwrites matching output files. Save any results you want to retain before another run.

A successful run reaches `v40 complete` without an assertion failure.

## Running the Experimental Reduction

From the repository root:

```matlab
addpath('hydrostatic_testing');

dataRoot = '/absolute/path/to/1 Sensor Module v1';

% Inspect workbook selection without loading recordings.
hydrostatic_testing_v19(dataRoot, Mode="audit");

% Reduce recordings, verify manuscript values, and generate Figures 5–6.
hydrostatic_testing_v19(dataRoot);

% Redraw figures from existing output tables.
hydrostatic_testing_v19(Mode="render");
```

Alternatively, set the `SENSOR_TRIAL_DATA_ROOT` environment variable and call `hydrostatic_testing_v19` without a data-path argument.

See [the hydrostatic testing README](hydrostatic_testing/README.md) for data layout, workbook options, trial-selection rules, and output definitions.

## Physical Model

Each sensing pathway consists of a membrane coupled to a sealed air cavity. External pressure deflects the membrane, changing cavity volume and internal pressure until equilibrium is reached.

The two cavities are solved independently. Their internal pressure difference gives the predicted differential sensor response.

### Default Parameters

| Parameter | Value |
|---|---:|
| Pressure-loaded radius | 5.5 mm |
| Unstretched membrane thickness | 0.508 mm |
| Effective Young's modulus | 0.6 MPa |
| Poisson's ratio | 0.49 |
| Initial cavity volume, per side | Approximately 691 mm³ |
| Initial absolute gas pressure | 101325 Pa |
| Water density | 1000 kg/m³ |
| Test depth | 0.3048 m |
| Total applied pressure difference | 800 Pa |
| Pressure offsets about hydrostatic loading | ±400 Pa |

Air compression is quasi-static and isothermal. Additional tubing, sensor, and fitting volumes were not measured and are not included in the cavity volume.

Internal calculations use SI units. Exported tables identify quantities expressed in other units, including thickness in millimeters and structural energy in microjoules.

## Model Cases

The script evaluates three cases:

1. **Full nominal-strain tension:** retains installation tension calculated from nominal installation strain.
2. **Effective-strain tension:** calculates center thickness using a volume-conserving clamping rule, then infers an equivalent strain and retained tension.
3. **Locally relaxed nonlinear response:** removes installation-induced tension while retaining resistance from bending and deformation-induced stretching.

The full nominal-strain and locally relaxed cases use a Poisson-ratio-based installed-thickness relation. The effective-strain case uses the mapping below.

### Installation Strain

Nominal engineering strain is defined as:

```text
epsilon_nom = (Dinst - Dcut) / Dcut
```

Here, `Dcut` is the unstretched membrane diameter and `Dinst` is the installed diameter.

The experimental model states use the manuscript's rounded nominal strains:

```text
0.111, 0.176, 0.250, 0.333
```

Replacing these values with exact diameter-ratio strains changes the model inputs.

### Prescribed Clamping Rule

The effective-strain case assumes incompressible installation stretching and complete inward redistribution of material displaced beneath the washers:

```text
h_nom = h0 / (1 + epsilon_nom)^2
h_center = h_nom + 0.5 * (h0 - h_nom)

A_total * h_nom = A_washer * h_washer + A_center * h_center

epsilon_eff = sqrt(h0 / h_center) - 1
```

The rule restores half the thickness lost during installation. This fraction is prescribed, not measured or fitted to the response data.

Across the four experimental installation states, the implied washer-region thickness reduction is approximately 3.9–13.0% relative to the stretched thickness.

Complete inward redistribution is also an assumption. Outward displacement would reduce center thickening for a given washer compression.

Effective strain is an inferred model coordinate, not a measured local strain. Converting it into retained tension requires a separate mechanical approximation. All three cases are plotted against the same effective-strain coordinate.

## Experimental Sensitivities

The model embeds experimental means, standard errors, and trial counts at full numerical precision. Rounded values are:

| Interface scale | Nominal strain | Effective strain | Response-plot trials | Normalized trials | Mean normalized sensitivity ± SEM |
|---|---:|---:|---:|---:|---:|
| 90% | 0.111 | 0.05113 | 11 | 10 | 1.152 ± 0.085 |
| 85% | 0.176 | 0.07736 | 6 | 6 | 1.327 ± 0.082 |
| 80% | 0.250 | 0.10432 | 6 | 3 | 1.270 ± 0.045 |
| 75% | 0.333 | 0.13127 | 10 | 7 | 1.225 ± 0.042 |

Trial sensitivity is the magnitude of the fitted voltage–pressure slope. Each sensitivity is divided by the mean sensitivity of retained bare-port trials measured on the same calendar day.

The response plots contain 33 membrane trials and 17 bare-port trials. The normalized analysis contains 26 membrane trials, with same-day references provided by 13 bare-port trials.

SEM is the sample standard deviation of the normalized trial sensitivities divided by the square root of the number of included trials. Uncertainty and correlations associated with shared reference means are not separately propagated.

### Relationship Between the Workflows

The hydrostatic workflow exports:

```text
hydrostatic_testing/results/both_good/experimental_summary.csv
```

The model reads its experimental values from the embedded `obs` table. It does not automatically import this CSV.

Before updating `obs`, compare interface scales, trial counts, means, and SEM with the experimental output. Investigate discrepancies before rerunning the model fits.

Figures 5–6 use between-trial SD of recording mean voltages. Figure 9 uses SEM of normalized trial sensitivities. The configuration-level fits drawn in Figure 6 do not determine Figure 9's experimental values.

## Bare-Port Normalization

Predicted pressure transmission is normalized as:

```text
R_relative = R_cavity / R_bare
```

Each model case has a separately fitted reference factor:

| Case | R_bare | Fit objective |
|---|---:|---|
| Full nominal-strain tension | Approximately 0.264 | Unweighted least squares at 85%, 80%, and 75% |
| Effective-strain tension | Approximately 0.364 | Unweighted least squares at 85%, 80%, and 75% |
| Locally relaxed nonlinear response | Approximately 0.727 | Unweighted least squares at all four scales |

The objectives fit sensitivity levels rather than slopes alone.

These factors are conditional on their respective models. They are not three independently measured properties of the bare-port pathway. Agreement with the data used for fitting does not provide independent model validation.

Changing `R_bare` rescales normalized sensitivity and its slope. It does not change membrane equilibrium, raw pressure transmission, or structural energy.

Inverse-SEM-squared weighted fits are also reported. Their factors are approximately 0.257, 0.356, and 0.733, respectively. These fits do not establish confidence intervals because shared-reference correlations are not included.

## Structural Energy

Structural energy is evaluated at the hydrostatic, high-side, and low-side equilibrium states.

This calculation uses the thickness and retained tension of the effective-strain case, together with deformation-induced stretching and a displaced-volume coefficient of:

```text
C_V = 1/2
```

The energy calculation is specific to this nonlinear formulation. It does not evaluate the structural energy of all three pressure-transmission cases.

Changing the clamping rule can change the calculated energies. Changing `R_bare` does not.

## Model Outputs

The script generates:

```text
results_v40/
    formatted/
        fig07_potential_energy_vs_tension.fig
        fig07_potential_energy_vs_tension.png
        fig07_potential_energy_vs_tension.pdf
        fig07_potential_energy_vs_tension.eps

        fig08_model_pressure_transmission.fig
        fig08_model_pressure_transmission.png
        fig08_model_pressure_transmission.pdf
        fig08_model_pressure_transmission.eps

        fig09_fitted_sensitivity_comparison.fig
        fig09_fitted_sensitivity_comparison.png
        fig09_fitted_sensitivity_comparison.pdf
        fig09_fitted_sensitivity_comparison.eps

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

Additional source figures and an energy-versus-strain diagnostic are also generated.

Formatted outputs include editable MATLAB figures, 600-dpi PNGs, and vector PDF/EPS files.

### Figure Definitions

- **Figure 7:** structural energy versus retained tension.
- **Figure 8:** raw pressure transmission, without fitted reference scaling.
- **Figure 9:** fitted bare-port-normalized sensitivity with experimental means and SEM.

Figure 9 displays sensitivities from 0.8 to 1.6. Full calculated curves remain available in the CSV outputs.

## Numerical Validation

Each model run checks:

- Convergence across the 243-state strain grid.
- Finite raw pressure transmission between zero and one.
- Volume conservation under the clamping rule.
- A strictly increasing effective-strain coordinate beginning at zero.
- Inclusion of the four experimental installation states.
- Agreement between each unweighted analytical normalization fit and a numerical parameter sweep.
- Agreement between bisection and under-relaxed fixed-point solutions for the two retained-tension cases.

The solver comparison includes zero strain and the four experimental states under hydrostatic and ±400 Pa loading, giving 30 comparisons.

For the default configuration, all 30 comparisons converge. The maximum cavity-pressure difference between solvers is approximately `7.88e-5 Pa`, below the `8e-5 Pa` assertion threshold.

### Solver Tolerances and Interpretation

- The retained-tension bisection solver uses a pressure-residual or pressure-bracket tolerance of `1e-5 Pa`.
- The nonlinear deflection solver stops when the pressure residual is below `1e-5 Pa` or the deflection bracket is below `1e-12 m`.
- Nonlinear iteration exhaustion returns `converged=false`.
- A deflection-bracket stopping condition does not directly bound the pressure residual.
- The fixed-point solver stops on the change between relaxed iterates. Inspect its reported equilibrium residual separately.
- The bisection solution is independent of the initial pressure estimate by construction. Agreement between direct and hydrostatic-start calculations does not independently establish physical loading-path independence.
- SEM-weighted fits are calculated analytically but are not checked against the numerical parameter sweep.
- The nonlinear case is not cross-checked with a second solver.

The model uses 1001-point radial integration and a small-argument plate approximation. Default validation does not include a quadrature-refinement study.

## Checking a Run

After execution:

1. Inspect `champion_validation.txt` and `solver_validation.txt`.
2. Confirm that `champion_fit_summary.csv` gives the expected rounded reference factors.
3. Compare `champion_experimental_comparison.csv` with the manuscript values.
4. Confirm that all 30 solver comparisons pass.
5. Inspect Figures 7–9 for labels, legends, uncertainty markers, and axis limits.

For the experimental workflow:

```matlab
addpath('hydrostatic_testing');
testExperimentalReduction;

reductionDir = fullfile('hydrostatic_testing','results','both_good');
verifyManuscriptSummary(reductionDir);
```

Review any skipped tests or diagnostic flags. Numerical checks establish consistency and repeatability within the implemented equations; they do not independently validate the physical assumptions.

## Interpretation and Limitations

None of the three model cases predicts amplification of the applied pressure difference.

The retained-tension cases capture the direction of the measured postpeak decrease but predict a larger decrease than observed. The locally relaxed case predicts an increasing response. The models do not capture the apparent intermediate-strain maximum.

Fitted normalization improves agreement in response magnitude without uniquely identifying membrane installation state or absolute bare-port transmission.

Additional limitations include:

- Quasi-static loading.
- No explicit membrane inertia, damping, external fluid inertia, or viscous fluid–structure interaction.
- No explicit material relaxation, hysteresis, clamp slip, or loading-history dependence.
- Unmeasured washer compression, local strain, and retained tension.
- A prescribed clamping rule that requires experimental assessment.
- No validated capillary-pressure correction through bare-port normalization.
- A nonlinear branch intended for compressive loading with nonnegative deflection and positive cavity volume.
- Numerical volume floors that do not establish physical validity outside the intended parameter range.

Results at other depths, loading conditions, or material parameters require additional numerical and experimental assessment.

## Data and Repository Practices

Raw recordings are stored separately from the code. Reproducing Figures 5–6 requires access to those recordings and the trial register; the model calculations and Figures 7–9 can run from the embedded summaries.

Keep trial-selection decisions, recording-order checks, diagnostics, and summary tables together when retaining experimental results.

Generated outputs and local configuration are excluded from version control by the repository's Git rules. Custom output locations may require additional exclusions. Git rules do not remove files that are already tracked; inspect the staged file list before committing.

## Citation

If you use this code, please cite the associated manuscript:

*Modeling and Experimental Validation of a Hydrodynamic Sensor Module for Enhanced Autonomous Underwater Vehicle Perception.*

Add the final journal citation and DOI when available.