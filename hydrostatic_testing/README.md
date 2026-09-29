# Hydrostatic Testing

MATLAB workflow for reducing longitudinal hydrostatic-test recordings and generating manuscript Figures 5–6.

The workflow selects trials from an experimental register, calculates voltage–pressure sensitivities, normalizes sensitivities using same-day bare-port references, and exports summary tables and figures.

## Requirements

- MATLAB R2025b, the release used for testing.
- Experimental recordings and the trial-register workbook.
- No additional MATLAB toolboxes are required.

Keep recordings outside the repository. The expected data layout is:

```text
<dataRoot>/
    SM v1 Test Cases.xlsx
    Trial 3/
        Data/
            <recording files>
    Trial 5/
        Data/
            <recording files>
    ...
```

The trial register contains eight columns, in this order:

1. Trial number
2. Date
3. Interface scale
4. Working fluid
5. Resolution
6. Qualitative reliability flag
7. Residual threshold flag
8. Notes

## Quick Start

From the repository root:

```matlab
addpath('hydrostatic_testing');

% Assign the external recording directory.
dataRoot = '/absolute/path/to/1 Sensor Module v1';

% Inspect trial selection without loading recordings.
hydrostatic_testing_v19(dataRoot, Mode="audit");

% Reduce recordings, verify manuscript values, and generate figures.
result = hydrostatic_testing_v19(dataRoot);
```

Alternatively, configure the `SENSOR_TRIAL_DATA_ROOT` environment variable and run:

```matlab
hydrostatic_testing_v19;
```

A full reduction reads every selected recording. Ensure cloud-hosted files are downloaded before running.

## Run Options

| Option | Default | Purpose |
|---|---|---|
| `Mode` | `"reduce"` | Reduce recordings, verify summaries, and render figures. |
| `Mode="audit"` | — | Read the workbook and export trial-selection decisions only. |
| `Mode="render"` | — | Verify existing summaries and redraw figures without reading recordings. |
| `OutputDir` | `results/both_good` beside the functions | Set the output directory. |
| `Workbook` | `SM v1 Test Cases.xlsx` inside `dataRoot` | Select another trial register. |
| `Sheet` | `1` | Select the worksheet containing the trial register. |

For example:

```matlab
result = hydrostatic_testing_v19(dataRoot, ...
    Workbook=fullfile(dataRoot,'SM v1 Test Cases.xlsx'), ...
    Sheet=1, ...
    OutputDir=fullfile(pwd,'hydrostatic_testing','results','comparison'));
```

Existing output files with matching names are overwritten. Use separate output directories when comparing runs.

## Trial Selection

The manuscript workflow includes trials that:

- Have `Good` qualitative reliability and residual threshold flags.
- Use an interface scale of 75%, 80%, 85%, 90%, or the 100% bare-port reference.
- Are not identified as mineral-oil trials.

The workbook flags determine inclusion. The code also recalculates each trial's voltage–pressure fit RMSE and reports values above 0.08 V. These diagnostics do not automatically override the workbook selection; investigate any disagreement before interpreting the results.

The lower-level `runExperimentalReduction` function also supports `Selection="longitudinal_good"`, which additionally accepts the `Bad (A)` residual flag. This alternative is not used for the manuscript results.

## Recording Reduction

Each selected trial must contain exactly 25 nonhidden recording files, corresponding to orientations from 0 to 360 degrees in 15-degree increments.

Files are assigned orientations in lexicographic filename order. Filename angle labels are checked and exported for inspection, but they do not change the assigned order.

For each recording, the workflow:

1. Reads the first 10,000 voltage samples from row 2.
2. Calculates the mean voltage and within-record sample standard deviation.
3. Assigns the longitudinal hydrostatic pressure difference:

   ```matlab
   pressure = 1000 * 9.81 * 0.040004 * sin(beta);
   ```

   Pressure is in pascals, and `beta` is in radians.

4. Fits recording mean voltage against pressure using unweighted linear least squares with an intercept.
5. Defines trial sensitivity as the magnitude of the fitted slope, in V/Pa.

Missing files, malformed recordings, and nonfinite voltage samples stop the reduction.

## Same-Day Normalization

Each trial sensitivity is divided by the mean sensitivity of selected bare-port trials measured on the same calendar day.

Trials without an eligible same-day reference remain available for the raw-voltage figures but receive no normalized sensitivity.

Normalized membrane sensitivities are grouped by interface scale. Each group's reported uncertainty is:

```text
SEM = sample standard deviation of normalized trial sensitivities / sqrt(n)
```

Here, `n` is the number of trials with finite normalized sensitivities. Groups with fewer than two such trials have an undefined SEM.

This SEM does not separately propagate uncertainty or correlations introduced by shared bare-port references.

## Figures and Uncertainty

### Figure 5: Voltage–Orientation Response

Markers show the mean recording voltage across selected trials at each orientation. Dashed bounds show the mean plus or minus one between-trial sample standard deviation.

### Figure 6: Voltage–Pressure Response

Markers show the same mean voltages plotted against pressure. Solid lines are linear fits to the configuration mean responses.

At repeated pressure values, dashed curves show the outer envelope of the pointwise mean ± SD bounds. All orientation-specific mean markers are retained.

Both figures use uncorrected recording mean voltages.

### Distinction from Figure 9

Three different uncertainty measures appear in the workflow:

| Quantity | Meaning |
|---|---|
| Within-record SD | Variation among voltage samples within one recording. |
| Figures 5–6 SD | Variation among trial recording means at a given orientation. |
| Figure 9 SEM | Uncertainty in the mean of individually normalized trial sensitivities. |

The configuration-level fits drawn in Figure 6 are not used to calculate Figure 9's experimental sensitivities.

## Output Files

Default output location:

```text
hydrostatic_testing/results/both_good/
```

| File | Contents |
|---|---|
| `trial_selection.csv` | Workbook metadata and inclusion decisions. |
| `trial_diagnostics.csv` | Trial slopes, intercepts, RMSE, endpoint voltage differences, and normalization. |
| `file_angle_manifest.csv` | Filenames, assigned angles, and filename-angle checks. |
| `calibration_points.csv` | Recording means, within-record SD, pressure, and offset-corrected voltage. |
| `experimental_summary.csv` | Trial counts, normalized means, and SEM by interface scale. |
| `figures/` | Figures 5–6 in editable FIG, PNG, and vector PDF formats. |

In `experimental_summary.csv`, `raw_n` counts selected trials included in the response plots. `normalized_n` counts trials with an eligible same-day reference.

## Redrawing Figures

To redraw the default results without loading recordings:

```matlab
hydrostatic_testing_v19(Mode="render");
```

For another output directory:

```matlab
hydrostatic_testing_v19(Mode="render", OutputDir=reductionDir);
```

To render tables directly without the manuscript-summary check:

```matlab
renderExperimentalFigures(reductionDir);
```

## Validation

Run the synthetic reduction tests with:

```matlab
testExperimentalReduction;
```

These tests cover recording ingestion, trial selection, filename-angle checks, sample SD, missing references, alternative selection, and rejection of incomplete trials. Check the console output for any explicitly skipped optional comparisons.

Check a reduction against the manuscript's Table II values with:

```matlab
verifyManuscriptSummary(reductionDir);
```

This check compares trial counts and reported mean/SEM values at manuscript precision. It does not establish full equality of individual trial results or independently validate the underlying measurements.

The manuscript analysis contains:

- 50 response-plot trials: 33 membrane and 17 bare-port.
- 26 membrane trials with eligible same-day references.
- 13 bare-port trials supplying references for those membrane trials.
- Seven membrane trials without same-day references: trials 3 and 5–10.

Trial 83 contains a filename-angle discrepancy: the file labeled `An 205` occupies the position assigned to 195 degrees. The workflow flags this discrepancy and preserves the recording-order assignment.

## Connection to the Strain Model

This workflow produces experimental Figures 5–6 and the normalized sensitivity summaries used by the model comparison.

`strainModel_v40.m` generates Figures 7–9 using embedded experimental summaries. Running the hydrostatic workflow does not automatically update those values.

Before transferring revised summaries into v40, compare interface scales, trial counts, means, and SEM, and investigate any differences.

## Repository Practices

Keep raw recordings, the trial register, and machine-specific configuration outside version control. The included Git exclusions cover the default results directory and designated local data and configuration paths.

Git exclusions do not remove files that are already tracked. Review the staged file list before committing, particularly when using a custom output directory.