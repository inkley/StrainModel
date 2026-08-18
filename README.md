# StrainModel

MATLAB model for strain-dependent pressure transmission through the membrane-covered working-fluid cavities of a hydrodynamic sensor module.

## Current Version

The authoritative implementation is [`strainModel_v37.m`](strainModel_v37.m). Version 37 supports the accompanying sensor-module manuscript and supersedes the earlier model iterations retained in the repository history.

The model is quasi-static and evaluates how membrane installation, membrane mechanics, hydrostatic loading, cavity-volume change, and isothermal compression of the working fluid affect differential-pressure transmission. It is intended to identify which observed trends are explained by the modeled mechanics and which remain sensitive to installation and boundary-condition uncertainty.

No model parameters are fitted to the measured strain-response data.

## Physical Model

Each side of the differential sensor is represented as an independently loaded membrane-covered cavity. The model:

1. Defines the installed membrane state from the measured nominal strain and an area-weighted effective strain over the pressure-loaded region.
2. Updates membrane thickness and bending rigidity for the selected strain state.
3. Establishes the common hydrostatic equilibrium at the calibration depth.
4. Applies the external differential pressure symmetrically about that equilibrium.
5. Solves the coupled membrane-deflection, cavity-volume, and working-fluid-pressure equilibrium on both sensing sides.
6. Computes raw cavity transmission and transmission relative to an assumed bare-port reference factor.

The primary calculations use:

- a fixed effective loaded radius of 5.5 mm;
- an unstretched latex thickness of 0.508 mm;
- an effective Young's modulus of 0.6 MPa;
- a Poisson ratio of 0.49;
- a per-side geometric cavity volume of approximately 691 mm^3;
- a calibration depth of 0.3048 m;
- an 800 Pa total external differential applied as +400 Pa and -400 Pa about the hydrostatic state; and
- an assumed bare-port transmission factor of 0.85.

The bare-port factor is an explicit comparison assumption. It is not fitted to the membrane-response measurements and is not presented as a validated capillary-pressure correction.

## Model Cases

Version 37 compares three prescribed membrane-mechanics cases:

1. **Full nominal-strain tension:** retains the ideal installed tension associated with the nominal installation strain and provides the lowest-compliance bound.
2. **Effective-strain locked tension:** retains the installed tension associated with the area-weighted effective strain over the pressure-loaded region and serves as the baseline physical candidate.
3. **Locally relaxed nonlinear response:** removes locally retained installation tension while preserving deformation-induced geometric stiffening and provides the highest-compliance bound.

The script also evaluates membrane structural potential energy at the hydrostatic, high-side, and low-side equilibrium states and verifies load-path independence for the active equilibrium formulation.

## Experimental Comparison

The embedded experimental summary contains four membrane installation states. Experimental sensitivity is normalized by the same-day bare-port response. Nominal strain is retained as measured installation geometry, while the primary model comparison is presented against area-weighted effective strain.

All active model cases remain passive when evaluated using raw pressure transmission. The modeled responses underpredict the measured bare-port-normalized sensitivity and do not reproduce its apparent intermediate-strain maximum.

## Run Modes

Set `runMode` near the beginning of `strainModel_v37.m`:

```matlab
runMode = 'publication';
```

Available modes are:

- `core` - baseline solve, compact summary, and pressure-state table.
- `meeting` - core outputs plus the current mechanics, energy, buckling, and candidate-comparison diagnostics.
- `publication` - manuscript figures and tables only, written to a clean publication-results folder. This is the default.
- `full` - all meeting products plus archived parameter sweeps and legacy diagnostics. This mode is substantially slower and produces many files.

Optional deployment-depth calculations are disabled by default. Set:

```matlab
p.includeKnoDeploymentDepths = true;
```

to include the anticipated 12-15 m Kilo Nalu Observatory depth range. These cases remain subject to the model-validity checks printed by the script.

## Running the Model

Open MATLAB in the repository directory and run:

```matlab
strainModel_v37
```

The default publication run writes its products to:

```text
results_v37_publication/
```

Other run modes use:

```text
results_v37/
```

## Primary Publication Outputs

The publication mode generates the manuscript-facing products, including:

- `01_potential_energy_vs_tension_v37`
- `03_candidate_model_comparison_v37`
- nominal-to-effective strain mapping tables
- experimental effective-strain tables
- model-comparison and potential-energy summary tables

Figures are exported in MATLAB `.fig`, PNG, EPS, and PDF formats when supported by the active export routine. The manuscript currently uses the high-resolution PNG files because they preserve the intended MATLAB text and legend rendering in the available LaTeX toolchain.

## Numerical and Physical Checks

Version 37 includes checks for:

- scalar equilibrium convergence;
- positive remaining cavity volume;
- passive raw pressure transmission;
- pressure-amplitude dependence;
- hydrostatic-depth dependence;
- consistency between direct and staged loading paths;
- structural and working-fluid energy balance; and
- sensitivity to selected membrane, geometry, and pneumatic-volume assumptions in the extended run modes.

## Repository Notes

- `strainModel_v37.m` is the authoritative model version used for the current manuscript.
- Earlier scripts are retained as development history and should not be interpreted as the current formulation.
- Internal variable names inherited from earlier iterations may use legacy terminology, but the v37 outputs and manuscript definitions govern interpretation.
- Unmeasured tubing, sensor, and fitting volumes are treated as model uncertainty rather than absorbed into a fitted scale factor.
- The model is quasi-static; dynamic membrane inertia, damping, external fluid inertia, and viscous fluid-structure interactions are outside its present scope.

## Citation

If you use this repository, please cite the associated hydrodynamic sensor-module manuscript once its final bibliographic information is available.
