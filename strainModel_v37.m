%% strainModel_v37.m
clear; clc; close all;

% v37 FINAL MODEL (37 is Luke Jackson's prisoner number in Cool Hand Luke)
% 1) Adds Mike's explicit load-path validation for incorporation of tension.
%    Direct loading from the sealed/rest reference to the final surface
%    pressure is compared with staged loading through hydrostatic equilibrium.
% 2) Audits both sides of the symmetric differential-pressure load using the
%    same trapped-gas inventory, geometry, and installation pre-tension.
% 3) Confirms that the active equilibrium formulation is path independent;
%    hydrostatic loading changes equilibrium deflection and cavity pressure,
%    but does not create a separate persistent pre-tension state.
% 4) Makes area-weighted effective strain the independent variable for both
%    the empirical dataset and model. Nominal installation strain is retained
%    only as measured-geometry provenance and is not used in the tension law.
% 5) Adds Mike's requested structural-potential-energy presentation as a
%    function of effective installed tension at hydrostatic and differential
%    high/low equilibrium states.
%
% v36 cleanup
% 1) Consolidates the many independent output switches into three run modes:
%    'core', 'meeting' (default), and 'full'.
% 2) Expensive pressure and parameter sweeps now run only when requested.
% 3) Retains the v35 governing equations and numerical parameters.
% 4) Adds the optional area-weighted volume-conserving strain diagnostic for
%    the observed 100%-interface buckling condition. It does not replace the
%    nominal-strain baseline.
%
% v35
% 1) Raises the nominal total differential-pressure perturbation from 5 Pa
%    to 400 Pa, applied symmetrically as +200 Pa and -200 Pa about the
%    hydrostatically loaded equilibrium.
% 2) Adds a differential-amplitude characterization over 5-800 Pa total
%    differential. The 800 Pa endpoint corresponds to +400/-400 Pa at the
%    two ports and therefore covers the maximum stated assembly condition.
% 3) Exports the amplitude characterization to CSV and plots the relaxed
%    nonlinear C_V=1/2 response at each reported strain. No amplitude is fit
%    to the experimental response.
% 4) Adds membrane structural-potential-energy versus maximum-deflection and
%    maximum-deflection versus signed port-pressure diagnostics. The original
%    locked-elastic baseline, locked nonlinear diagnostic, and locally relaxed
%    C_V=1/2 nonlinear bound are all shown separately.
% 5) Extends the meeting pressure-state table with each cavity's signed
%    pressure change from hydrostatic equilibrium, individual transmission,
%    high/low symmetry ratio, and pressure-change imbalance.
%
% Retained from v34:
% 1) Adds a conservative reduced-order nonlinear membrane diagnostic to test
%    whether hydrostatically stored membrane energy can be released during
%    differential unloading and materially enhance pressure transmission.
% 2) The diagnostic uses q = A*w + B*w^3, explicit structural and isothermal
%    gas energies, and symmetric loaded/unloaded equilibrium states. It is
%    diagnostic only and does not replace the verified v33 baseline.
% 3) Adds a parameter-free local-prestress relaxation/redistribution bound:
%    installed residual tension is set to zero while hydrostatic deformation-
%    induced geometric tension remains active. This is an upper-compliance
%    limiting case, not a fitted or selected baseline.
%
% Retained from v33:
% 1) Adds a meeting-ready pressure-state audit at zero strain and each
%    experimental strain. The audit reports cavity pressure after the
%    hydrostatic preload but before the differential, then reports both
%    cavity pressures after the symmetric differential is applied.
% 2) Reports hydrostatic pressure loss, achieved internal differential, and
%    differential transmission, and exports the same values to CSV.
% 3) No governing equations, parameters, or baseline physics are changed
%    relative to v32.
% 4) The mechanics diagnostic additionally reports normalized effective
%    transverse stiffness, Ktrans = q/wmax, from the combined bending- and
%    pre-tension-dependent center-deflection compliance.
%
% FUTURE MODEL EXTENSION (not active in v33): allow membrane tension and
% stored strain energy to evolve between the hydrostatic and differential
% loading states. This requires a nonlinear, state-dependent energy model
% rather than absorbing tension into Kbend.
% 5) Optionally extends the hydrostatic-depth audit to the anticipated 12-15 m
%    Kilo Nalu Observatory deployment range. Set
%    p.includeKnoDeploymentDepths = true to enable these cases. Deep-water
%    results remain subject to the printed validity audit.
%
% Retained from v32:
% 1) Installation pre-tension is restored to the baseline tensioned-plate
%    model. The ideal equibiaxial tensile resultant is multiplied by the
%    explicit realization factor kT0 before entering the plate equation as
%    Tpre_eff. The baseline uses kT0 = 1; sensitivity cases can still vary it.
% 2) The differential surface-pressure perturbation is applied symmetrically
%    about the hydrostatically loaded equilibrium: +deltaP/2 on one funnel
%    and -deltaP/2 on the other.
%
% Retained from v31
% This manuscript-refinement release removes quantitative model/experiment
% overlays from parameter-sensitivity plots, adds w/h to the validity audit,
% and documents the CAD-derived geometric cavity volume. It introduces no
% new membrane physics relative to v29.
%
% Retained from v29:
% This housekeeping release makes the model/experiment reference distinction
% explicit and audits passive pressure transmission. It introduces no new
% membrane physics relative to v28.
%
% New in v29:
% 1) Raw R_cav, zero-strain-relative modeled transmission, and experimental
%    bare-port-normalized sensitivity are plotted and tabulated separately.
% 2) Passivity is checked over the full strain range at every calibration-
%    depth case: 0 <= R_cav <= 1.
% 3) Active output terminology reserves "bare-port-normalized sensitivity
%    enhancement" for the experiment and uses "relative transmission
%    amplification" or "strain-dependent transmission enhancement" for the
%    model. No bare-port/zero-strain-membrane equivalence is assumed.
%
% Retained from v28:
% This version refines the two-stage hydrostatic preload formulation using
% the estimated calibration depth and adds auditable normalization,
% statistical, loaded-radius, and pneumatic-volume diagnostics.
%
% New in v28:
% 1) Baseline calibration depth is 0.3048 m (12 in), with 0.20-0.40 m
%    model-sensitivity cases plus the zero-depth regression case.
% 2) Experimental/model comparisons are shown using both bare-port and
%    90%-membrane references. No fitted scale bridges the two references.
% 3) Aggregate-data peak diagnostics compare constant, linear, and quadratic
%    weighted fits and report pairwise combined-uncertainty separations.
% 4) A bounded +/-5% effective loaded-radius diagnostic tests whether a
%    physically modest, state-dependent radius change can match each point.
% 5) Pneumatic volume is separated into cavity, tubing, sensor, and fitting
%    contributions, with an added-volume sensitivity diagnostic.
%
% Retained from v27:
% The model uses an explicit two-stage hydrostatic preload formulation:
% (1) solve the sealed cavity/membrane equilibrium at the prescribed external
%     static pressure while retaining the original trapped-gas inventory;
% (2) apply the small differential perturbation about that loaded equilibrium.
% This distinguishes gas fill pressure from external pressure at depth.
%
% The analytical tensioned-plate compliance is used by the v32 baseline.
% Earlier zero-local-tension behavior remains available only as a diagnostic
% through p.preTensionMode = 'none' or p.kT0 = 0.
%
% Introduced in v27:
% 1) Separate trapped-gas reference pressure P_gas0 from external static
%    pressure P_static.  The gas invariant is P_gas0*V00.
% 2) Hydrostatic equilibrium is solved before the small-signal response.
% 3) Static membrane deflection, cavity compression, membrane strain energy,
%    and reversible isothermal gas-compression work are reported by depth.
% 4) At zero depth (P_static = P_gas0), v27 reduces to the v26 formulation.
%
% Retained from v26:
% 1) Exact Bessel-function compliance available for the generalized
%    tensioned-plate formulation.
% 2) Baseline formulation uses pre-strain through h(eps0) and Kbend(eps0),
%    without applying installation pre-tension as a local tensile-restraint
%    term in the plate compliance.
% 3) Solver mode toggle:
%       p.solverMode = 'root';        % scalar pressure-equilibrium solve
%       p.solverMode = 'fixedPoint';  % original relaxed fixed-point solver
% 4) Volume-cap toggle:
%       p.useVolumeCap = true;        % tanh-capped dV regularization
%       p.useVolumeCap = false;       % uncapped dV = dV_raw
% 5) Diagnostic model comparisons retained for direct, capped, and
%    spatially attenuated local pre-tension cases.
% 6) Calibration pressure-context and static-depth diagnostics retained.
%
% v32 baseline assumptions / model formulations:
% 1) axisymmetric membrane/interface response
% 2) small-signal linear circular-plate compliance about the installed state
% 3) installation pre-strain affects membrane thickness h(eps0)
% 4) bending stiffness Kbend(eps0) is evaluated from the strain-updated thickness
% 5) the ideal equibiaxial installation tensile resultant enters the local
%    plate compliance through the explicit realization factor kT0
% 6) analytical radial deflection profile used consistently for center
%    deflection and numerical cavity-volume integration
% 7) effective circular loaded region with fixed loaded radius
% 8) cavity response solved using isothermal pressure-volume compression
% 9) unresolved installation, seating, and boundary-condition effects are
%    interpreted through effective stiffness behavior

%% -----------------------------
% Results folder
% ------------------------------
versionTag = 'v37';
resultsDir = fullfile(pwd, ['results_', versionTag]);
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

%% -----------------------------
% Run mode and optional outputs
% ------------------------------
% 'core'        : baseline solve, compact summary, and pressure-state table
% 'meeting'     : core + current mechanics, energy, and buckling diagnostics
% 'publication' : final manuscript figures/tables only, in a clean folder
% 'full'        : meeting + archived parameter sweeps and legacy diagnostics
runMode = 'publication';

run = struct();
run.mainFigure = false;
run.transmissionPlot = false;
run.transmissionAuditFigure = false;
run.mechanicsPlot = false;
run.geometryPlot = false;
run.dualNormalizationPlot = false;
run.staticDepthPlots = false;
run.modelClosurePlot = false;
run.potentialEnergy = false;
run.extendedEnergyGauntlet = false;
run.bucklingWeightedStrain = false;
run.candidateModelComparison = false;
run.pressureLevelSweep = false;
run.parameterSweeps = false;
run.legacyDiagnostics = false;
run.printConfiguration = false;
run.printModelSummary = true;
run.printHydrostaticValidation = false;
run.printLoadPathValidation = true;
run.printPeakStatistics = false;
run.printTransmissionAudit = true;
run.printMikePressureReport = true;
run.printBaselineNormalization = true;
run.publicationOnly = false;

switch lower(runMode)
    case 'core'
        % Defaults above intentionally produce a short, auditable run.

    case 'meeting'
        run.transmissionAuditFigure = true;
        run.mechanicsPlot = true;
        run.potentialEnergy = true;
        run.bucklingWeightedStrain = true;
        run.candidateModelComparison = true;
        run.printHydrostaticValidation = true;

    case 'publication'
        run.transmissionAuditFigure = true;
        run.potentialEnergy = true;
        run.candidateModelComparison = true;
        run.publicationOnly = true;
        run.printModelSummary = false;
        run.printLoadPathValidation = false;
        run.printTransmissionAudit = false;
        run.printMikePressureReport = false;
        run.printBaselineNormalization = false;

        % Keep publication products isolated from meeting/debug artifacts.
        resultsDir = fullfile(pwd, ['results_', versionTag, '_publication']);
        if ~exist(resultsDir, 'dir')
            mkdir(resultsDir);
        end

    case 'full'
        names = fieldnames(run);
        for iRun = 1:numel(names)
            run.(names{iRun}) = true;
        end

    otherwise
        error(['Unknown runMode: %s. Use core, meeting, publication, ', ...
            'or full.'], runMode);
end

% Legacy aliases keep the well-tested v35 plotting functions isolated from
% the cleaner run-mode interface above.
makeMainFigure = run.mainFigure;
makeTransmissionPlot = run.transmissionPlot;
makeTransmissionAuditFigure = run.transmissionAuditFigure;
makeMechanicsPlot = run.mechanicsPlot;
makeGeometryPlot = run.geometryPlot;
makeDualNormalizationPlot = run.dualNormalizationPlot;
makeStaticDepthDiagnosticPlots = run.staticDepthPlots;
makeModelClosureDiagnosticPlot = run.modelClosurePlot;
runPotentialEnergyReleaseTests = run.potentialEnergy;
runBucklingWeightedStrainTest = run.bucklingWeightedStrain;
printConfigSummary = run.printConfiguration;
printModelPointSummary = run.printModelSummary;
printHydrostaticValidation = run.printHydrostaticValidation;
printLoadPathValidation = run.printLoadPathValidation;
printPeakDiagnostics = run.printPeakStatistics;
printTransmissionAudit = run.printTransmissionAudit;
printMikePressureReport = run.printMikePressureReport;

makeMaterialFitOverlayPlots = run.legacyDiagnostics;
makeKT0FitOverlayPlot = run.legacyDiagnostics;
makeKT0IntermediateDiagnostics = run.legacyDiagnostics;
makeKT0RawDebugPlots = run.legacyDiagnostics;
makeSolverComparisonPlots = run.legacyDiagnostics;
makeThicknessDiagnosticPlots = run.legacyDiagnostics;
makePreTensionThicknessDiagnosticPlots = run.legacyDiagnostics;
makeKT0ThicknessReferenceFitSweep = run.legacyDiagnostics;
makeLoadedRadiusDiagnosticPlot = run.legacyDiagnostics;
makePneumaticVolumeDiagnosticPlot = run.legacyDiagnostics;
makeSpatialPreTensionDiagnosticPlots = run.legacyDiagnostics;
makeLocalTensionScaleDiagnosticPlots = run.legacyDiagnostics;
makeTensionCapDiagnosticPlots = run.legacyDiagnostics;
makeModelComparisonDiagnosticPlots = run.legacyDiagnostics;
makePressureSweepPlot = run.pressureLevelSweep;
printPressureSweepSummary = run.pressureLevelSweep;
printSweepSummary = run.parameterSweeps;
printCalibrationPressureContext = run.legacyDiagnostics;

%% -----------------------------
% Experimental data / interface states
% ------------------------------
d0 = 21.978; % mm, funnel diameter at zero installation strain
effectiveLoadedRadius_mm = 5.5; % mm, pressure-loaded center region

scaleDiameters_raw = [16.484, 17.582, 18.681, 19.780, 21.978]; % mm
nominalStrains_raw = (d0 - scaleDiameters_raw) ./ scaleDiameters_raw; % (-)

[nominalStrains, sortIdx] = sort(nominalStrains_raw);
scaleDiameters = scaleDiameters_raw(sortIdx);

% Convert the measured global installation strain to the effective strain
% carried by the annulus outside the pressure-loaded center. Incompressible
% thickness stretch is area averaged, then converted back to equibiaxial
% in-plane stretch. This mapping is applied exactly once.
effectiveAreaFraction = (effectiveLoadedRadius_mm / (0.5*d0))^2;
lambdaZCenter = (1 + nominalStrains).^(-2);
lambdaZEffective = (1-effectiveAreaFraction) + ...
    effectiveAreaFraction.*lambdaZCenter;
effectiveStrains = lambdaZEffective.^(-1/2) - 1;

eps_data_nominal = [0.111, 0.176, 0.250, 0.333];
lambdaZDataCenter = (1 + eps_data_nominal).^(-2);
lambdaZDataEffective = (1-effectiveAreaFraction) + ...
    effectiveAreaFraction.*lambdaZDataCenter;
eps_data = lambdaZDataEffective.^(-1/2) - 1;
y_data   = [1.167, 1.328, 1.270, 1.224];
y_err    = [0.093, 0.082, 0.045, 0.042];

% Export the coordinate change used to present both geometry states and the
% empirical response. These tables make the nominal-to-effective mapping
% auditable without allowing nominal strain into the baseline mechanics.
geometryStrainMap = table(scaleDiameters(:), nominalStrains(:), ...
    effectiveStrains(:), 'VariableNames', ...
    {'scaleDiameter_mm','nominalStrain','effectiveStrain'});
writetable(geometryStrainMap, fullfile(resultsDir, ...
    ['geometry_strain_mapping_',versionTag,'.csv']));
experimentalStrainMap = table(eps_data_nominal(:), eps_data(:), ...
    y_data(:), y_err(:), 'VariableNames', ...
    {'nominalStrain','effectiveStrain','response','uncertainty'});
writetable(experimentalStrainMap, fullfile(resultsDir, ...
    ['experimental_effective_strain_',versionTag,'.csv']));

% The experimental values above are normalized to the bare-port condition.
% A second comparison uses the 90% membrane (eps = 0.111) as its reference,
% matching the model's membrane-to-membrane comparison without assuming that
% a zero-strain membrane is equivalent to the bare port.
membraneRefIdx = 1;
y_data_mem90 = y_data ./ y_data(membraneRefIdx);
y_err_mem90 = y_data_mem90 .* sqrt((y_err ./ y_data).^2 + ...
    (y_err(membraneRefIdx) / y_data(membraneRefIdx))^2);
y_err_mem90(membraneRefIdx) = 0.0;

%% -----------------------------
% Plot styling
% ------------------------------
modelColor      = [0.00, 0.20, 0.60];
physColor       = [0.20, 0.55, 0.20];
geomColor       = [0.50, 0.10, 0.55];
bareColor       = [0.35, 0.35, 0.35];
c75             = [0.00, 0.4470, 0.7410];
c80             = [0.8500, 0.3250, 0.0980];
c85             = [0.9290, 0.6940, 0.1250];
c90             = [0.4940, 0.1840, 0.5560];

markerSize      = 9;
modelLineWidth  = 3.0;
errorLineWidth  = 1.8;

%% -----------------------------
% Model parameters
% ------------------------------
p = struct();
p.publicationMode = run.publicationOnly;

% Reference / perturbation pressures
p.P_ref = 101325; % Pa
% Nominal total port-to-port differential. Because the model applies the
% load symmetrically, 400 Pa means +200 Pa on one port and -200 Pa on the
% other. The amplitude audit also includes 800 Pa total (+400/-400 Pa).
p.dP0   = 800.0;  % Pa total differential

% Hydrostatic preload state.  P_gas0 and V00 define the gas inventory at
% sealing. P_static is the common external pressure before the differential
% perturbation is applied. baselineDepth_m = 0 reproduces v26.
p.P_gas0 = p.P_ref;  % Pa, cavity gas pressure when sealed at volume V00
% Polytropic gas exponent. n=1 is the quasi-static isothermal baseline;
% n=1.4 is retained only as a rapid/adiabatic upper-bound diagnostic.
p.gasExponent = 1.0;
p.rho_water = 1025;  % kg/m^3
p.g_water = 9.81;    % m/s^2
p.baselineDepth_m = 0.3048; % m, estimated membrane-center depth (12 in)
p.depthSensitivityRange_m = [0.20, 0.40]; % m, model-sensitivity cases
p.includeKnoDeploymentDepths = false; % true adds the optional KNO cases below
p.knoDeploymentDepthRange_m = [12.0, 15.0]; % m, anticipated KNO range
p.depth_sweep_m = [0, 0.20, 0.25, 0.3048, 0.35, 0.40];
if p.includeKnoDeploymentDepths
    p.deploymentDepthRange_m = p.knoDeploymentDepthRange_m;
    p.depth_sweep_m = [p.depth_sweep_m, p.deploymentDepthRange_m];
else
    p.deploymentDepthRange_m = [];
end
p.P_static = p.P_ref + p.rho_water * p.g_water * p.baselineDepth_m;

% Measured geometry states
p.d0_mm             = d0;
p.scaleDiameters_mm = scaleDiameters;
p.nominalEngStrains = nominalStrains;
p.engStrains        = effectiveStrains;
p.experimentalNominalStrains = eps_data_nominal;
p.experimentalEffectiveStrains = eps_data;

% Pre-tension model options: 'ideal', 'activated', 'none'. The v32 baseline
% restores the installed pre-tension directly; the alternatives are retained
% only for sensitivity checks.
p.preTensionMode = 'ideal';

% Constitutive options for the installed tensile resultant:
% 'linearEquibiaxial' : plane-stress linear-elastic approximation
% 'neoHookean'        : incompressible neo-Hookean equibiaxial membrane
p.preTensionLaw = 'linearEquibiaxial';

% Pre-tension spatial coupling options:
% 'direct'        : current behavior, uses p.kT0 directly
% 'radiusPower'   : global engagement attenuated by (a_load/a_geom)^n
% 'tensionCap'    : computes ideal local tension, then caps it at Tpre_cap_Npm
p.preTensionCouplingMode = 'direct';
p.kT0_global = 0.75;          % realistic global installation engagement
p.radiusPowerExponent = 6;    % spatial attenuation exponent

% Pre-tension cap options
p.Tpre_cap_Npm = Inf;   % N/m; only used for 'tensionCap'

%% -----------------------------
% Cavity geometry, per sensing side
% ------------------------------
h1 = 2.506; % mm
h2 = 6.617; % mm
r  = effectiveLoadedRadius_mm; % mm

p.r_forced_m = r * 1e-3;
p.A_forced   = pi * p.r_forced_m^2;

% Effective strain is now supplied directly to all baseline model equations.
% The areaWeightedVolume option remains only for the nominal-to-effective
% traceability audit and must not be applied to baseline effective strain.
p.strainMappingMode = 'effectiveDirect';
p.strainWeightingOuterRadius_m = 0.5*p.d0_mm*1e-3;
p.strainWeightingCenterRadius_m = p.r_forced_m;

% Effective loaded-radius formulation. The baseline remains fixed.
% 'fixed'       : unsupported pressure-loaded throat radius
% 'stateTable'  : prescribed diagnostic scale applied to the throat radius
p.loadedRadiusMode = 'fixed';
p.radiusStateStrains = [0, eps_data, max(effectiveStrains)];
p.radiusStateScales = ones(size(p.radiusStateStrains));
p.radiusDiagnosticBounds = [0.95, 1.05];

p.V11 = pi * r^2 * h1;
p.V22 = pi * h2^2 * r - pi * h2^3 / 3;
p.V_cavity_geom = (p.V11 + p.V22) * 1e-9; % m^3, per-side cavity volume from funnel CAD solid model
p.V_tubing = 0.0;   % m^3, placeholder pending measurement
p.V_sensor = 0.0;   % m^3, placeholder pending manufacturer/measurement data
p.V_fittings = 0.0; % m^3, placeholder pending measurement
p.V00 = p.V_cavity_geom + p.V_tubing + p.V_sensor + p.V_fittings;

%% -----------------------------
% Plate / membrane model
% ------------------------------
p.t_plate0 = 0.020 * 0.0254; % m, 0.020 in latex membrane
p.nu_plate = 0.49;
p.E_plate0 = 6.0e5;          % Pa

% Fraction of the ideal equibiaxial installation tensile resultant entering
% the local tensioned-plate equation. The physical baseline uses the full
% ideal resultant; this is an assumption, not a fitted parameter.
p.kT0 = 1.0;
p.eps_char = 0.160;

% Assumed bare-port pressure-transmission ratio. This represents
% attenuation caused by capillary resistance and trapped-air compression
% in the submerged vinyl-tubing pathway.
p.R_bare = 0.85;

% Optional strain-dependent modulus.
p.useStrainDependentE = false;
p.c1 = 0.0;
p.c2 = 0.0;

% Thickness update options:
% 'poisson'        : current Poisson-ratio-dependent update
% 'incompressible' : t = t0/lambda^2
% 'constant'       : t = t0
p.thicknessMode = 'poisson';

% Pre-tension thickness reference options:
% 'current'   : Tpre uses strain-updated thickness, t(eps)
% 'reference' : Tpre uses original thickness, t0
p.preTensionThicknessMode = 'current';

%% -----------------------------
% Solver and regularization settings
% ------------------------------
p.tolP    = 1e-5;
p.maxIter = 1000;
p.relax   = 0.08;

% Solver options: 'root' or 'fixedPoint'
p.solverMode = 'root';

% Volume regularization options.
% true  -> dV = dV_cap*tanh(dV_raw/dV_cap)
% false -> dV = dV_raw
p.useVolumeCap    = false;
p.dV_cap_fraction = 0.10;
p.minVolumeFraction = 1e-6;

%% -----------------------------
% Strain range and parameter sweeps
% ------------------------------
eps_plot = linspace(0, 1.08*max(effectiveStrains), 150);

% Verify trapezoidal integration of the zero-tension radial profile.
rho_check = linspace(0, 1, 1001);
shape_check = (1 - rho_check.^2).^2;

volume_factor_numeric = ...
    2 * trapz(rho_check, shape_check .* rho_check);

volume_factor_exact = 1/3;

volume_factor_error = abs( ...
    volume_factor_numeric - volume_factor_exact);

if volume_factor_error > 1e-6
    error(['Zero-tension volume integration check failed: ', ...
           'numeric factor = %.12f, expected factor = %.12f.'], ...
        volume_factor_numeric, volume_factor_exact);
end

kT0_sweep       = [0, 0.0075, 0.02, 0.10, 0.25, 0.50, 1.00];
kT0_debug_sweep = [0, 0.10, 0.25, 0.50, 0.75, 1.00];

kT0_grid = [0, 0.0025, 0.005, 0.0075, 0.010, 0.015, 0.020, 0.025, 0.050, 0.075, 0.100];

E_sweep     = [1.0, 1.25, 1.5, 2.0] * p.E_plate0; % Pa
nu_sweep    = [0.45, 0.47, 0.49];                 % (-)
t_sweep_mm  = [0.508, 0.65, 0.75, 0.90, 1.00];   % mm
t_sweep_m   = t_sweep_mm * 1e-3;                  % m

%% -----------------------------
% Pressure-level sweep settings
% ------------------------------
eps_pressure_cases = eps_data;
P3_offset_sweep = [-400, -200, 0, 200, 400]; % Pa relative to p.P_static
P3_sweep = p.P_static + P3_offset_sweep;
p.dP_local = p.dP0;

%% -----------------------------
% Static pressure / depth diagnostic settings
% ------------------------------
rho_water       = p.rho_water;                  % kg/m^3, seawater
g_water         = p.g_water;                    % m/s^2
depth_sweep_m   = p.depth_sweep_m;               % m
Pstatic_sweep   = p.P_ref + rho_water * g_water * depth_sweep_m; % Pa

%% -----------------------------
% Calibration pressure context
% ------------------------------
% Set this to the largest expected static water-height offset during calibration.
% This is not the oscillatory/small-signal pressure; it is the static bias.
delta_h_cal_m = 0.010; % m, example = 10 mm water head

Pstatic_cal_offset = rho_water * g_water * delta_h_cal_m;
depth_equiv_cal_m  = Pstatic_cal_offset / (rho_water * g_water);

if printCalibrationPressureContext
    fprintf('\n--- calibration pressure context ---\n');
    fprintf('Assumed calibration static head offset: %.3f mm\n', 1e3 * delta_h_cal_m);
    fprintf('Equivalent static pressure bias: %.4f Pa\n', Pstatic_cal_offset);
    fprintf('Equivalent water depth: %.5f m\n', depth_equiv_cal_m);
    fprintf('Baseline membrane-center depth: %.4f m.\n\n', p.baselineDepth_m);
end

%% -----------------------------
% Print configuration summary
% ------------------------------
if printConfigSummary
    printConfigurationSummary(p, scaleDiameters, nominalStrains, effectiveStrains);
end

%% -----------------------------
% Evaluate baseline model over strain range
% ------------------------------
[T_cav, Dplate_curve, E_curve, Tpre_eff_curve, phi_pre_curve, ...
 a_load_curve, a_geom_curve, D_installed_curve_mm, eps_pre_curve, ...
 converged_curve, iter_curve] = runModelOverStrainRange(eps_plot, p);

Tref = interp1(eps_plot, T_cav, 0.00, 'linear');
if abs(Tref) < 1e-12 || ~isfinite(Tref)
    error(['Zero-strain membrane transmission is too small or invalid. ', ...
           'Inspect raw T_cav first.']);
end

% Secondary diagnostic: transmission relative to the zero-strain membrane.
T_cav_rel0 = T_cav ./ Tref;

% Primary comparison: transmission relative to the assumed bare port.
if ~isfinite(p.R_bare) || p.R_bare <= 0 || p.R_bare > 1
    error('p.R_bare must be finite and satisfy 0 < R_bare <= 1.');
end

R_cav_rel_bare = T_cav ./ p.R_bare;
R_cav_rel_bare(converged_curve < 0.5) = NaN;

S_model = R_cav_rel_bare;
S_interface_model = R_cav_rel_bare;

% Meeting-requested pressure states: common hydrostatic equilibrium followed
% by the two cavity equilibria under the symmetric differential loading.
if printMikePressureReport
    mikePressureStrains = [0, eps_data];
    compileMikePressureReport(mikePressureStrains, p, resultsDir, versionTag);
end

% Compact KNO report is evaluated independently of the full depth-curve
% diagnostic so the deployment cases remain quick to inspect and export.
if p.includeKnoDeploymentDepths
    compileKnoDepthReport([0, eps_data], p, resultsDir, versionTag);
end

if runPotentialEnergyReleaseTests
    runPotentialEnergyReleaseDiagnostic([0, eps_data], p, resultsDir, ...
        versionTag, run.extendedEnergyGauntlet, run.publicationOnly);
end

if runBucklingWeightedStrainTest
    runBucklingWeightedStrainDiagnostic( ...
        [0, p.experimentalNominalStrains], p, resultsDir, versionTag);
end

if run.candidateModelComparison
    plotCandidateModelComparison(eps_plot, eps_data, y_data, y_err, ...
        p, resultsDir, versionTag);
end

if run.printBaselineNormalization
    fprintf('\n--- assumed bare-port normalization ---\n');
    fprintf('Assumed R_bare = %.3f\n', p.R_bare);
    fprintf('eps       R_cav     R_cav/R_bare\n');

    for i = 1:numel(eps_data)
        R_cav_i = interp1(eps_plot, T_cav, eps_data(i), 'linear');
        R_rel_i = interp1(eps_plot, R_cav_rel_bare, eps_data(i), 'linear');

        fprintf('%6.3f    %7.4f       %7.4f\n', ...
            eps_data(i), R_cav_i, R_rel_i);
    end
end

Tref_mem90 = interp1(eps_plot, T_cav, eps_data(membraneRefIdx), 'linear');
if abs(Tref_mem90) < 1e-12 || ~isfinite(Tref_mem90)
    error('90%%-membrane reference transmission is too small or invalid.');
end
T_cav_rel_mem90 = T_cav ./ Tref_mem90;

%% -----------------------------
% Print model summary at experimental points
% ------------------------------
[resid, wres, rmse_base, wrmse_base] = computePointResiduals(eps_plot, S_interface_model, eps_data, y_data, y_err);

if printModelPointSummary
    printModelPointSummaryTable(eps_plot, eps_data, y_data, y_err, S_interface_model, T_cav_rel0, ...
        a_load_curve, a_geom_curve, Dplate_curve, Tpre_eff_curve, phi_pre_curve, resid, wres);
    fprintf('\n RMSE  = %.4f\n', rmse_base);
    fprintf('WRMSE = %.4f\n', wrmse_base);
    fprintf('converged points: %d / %d\n', sum(converged_curve > 0.5), numel(converged_curve));
    fprintf('max solver iterations: %d\n', max(iter_curve));
    fprintf('mean solver iterations: %.1f\n', mean(iter_curve));
end

if printHydrostaticValidation
    runHydrostaticPreloadValidation(p);
end

if printLoadPathValidation
    runLoadPathIndependenceValidation([0, eps_data], p, resultsDir, versionTag);
end

if printTransmissionAudit
    printTransmissionReferenceAudit( ...
        eps_plot, T_cav, R_cav_rel_bare, ...
        eps_data, y_data, y_err, p);
end

if printPeakDiagnostics
    printAggregatePeakDiagnostics(eps_data, y_data, y_err, y_data_mem90);
end

%% -----------------------------
% Optional pressure-level dependence check
% ------------------------------
if run.pressureLevelSweep
    [pressureGain, pressureGain_relLocal, pressureGain_relZeroStrain, ...
     S_pressure, pressureConverged] = runPressureLevelSweep( ...
        eps_pressure_cases, P3_sweep, P3_offset_sweep, p);

if makePressureSweepPlot
    figure; hold on;
    for ie = 1:numel(eps_pressure_cases)
        plot(P3_offset_sweep, pressureGain_relLocal(ie, :), '-o', ...
            'LineWidth', 2.2, 'MarkerSize', 7, ...
            'DisplayName', sprintf('\\epsilon = %.3f', eps_pressure_cases(ie)));
    end
    yline(1.0, ':k', 'HandleVisibility', 'off');
    xlabel('P_3 - P_{static} (Pa)');
    ylabel('Relative Local Gain (-)');
    title('Pressure-Level Dependence at Fixed Effective Strain');
    legend('Location', 'northwest');
    grid on; box on; formatAxes(gca);
    saveCurrentFigure(resultsDir, ['07_pressure_level_dependence_', versionTag]);
end

if printPressureSweepSummary
    printPressureSweepSummaryTable(eps_pressure_cases, P3_offset_sweep, pressureGain, ...
        pressureGain_relLocal, pressureGain_relZeroStrain, S_pressure, pressureConverged);
end
end

%% -----------------------------
% Main response figures
% ------------------------------
if makeMainFigure
    figure; hold on;
    hBareMean = plotBareReference(bareColor);
    [h75, h80, h85, h90] = plotExperimentalData(eps_data, y_data, y_err, c75, c80, c85, c90, markerSize, errorLineWidth);
    hModel = plot(eps_plot, S_model, '-', 'Color', modelColor, 'LineWidth', modelLineWidth, ...
        'DisplayName', 'Modeled Relative Transmission, h_0 = 0.508 mm');
    xlabel('Effective Engineering Strain, \epsilon_{eff} (-)');
    ylabel('Normalized Response (-)');
    title('Experimental Sensitivity and Modeled Relative Transmission');
    xlim([eps_plot(1) eps_plot(end)]); ylim([0.85 1.5]);
    grid on; box on; formatAxes(gca);
    legend([h75, h80, h85, h90, hModel, hBareMean], ...
        {'75%', '80%', '85%', '90%', 'Modeled Relative Transmission', 'Experimental Bare-Port Reference'}, ...
        'Location', 'northwest');
    saveCurrentFigure(resultsDir, ['01_main_normalized_sensitivity_', versionTag]);
end

if makeTransmissionPlot
    figure; hold on;
    plot(eps_plot, T_cav_rel0, '-', 'Color', physColor, 'LineWidth', 2.5);
    yline(1.0, '-', 'Color', bareColor, 'LineWidth', 1.8);
    xlabel('Effective Engineering Strain, \epsilon_{eff} (-)');
    ylabel('Zero-Strain-Relative Transmission (-)');
    title('Modeled Relative Transmission Enhancement');
    legend('R_{cav}(\epsilon)/R_{cav}(0)', 'Unity Relative Reference', 'Location', 'northwest');
    grid on; box on; formatAxes(gca);
    xlim([eps_plot(1) eps_plot(end)]);
    saveCurrentFigure(resultsDir, ['02_transmission_rel0_', versionTag]);
end


if makeTransmissionAuditFigure
    plotTransmissionReferenceAudit( ...
        eps_plot, T_cav, R_cav_rel_bare, ...
        eps_data, y_data, y_err, p.R_bare, ...
        modelColor, physColor, bareColor, ...
        resultsDir, versionTag, run.publicationOnly);
end


if makeDualNormalizationPlot
    plotDualNormalizationComparison(eps_plot, T_cav_rel0, T_cav_rel_mem90, ...
        eps_data, y_data, y_err, y_data_mem90, y_err_mem90, ...
        modelColor, physColor, bareColor, resultsDir, versionTag);
end

if makeMechanicsPlot
    plotMechanicsDiagnostic(eps_plot, Dplate_curve, Tpre_eff_curve, ...
        a_load_curve, p, resultsDir, versionTag);
end

if makeThicknessDiagnosticPlots
    plotThicknessUpdateDiagnostics(eps_plot, p, resultsDir, versionTag);
end

if makeSpatialPreTensionDiagnosticPlots
    plotSpatialPreTensionDiagnostic(eps_plot, eps_data, y_data, y_err, p, ...
        c75, c80, c85, c90, markerSize, errorLineWidth, bareColor, resultsDir, versionTag);
end

if makeLocalTensionScaleDiagnosticPlots
    plotLocalTensionScaleDiagnostic(eps_plot, eps_data, y_data, y_err, p, ...
        c75, c80, c85, c90, markerSize, errorLineWidth, ...
        bareColor, resultsDir, versionTag);
end

if makeTensionCapDiagnosticPlots
    plotTensionCapDiagnostic(eps_plot, eps_data, y_data, y_err, p, ...
        c75, c80, c85, c90, markerSize, errorLineWidth, ...
        bareColor, resultsDir, versionTag);
end

if makeModelComparisonDiagnosticPlots
    plotModelComparisonDiagnostic(eps_plot, eps_data, y_data, y_err, p, ...
        c75, c80, c85, c90, markerSize, errorLineWidth, ...
        bareColor, resultsDir, versionTag);
end

if makeModelClosureDiagnosticPlot
    plotModelClosureDiagnostic(eps_plot, eps_data, y_data, y_err, ...
        p, bareColor, resultsDir, versionTag);
end

if makeGeometryPlot
    figure; hold on;
    plot(eps_plot, 1e3 * a_geom_curve, '-', 'Color', geomColor, 'LineWidth', 2.5, 'DisplayName', 'a_{installed}(\epsilon)');
    plot(eps_plot, 1e3 * a_load_curve, '--', 'Color', modelColor, 'LineWidth', 2.5, 'DisplayName', 'a_{load}(\epsilon)');
    yline(1e3 * p.r_forced_m, ':', 'Color', bareColor, 'LineWidth', 1.8, 'DisplayName', 'r_{forced}');
    xlabel('Effective Engineering Strain, \epsilon_{eff} (-)');
    ylabel('Radius (mm)');
    title('Geometry Diagnostics: Fixed Loaded Radius');
    legend('Location', 'best');
    grid on; box on; formatAxes(gca);
    saveCurrentFigure(resultsDir, ['04_geometry_loaded_radius_', versionTag]);
end

%% -----------------------------
% Optional fit / sweep calculations
% ------------------------------
if run.parameterSweeps
    [kT0_rmse, kT0_wrmse] = runFocusedKTSweep( ...
        eps_plot, eps_data, y_data, y_err, p, kT0_sweep);
end

if makeMaterialFitOverlayPlots
    plotMaterialFitOverlays(eps_plot, eps_data, y_data, y_err, p, E_sweep, nu_sweep, t_sweep_m, t_sweep_mm, ...
        c75, c80, c85, c90, markerSize, errorLineWidth, bareColor, resultsDir, versionTag);
end

if makeKT0FitOverlayPlot
    plotKT0FitOverlay(eps_plot, eps_data, y_data, y_err, p, kT0_sweep, ...
        c75, c80, c85, c90, markerSize, errorLineWidth, bareColor, resultsDir, versionTag);
end

if makeKT0IntermediateDiagnostics
    plotKT0IntermediateDiagnostics(eps_plot, p, kT0_sweep, resultsDir, versionTag);
end

if makePreTensionThicknessDiagnosticPlots
    plotPreTensionThicknessDiagnostics(eps_plot, p, resultsDir, versionTag);
end

if makeStaticDepthDiagnosticPlots
    plotStaticDepthDiagnostics(eps_plot, p, resultsDir, versionTag);
end


if makeLoadedRadiusDiagnosticPlot
    plotLoadedRadiusDiagnostic(eps_data, y_data_mem90, y_err_mem90, p, ...
        modelColor, resultsDir, versionTag);
end

if makePneumaticVolumeDiagnosticPlot
    plotPneumaticVolumeDiagnostic(eps_plot, eps_data, y_data_mem90, y_err_mem90, ...
        p, resultsDir, versionTag);
end

if makeKT0ThicknessReferenceFitSweep
    plotKT0ThicknessReferenceFitSweep(eps_plot, eps_data, y_data, y_err, p, ...
        kT0_grid, c75, c80, c85, c90, markerSize, errorLineWidth, ...
        bareColor, resultsDir, versionTag);
end

if makeKT0RawDebugPlots
    plotKT0RawImplementationDiagnostics(eps_plot, p, kT0_debug_sweep, resultsDir, versionTag);
end

if makeSolverComparisonPlots
    plotSolverAndVolumeCapComparison(eps_plot, p, kT0_debug_sweep, resultsDir, versionTag);
end

if printSweepSummary
    printSensitivitySweepSummary( ...
        kT0_sweep, kT0_rmse, kT0_wrmse);
end

%% -----------------------------
% Local functions
% ------------------------------
function printConfigurationSummary(p, scaleDiameters, nominalStrains, effectiveStrains)
    fprintf('Solver mode: %s\n', p.solverMode);
    fprintf('Trapped-gas fill pressure P_gas0: %.3f Pa\n', p.P_gas0);
    fprintf('External static pressure P_static: %.3f Pa\n', p.P_static);
    fprintf('Baseline hydrostatic depth: %.3f m\n', p.baselineDepth_m);
    fprintf('Use volume cap: %d\n', p.useVolumeCap);
    fprintf('Volume-change saturation fraction: %.3f\n', p.dV_cap_fraction);
    fprintf('Minimum volume fraction safeguard: %.3e\n', p.minVolumeFraction);
    fprintf('Nominal cavity volume V00: %.3e m^3\n', p.V00);
    fprintf('Geometric cavity volume: %.3e m^3\n', p.V_cavity_geom);
    fprintf('Tubing/sensor/fitting added volume: %.3e m^3\n', ...
        p.V_tubing + p.V_sensor + p.V_fittings);
    fprintf('Thickness mode: %s\n', p.thicknessMode);
    fprintf('Pre-tension thickness mode: %s\n', p.preTensionThicknessMode);
    fprintf('Plate thickness: %.3e m\n', p.t_plate0);
    fprintf('Baseline modulus: %.3e Pa\n', p.E_plate0);
    fprintf('Poisson ratio: %.3f\n', p.nu_plate);
    fprintf('Forced loading radius: %.3e m\n', p.r_forced_m);
    fprintf('Pressure tolerance: %.3e Pa\n', p.tolP);
    fprintf('Max iterations: %d\n', p.maxIter);
    fprintf('Pressure relaxation factor: %.3f\n', p.relax);
    fprintf('Pre-tension mode: %s\n', p.preTensionMode);
    fprintf('Pre-tension scale factor kT0: %.4f\n', p.kT0);
    fprintf('Assumed bare-port transmission R_bare: %.3f\n', p.R_bare);
    fprintf('Use strain-dependent E: %d\n', p.useStrainDependentE);
    fprintf('Pressure sweep local perturbation dP_local: %.3f Pa\n\n', p.dP_local);
    fprintf('Sampling interface diameters (mm): ');
    fprintf('%.3f ', scaleDiameters);
    fprintf('\n');
    fprintf('Nominal installation strain (geometry provenance only): ');
    fprintf('%.4f ', nominalStrains);
    fprintf('\n');
    fprintf('Effective engineering strain used by data and model: ');
    fprintf('%.4f ', effectiveStrains);
    fprintf('\n\n');
end

function runHydrostaticPreloadValidation(p)
    depth_cases = [0, p.depthSensitivityRange_m(1), p.baselineDepth_m, ...
        p.depthSensitivityRange_m(2), 1.0, p.deploymentDepthRange_m];
    eps_cases = [0, p.experimentalEffectiveStrains(2), ...
        p.experimentalEffectiveStrains(4)];

    fprintf('\n--- v31 hydrostatic-preload validation ---\n');
    fprintf('depth(m)  eps   rootConv  fixedConv  root/fp dPi(Pa)   gas residual   dV/V0    w/a    w/h    gain\n');

    for depth_i = depth_cases
        pRoot = p;
        pRoot.solverMode = 'root';
        pRoot.P_static = pRoot.P_gas0 + pRoot.rho_water * pRoot.g_water * depth_i;

        pFP = pRoot;
        pFP.solverMode = 'fixedPoint';

        for eps_i = eps_cases
            [Dplate, ~, Tpre_eff, ~, a_load, ~, ~, ~] = evaluateMembraneState(eps_i, pRoot);
            [~, ~, ~, ~, ~, t_installed] = evaluatePlateGeometry(eps_i, pRoot);
            initialState.Pi = pRoot.P_gas0;

            [Pi_root, Vi_root, dV_root, w_root, ~, conv_root, staticState] = ...
                solveCavityPressure(pRoot.P_static, a_load, Dplate, Tpre_eff, pRoot, initialState);
            [Pi_fp, ~, ~, ~, ~, conv_fp] = ...
                solveCavityPressure(pFP.P_static, a_load, Dplate, Tpre_eff, pFP, initialState);
            [gain_i, conv_gain] = evaluateLocalP3GainFromState( ...
                pRoot.P_static, a_load, Dplate, Tpre_eff, pRoot, staticState);

            gasResidual = abs(Pi_root * Vi_root - pRoot.P_gas0 * pRoot.V00) / ...
                (pRoot.P_gas0 * pRoot.V00);
            validRoot = conv_root && conv_gain;
            if conv_root && conv_fp
                solverDelta = abs(Pi_root - Pi_fp);
            else
                solverDelta = NaN;
            end

            fprintf('%8.4f %6.3f      %1.0f         %1.0f       %11.4e   %11.4e  %7.4f  %6.3f  %6.3f  %7.4f\n', ...
                depth_i, eps_i, validRoot, conv_fp, solverDelta, gasResidual, ...
                dV_root / pRoot.V00, abs(w_root) / a_load, abs(w_root) / t_installed, gain_i);
        end
    end

    fprintf('Validity note: w/a above about 0.1, w/h of order unity or larger, or\n');
    fprintf('large dV/V0 indicates that the linear small-deflection plate formulation is\n');
    fprintf('being extrapolated; depth results are mechanistic, not predictive.\n');
end

function pathTable = runLoadPathIndependenceValidation(eps_cases, p, resultsDir, versionTag)
    % Mike's load-path check:
    %   direct: sealed/rest reference -> final (hydrostatic +/- DeltaP/2)
    %   staged: sealed/rest reference -> hydrostatic -> final
    % The root solver determines equilibrium from the final boundary pressure
    % and conserved gas inventory, so the hydrostatic state is an initial
    % iterate only. With fixed constitutive properties, both paths must agree.

    nCases = numel(eps_cases);
    strain = eps_cases(:);
    Tpre_Npm = NaN(nCases,1);
    max_dPi_Pa = NaN(nCases,1);
    max_dV_m3 = NaN(nCases,1);
    max_dw_m = NaN(nCases,1);
    directConverged = false(nCases,1);
    stagedConverged = false(nCases,1);
    pathIndependent = false(nCases,1);

    fprintf('\n--- v37 load-path / tension-incorporation validation (Mike) ---\n');
    fprintf(['Compare direct rest-to-final equilibrium with rest-to-HS-to-final ', ...
        'equilibrium at Pstatic +/- DeltaP/2.\n']);
    fprintf(['Installation Tpre is fixed by prescribed pre-strain; hydrostatic ', ...
        'loading does not overwrite it.\n']);
    fprintf('eps     Tpre(N/m)    max|dPi|(Pa)    max|dV|(m^3)    max|dw|(m)    same\n');

    restState.Pi = p.P_gas0;
    for i = 1:nCases
        eps_i = eps_cases(i);
        [Dplate, ~, Tpre_eff, ~, a_load, ~, ~, ~] = evaluateMembraneState(eps_i, p);
        Tpre_Npm(i) = Tpre_eff;

        [~, ~, ~, ~, ~, convHS, hsState] = ...
            solveCavityPressure(p.P_static, a_load, Dplate, Tpre_eff, p, restState);

        finalPressures = p.P_static + [-p.dP0/2, p.dP0/2];
        dPi = NaN(1,2); dV = NaN(1,2); dw = NaN(1,2);
        convDirect = false(1,2); convStaged = false(1,2);
        for j = 1:2
            [PiDirect, ~, dVDirect, wDirect, ~, convDirect(j)] = ...
                solveCavityPressure(finalPressures(j), a_load, Dplate, Tpre_eff, p, restState);
            [PiStaged, ~, dVStaged, wStaged, ~, convStaged(j)] = ...
                solveCavityPressure(finalPressures(j), a_load, Dplate, Tpre_eff, p, hsState);
            dPi(j) = PiDirect - PiStaged;
            dV(j) = dVDirect - dVStaged;
            dw(j) = wDirect - wStaged;
        end

        max_dPi_Pa(i) = max(abs(dPi));
        max_dV_m3(i) = max(abs(dV));
        max_dw_m(i) = max(abs(dw));
        directConverged(i) = all(convDirect);
        stagedConverged(i) = convHS && all(convStaged);
        pathIndependent(i) = directConverged(i) && stagedConverged(i) && ...
            max_dPi_Pa(i) <= 10*p.tolP && ...
            max_dV_m3(i) <= 1e-12*max(p.V00, realmin) && ...
            max_dw_m(i) <= 1e-12*max(a_load, realmin);

        fprintf('%5.3f   %10.4f    %12.4e    %12.4e   %12.4e     %d\n', ...
            eps_i, Tpre_eff, max_dPi_Pa(i), max_dV_m3(i), max_dw_m(i), ...
            pathIndependent(i));
    end

    pathTable = table(strain, Tpre_Npm, max_dPi_Pa, max_dV_m3, max_dw_m, ...
        directConverged, stagedConverged, pathIndependent);
    writetable(pathTable, fullfile(resultsDir, ...
        ['load_path_validation_', versionTag, '.csv']));

    if ~all(pathIndependent)
        error(['Load-path validation failed. The final equilibrium depends on ', ...
            'whether the hydrostatic state was supplied as solver history.']);
    end
    fprintf(['PASS: direct and staged routes give the same final equilibrium. ', ...
        'Any future history-dependent tension law must intentionally revise this audit.\n']);
end

function printTransmissionReferenceAudit( ...
    eps_plot, rawTransmission, relBare, ...
    eps_data, y_data, y_err, p)

    fprintf('\n--- transmission reference and passivity audit ---\n');
    fprintf('Assumed bare-port transmission R_bare = %.3f.\n', ...
        p.R_bare);
    fprintf(['Modeled and experimental responses are presented on a ', ...
             'common bare-port-relative basis.\n\n']);

    fprintf(['eps_eff raw R_cav   model/R_bare', ...
             '   experiment/bare port   uncertainty\n']);

    for i = 1:numel(eps_data)
        raw_i = interp1( ...
            eps_plot, rawTransmission, eps_data(i), 'linear');

        relBare_i = interp1( ...
            eps_plot, relBare, eps_data(i), 'linear');

        fprintf(['%.3f    %8.5f      %8.5f', ...
                 '             %8.5f          %8.5f\n'], ...
            eps_data(i), raw_i, relBare_i, ...
            y_data(i), y_err(i));
    end

    tolerance = 1e-8;
    if p.includeKnoDeploymentDepths
        depthAuditLabel = 'calibration and deployment-depth';
    else
        depthAuditLabel = 'calibration-depth';
    end
    fprintf('\npassivity by %s case over eps = [%.3f, %.3f]\n', ...
        depthAuditLabel, min(eps_plot), max(eps_plot));
    fprintf('depth(m)   converged   min R_cav   max R_cav   violations   status\n');

    allCasesPassive = true;
    for depth_i = p.depth_sweep_m
        pDepth = p;
        pDepth.baselineDepth_m = depth_i;
        pDepth.P_static = pDepth.P_gas0 + pDepth.rho_water * pDepth.g_water * depth_i;
        [T_depth, ~, ~, ~, ~, ~, ~, ~, ~, converged] = ...
            runModelOverStrainRange(eps_plot, pDepth);

        finiteMask = isfinite(T_depth);
        if any(finiteMask)
            minTransmission = min(T_depth(finiteMask));
            maxTransmission = max(T_depth(finiteMask));
        else
            minTransmission = NaN;
            maxTransmission = NaN;
        end
        violationMask = ~finiteMask | T_depth < -tolerance | T_depth > 1 + tolerance;
        nViolations = nnz(violationMask);
        allConverged = all(converged > 0.5);
        passiveCase = allConverged && nViolations == 0;
        allCasesPassive = allCasesPassive && passiveCase;

        if passiveCase
            statusText = 'PASS';
        else
            statusText = 'CHECK';
        end
        fprintf('%8.4f      %1.0f       %9.6f   %9.6f      %4d       %s\n', ...
            depth_i, allConverged, minTransmission, maxTransmission, ...
            nViolations, statusText);
    end

    if allCasesPassive
        fprintf('PASS: 0 <= R_cav <= 1 for every converged baseline strain/depth case.\n');
    else
        warning('Transmission passivity audit found a failed or nonconverged case.');
    end
end

function plotTransmissionReferenceAudit( ...
    eps_plot, rawTransmission, relBare, ...
    eps_data, y_data, y_err, R_bare, ...
    modelColor, physColor, bareColor, resultsDir, versionTag, publicationOnly)

    figure('Position', [100 100 1150 520], ...
        'Tag', 'transmissionAudit');

    tiledlayout(1, 2, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    % -------------------------------------------------------------
    % Panel 1: raw modeled transmission
    % -------------------------------------------------------------
    nexttile;
    hold on;

    plot(eps_plot, rawTransmission, '-', ...
        'Color', modelColor, ...
        'LineWidth', 2.5, ...
        'DisplayName', 'Modeled R_{cav}');

    yline(R_bare, '-.', ...
        'Color', physColor, ...
        'LineWidth', 2.0, ...
        'DisplayName', sprintf('Assumed R_{bare} = %.2f', R_bare));

    yline(1.0, '--', ...
        'Color', bareColor, ...
        'LineWidth', 2.0, ...
        'DisplayName', 'Passivity Ceiling');

    xlabel('Effective Engineering Strain, \epsilon_{eff} (-)');
    ylabel('Raw Transmission');
    title('Raw Pressure Transmission');

    xlim([eps_plot(1) eps_plot(end)]);
    rawLower = max(0, 0.95 * min(rawTransmission(isfinite(rawTransmission))));
    ylim([rawLower 1.005]);

    legend('Location', 'southeast');
    grid on;
    box on;
    formatAxes(gca);

    % -------------------------------------------------------------
    % Panel 2: common bare-port normalization
    % -------------------------------------------------------------
    nexttile;
    hold on;

    plot(eps_plot, relBare, '-', ...
        'Color', physColor, ...
        'LineWidth', 2.5, ...
        'DisplayName', 'Model: R_{cav}/R_{bare}');

    errorbar(eps_data, y_data, y_err, 'ks', ...
        'MarkerFaceColor', 'k', ...
        'MarkerSize', 7, ...
        'LineWidth', 1.5, ...
        'CapSize', 7, ...
        'DisplayName', 'Experiment / Bare Port');

    yline(1.0, '--', ...
        'Color', bareColor, ...
        'LineWidth', 2.0, ...
        'DisplayName', 'Bare-Port Reference');

    xlabel('Effective Engineering Strain, \epsilon_{eff} (-)');
    ylabel('Bare-Port-Normalized Response');
    title('Model and Experimental Comparison');

    xlim([eps_plot(1) eps_plot(end)]);
    comparisonLower = max(0, 0.95 * min(relBare(isfinite(relBare))));
    ylim([comparisonLower 1.45]);

    legend('Location', 'best');
    grid on;
    box on;
    formatAxes(gca);

    if publicationOnly
        outName = ['02_transmission_and_experiment_',versionTag];
    else
        outName = ['09_transmission_reference_audit_',versionTag];
    end
    saveCurrentFigure(resultsDir,outName);
end

function printAggregatePeakDiagnostics(eps_data, y_data, y_err, y_data_mem90)
    fprintf('\n--- aggregate-data peak diagnostics ---\n');
    fprintf('These checks use four reported means only; they do not replace a trial-level analysis.\n');
    peakIdx = 2;
    for i = 1:numel(y_data)
        if i == peakIdx
            continue;
        end
        zsep = (y_data(peakIdx) - y_data(i)) / ...
            sqrt(y_err(peakIdx)^2 + y_err(i)^2);
        fprintf('eps %.3f maximum vs eps %.3f: combined-uncertainty separation = %.3f\n', ...
            eps_data(peakIdx), eps_data(i), zsep);
    end

    W = diag(1 ./ y_err.^2);
    modelNames = {'constant/plateau', 'linear', 'quadratic'};
    fprintf('\nmodel                 k       chi2        BIC       weighted RMSE\n');
    for degree = 0:2
        X = ones(numel(eps_data), degree + 1);
        for j = 1:degree
            X(:, j + 1) = eps_data(:).^j;
        end
        beta = (X' * W * X) \ (X' * W * y_data(:));
        yhat = X * beta;
        standardizedResidual = (y_data(:) - yhat) ./ y_err(:);
        chi2 = sum(standardizedResidual.^2);
        k = degree + 1;
        bic = chi2 + k * log(numel(eps_data));
        wrmse = sqrt(mean(standardizedResidual.^2));
        fprintf('%-20s %2d   %10.4f   %10.4f   %12.4f\n', ...
            modelNames{degree + 1}, k, chi2, bic, wrmse);

        if degree == 1
            fprintf('  linear slope = %.4f normalized-sensitivity/strain\n', beta(2));
        elseif degree == 2 && abs(beta(3)) > eps
            vertex = -beta(2) / (2 * beta(3));
            fprintf('  quadratic vertex = %.4f strain (descriptive only)\n', vertex);
        end
    end
    fprintf('90%%-membrane-normalized means: ');
    fprintf('%.4f ', y_data_mem90);
    fprintf('\nInterpret BIC cautiously because n = 4; retain "apparent maximum" unless trial-level data support a peak.\n');
end

function plotDualNormalizationComparison(eps_plot, modelRelZero, modelRelMem90, ...
    eps_data, y_data, y_err, y_data_mem90, y_err_mem90, ...
    modelColor, physColor, bareColor, resultsDir, versionTag)

    figure;
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile; hold on;
    errorbar(eps_data, y_data, y_err, 'ks', 'MarkerFaceColor', 'k', ...
        'LineWidth', 1.5, 'CapSize', 7, 'DisplayName', 'Experiment / bare port');
    plot(eps_plot, modelRelZero, '-', 'Color', modelColor, 'LineWidth', 2.5, ...
        'DisplayName', 'Model / zero-strain membrane');
    yline(1, ':', 'Color', bareColor, 'HandleVisibility', 'off');
    xlabel('Effective Engineering Strain, \epsilon_{eff} (-)');
    ylabel('Normalized Response (-)');
    title('Bare-Port vs. Zero-Strain References');
    legend('Location', 'best'); grid on; box on; formatAxes(gca);

    nexttile; hold on;
    errorbar(eps_data, y_data_mem90, y_err_mem90, 'ks', 'MarkerFaceColor', 'k', ...
        'LineWidth', 1.5, 'CapSize', 7, 'DisplayName', 'Experiment / 90% membrane');
    plot(eps_plot, modelRelMem90, '-', 'Color', physColor, 'LineWidth', 2.5, ...
        'DisplayName', 'Model / 90% membrane');
    yline(1, ':', 'Color', bareColor, 'HandleVisibility', 'off');
    xlabel('Effective Engineering Strain, \epsilon_{eff} (-)');
    ylabel('90%-Membrane-Relative Response (-)');
    title('90% Membrane Reference');
    legend('Location', 'best'); grid on; box on; formatAxes(gca);

    saveCurrentFigure(resultsDir, ['05_dual_normalization_', versionTag]);
end

function T = evaluateTransmissionAtState(eps_i, p)
    state.Pi = p.P_gas0;
    T = evaluateCavityTransmission(eps_i, p, state);
end

function plotLoadedRadiusDiagnostic(eps_data, targetMem90, targetErrMem90, p, ...
    modelColor, resultsDir, versionTag)

    nominalRadius = p.r_forced_m;
    refEps = eps_data(1);
    Tref = evaluateTransmissionAtState(refEps, p);
    scaleGrid = linspace(p.radiusDiagnosticBounds(1), p.radiusDiagnosticBounds(2), 101);
    extendedGrid = linspace(0.50, 2.00, 301);

    n = numel(eps_data);
    fixedResponse = NaN(1, n);
    boundedBestResponse = NaN(1, n);
    boundedBestScale = NaN(1, n);
    boundedMin = NaN(1, n);
    boundedMax = NaN(1, n);

    fprintf('\n--- bounded loaded-radius diagnostic (90%% membrane reference) ---\n');
    fprintf('eps    target   fixed    bestScale   bestBounded   boundedRange       extendedFeasible\n');

    for i = 1:n
        fixedResponse(i) = evaluateTransmissionAtState(eps_data(i), p) / Tref;
        if i == 1
            boundedBestScale(i) = 1.0;
            boundedBestResponse(i) = 1.0;
            boundedMin(i) = 1.0;
            boundedMax(i) = 1.0;
            extendedFeasible = true;
        else
            responseGrid = NaN(size(scaleGrid));
            for j = 1:numel(scaleGrid)
                pPoint = p;
                pPoint.r_forced_m = nominalRadius * scaleGrid(j);
                responseGrid(j) = evaluateTransmissionAtState(eps_data(i), pPoint) / Tref;
            end
            [~, bestIdx] = min(abs(responseGrid - targetMem90(i)));
            boundedBestScale(i) = scaleGrid(bestIdx);
            boundedBestResponse(i) = responseGrid(bestIdx);
            boundedMin(i) = min(responseGrid);
            boundedMax(i) = max(responseGrid);

            extendedResponse = NaN(size(extendedGrid));
            for j = 1:numel(extendedGrid)
                pPoint = p;
                pPoint.r_forced_m = nominalRadius * extendedGrid(j);
                extendedResponse(j) = evaluateTransmissionAtState(eps_data(i), pPoint) / Tref;
            end
            extendedFeasible = targetMem90(i) >= min(extendedResponse) && ...
                targetMem90(i) <= max(extendedResponse);
        end

        fprintf('%.3f  %7.4f  %7.4f    %7.4f      %7.4f     [%7.4f,%7.4f]        %d\n', ...
            eps_data(i), targetMem90(i), fixedResponse(i), boundedBestScale(i), ...
            boundedBestResponse(i), boundedMin(i), boundedMax(i), extendedFeasible);
    end

    figure;
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    nexttile; hold on;
    errorbar(eps_data, targetMem90, targetErrMem90, 'ks', 'MarkerFaceColor', 'k', ...
        'LineWidth', 1.5, 'CapSize', 7, 'DisplayName', 'Experiment / 90% membrane');
    plot(eps_data, fixedResponse, '-o', 'Color', modelColor, 'LineWidth', 2.2, ...
        'DisplayName', 'Fixed radius');
    plot(eps_data, boundedBestResponse, '--d', 'LineWidth', 2.2, ...
        'DisplayName', 'Best pointwise radius within +/-5%');
    xlabel('Effective Engineering Strain, \epsilon_{eff} (-)');
    ylabel('90%-Membrane-Relative Response (-)');
    title('Bounded Loaded-Radius Feasibility');
    legend('Location', 'best'); grid on; box on; formatAxes(gca);

    nexttile; hold on;
    plot(eps_data, boundedBestScale, '-o', 'LineWidth', 2.2);
    yline(p.radiusDiagnosticBounds(1), ':k');
    yline(p.radiusDiagnosticBounds(2), ':k');
    xlabel('Effective Engineering Strain, \epsilon_{eff} (-)');
    ylabel('Best Effective Radius Scale (-)');
    title('Pointwise Diagnostic - Not a Fitted Law');
    grid on; box on; formatAxes(gca);

    saveCurrentFigure(resultsDir, ['06_loaded_radius_bounded_', versionTag]);
end

function plotPneumaticVolumeDiagnostic(eps_plot, eps_data, targetMem90, targetErrMem90, ...
    p, resultsDir, versionTag)

    addedFractions = [0, 0.25, 0.50, 1.00];
    figure; hold on;
    errorbar(eps_data, targetMem90, targetErrMem90, 'ks', 'MarkerFaceColor', 'k', ...
        'LineWidth', 1.5, 'CapSize', 7, 'DisplayName', 'Experiment / 90% membrane');

    fprintf('\n--- added pneumatic-volume diagnostic ---\n');
    fprintf('Added volume is diagnostic until tubing, sensor, and fitting volumes are measured.\n');
    fprintf('added/Vcavity   totalVolume(mm^3)   response at experimental effective strains\n');

    for i = 1:numel(addedFractions)
        pVol = p;
        pVol.V00 = p.V_cavity_geom * (1 + addedFractions(i));
        [T, ~, ~, ~, ~, ~, ~, ~, ~, converged] = runModelOverStrainRange(eps_plot, pVol);
        ref = interp1(eps_plot, T, eps_data(1), 'linear');
        response = T ./ ref;
        response(converged < 0.5) = NaN;
        responsePts = interp1(eps_plot, response, eps_data, 'linear');

        plot(eps_plot, response, '-', 'LineWidth', 2.1, ...
            'DisplayName', sprintf('V_{added}/V_{cavity}=%.2f', addedFractions(i)));
        fprintf('%13.2f   %16.2f      ', addedFractions(i), 1e9 * pVol.V00);
        fprintf('%8.4f ', responsePts);
        fprintf('\n');
    end

    xlabel('Effective Engineering Strain, \epsilon_{eff} (-)');
    ylabel('90%-Membrane-Relative Response (-)');
    title('Sensitivity to Unmeasured Pneumatic Volume');
    legend('Location', 'best'); grid on; box on; formatAxes(gca);
    saveCurrentFigure(resultsDir, ['08_pneumatic_volume_', versionTag]);
end

function [T_cav, Dplate_curve, E_curve, Tpre_eff_curve, phi_pre_curve, ...
          a_load_curve, a_geom_curve, D_installed_curve_mm, eps_pre_curve, ...
          converged_curve, iter_curve] = runModelOverStrainRange(eps_plot, p)

    T_cav                = zeros(size(eps_plot));
    Dplate_curve         = zeros(size(eps_plot));
    E_curve              = zeros(size(eps_plot));
    Tpre_eff_curve       = zeros(size(eps_plot));
    phi_pre_curve        = zeros(size(eps_plot));
    a_load_curve         = zeros(size(eps_plot));
    a_geom_curve         = zeros(size(eps_plot));
    D_installed_curve_mm = zeros(size(eps_plot));
    eps_pre_curve        = zeros(size(eps_plot));
    converged_curve      = zeros(size(eps_plot));
    iter_curve           = zeros(size(eps_plot));

    solverState.Pi = p.P_gas0;

    for k = 1:numel(eps_plot)
        eps_query = eps_plot(k);
        [T_cav(k), Dplate_curve(k), E_curve(k), Tpre_eff_curve(k), ...
         phi_pre_curve(k), a_load_curve(k), a_geom_curve(k), ...
         D_installed_curve_mm(k), eps_pre_curve(k), converged_curve(k), ...
         iter_curve(k), solverState] = evaluateCavityTransmission(eps_query, p, solverState);
    end
end

function [R_rel_bare, T_cav_rel0] = runBareReferencedModelResponse(eps_plot, p)
    [T_cav, ~, ~, ~, ~, ~, ~, ~, ~, converged, ~] = ...
        runModelOverStrainRange(eps_plot, p);

    Tref = interp1(eps_plot, T_cav, 0.00, 'linear');
    if abs(Tref) < 1e-12 || ~isfinite(Tref)
        error(['Zero-strain membrane transmission is too small or invalid. ', ...
               'Inspect raw T_cav first.']);
    end

    T_cav_rel0 = T_cav ./ Tref;
    R_rel_bare = T_cav ./ p.R_bare;

    invalid = ~isfinite(T_cav) | converged < 0.5;
    T_cav_rel0(invalid) = NaN;
    R_rel_bare(invalid) = NaN;
end

function [rmse, wrmse] = computeModelErrors(eps_plot, S_model, eps_data, y_data, y_err)
    [~, ~, rmse, wrmse] = computePointResiduals(eps_plot, S_model, eps_data, y_data, y_err);
end

function [resid, wres, rmse, wrmse] = computePointResiduals(eps_plot, S_model, eps_data, y_data, y_err)
    resid = NaN(size(eps_data));
    wres  = NaN(size(eps_data));
    for i = 1:numel(eps_data)
        Si_model = interp1(eps_plot, S_model, eps_data(i), 'linear');
        resid(i) = Si_model - y_data(i);
        wres(i)  = resid(i) / y_err(i);
    end
    rmse  = sqrt(mean(resid .^ 2, 'omitnan'));
    wrmse = sqrt(mean(wres .^ 2, 'omitnan'));
end

function summaryTable = runPotentialEnergyReleaseDiagnostic( ...
        eps_cases, p, resultsDir, versionTag, runExtendedGauntlet, publicationOnly)
    % Conservative reduced-order large-deflection diagnostic. The pressure-
    % deflection law combines clamped-plate bending, installed pre-tension,
    % and the standard circular-membrane cubic stretching contribution:
    %   q = A*w + B*w^3
    %   A = 64*Kbend/a^4 + 4*Tpre/a^2
    %   B = kGeo*8*E*h/[3*(1-nu)*a^4]
    % With dV = Cv*pi*a^2*w, structural energy is defined so that
    % dUstruct/dw = q*d(dV)/dw. This makes loading and unloading conservative.

    Cv = 1/3;                   % clamped-plate volume/center-deflection factor
    kGeo = 1.0;                 % nominal cubic geometric-stretching coefficient
    dP = p.dP0;
    nCases = numel(eps_cases);

    strain = eps_cases(:);
    installed_Tpre_Npm = NaN(nCases,1);
    fixed_tension_gain = NaN(nCases,1);
    nonlinear_gain = NaN(nCases,1);
    Pi_HS_kPa = NaN(nCases,1);
    Pi_plus_kPa = NaN(nCases,1);
    Pi_minus_kPa = NaN(nCases,1);
    achieved_dPi_Pa = NaN(nCases,1);
    w_HS_mm = NaN(nCases,1);
    w_plus_mm = NaN(nCases,1);
    w_minus_mm = NaN(nCases,1);
    Tgeo_tangent_HS_Npm = NaN(nCases,1);
    Tgeo_tangent_plus_Npm = NaN(nCases,1);
    Tgeo_tangent_minus_Npm = NaN(nCases,1);
    Ustruct_HS_uJ = NaN(nCases,1);
    Ustruct_plus_uJ = NaN(nCases,1);
    Ustruct_minus_uJ = NaN(nCases,1);
    Ustruct_released_low_uJ = NaN(nCases,1);
    Ustruct_added_high_uJ = NaN(nCases,1);
    net_structural_change_uJ = NaN(nCases,1);
    net_total_stored_change_uJ = NaN(nCases,1);
    energy_balance_residual_nJ = NaN(nCases,1);
    equilibrium_residual_Pa = NaN(nCases,1);
    passive = false(nCases,1);
    converged = false(nCases,1);

    fprintf('\n--- potential-energy release diagnostic (conservative nonlinear membrane) ---\n');
    fprintf(['Cv=%.4f, kGeo=%.3f, depth=%.4f m, DeltaP=%.3f Pa.\n', ...
        'Stored energy may decrease on the unloaded side; passivity tests whether ', ...
        'that release produces dPi/dPext > 1.\n'], Cv, kGeo, p.baselineDepth_m, dP);
    fprintf(['eps    fixedGain  nonlinearGain  releaseLow(uJ)  addHigh(uJ)  ', ...
        'netStruct(uJ)  netStored(uJ)  TgeoHS(N/m)  passive\n']);

    for i = 1:nCases
        eps_i = eps_cases(i);
        [Kbend, E, Tpre, ~, a, ~, ~, ~] = evaluateMembraneState(eps_i, p);
        installed_Tpre_Npm(i) = Tpre;
        [~, ~, ~, ~, ~, h] = evaluatePlateGeometry(eps_i, p);

        [fixed_tension_gain(i), ~] = evaluateCavityTransmission(eps_i, p);
        hs = solveConservativeNonlinearState(p.P_static, a, Kbend, Tpre, E, h, Cv, kGeo, p);
        plusState = solveConservativeNonlinearState(p.P_static + dP/2, a, Kbend, Tpre, E, h, Cv, kGeo, p);
        minusState = solveConservativeNonlinearState(p.P_static - dP/2, a, Kbend, Tpre, E, h, Cv, kGeo, p);

        nonlinear_gain(i) = (plusState.Pi - minusState.Pi) / dP;
        Pi_HS_kPa(i) = hs.Pi / 1e3;
        Pi_plus_kPa(i) = plusState.Pi / 1e3;
        Pi_minus_kPa(i) = minusState.Pi / 1e3;
        achieved_dPi_Pa(i) = plusState.Pi - minusState.Pi;
        w_HS_mm(i) = 1e3 * hs.w;
        w_plus_mm(i) = 1e3 * plusState.w;
        w_minus_mm(i) = 1e3 * minusState.w;
        Tgeo_tangent_HS_Npm(i) = hs.TgeoTangent;
        Tgeo_tangent_plus_Npm(i) = plusState.TgeoTangent;
        Tgeo_tangent_minus_Npm(i) = minusState.TgeoTangent;
        Ustruct_HS_uJ(i) = 1e6 * hs.Ustruct;
        Ustruct_plus_uJ(i) = 1e6 * plusState.Ustruct;
        Ustruct_minus_uJ(i) = 1e6 * minusState.Ustruct;
        Ustruct_released_low_uJ(i) = 1e6 * (hs.Ustruct - minusState.Ustruct);
        Ustruct_added_high_uJ(i) = 1e6 * (plusState.Ustruct - hs.Ustruct);
        net_structural_change_uJ(i) = 1e6 * ...
            (plusState.Ustruct + minusState.Ustruct - 2*hs.Ustruct);
        net_total_stored_change_uJ(i) = 1e6 * ...
            (plusState.Ustruct + plusState.Wgas + minusState.Ustruct + minusState.Wgas ...
            - 2*(hs.Ustruct + hs.Wgas));
        workHigh = integrateExternalWork(hs, plusState, p);
        workLow = integrateExternalWork(hs, minusState, p);
        balanceHigh = (plusState.Ustruct + plusState.Wgas) ...
            - (hs.Ustruct + hs.Wgas) - workHigh;
        balanceLow = (minusState.Ustruct + minusState.Wgas) ...
            - (hs.Ustruct + hs.Wgas) - workLow;
        energy_balance_residual_nJ(i) = 1e9 * max(abs([balanceHigh,balanceLow]));
        equilibrium_residual_Pa(i) = max(abs([hs.residual, plusState.residual, minusState.residual]));
        passive(i) = nonlinear_gain(i) >= -1e-8 && nonlinear_gain(i) <= 1 + 1e-8;
        converged(i) = hs.converged && plusState.converged && minusState.converged;

        fprintf('%5.3f    %8.4f      %8.4f       %10.4f    %10.4f    %10.4e    %10.4e     %10.3f      %d\n', ...
            eps_i, fixed_tension_gain(i), nonlinear_gain(i), ...
            Ustruct_released_low_uJ(i), Ustruct_added_high_uJ(i), ...
            net_structural_change_uJ(i), net_total_stored_change_uJ(i), ...
            Tgeo_tangent_HS_Npm(i), passive(i));
    end

    summaryTable = table(strain, installed_Tpre_Npm, ...
        fixed_tension_gain, nonlinear_gain, ...
        Pi_HS_kPa, Pi_plus_kPa, Pi_minus_kPa, achieved_dPi_Pa, ...
        w_HS_mm, w_plus_mm, w_minus_mm, Tgeo_tangent_HS_Npm, ...
        Tgeo_tangent_plus_Npm, Tgeo_tangent_minus_Npm, Ustruct_HS_uJ, ...
        Ustruct_plus_uJ, Ustruct_minus_uJ, ...
        Ustruct_released_low_uJ, Ustruct_added_high_uJ, ...
        net_structural_change_uJ, net_total_stored_change_uJ, energy_balance_residual_nJ, ...
        equilibrium_residual_Pa, passive, converged);
    writetable(summaryTable, fullfile(resultsDir, ...
        ['potential_energy_summary_', versionTag, '.csv']));

    if ~publicationOnly
        figure; hold on;
        plot(eps_cases, fixed_tension_gain, '-o', 'LineWidth', 2.5, ...
            'DisplayName', 'v33 fixed-tension baseline');
        plot(eps_cases, nonlinear_gain, '-s', 'LineWidth', 2.5, ...
            'DisplayName', 'Conservative state-dependent tension');
        yline(1, ':k', 'DisplayName', 'Unity transmission');
        xlabel('Effective Engineering Strain, \epsilon_{eff} (-)');
        ylabel('Raw Differential Transmission, \Delta P_i/\Delta P_{ext} (-)');
        title('Potential-Energy Release Test');
        subtitle('Nonlinear geometric stretching allows tension and stored energy to change with load');
        grid on; box on; formatAxes(gca); legend('Location', 'best');
        saveCurrentFigure(resultsDir, ['80_potential_energy_transmission_', versionTag]);

        figure; hold on;
        plot(eps_cases, Ustruct_released_low_uJ, '-o', 'LineWidth', 2.5, ...
            'DisplayName', 'Released on unloaded side');
        plot(eps_cases, Ustruct_added_high_uJ, '-s', 'LineWidth', 2.5, ...
            'DisplayName', 'Added on loaded side');
        plot(eps_cases, net_structural_change_uJ, '-^', 'LineWidth', 2.5, ...
            'DisplayName', 'Net two-side structural change');
        xlabel('Effective Engineering Strain, \epsilon_{eff} (-)');
        ylabel('Structural Energy Change (\muJ)');
        title('Membrane-Energy Redistribution Under Symmetric Differential Loading');
        grid on; box on; formatAxes(gca); legend('Location', 'best');
        saveCurrentFigure(resultsDir, ['81_potential_energy_redistribution_', versionTag]);
    end

    figure; hold on;
    plot(installed_Tpre_Npm, Ustruct_HS_uJ, '-o', ...
        'LineWidth', 2.5, ...
        'DisplayName', 'Hydrostatic equilibrium');
    plot(installed_Tpre_Npm, Ustruct_plus_uJ, '-s', ...
        'LineWidth', 2.5, ...
        'DisplayName', 'High side: $P_{\mathrm{static}}+\Delta P/2$');
    plot(installed_Tpre_Npm, Ustruct_minus_uJ, '-d', ...
        'LineWidth', 2.5, ...
        'DisplayName', 'Low side: $P_{\mathrm{static}}-\Delta P/2$');
    xlabel('Effective Installed Tension, $T_0$ (N/m)', ...
    'Interpreter', 'latex');
    ylabel('Structural Potential Energy, $U_{\mathrm{struct}}$ ($\mu$J)', ...
    'Interpreter', 'latex');
    title('Potential Energy vs. Effective Installed Tension');
    grid on; box on; formatAxes(gca); legend('Location', 'best', 'Interpreter', 'latex');
    if publicationOnly
        outName = ['01_potential_energy_vs_tension_',versionTag];
    else
        outName = ['88_potential_energy_vs_tension_',versionTag];
    end
    saveCurrentFigure(resultsDir,outName);

    if runExtendedGauntlet
        runPotentialEnergyCoefficientAndAmplitudeAudit( ...
            eps_cases, p, Cv, resultsDir, versionTag);
        runLocalPrestressRelaxationBound( ...
            eps_cases, fixed_tension_gain, p, Cv, resultsDir, versionTag);
    end
    if ~publicationOnly
        plotPotentialEnergyAndDeflectionPaths(eps_cases, p, resultsDir, versionTag);
    end
end

function state = solveConservativeNonlinearState(Psurf, a, Kbend, Tpre, E, h, Cv, kGeo, p)
    areaFactor = Cv * pi * a^2;
    A = 64 * Kbend / a^4 + 4 * Tpre / a^2;
    B = kGeo * 8 * E * h / (3 * (1 - p.nu_plate) * a^4);
    wLimit = 0.999 * p.V00 / areaFactor;
    residual = @(w) nonlinearEnergyResidual(w, Psurf, A, B, areaFactor, p);

    wLo = 0;
    wHi = min(wLimit, max(1e-9, (max(Psurf - p.P_gas0, 0) / max(B,1e-30))^(1/3) * 4));
    while residual(wHi) > 0 && wHi < wLimit
        wHi = min(2*wHi, wLimit);
    end

    fLo = residual(wLo); fHi = residual(wHi);
    converged = isfinite(fLo) && isfinite(fHi) && fLo >= 0 && fHi <= 0;
    w = NaN;
    if converged
        for iter = 1:p.maxIter
            wMid = 0.5 * (wLo + wHi);
            fMid = residual(wMid);
            if abs(fMid) < p.tolP || abs(wHi-wLo) < 1e-12
                w = wMid;
                break;
            elseif fMid > 0
                wLo = wMid;
            else
                wHi = wMid;
            end
        end
        if ~isfinite(w), w = 0.5*(wLo+wHi); end
    end

    dV = areaFactor * w;
    V = p.V00 - dV;
    Pi = p.P_gas0 * (p.V00 / V)^p.gasExponent;
    q = A*w + B*w^3;
    Ustruct = areaFactor * (0.5*A*w^2 + 0.25*B*w^4);
    if abs(p.gasExponent - 1) < 1e-12
        Wgas = p.P_gas0 * p.V00 * log(p.V00 / V);
    else
        Wgas = p.P_gas0 * p.V00 / (p.gasExponent - 1) * ...
            ((p.V00/V)^(p.gasExponent-1) - 1);
    end
    state = struct('w',w,'dV',dV,'V',V,'Pi',Pi,'q',q,'A',A,'B',B, ...
        'areaFactor',areaFactor, ...
        'Ustruct',Ustruct,'Wgas',Wgas,'TgeoTangent',0.75*B*a^2*w^2, ...
        'residual',Psurf-Pi-q,'converged',converged && isfinite(w));
end

function Wext = integrateExternalWork(stateA, stateB, p)
    % Independent numerical path integral: Wext = integral(Psurf dDeltaV).
    % Agreement with Delta(Ustruct+Wgas) audits conservative energy closure.
    wPath = linspace(stateA.w, stateB.w, 2001);
    VPath = p.V00 - stateA.areaFactor*wPath;
    PiPath = p.P_gas0 * (p.V00 ./ VPath).^p.gasExponent;
    qPath = stateA.A*wPath + stateA.B*wPath.^3;
    PsurfPath = PiPath + qPath;
    Wext = trapz(wPath, PsurfPath * stateA.areaFactor);
end

function F = nonlinearEnergyResidual(w, Psurf, A, B, areaFactor, p)
    V = p.V00 - areaFactor*w;
    if V <= 0
        F = -Inf;
        return;
    end
    Pi = p.P_gas0 * (p.V00/V)^p.gasExponent;
    F = Psurf - Pi - A*w - B*w^3;
end

function pathTable = plotPotentialEnergyAndDeflectionPaths(eps_cases, p, resultsDir, versionTag)
    % Trace equilibrium paths across the full symmetric assembly excursion.
    % A signed port offset of +/-400 Pa corresponds to an 800 Pa total
    % port-to-port differential. Ustruct is membrane structural potential
    % energy referenced to the undeformed state, not gas-compression energy.
    portOffsetCases = linspace(-p.dP0/2, p.dP0/2, 81);
    modelNames = ["Locked elastic (v33)", ...
                  "Locked nonlinear, C_V=1/3", ...
                  "Locally relaxed nonlinear, C_V=1/2"];
    modelTypes = ["lockedElastic", "lockedNonlinear", "relaxedNonlinear"];
    CvCases = [NaN, 1/3, 1/2];

    nRows = numel(modelNames)*numel(eps_cases)*numel(portOffsetCases);
    model_case = strings(nRows,1);
    strain = NaN(nRows,1);
    signed_port_offset_Pa = NaN(nRows,1);
    equivalent_total_differential_Pa = NaN(nRows,1);
    wmax_mm = NaN(nRows,1);
    Ustruct_uJ = NaN(nRows,1);
    Ulinear_uJ = NaN(nRows,1);
    Ugeometric_uJ = NaN(nRows,1);
    cubic_to_linear_load_ratio = NaN(nRows,1);
    Wgas_uJ = NaN(nRows,1);
    P_cavity_kPa_abs = NaN(nRows,1);
    converged = false(nRows,1);
    row = 0;

    for im = 1:numel(modelNames)
        for eps_i = eps_cases
            [Kbend,E,Tpre,~,a,~,~,~] = evaluateMembraneState(eps_i,p);
            [~,~,~,~,~,h] = evaluatePlateGeometry(eps_i,p);
            if modelTypes(im)=="relaxedNonlinear"
                Tpre = 0;
            end
            for offset_i = portOffsetCases
                if modelTypes(im)=="lockedElastic"
                    initialState.Pi = p.P_gas0;
                    [Pi,Vi,~,w,~,conv,~,dVraw,q] = solveCavityPressure( ...
                        p.P_static+offset_i,a,Kbend,Tpre,p,initialState);
                    Ustruct = 0.5*q*dVraw;
                    Ulinear = Ustruct;
                    Ugeometric = 0;
                    cubicRatio = 0;
                    if abs(p.gasExponent-1)<1e-12
                        Wgas = p.P_gas0*p.V00*log(p.V00/Vi);
                    else
                        Wgas = p.P_gas0*p.V00/(p.gasExponent-1)* ...
                            ((p.V00/Vi)^(p.gasExponent-1)-1);
                    end
                else
                    state = solveConservativeNonlinearState( ...
                        p.P_static+offset_i,a,Kbend,Tpre,E,h,CvCases(im),1,p);
                    Pi = state.Pi; w = state.w; Ustruct = state.Ustruct;
                    Ulinear = state.areaFactor*0.5*state.A*w^2;
                    Ugeometric = state.areaFactor*0.25*state.B*w^4;
                    cubicRatio = state.B*w^2/max(state.A,1e-30);
                    Wgas = state.Wgas; conv = state.converged;
                end
                row = row+1;
                model_case(row) = modelNames(im);
                strain(row) = eps_i;
                signed_port_offset_Pa(row) = offset_i;
                equivalent_total_differential_Pa(row) = 2*abs(offset_i);
                wmax_mm(row) = 1e3*w;
                Ustruct_uJ(row) = 1e6*Ustruct;
                Ulinear_uJ(row) = 1e6*Ulinear;
                Ugeometric_uJ(row) = 1e6*Ugeometric;
                cubic_to_linear_load_ratio(row) = cubicRatio;
                Wgas_uJ(row) = 1e6*Wgas;
                P_cavity_kPa_abs(row) = Pi/1e3;
                converged(row) = conv;
            end
        end
    end

    pathTable = table(model_case,strain,signed_port_offset_Pa, ...
        equivalent_total_differential_Pa,wmax_mm,Ustruct_uJ,Ulinear_uJ, ...
        Ugeometric_uJ,cubic_to_linear_load_ratio,Wgas_uJ,P_cavity_kPa_abs,converged);
    writetable(pathTable,fullfile(resultsDir, ...
        ['potential_energy_deflection_paths_',versionTag,'.csv']));

    figure('Position',[40 80 1900 680]);
    energyLayout = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
    for im = 1:numel(modelNames)
        ax = nexttile; hold(ax,'on');
        for eps_i = eps_cases
            mask = model_case==modelNames(im) & abs(strain-eps_i)<1e-12;
            hLine = plot(ax,wmax_mm(mask),Ustruct_uJ(mask),'-','LineWidth',2.3, ...
                'DisplayName',sprintf('strain = %.3f',eps_i));
            hsMask = mask & abs(signed_port_offset_Pa)<1e-12;
            plot(ax,wmax_mm(hsMask),Ustruct_uJ(hsMask),'o', ...
                'MarkerSize',7,'MarkerFaceColor',hLine.Color, ...
                'MarkerEdgeColor',hLine.Color,'HandleVisibility','off');
        end
        xlabel(ax,'Maximum Center Deflection, w_{max} (mm)');
        ylabel(ax,'Membrane Structural Potential Energy (\muJ)');
        title(ax,modelNames(im));
        grid(ax,'on'); box(ax,'on'); formatAxes(ax);
        legend(ax,'Location','best');
    end
    title(energyLayout,'Membrane Potential Energy Along the Differential-Loading Path');
    subtitle(energyLayout,'Filled markers identify the common hydrostatic equilibrium (zero port offset)');
    saveCurrentFigure(resultsDir,['84_potential_energy_vs_wmax_',versionTag]);

    figure('Position',[40 80 1900 680]);
    deflectionLayout = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
    for im = 1:numel(modelNames)
        ax = nexttile; hold(ax,'on');
        for eps_i = eps_cases
            mask = model_case==modelNames(im) & abs(strain-eps_i)<1e-12;
            plot(ax,signed_port_offset_Pa(mask),wmax_mm(mask),'-','LineWidth',2.3, ...
                'DisplayName',sprintf('strain = %.3f',eps_i));
        end
        xline(ax,0,':k','Hydrostatic state','LineWidth',1.5, ...
            'HandleVisibility','off','LabelVerticalAlignment','bottom');
        xlabel(ax,'Signed Port Pressure Perturbation, \deltaP_{port} (Pa)');
        ylabel(ax,'Maximum Center Deflection, w_{max} (mm)');
        title(ax,modelNames(im));
        grid(ax,'on'); box(ax,'on'); formatAxes(ax);
        legend(ax,'Location','best');
    end
    title(deflectionLayout,'Maximum Deflection Versus Differential Port Loading');
    subtitle(deflectionLayout, ...
        'For symmetric loading, total port-to-port differential = 2|\deltaP_{port}|');
    saveCurrentFigure(resultsDir,['85_wmax_vs_differential_pressure_',versionTag]);

    fprintf('\n--- potential-energy/deflection path diagnostics ---\n');
    fprintf(['Saved locked-elastic, locked-nonlinear, and relaxed-nonlinear paths over signed port offsets ', ...
        '[%.1f, %.1f] Pa (%d converged of %d).\n'], ...
        min(portOffsetCases),max(portOffsetCases),nnz(converged),numel(converged));
    for im = 1:numel(modelNames)
        mask = model_case==modelNames(im);
        fprintf('%-38s cubic/linear load ratio range: %.4g to %.4g\n', ...
            modelNames(im),min(cubic_to_linear_load_ratio(mask)), ...
            max(cubic_to_linear_load_ratio(mask)));
    end
end

function auditTable = runPotentialEnergyCoefficientAndAmplitudeAudit(eps_cases, p, CvNominal, resultsDir, versionTag)
    kGeoCases = [0, 0.25, 0.5, 1, 2];
    CvCases = [1/3, 1/2];
    amplitudeCases = [5, 10, 25, 50, 100, 200, 400, 800];
    nRows = numel(eps_cases)*numel(kGeoCases)*numel(CvCases)*numel(amplitudeCases);
    strain = NaN(nRows,1); kGeo = NaN(nRows,1); Cv = NaN(nRows,1);
    differential_Pa = NaN(nRows,1); raw_gain = NaN(nRows,1);
    released_low_uJ = NaN(nRows,1); passive = false(nRows,1); converged = false(nRows,1);
    row = 0;
    for Cv_i = CvCases
        for kGeo_i = kGeoCases
            for dP_i = amplitudeCases
                for eps_i = eps_cases
                    row = row + 1;
                    [Kbend, E, Tpre, ~, a, ~, ~, ~] = evaluateMembraneState(eps_i, p);
                    [~, ~, ~, ~, ~, h] = evaluatePlateGeometry(eps_i, p);
                    hs = solveConservativeNonlinearState(p.P_static, a, Kbend, Tpre, E, h, Cv_i, kGeo_i, p);
                    plusState = solveConservativeNonlinearState(p.P_static+dP_i/2, a, Kbend, Tpre, E, h, Cv_i, kGeo_i, p);
                    minusState = solveConservativeNonlinearState(p.P_static-dP_i/2, a, Kbend, Tpre, E, h, Cv_i, kGeo_i, p);
                    gain_i = (plusState.Pi-minusState.Pi)/dP_i;
                    strain(row)=eps_i; kGeo(row)=kGeo_i; Cv(row)=Cv_i;
                    differential_Pa(row)=dP_i; raw_gain(row)=gain_i;
                    released_low_uJ(row)=1e6*(hs.Ustruct-minusState.Ustruct);
                    passive(row)=gain_i>=-1e-8 && gain_i<=1+1e-8;
                    converged(row)=hs.converged&&plusState.converged&&minusState.converged;
                end
            end
        end
    end
    auditTable = table(strain,kGeo,Cv,differential_Pa,raw_gain,released_low_uJ,passive,converged);
    writetable(auditTable, fullfile(resultsDir, ['potential_energy_gauntlet_',versionTag,'.csv']));
    fprintf('\n--- potential-energy coefficient/amplitude gauntlet ---\n');
    fprintf('Cases: %d; converged: %d; passive: %d; max gain: %.6f; min gain: %.6f\n', ...
        height(auditTable), nnz(converged), nnz(passive), max(raw_gain), min(raw_gain));
    nominalMask = abs(Cv-CvNominal)<1e-12 & abs(kGeo-1)<1e-12 & abs(differential_Pa-p.dP0)<1e-12;
    fprintf('Nominal %.0f-Pa gains: ',p.dP0); fprintf('%.6f ',raw_gain(nominalMask)); fprintf('\n');
end

function boundTable = runLocalPrestressRelaxationBound(eps_cases, fixedGain, p, Cv, resultsDir, versionTag)
    % Parameter-free endpoints only:
    % locked elastic = existing full installed tension;
    % relaxed nonlinear = zero residual installed tension plus B*w^3;
    % relaxed linear = zero residual installed tension and B=0.
    nCases = numel(eps_cases);
    relaxed_nonlinear_gain = NaN(nCases,1);
    relaxed_nonlinear_membrane_shape_gain = NaN(nCases,1);
    relaxed_linear_gain = NaN(nCases,1);
    relaxed_nonlinear_over_Rbare = NaN(nCases,1);
    relaxed_nonlinear_membrane_shape_over_Rbare = NaN(nCases,1);
    relaxed_linear_over_Rbare = NaN(nCases,1);
    Tgeo_HS_Npm = NaN(nCases,1);
    wHS_over_h_relaxed_nonlinear = NaN(nCases,1);
    wHS_over_h_relaxed_linear = NaN(nCases,1);
    passive_nonlinear = false(nCases,1);
    passive_linear = false(nCases,1);
    converged = false(nCases,1);

    fprintf('\n--- local installed-prestress relaxation/redistribution bound ---\n');
    fprintf(['No empirical scale is used. Residual installed tension is either fully ', ...
        'locked (v33) or zero locally; geometric tension remains in the nonlinear bound.\n']);
    fprintf(['eps    lockedGain  relaxedNL(Cv=1/3)  relaxedNL(Cv=1/2)  relaxedLinear  ', ...
        'NLhalf/Rbare  Linear/Rbare  TgeoHS(N/m)  wNL/h  wLinear/h\n']);

    for i = 1:nCases
        eps_i = eps_cases(i);
        [Kbend, E, ~, ~, a, ~, ~, ~] = evaluateMembraneState(eps_i, p);
        [~, ~, ~, ~, ~, h] = evaluatePlateGeometry(eps_i, p);
        dP = p.dP0;

        hsNL = solveConservativeNonlinearState(p.P_static, a, Kbend, 0, E, h, Cv, 1, p);
        plusNL = solveConservativeNonlinearState(p.P_static+dP/2, a, Kbend, 0, E, h, Cv, 1, p);
        minusNL = solveConservativeNonlinearState(p.P_static-dP/2, a, Kbend, 0, E, h, Cv, 1, p);
        hsNLhalf = solveConservativeNonlinearState(p.P_static, a, Kbend, 0, E, h, 1/2, 1, p);
        plusNLhalf = solveConservativeNonlinearState(p.P_static+dP/2, a, Kbend, 0, E, h, 1/2, 1, p);
        minusNLhalf = solveConservativeNonlinearState(p.P_static-dP/2, a, Kbend, 0, E, h, 1/2, 1, p);
        hsLinear = solveConservativeNonlinearState(p.P_static, a, Kbend, 0, E, h, Cv, 0, p);
        plusLinear = solveConservativeNonlinearState(p.P_static+dP/2, a, Kbend, 0, E, h, Cv, 0, p);
        minusLinear = solveConservativeNonlinearState(p.P_static-dP/2, a, Kbend, 0, E, h, Cv, 0, p);

        relaxed_nonlinear_gain(i) = (plusNL.Pi-minusNL.Pi)/dP;
        relaxed_nonlinear_membrane_shape_gain(i) = (plusNLhalf.Pi-minusNLhalf.Pi)/dP;
        relaxed_linear_gain(i) = (plusLinear.Pi-minusLinear.Pi)/dP;
        relaxed_nonlinear_over_Rbare(i) = relaxed_nonlinear_gain(i)/p.R_bare;
        relaxed_nonlinear_membrane_shape_over_Rbare(i) = relaxed_nonlinear_membrane_shape_gain(i)/p.R_bare;
        relaxed_linear_over_Rbare(i) = relaxed_linear_gain(i)/p.R_bare;
        Tgeo_HS_Npm(i) = hsNL.TgeoTangent;
        wHS_over_h_relaxed_nonlinear(i) = abs(hsNLhalf.w)/h;
        wHS_over_h_relaxed_linear(i) = abs(hsLinear.w)/h;
        passive_nonlinear(i) = relaxed_nonlinear_gain(i)>=0 && relaxed_nonlinear_gain(i)<=1;
        passive_linear(i) = relaxed_linear_gain(i)>=0 && relaxed_linear_gain(i)<=1;
        converged(i) = hsNL.converged && plusNL.converged && minusNL.converged ...
            && hsNLhalf.converged && plusNLhalf.converged && minusNLhalf.converged ...
            && hsLinear.converged ...
            && plusLinear.converged && minusLinear.converged;

        fprintf('%5.3f    %8.4f         %8.4f            %8.4f        %8.4f        %8.4f       %8.4f     %9.3f   %6.2f    %7.2f\n', ...
            eps_i, fixedGain(i), relaxed_nonlinear_gain(i), relaxed_nonlinear_membrane_shape_gain(i), ...
            relaxed_linear_gain(i), relaxed_nonlinear_membrane_shape_over_Rbare(i), ...
            relaxed_linear_over_Rbare(i), Tgeo_HS_Npm(i), ...
            wHS_over_h_relaxed_nonlinear(i), wHS_over_h_relaxed_linear(i));
    end

    strain = eps_cases(:);
    locked_elastic_gain = fixedGain(:);
    boundTable = table(strain,locked_elastic_gain,relaxed_nonlinear_gain, ...
        relaxed_nonlinear_membrane_shape_gain,relaxed_linear_gain, ...
        relaxed_nonlinear_over_Rbare,relaxed_nonlinear_membrane_shape_over_Rbare, ...
        relaxed_linear_over_Rbare,Tgeo_HS_Npm,wHS_over_h_relaxed_nonlinear, ...
        wHS_over_h_relaxed_linear,passive_nonlinear,passive_linear,converged);
    writetable(boundTable,fullfile(resultsDir,['local_prestress_relaxation_bound_',versionTag,'.csv']));

    figure('Position',[100 100 1200 800]); hold on;
    plot(eps_cases,fixedGain,'-o','LineWidth',2.5,'DisplayName','Locked elastic installation (v33)');
    plot(eps_cases,relaxed_nonlinear_gain,'-s','LineWidth',2.5, ...
        'DisplayName','Relaxed + geometric tension, C_V=1/3');
    plot(eps_cases,relaxed_nonlinear_membrane_shape_gain,'-d','LineWidth',2.5, ...
        'DisplayName','Relaxed + geometric tension, C_V=1/2');
    plot(eps_cases,relaxed_linear_gain,'-^','LineWidth',2.5, ...
        'DisplayName','Locally relaxed linear-plate endpoint');
    yline(p.R_bare,'--','LineWidth',1.8,'DisplayName',sprintf('Assumed bare path R_{bare}=%.2f',p.R_bare));
    yline(1,':k','LineWidth',1.5,'DisplayName','Unity raw transmission');
    xlabel('Effective Engineering Strain, \epsilon_{eff} (-)');
    ylabel('Raw Differential Transmission (-)');
    title('Parameter-Free Local Prestress-Relaxation Bounds');
    subtitle('Relaxed cases remove installed residual tension without fitting; nonlinear case retains geometric stiffening');
    grid on; box on; formatAxes(gca); legend('Location','best');
    saveCurrentFigure(resultsDir,['82_local_prestress_relaxation_bound_',versionTag]);

    % Finite-amplitude passivity audit for the physically stronger relaxed-
    % nonlinear case, evaluated across both common shape-factor endpoints.
    amplitudeCases = [5,10,25,50,100,200,400,800];
    CvCases = [1/3,1/2];
    nAuditRows = numel(CvCases)*numel(amplitudeCases)*numel(eps_cases);
    auditStrain = NaN(nAuditRows,1);
    auditCv = NaN(nAuditRows,1);
    totalDifferential_Pa = NaN(nAuditRows,1);
    symmetricPortOffset_Pa = NaN(nAuditRows,1);
    auditRawGain = NaN(nAuditRows,1);
    auditOverRbare = NaN(nAuditRows,1);
    auditPassive = false(nAuditRows,1);
    auditConverged = false(nAuditRows,1);
    auditRow = 0;
    maxGain = -Inf; minGain = Inf; allPassive = true; allConverged = true;
    for Cv_i = CvCases
        for dP_i = amplitudeCases
            for eps_i = eps_cases
                [Kbend,E,~,~,a,~,~,~] = evaluateMembraneState(eps_i,p);
                [~,~,~,~,~,h] = evaluatePlateGeometry(eps_i,p);
                plusState = solveConservativeNonlinearState(p.P_static+dP_i/2,a,Kbend,0,E,h,Cv_i,1,p);
                minusState = solveConservativeNonlinearState(p.P_static-dP_i/2,a,Kbend,0,E,h,Cv_i,1,p);
                gain_i = (plusState.Pi-minusState.Pi)/dP_i;
                auditRow = auditRow + 1;
                auditStrain(auditRow) = eps_i;
                auditCv(auditRow) = Cv_i;
                totalDifferential_Pa(auditRow) = dP_i;
                symmetricPortOffset_Pa(auditRow) = dP_i/2;
                auditRawGain(auditRow) = gain_i;
                auditOverRbare(auditRow) = gain_i/p.R_bare;
                auditPassive(auditRow) = gain_i>=0 && gain_i<=1;
                auditConverged(auditRow) = plusState.converged && minusState.converged;
                maxGain=max(maxGain,gain_i); minGain=min(minGain,gain_i);
                allPassive=allPassive&&auditPassive(auditRow);
                allConverged=allConverged&&auditConverged(auditRow);
            end
        end
    end
    amplitudeTable = table(auditStrain,auditCv,totalDifferential_Pa, ...
        symmetricPortOffset_Pa,auditRawGain,auditOverRbare,auditPassive,auditConverged, ...
        'VariableNames',{'strain','Cv','total_differential_Pa','port_offset_plus_minus_Pa', ...
        'raw_gain','gain_over_Rbare','passive','converged'});
    writetable(amplitudeTable,fullfile(resultsDir, ...
        ['differential_amplitude_characterization_',versionTag,'.csv']));

    figure; hold on;
    for eps_i = eps_cases
        mask = abs(auditCv-0.5)<1e-12 & abs(auditStrain-eps_i)<1e-12;
        plot(totalDifferential_Pa(mask),auditRawGain(mask),'-o','LineWidth',2.2, ...
            'DisplayName',sprintf('strain = %.3f',eps_i));
    end
    yline(p.R_bare,'--','LineWidth',1.8, ...
        'DisplayName',sprintf('Assumed bare path R_{bare}=%.2f',p.R_bare));
    yline(1,':k','LineWidth',1.5,'DisplayName','Unity raw transmission');
    xline(400,':','LineWidth',1.5,'DisplayName','Nominal 400 Pa total');
    xline(800,'-.','LineWidth',1.5,'DisplayName','+400/-400 Pa ports');
    set(gca,'XScale','log');
    xlabel('Total Port-to-Port Differential Pressure (Pa)');
    ylabel('Raw Differential Transmission (-)');
    title('Differential-Amplitude Dependence of Relaxed Nonlinear Bound');
    subtitle('C_V=1/2; total differential is applied symmetrically as +DeltaP/2 and -DeltaP/2');
    grid on; box on; formatAxes(gca); legend('Location','best');
    saveCurrentFigure(resultsDir,['83_differential_amplitude_characterization_',versionTag]);

    fprintf('Relaxed-nonlinear amplitude/shape audit: minGain=%.6f maxGain=%.6f converged=%d passive=%d\n', ...
        minGain,maxGain,allConverged,allPassive);
end

function knoTable = compileKnoDepthReport(eps_cases, p, resultsDir, versionTag)
    deploymentDepths = p.deploymentDepthRange_m;
    nRows = numel(deploymentDepths) * numel(eps_cases);
    depth_m = NaN(nRows,1); strain = NaN(nRows,1);
    P_external_abs_Pa = NaN(nRows,1); P_cavity_abs_Pa = NaN(nRows,1);
    raw_gain = NaN(nRows,1); wstatic_mm = NaN(nRows,1);
    dVstatic_over_V0 = NaN(nRows,1); w_over_a = NaN(nRows,1);
    w_over_h = NaN(nRows,1); within_small_deflection = false(nRows,1);
    converged = false(nRows,1);
    row = 0;

    fprintf('\n--- KNO deployment-depth report ---\n');
    fprintf('depth(m)  eps    Pext(abs Pa)  Pi(abs Pa)   raw gain   dV/V0    w/a    w/h   valid\n');
    for depth_i = deploymentDepths
        pDepth = p;
        pDepth.baselineDepth_m = depth_i;
        pDepth.P_static = pDepth.P_gas0 + pDepth.rho_water * pDepth.g_water * depth_i;
        for eps_i = eps_cases
            row = row + 1;
            [Kbend_i, ~, Tpre_i, ~, a_i, ~, ~, ~] = evaluateMembraneState(eps_i, pDepth);
            [~, ~, ~, ~, ~, h_i] = evaluatePlateGeometry(eps_i, pDepth);
            initialState.Pi = pDepth.P_gas0;
            [Pi_i, ~, dV_i, w_i, ~, convEq, staticState] = ...
                solveCavityPressure(pDepth.P_static, a_i, Kbend_i, Tpre_i, pDepth, initialState);
            [gain_i, convGain] = evaluateLocalP3GainFromState( ...
                pDepth.P_static, a_i, Kbend_i, Tpre_i, pDepth, staticState);

            depth_m(row) = depth_i;
            strain(row) = eps_i;
            P_external_abs_Pa(row) = pDepth.P_static;
            P_cavity_abs_Pa(row) = Pi_i;
            raw_gain(row) = gain_i;
            wstatic_mm(row) = 1e3 * w_i;
            dVstatic_over_V0(row) = dV_i / pDepth.V00;
            w_over_a(row) = abs(w_i) / a_i;
            w_over_h(row) = abs(w_i) / h_i;
            within_small_deflection(row) = w_over_a(row) <= 0.1 && w_over_h(row) <= 1.0;
            converged(row) = convEq && convGain;

            fprintf('%8.1f %6.3f  %12.2f  %10.2f   %8.4f  %7.4f  %6.3f %6.2f     %d\n', ...
                depth_i, eps_i, pDepth.P_static, Pi_i, gain_i, ...
                dVstatic_over_V0(row), w_over_a(row), w_over_h(row), ...
                within_small_deflection(row));
        end
    end

    knoTable = table(depth_m, strain, P_external_abs_Pa, P_cavity_abs_Pa, ...
        raw_gain, wstatic_mm, dVstatic_over_V0, w_over_a, w_over_h, ...
        within_small_deflection, converged);
    csvPath = fullfile(resultsDir, ['kno_depth_states_', versionTag, '.csv']);
    writetable(knoTable, csvPath);
    fprintf('Saved KNO deployment-depth table: %s\n', csvPath);
    fprintf(['Validity note: converged deep-water solutions are mechanistic ', ...
        'extrapolations when the valid flag is 0.\n']);
end

function pressureTable = compileMikePressureReport(eps_cases, p, resultsDir, versionTag)
    nCases = numel(eps_cases);
    P_external_HS = repmat(p.P_static, nCases, 1);
    P_cavity_HS = NaN(nCases, 1);
    HS_loss = NaN(nCases, 1);
    P_external_plus = repmat(p.P_static + p.dP0/2, nCases, 1);
    P_external_minus = repmat(p.P_static - p.dP0/2, nCases, 1);
    P_cavity_plus = NaN(nCases, 1);
    P_cavity_minus = NaN(nCases, 1);
    cavity_high_change = NaN(nCases, 1);
    cavity_low_change_signed = NaN(nCases, 1);
    cavity_low_change_magnitude = NaN(nCases, 1);
    high_cavity_transmission = NaN(nCases, 1);
    low_cavity_transmission = NaN(nCases, 1);
    high_to_low_change_ratio = NaN(nCases, 1);
    high_low_change_imbalance = NaN(nCases, 1);
    cavity_differential = NaN(nCases, 1);
    differential_transmission = NaN(nCases, 1);
    converged = false(nCases, 1);

    for i = 1:nCases
        [Dplate, ~, Tpre_eff, ~, a_load, ~, ~, ~] = ...
            evaluateMembraneState(eps_cases(i), p);
        initialState.Pi = p.P_gas0;
        [P_cavity_HS(i), ~, ~, ~, ~, convHS, staticState] = ...
            solveCavityPressure(p.P_static, a_load, Dplate, Tpre_eff, p, initialState);
        [P_cavity_plus(i), ~, ~, ~, ~, convPlus] = ...
            solveCavityPressure(P_external_plus(i), a_load, Dplate, Tpre_eff, p, staticState);
        [P_cavity_minus(i), ~, ~, ~, ~, convMinus] = ...
            solveCavityPressure(P_external_minus(i), a_load, Dplate, Tpre_eff, p, staticState);

        HS_loss(i) = P_external_HS(i) - P_cavity_HS(i);
        cavity_high_change(i) = P_cavity_plus(i) - P_cavity_HS(i);
        cavity_low_change_signed(i) = P_cavity_minus(i) - P_cavity_HS(i);
        cavity_low_change_magnitude(i) = -cavity_low_change_signed(i);
        high_cavity_transmission(i) = cavity_high_change(i) / (p.dP0/2);
        low_cavity_transmission(i) = cavity_low_change_magnitude(i) / (p.dP0/2);
        high_to_low_change_ratio(i) = cavity_high_change(i) / ...
            max(cavity_low_change_magnitude(i), eps);
        high_low_change_imbalance(i) = cavity_high_change(i) - ...
            cavity_low_change_magnitude(i);
        cavity_differential(i) = P_cavity_plus(i) - P_cavity_minus(i);
        differential_transmission(i) = cavity_differential(i) / p.dP0;
        converged(i) = convHS && convPlus && convMinus;
    end

    % Use kPa for absolute/common-mode pressures and hydrostatic shortfall;
    % retain Pa for the small applied and achieved differential pressures.
    P_external_HS_kPa = P_external_HS / 1e3;
    P_cavity_HS_kPa = P_cavity_HS / 1e3;
    HS_loss_kPa = HS_loss / 1e3;
    P_external_plus_kPa = P_external_plus / 1e3;
    P_external_minus_kPa = P_external_minus / 1e3;
    P_cavity_plus_kPa = P_cavity_plus / 1e3;
    P_cavity_minus_kPa = P_cavity_minus / 1e3;

    pressureTable = table(eps_cases(:), P_external_HS_kPa, P_cavity_HS_kPa, HS_loss_kPa, ...
        P_external_plus_kPa, P_external_minus_kPa, P_cavity_plus_kPa, P_cavity_minus_kPa, ...
        cavity_high_change, cavity_low_change_signed, cavity_low_change_magnitude, ...
        high_cavity_transmission, low_cavity_transmission, high_to_low_change_ratio, ...
        high_low_change_imbalance, cavity_differential, differential_transmission, converged, ...
        'VariableNames', {'strain', 'P_external_HS_kPa', 'P_cavity_HS_kPa', ...
        'HS_loss_kPa', 'P_external_plus_kPa', 'P_external_minus_kPa', ...
        'P_cavity_plus_kPa', 'P_cavity_minus_kPa', 'cavity_high_change_Pa', ...
        'cavity_low_change_signed_Pa', 'cavity_low_change_magnitude_Pa', ...
        'high_cavity_transmission', 'low_cavity_transmission', ...
        'high_to_low_change_ratio', 'high_low_change_imbalance_Pa', ...
        'cavity_differential_Pa', 'differential_transmission', 'converged'});

    fprintf('\n--- Mike meeting pressure-state report ---\n');
    fprintf('Depth = %.4f m; applied external hydrostatic pressure = %.6f kPa absolute.\n', ...
        p.baselineDepth_m, p.P_static / 1e3);
    fprintf('Applied differential = %.6f Pa, imposed symmetrically as +/- %.6f Pa.\n', ...
        p.dP0, p.dP0/2);
    fprintf(['eps      Pi,HS(kPa)  HS loss(kPa)  Pi,plus(kPa)  Pi,minus(kPa) ', ...
        'riseHigh(Pa) dropLow(Pa) gainHigh gainLow H/L ratio imbalance(Pa) ', ...
        'net dPi(Pa) netGain conv\n']);
    for i = 1:nCases
        fprintf(['%6.3f  %12.6f  %12.6f  %13.6f  %14.6f  %12.6f ', ...
            '%11.6f  %8.6f %7.6f %9.6f %13.6f %11.6f %8.6f   %d\n'], ...
            eps_cases(i), P_cavity_HS_kPa(i), HS_loss_kPa(i), P_cavity_plus_kPa(i), ...
            P_cavity_minus_kPa(i), cavity_high_change(i), cavity_low_change_magnitude(i), ...
            high_cavity_transmission(i), low_cavity_transmission(i), ...
            high_to_low_change_ratio(i), high_low_change_imbalance(i), ...
            cavity_differential(i), differential_transmission(i), converged(i));
    end

    csvPath = fullfile(resultsDir, ['mike_pressure_states_', versionTag, '.csv']);
    writetable(pressureTable, csvPath);
    fprintf('Saved meeting pressure table: %s\n', csvPath);
end

function printModelPointSummaryTable(eps_plot, eps_data, y_data, y_err, S_interface_model, T_cav_rel0, ...
    a_load_curve, a_geom_curve, Dplate_curve, Tpre_eff_curve, phi_pre_curve, resid, wres)

    fprintf('\n--- model summary at experimental strain values ---\n');
    fprintf(' eps_eff    data      model    T_rel0      a_load(mm)   a_geom(mm)   Kbend(Nm)   Tpre_eff(N/m)   phi_pre   resid    w_resid\n');

    for i = 1:numel(eps_data)
        Si_model  = interp1(eps_plot, S_interface_model, eps_data(i), 'linear');
        Ti_rel0   = interp1(eps_plot, T_cav_rel0, eps_data(i), 'linear');
        ai_load   = interp1(eps_plot, 1e3 * a_load_curve, eps_data(i), 'linear');
        ai_geom   = interp1(eps_plot, 1e3 * a_geom_curve, eps_data(i), 'linear');
        Di        = interp1(eps_plot, Dplate_curve, eps_data(i), 'linear');
        Tpre_effi = interp1(eps_plot, Tpre_eff_curve, eps_data(i), 'linear');
        phii      = interp1(eps_plot, phi_pre_curve, eps_data(i), 'linear');
        fprintf('%7.3f   %7.3f   %7.3f   %7.3f   %10.3f   %10.3f   %8.3e   %8.3e   %7.3f   %7.3f   %7.3f\n', ...
            eps_data(i), y_data(i), Si_model, Ti_rel0, ai_load, ai_geom, Di, Tpre_effi, phii, resid(i), wres(i));
    end
    fprintf('Measurement uncertainty used in weighted residuals: ');
    fprintf('%.3f ', y_err);
    fprintf('\n');
end

function diagnosticTable = runBucklingWeightedStrainDiagnostic(eps_cases, p, resultsDir, versionTag)
    pNominal = p;
    pNominal.strainMappingMode = 'nominal';
    pNominal.engStrains = p.nominalEngStrains;
    pWeighted = p;
    pWeighted.strainMappingMode = 'areaWeightedVolume';
    pWeighted.engStrains = p.nominalEngStrains;
    pWeighted.thicknessMode = 'incompressible';

    A_total = pi*p.strainWeightingOuterRadius_m^2;
    A_center = pi*p.strainWeightingCenterRadius_m^2;
    centerFraction = A_center/A_total;
    washerFraction = 1-centerFraction;
    n = numel(eps_cases);

    nominal_strain = eps_cases(:);
    effective_weighted_strain = NaN(n,1);
    lambda_z_center = NaN(n,1);
    lambda_z_weighted = NaN(n,1);
    baseline_thickness_mm = NaN(n,1);
    weighted_effective_thickness_mm = NaN(n,1);
    baseline_Tpre_Npm = NaN(n,1);
    weighted_Tpre_Npm = NaN(n,1);
    baseline_Kbend_Nm = NaN(n,1);
    weighted_Kbend_Nm = NaN(n,1);
    baseline_raw_gain = NaN(n,1);
    weighted_raw_gain = NaN(n,1);
    baseline_converged = false(n,1);
    weighted_converged = false(n,1);

    for i = 1:n
        eps_i = eps_cases(i);
        lambda_z_center(i) = (1+max(eps_i,0))^(-2);
        lambda_z_weighted(i) = washerFraction + centerFraction*lambda_z_center(i);
        effective_weighted_strain(i) = lambda_z_weighted(i)^(-1/2)-1;

        [~,~,~,~,~,hBase] = evaluatePlateGeometry(eps_i,pNominal);
        [~,~,~,~,~,hWeighted] = evaluatePlateGeometry(eps_i,pWeighted);
        baseline_thickness_mm(i) = 1e3*hBase;
        weighted_effective_thickness_mm(i) = 1e3*hWeighted;
        [baseline_Kbend_Nm(i),~,baseline_Tpre_Npm(i)] = evaluateMembraneState(eps_i,pNominal);
        [weighted_Kbend_Nm(i),~,weighted_Tpre_Npm(i)] = evaluateMembraneState(eps_i,pWeighted);
        [baseline_raw_gain(i),~,~,~,~,~,~,~,~,baseline_converged(i)] = ...
            evaluateCavityTransmission(eps_i,pNominal);
        [weighted_raw_gain(i),~,~,~,~,~,~,~,~,weighted_converged(i)] = ...
            evaluateCavityTransmission(eps_i,pWeighted);
    end

    center_area_fraction = repmat(centerFraction,n,1);
    washer_area_fraction = repmat(washerFraction,n,1);
    baseline_over_Rbare = baseline_raw_gain/p.R_bare;
    weighted_over_Rbare = weighted_raw_gain/p.R_bare;
    diagnosticTable = table(nominal_strain,effective_weighted_strain, ...
        center_area_fraction,washer_area_fraction,lambda_z_center,lambda_z_weighted, ...
        baseline_thickness_mm,weighted_effective_thickness_mm,baseline_Tpre_Npm, ...
        weighted_Tpre_Npm,baseline_Kbend_Nm,weighted_Kbend_Nm,baseline_raw_gain, ...
        weighted_raw_gain,baseline_over_Rbare,weighted_over_Rbare, ...
        baseline_converged,weighted_converged);
    writetable(diagnosticTable,fullfile(resultsDir, ...
        ['buckling_weighted_strain_diagnostic_',versionTag,'.csv']));

    fprintf('\n--- area-weighted volume-conserving strain diagnostic ---\n');
    fprintf('Area fractions: washer = %.4f, exposed center = %.4f\n', ...
        washerFraction,centerFraction);
    fprintf('eps_nom   eps_eff   Tpre_base   Tpre_weight   gain_base   gain_weight\n');
    for i = 1:n
        fprintf('%7.3f   %7.4f    %9.3f    %11.3f    %9.4f    %11.4f\n', ...
            nominal_strain(i),effective_weighted_strain(i),baseline_Tpre_Npm(i), ...
            weighted_Tpre_Npm(i),baseline_raw_gain(i),weighted_raw_gain(i));
    end

    figure('Position',[80 80 2100 760],'Tag','bucklingWeightedDiagnostic');
    tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
    nexttile; hold on;
    plot(nominal_strain,nominal_strain,'--k','DisplayName','Nominal mapping');
    plot(nominal_strain,effective_weighted_strain,'-o','DisplayName','Area-weighted mapping');
    xlabel('Nominal Engineering Strain (-)'); ylabel('Strain Used in Model (-)');
    title('Equivalent Strain'); legend('Location','northwest');

    nexttile; hold on;
    plot(nominal_strain,baseline_thickness_mm,'-o','DisplayName','Current baseline');
    plot(nominal_strain,weighted_effective_thickness_mm,'-d','DisplayName','Area-weighted effective');
    xlabel('Nominal Engineering Strain (-)'); ylabel('Thickness (mm)');
    title('Thickness Interpretation'); legend('Location','southwest');

    nexttile; hold on;
    plot(nominal_strain,baseline_raw_gain,'-o','DisplayName','Nominal strain, locked tension');
    plot(nominal_strain,weighted_raw_gain,'-d','DisplayName','Weighted strain, locked tension');
    yline(p.R_bare,'--','DisplayName',sprintf('R_{bare}=%.2f',p.R_bare));
    yline(1,':k','DisplayName','Unity raw transmission');
    xlabel('Nominal Engineering Strain (-)'); ylabel('Raw Differential Transmission (-)');
    title('Locked-Tension Consequence'); legend('Location','best');
    sgtitle({'Area-Weighted Volume-Conservation Diagnostic', ...
        'Diagnostic only: no measured buckle volume or fitted strain offset'});
    saveCurrentFigure(resultsDir,['86_buckling_weighted_strain_',versionTag]);
end

function comparisonTable = plotCandidateModelComparison( ...
        eps_plot, eps_data, y_data, y_err, p, resultsDir, versionTag)
    % Direct comparison on the experiment's bare-port-normalized basis.
    % These are physics cases/bounds, not fitted curves.
    % Common presentation coordinate: effective strain. The conservative
    % full-tension bound is evaluated with diameter-derived nominal strain,
    % then re-parameterized onto the same effective-strain axis.
    [effective_locked_over_Rbare, ~] = ...
        runBareReferencedModelResponse(eps_plot, p);

    centerAreaFraction = (p.strainWeightingCenterRadius_m / ...
        p.strainWeightingOuterRadius_m)^2;
    lambdaZEffective = (1 + eps_plot).^(-2);
    eps_nominal_plot = sqrt(centerAreaFraction ./ ...
        (lambdaZEffective - (1-centerAreaFraction))) - 1;

    pNominal = p;
    pNominal.strainMappingMode = 'nominal';
    pNominal.engStrains = p.nominalEngStrains;
    [full_nominal_locked_over_Rbare, ~] = ...
        runBareReferencedModelResponse(eps_nominal_plot, pNominal);

    relaxed_Cv_half_over_Rbare = NaN(size(eps_plot));
    relaxed_converged = false(size(eps_plot));
    for i = 1:numel(eps_plot)
        eps_nominal_i = eps_nominal_plot(i);
        [Kbend,E,~,~,a,~,~,~] = ...
            evaluateMembraneState(eps_nominal_i,pNominal);
        [~,~,~,~,~,h] = evaluatePlateGeometry(eps_nominal_i,pNominal);
        plusState = solveConservativeNonlinearState( ...
            pNominal.P_static+pNominal.dP0/2,a,Kbend,0,E,h,1/2,1,pNominal);
        minusState = solveConservativeNonlinearState( ...
            pNominal.P_static-pNominal.dP0/2,a,Kbend,0,E,h,1/2,1,pNominal);
        rawGain = (plusState.Pi-minusState.Pi)/pNominal.dP0;
        relaxed_Cv_half_over_Rbare(i) = rawGain/pNominal.R_bare;
        relaxed_converged(i) = plusState.converged && minusState.converged;
    end
    relaxed_Cv_half_over_Rbare(~relaxed_converged) = NaN;

    strain = eps_plot(:);
    nominal_strain = eps_nominal_plot(:);
    full_nominal_locked_tension_over_Rbare = ...
        full_nominal_locked_over_Rbare(:);
    effective_locked_tension_over_Rbare = effective_locked_over_Rbare(:);
    relaxed_geometric_Cv_half_over_Rbare = relaxed_Cv_half_over_Rbare(:);
    relaxed_Cv_half_converged = relaxed_converged(:);
    comparisonTable = table(strain,nominal_strain, ...
        full_nominal_locked_tension_over_Rbare, ...
        effective_locked_tension_over_Rbare, ...
        relaxed_geometric_Cv_half_over_Rbare, ...
        relaxed_Cv_half_converged);
    writetable(comparisonTable,fullfile(resultsDir, ...
        ['candidate_model_comparison_',versionTag,'.csv']));

    figure; hold on;
    hData = errorbar(eps_data,y_data,y_err,'ko','LineWidth',1.8, ...
        'MarkerFaceColor','k','MarkerSize',8, ...
        'DisplayName','Experiment (bare-port normalized)');
    hFull = plot(eps_plot,full_nominal_locked_over_Rbare,'-','LineWidth',2.7, ...
        'DisplayName','Full nominal-strain tension (pessimistic bound)');
    hEffective = plot(eps_plot,effective_locked_over_Rbare,'-','LineWidth',2.7, ...
        'DisplayName','Effective-strain locked tension (candidate)');
    hRelaxed = plot(eps_plot,relaxed_Cv_half_over_Rbare,'-','LineWidth',2.7, ...
        'DisplayName','Locally relaxed + geometric tension, $C_V=1/2$');
    hUnity = yline(1,':k','LineWidth',1.5, ...
        'DisplayName','Bare-port reference');
    xlabel('Effective Engineering Strain, \epsilon_{eff} (-)');
    ylabel('Bare-Port-Normalized Response (-)');
    title('Tension-Model Bounds and Empirical Response');
    subtitle('Physics cases and bounds; no parameters fitted to the response data');
    xlim([eps_plot(1),eps_plot(end)]);
    ylim([0.25,1.48]);
    grid on; box on; formatAxes(gca);
    legend([hData,hFull,hEffective,hRelaxed,hUnity], ...
        'Location','southoutside', ...
        'NumColumns',2, ...
        'Interpreter','latex');
    if isfield(p,'publicationMode') && p.publicationMode
        outName = ['03_candidate_model_comparison_',versionTag];
    else
        outName = ['87_candidate_model_comparison_',versionTag];
    end
    saveCurrentFigure(resultsDir,outName);

    fprintf('\n--- candidate model comparison at experimental strains ---\n');
    fprintf(['eps_eff eps_nom experiment fullT/Rbare effectiveT/Rbare ', ...
        'relaxed Cv=1/2/Rbare\n']);
    for i = 1:numel(eps_data)
        fprintf('%7.4f %7.4f   %7.4f      %7.4f         %7.4f              %7.4f\n', ...
            eps_data(i),p.experimentalNominalStrains(i),y_data(i), ...
            interp1(eps_plot,full_nominal_locked_over_Rbare,eps_data(i),'linear'), ...
            interp1(eps_plot,effective_locked_over_Rbare,eps_data(i),'linear'), ...
            interp1(eps_plot,relaxed_Cv_half_over_Rbare,eps_data(i),'linear'));
    end
end

function [Tcav, Dplate, E, Tpre_eff, phi_pre, a_load, a_geom, ...
          D_installed_mm, eps_pre, converged, iter, solverStateOut] = evaluateCavityTransmission(eps_query, p, solverStateIn)

    [Dplate, E, Tpre_eff, phi_pre, a_load, a_geom, D_installed_mm, eps_pre] = evaluateMembraneState(eps_query, p);

    if nargin < 3 || isempty(solverStateIn)
        solverStateIn.Pi = p.P_gas0;
    end

    % Stage 1: establish the common-mode hydrostatic equilibrium of the
    % sealed cavity before applying the differential perturbation.
    [~, ~, ~, ~, iterStatic, convStatic, staticState] = ...
        solveCavityPressure(p.P_static, a_load, Dplate, Tpre_eff, p, solverStateIn);

    % Stage 2: perturb the two external surfaces symmetrically about the
    % hydrostatically loaded state. The trapped-gas inventory is unchanged.
    dPsens_plus  = sensorDifferentialResponse(+p.dP0, a_load, Dplate, Tpre_eff, p, staticState);
    dPsens_minus = sensorDifferentialResponse(-p.dP0, a_load, Dplate, Tpre_eff, p, staticState);
    Tcav = (dPsens_plus - dPsens_minus) / (2 * p.dP0);

    [~, ~, ~, ~, iterPerturb, convPerturb, solverStateOut] = ...
        solveCavityPressure(p.P_static + p.dP0/2, a_load, Dplate, Tpre_eff, p, staticState);
    converged = convStatic && convPerturb && isfinite(Tcav);
    iter = iterStatic + iterPerturb;
end

function [eps_pre, lambda_pre, D_installed_mm] = evaluateInstalledPrestretchState(eps_query, p)
    D_installed_mm = interp1(p.engStrains, p.scaleDiameters_mm, eps_query, 'linear', 'extrap');
    eps_input = max(eps_query, 0);
    if ~isfield(p,'strainMappingMode')
        p.strainMappingMode = 'nominal';
    end

    switch lower(p.strainMappingMode)
        case {'effectivedirect','nominal'}
            % 'nominal' is retained as a legacy alias for archived diagnostics.
            % The v37 baseline passes effective strain via 'effectiveDirect'.
            eps_pre = eps_input;
            lambda_pre = 1+eps_pre;

        case 'areaweightedvolume'
            A_total = pi*p.strainWeightingOuterRadius_m^2;
            A_center = pi*p.strainWeightingCenterRadius_m^2;
            if A_total <= 0 || A_center < 0 || A_center > A_total
                error('Invalid areas for weighted strain mapping.');
            end
            lambda_z_center = (1+eps_input)^(-2);
            lambda_z_effective = ((A_total-A_center) + ...
                A_center*lambda_z_center)/A_total;
            lambda_pre = lambda_z_effective^(-1/2);
            eps_pre = lambda_pre-1;

        otherwise
            error('Unknown strainMappingMode: %s',p.strainMappingMode);
    end
end

function [a_load, a_geom, D_installed_mm, eps_pre, lambda_pre, t_plate] = evaluatePlateGeometry(eps_query, p)
    [eps_pre, lambda_pre, D_installed_mm] = evaluateInstalledPrestretchState(eps_query, p);
    a_geom = 0.5 * D_installed_mm * 1e-3;
    switch lower(p.loadedRadiusMode)
        case 'fixed'
            radiusScale = 1.0;
        case 'statetable'
            radiusScale = interp1(p.radiusStateStrains, p.radiusStateScales, ...
                eps_query, 'linear', 'extrap');
        otherwise
            error('Unknown loadedRadiusMode: %s', p.loadedRadiusMode);
    end
    a_load = p.r_forced_m * radiusScale;
    switch lower(p.thicknessMode)
    case 'poisson'
        t_plate = p.t_plate0 * lambda_pre ^ (-2 * p.nu_plate / (1 - p.nu_plate));
    case 'incompressible'
        t_plate = p.t_plate0 / lambda_pre^2;
    case 'constant'
        t_plate = p.t_plate0;
    otherwise
        error('Unknown thicknessMode: %s', p.thicknessMode);
    end
end

function [Dplate, E, Tpre_eff, phi_pre, a_load, a_geom, D_installed_mm, eps_pre] = evaluateMembraneState(eps_query, p)
    [a_load, a_geom, D_installed_mm, eps_pre, ~, t_plate] = evaluatePlateGeometry(eps_query, p);

    if p.useStrainDependentE
        E = p.E_plate0 * (1 + p.c1 * eps_pre + p.c2 * eps_pre ^ 2);
    else
        E = p.E_plate0;
    end

    Dplate = E * t_plate ^ 3 / (12 * (1 - p.nu_plate ^ 2));

    if ~isfield(p, 'preTensionThicknessMode')
        p.preTensionThicknessMode = 'current';
    end

    switch lower(p.preTensionThicknessMode)
        case 'current'
            t_for_Tpre = t_plate;
        case 'reference'
            t_for_Tpre = p.t_plate0;
        otherwise
            error('Unknown preTensionThicknessMode: %s', p.preTensionThicknessMode);
    end

    if ~isfield(p, 'preTensionCouplingMode')
    p.preTensionCouplingMode = 'direct';
    end
    
    switch lower(p.preTensionCouplingMode)
        case 'direct'
            kT_effective = p.kT0;
    
        case 'radiuspower'
            radiusRatio = min(max(a_load / max(a_geom, 1e-12), 0), 1);
            kT_effective = p.kT0_global * radiusRatio ^ p.radiusPowerExponent;
    
        case 'tensioncap'
            kT_effective = 1.0;  % compute full ideal tension first, then cap it
    
        otherwise
            error('Unknown preTensionCouplingMode: %s', p.preTensionCouplingMode);
    end

    if ~isfield(p, 'preTensionLaw')
        p.preTensionLaw = 'linearEquibiaxial';
    end

    switch lower(p.preTensionLaw)
        case 'linearequibiaxial'
            % Plane-stress equibiaxial linear-elastic approximation.
            Tpre_ideal = E * t_for_Tpre * eps_pre / (1 - p.nu_plate);

        case 'neohookean'
            % Incompressible neo-Hookean equibiaxial membrane. With
            % lambda_1=lambda_2=lambda and lambda_3=lambda^-2, the Cauchy
            % stress is mu(lambda^2-lambda^-4). Multiplication by the
            % current thickness t0/lambda^2 gives the force resultant below.
            lambda_pre = 1 + eps_pre;
            mu = E / (2 * (1 + p.nu_plate));
            Tpre_ideal = mu * p.t_plate0 * (1 - lambda_pre ^ (-6));

        otherwise
            error('Unknown preTensionLaw: %s', p.preTensionLaw);
    end

    Tpre_uncapped = kT_effective * Tpre_ideal;
    
    if strcmpi(p.preTensionCouplingMode, 'tensionCap')
        if ~isfield(p, 'Tpre_cap_Npm')
            p.Tpre_cap_Npm = Inf;
        end
        Tpre_uncapped = min(Tpre_uncapped, p.Tpre_cap_Npm);
    end
    
    switch lower(p.preTensionMode)
        case 'ideal'
            phi_pre = 1.0;
            Tpre_eff = Tpre_uncapped;
        case 'activated'
            phi_pre = 1 - exp(-(eps_pre / max(p.eps_char, 1e-12)) ^ 2);
            Tpre_eff = phi_pre * Tpre_uncapped;
        case 'none'
            phi_pre = 0.0;
            Tpre_eff = 0.0;
        otherwise
            error('Unknown preTensionMode: %s', p.preTensionMode);
    end
end

function dPsens = sensorDifferentialResponse(dPsurf, a_plate, Dplate, Tpre_eff, p, solverStateIn)
    P3 = p.P_static + dPsurf / 2;
    P4 = p.P_static - dPsurf / 2;
    [P1, ~, ~, ~, ~, conv1, ~] = solveCavityPressure(P3, a_plate, Dplate, Tpre_eff, p, solverStateIn);
    [P2, ~, ~, ~, ~, conv2, ~] = solveCavityPressure(P4, a_plate, Dplate, Tpre_eff, p, solverStateIn);
    if ~(conv1 && conv2)
        dPsens = NaN;
    else
        dPsens = P1 - P2;
    end
end

function [pressureGain, pressureGain_relLocal, pressureGain_relZeroStrain, S_pressure, pressureConverged] = ...
    runPressureLevelSweep(eps_pressure_cases, P3_sweep, P3_offset_sweep, p)

    pressureGain = zeros(numel(eps_pressure_cases), numel(P3_sweep));
    pressureGain_relLocal = zeros(size(pressureGain));
    pressureConverged = zeros(size(pressureGain));

    for ie = 1:numel(eps_pressure_cases)
        eps_i = eps_pressure_cases(ie);
        [Dplate_i, ~, Tpre_eff_i, ~, a_load_i, ~, ~, ~] = evaluateMembraneState(eps_i, p);

        for ip = 1:numel(P3_sweep)
            [pressureGain(ie, ip), pressureConverged(ie, ip)] = evaluateLocalP3Gain(P3_sweep(ip), a_load_i, Dplate_i, Tpre_eff_i, p);
        end

        zeroOffsetIdx = find(abs(P3_offset_sweep) < 1e-12, 1);
        if isempty(zeroOffsetIdx)
            [~, zeroOffsetIdx] = min(abs(P3_offset_sweep));
        end

        localRef = pressureGain(ie, zeroOffsetIdx);
        if abs(localRef) < 1e-12 || ~isfinite(localRef)
            warning('Local pressure-gain reference is near zero or invalid at eps = %.3f.', eps_i);
            localRef = NaN;
        end
        pressureGain_relLocal(ie, :) = pressureGain(ie, :) ./ localRef;
    end

    [Dplate_0, ~, Tpre_eff_0, ~, a_load_0, ~, ~, ~] = evaluateMembraneState(0.0, p);
    [G0_pressure, ~] = evaluateLocalP3Gain(p.P_static, a_load_0, Dplate_0, Tpre_eff_0, p);
    if abs(G0_pressure) < 1e-12 || ~isfinite(G0_pressure)
        error('Zero-strain pressure-gain reference is too small or invalid. Inspect pressureGain before normalization.');
    end

    pressureGain_relZeroStrain = pressureGain ./ G0_pressure;
    S_pressure = pressureGain ./ p.R_bare;
end

function printPressureSweepSummaryTable(eps_pressure_cases, P3_offset_sweep, pressureGain, pressureGain_relLocal, pressureGain_relZeroStrain, S_pressure, pressureConverged)
    fprintf('\n--- pressure-level dependence check ---\n');
    fprintf('Rows are effective-strain cases. Columns are P3 offsets from P_ref.\n\n');
    fprintf('P3 offsets (Pa): '); fprintf('%10.1f ', P3_offset_sweep); fprintf('\n');
    printPressureMatrix('Raw local gain dP1/dP3:', eps_pressure_cases, pressureGain);
    printPressureMatrix('Gain normalized within each pre-strain case to P3 offset = 0:', eps_pressure_cases, pressureGain_relLocal);
    printPressureMatrix('Gain normalized by zero-strain, zero-offset gain:', eps_pressure_cases, pressureGain_relZeroStrain);
    printPressureMatrix('Interface-relative sensitivity prediction across pressure levels:', eps_pressure_cases, S_pressure);
    fprintf('\nPressure sweep converged points: %d / %d\n', nnz(pressureConverged(:) > 0.5), numel(pressureConverged));
end

function printPressureMatrix(label, eps_cases, M)
    fprintf('\n%s\n', label);
    for ie = 1:numel(eps_cases)
        fprintf('eps = %.3f: ', eps_cases(ie));
        fprintf('%10.5f ', M(ie, :));
        fprintf('\n');
    end
end

function [Glocal, convergedBoth] = evaluateLocalP3Gain(P3_center, a_plate, Dplate, Tpre_eff, p)
    dP = p.dP_local;
    solverState.Pi = p.P_gas0;
    P3_plus  = P3_center + dP / 2;
    P3_minus = P3_center - dP / 2;

    [P1_plus, ~, ~, ~, ~, convPlus, ~] = solveCavityPressure(P3_plus, a_plate, Dplate, Tpre_eff, p, solverState);
    [P1_minus, ~, ~, ~, ~, convMinus, ~] = solveCavityPressure(P3_minus, a_plate, Dplate, Tpre_eff, p, solverState);

    Glocal = (P1_plus - P1_minus) / dP;
    convergedBoth = convPlus && convMinus && isfinite(Glocal);
end

function [Pi, Vi, dV, wmax, iter, converged, solverStateOut, dV_raw, q_last] = solveCavityPressure(Psurf, a_plate, Dplate, Tpre_eff, p, solverStateIn)
    if ~isfield(p, 'solverMode')
        p.solverMode = 'fixedPoint';
    end

    switch lower(p.solverMode)
        case 'fixedpoint'
            [Pi, Vi, dV, wmax, iter, converged, solverStateOut, dV_raw, q_last] = ...
                solveCavityPressure_fixedPoint(Psurf, a_plate, Dplate, Tpre_eff, p, solverStateIn);
        case 'root'
            [Pi, Vi, dV, wmax, iter, converged, solverStateOut, dV_raw, q_last] = ...
                solveCavityPressure_root(Psurf, a_plate, Dplate, Tpre_eff, p, solverStateIn);
        otherwise
            error('Unknown solverMode: %s', p.solverMode);
    end
end

function [Pi, Vi, dV, wmax, iter, converged, solverStateOut, dV_raw, q_last] = solveCavityPressure_fixedPoint(Psurf, a_plate, Dplate, Tpre_eff, p, solverStateIn)
    if nargin < 6 || isempty(solverStateIn) || ~isfield(solverStateIn, 'Pi') || ~isfinite(solverStateIn.Pi)
        Pi = p.P_gas0;
    else
        Pi = solverStateIn.Pi;
    end

    converged = false;
    Vi = p.V00;
    dV = 0;
    dV_raw = 0;
    wmax = 0;
    q_last = 0;

    for iter = 1:p.maxIter
        q_membrane = Psurf - Pi;
        [Vi_new, dV_new, wmax_new, dV_raw_new] = loadedCavityState_fromMembraneLoad(q_membrane, p.V00, a_plate, Dplate, Tpre_eff, p);
        Pi_new_raw = p.P_gas0 * (p.V00 / Vi_new) ^ p.gasExponent;
        Pi_new = (1 - p.relax) * Pi + p.relax * Pi_new_raw;

        Vi = Vi_new;
        dV = dV_new;
        dV_raw = dV_raw_new;
        wmax = wmax_new;
        q_last = q_membrane;

        if abs(Pi_new - Pi) < p.tolP
            Pi = Pi_new;
            converged = true;
            solverStateOut.Pi = Pi;
            return;
        end
        Pi = Pi_new;
    end

    solverStateOut.Pi = Pi;
end

function [Pi, Vi, dV, wmax, iter, converged, solverStateOut, dV_raw, q_last] = solveCavityPressure_root(Psurf, a_plate, Dplate, Tpre_eff, p, ~)
    residual = @(Pi_trial) cavityPressureResidual(Pi_trial, Psurf, a_plate, Dplate, Tpre_eff, p);

    % For the current sign convention, Pi is physically between the trapped-
    % gas fill pressure and the applied external surface pressure.
    % Bisection is intentionally used instead of fzero so iteration count is
    % deterministic and convergence behavior is easy to audit.
    Pi_lo = min(p.P_gas0, Psurf);
    Pi_hi = max(p.P_gas0, Psurf);
    f_lo = residual(Pi_lo);
    f_hi = residual(Pi_hi);

    converged = false;

    if ~isfinite(f_lo) || ~isfinite(f_hi)
        [Pi, Vi, dV, wmax, iter, solverStateOut, dV_raw, q_last] = returnFailedRoot(p);
        return;
    end

    if abs(f_lo) < p.tolP
        Pi = Pi_lo;
        converged = true;
        iter = 0;
    elseif abs(f_hi) < p.tolP
        Pi = Pi_hi;
        converged = true;
        iter = 0;
    elseif sign(f_lo) == sign(f_hi)
        [Pi, Vi, dV, wmax, iter, solverStateOut, dV_raw, q_last] = returnFailedRoot(p);
        return;
    else
        for iter = 1:p.maxIter
            Pi_mid = 0.5 * (Pi_lo + Pi_hi);
            f_mid = residual(Pi_mid);

            if ~isfinite(f_mid)
                break;
            end

            if abs(f_mid) < p.tolP || abs(Pi_hi - Pi_lo) < p.tolP
                Pi = Pi_mid;
                converged = true;
                break;
            end

            if sign(f_mid) == sign(f_lo)
                Pi_lo = Pi_mid;
                f_lo = f_mid;
            else
                Pi_hi = Pi_mid;
            end
        end

        if ~converged
            Pi = 0.5 * (Pi_lo + Pi_hi);
        end
    end

    q_last = Psurf - Pi;
    [Vi, dV, wmax, dV_raw] = loadedCavityState_fromMembraneLoad(q_last, p.V00, a_plate, Dplate, Tpre_eff, p);
    solverStateOut.Pi = Pi;
end

function [Pi, Vi, dV, wmax, iter, solverStateOut, dV_raw, q_last] = returnFailedRoot(p)
    Pi = NaN;
    Vi = NaN;
    dV = NaN;
    wmax = NaN;
    iter = 0;
    dV_raw = NaN;
    q_last = NaN;
    solverStateOut.Pi = p.P_gas0;
end

function F = cavityPressureResidual(Pi_trial, Psurf, a_plate, Dplate, Tpre_eff, p)
    q_membrane = Psurf - Pi_trial;
    [Vi_trial, ~, ~, ~] = loadedCavityState_fromMembraneLoad(q_membrane, p.V00, a_plate, Dplate, Tpre_eff, p);
    Pi_from_volume = p.P_gas0 * (p.V00 / Vi_trial) ^ p.gasExponent;
    F = Pi_trial - Pi_from_volume;
end

function [Vi, dV, wmax, dV_raw] = loadedCavityState_fromMembraneLoad(q_membrane, V00, a, Dplate, Tpre_eff, p)
    Dplate_safe = max(Dplate, 1e-30);
    T_safe = max(Tpre_eff, 0);

    % Evaluate the radial deflection profile from the same analytical
    % plate solution used to determine the center deflection.
    rho_grid = linspace(0, 1, 1001);
    r_grid   = a * rho_grid;
    
    if T_safe <= 0
        % Zero-tension limit: uniformly loaded, clamped circular plate.
        w_profile = (q_membrane * a^4 / (64 * Dplate_safe)) .* ...
                    (1 - rho_grid.^2).^2;
    else
        alpha = a * sqrt(T_safe / Dplate_safe);
    
        if alpha < 1e-2
            % Avoid cancellation in the Bessel expression near alpha = 0.
            w_profile = (q_membrane * a^4 / (64 * Dplate_safe)) .* ...
                        (1 - rho_grid.^2).^2;
        else
            w_profile = (q_membrane * a^4 / Dplate_safe) .* ...
                ((1 - rho_grid.^2) ./ (4 * alpha^2) + ...
                (besseli(0, alpha .* rho_grid) - besseli(0, alpha)) ./ ...
                (2 * alpha^3 * besseli(1, alpha)));
        end
    end
    
    wmax   = w_profile(1);
    dV_raw = 2 * pi * trapz(r_grid, w_profile .* r_grid);

    if isfield(p, 'useVolumeCap') && p.useVolumeCap
        dV_cap = p.dV_cap_fraction * V00;
        dV = dV_cap * tanh(dV_raw / dV_cap);
    else
        dV = dV_raw;
    end

    minVolume = max(p.minVolumeFraction * V00, realmin);
    Vi = max(V00 - dV, minVolume);
end

function hBareMean = plotBareReference(bareColor)
    hBareMean = yline(1.0, '-', 'Color', bareColor, 'LineWidth', 1.8, 'DisplayName', 'Bare-Port Reference');
    yline(1.04, '--', 'Color', bareColor, 'LineWidth', 1.2, 'HandleVisibility', 'off');
    yline(0.96, '--', 'Color', bareColor, 'LineWidth', 1.2, 'HandleVisibility', 'off');
end

function [h75, h80, h85, h90] = plotExperimentalData(eps_data, y_data, y_err, c75, c80, c85, c90, markerSize, errorLineWidth)
    h75 = errorbar(eps_data(4), y_data(4), y_err(4), 's', 'Color', c75, 'MarkerFaceColor', c75, 'MarkerEdgeColor', c75, 'MarkerSize', markerSize, 'LineWidth', errorLineWidth, 'CapSize', 8, 'DisplayName', '75%');
    h80 = errorbar(eps_data(3), y_data(3), y_err(3), 's', 'Color', c80, 'MarkerFaceColor', c80, 'MarkerEdgeColor', c80, 'MarkerSize', markerSize, 'LineWidth', errorLineWidth, 'CapSize', 8, 'DisplayName', '80%');
    h85 = errorbar(eps_data(2), y_data(2), y_err(2), 's', 'Color', c85, 'MarkerFaceColor', c85, 'MarkerEdgeColor', c85, 'MarkerSize', markerSize, 'LineWidth', errorLineWidth, 'CapSize', 8, 'DisplayName', '85%');
    h90 = errorbar(eps_data(1), y_data(1), y_err(1), 's', 'Color', c90, 'MarkerFaceColor', c90, 'MarkerEdgeColor', c90, 'MarkerSize', markerSize, 'LineWidth', errorLineWidth, 'CapSize', 8, 'DisplayName', '90%');
end

function plotMaterialFitOverlays(eps_plot, eps_data, y_data, y_err, p, E_sweep, nu_sweep, t_sweep_m, t_sweep_mm, c75, c80, c85, c90, markerSize, errorLineWidth, bareColor, resultsDir, versionTag) %#ok<INUSD>
    figure; hold on;
    for iE = 1:numel(E_sweep)
        pSweep = p;
        pSweep.E_plate0 = E_sweep(iE);
        [S_E, ~] = runBareReferencedModelResponse(eps_plot, pSweep);
        plot(eps_plot, S_E, '-', 'LineWidth', 2.2, 'DisplayName', sprintf('E = %.1e Pa', E_sweep(iE)));
    end
    xlabel('Effective Engineering Strain, \epsilon_{eff} (-)');
    ylabel('R_{cav}(\epsilon)/R_{bare} (-)');
    title('Modeled Young''s-Modulus Sensitivity');
    xlim([eps_plot(1) eps_plot(end)]); grid on; box on; formatAxes(gca); legend('Location', 'northwest');
    saveCurrentFigure(resultsDir, ['14_E_sweep_model_only_', versionTag]);

    figure; hold on;
    for inu = 1:numel(nu_sweep)
        pSweep = p;
        pSweep.nu_plate = nu_sweep(inu);
        [S_nu, ~] = runBareReferencedModelResponse(eps_plot, pSweep);
        plot(eps_plot, S_nu, '-', 'LineWidth', 2.2, 'DisplayName', sprintf('\\nu = %.2f', nu_sweep(inu)));
    end
    xlabel('Effective Engineering Strain, \epsilon_{eff} (-)');
    ylabel('R_{cav}(\epsilon)/R_{bare} (-)');
    title('Modeled Poisson-Ratio Sensitivity');
    xlim([eps_plot(1) eps_plot(end)]); grid on; box on; formatAxes(gca); legend('Location', 'northwest');
    saveCurrentFigure(resultsDir, ['15_nu_sweep_model_only_', versionTag]);

    figure; hold on;
    for it = 1:numel(t_sweep_m)
        pSweep = p;
        pSweep.t_plate0 = t_sweep_m(it);
        [S_t, ~] = runBareReferencedModelResponse(eps_plot, pSweep);
        plot(eps_plot, S_t, '-', 'LineWidth', 2.2, 'DisplayName', sprintf('h_0 = %.3f mm', t_sweep_mm(it)));
    end
    xlabel('Effective Engineering Strain, \epsilon_{eff} (-)');
    ylabel('R_{cav}(\epsilon)/R_{bare} (-)');
    title('Modeled Membrane-Thickness Sensitivity');
    xlim([eps_plot(1) eps_plot(end)]); grid on; box on; formatAxes(gca); legend('Location', 'southeast');
    saveCurrentFigure(resultsDir, ['16_thickness_sweep_model_only_', versionTag]);
end

function plotKT0FitOverlay(eps_plot, eps_data, y_data, y_err, p, kT0_sweep, c75, c80, c85, c90, markerSize, errorLineWidth, bareColor, resultsDir, versionTag)
    figure; hold on;
    plotExperimentalData(eps_data, y_data, y_err, c75, c80, c85, c90, markerSize, errorLineWidth);

    for ik = 1:numel(kT0_sweep)
        pSweep = p;
        pSweep.preTensionMode = 'ideal';
        pSweep.kT0 = kT0_sweep(ik);
        [S_kT0, ~] = runBareReferencedModelResponse(eps_plot, pSweep);
        invalid = ~isfinite(S_kT0) | S_kT0 < 0.5 | S_kT0 > 2.0;
        S_kT0(invalid) = NaN;
        plot(eps_plot, S_kT0, '-', 'LineWidth', 2.2, 'DisplayName', sprintf('k_{T0} = %.4f', kT0_sweep(ik)));
    end

    addFitOverlayLabels('Effect of Effective Pre-Tension Realization on Model Fit', bareColor, [0.85 1.70]);
    saveCurrentFigure(resultsDir, ['17_kT0_sweep_fit_overlay_', versionTag]);
end

function addFitOverlayLabels(titleText, bareColor, yLimits)
    yline(1.0, '-', 'Color', bareColor, 'LineWidth', 1.8, 'DisplayName', 'Bare-Port Reference');
    xlabel('Effective Engineering Strain, \epsilon_{eff} (-)');
    ylabel('Normalized Sensitivity, S/S_0 (-)');
    title(titleText);
    xlim([eps_plot(1) eps_plot(end)]); %#ok<FNCOLND>
    ylim(yLimits);
    legend('Location', 'northwest');
    grid on; box on; formatAxes(gca);
end

function plotMechanicsDiagnostic(eps_plot, Dplate_curve, Tpre_eff_curve, ...
    a_load_curve, p, resultsDir, versionTag)

    % Recompute strain-dependent thickness for the same baseline model state.
    % Kbend is flexural rigidity only; it is not the membrane's total resistance
    % to transverse deflection. The latter also includes Tpre_eff, which grows
    % with pre-strain and dominates the strained baseline response.
    
    % Since Dplate_curve is already computed, normalize it directly.
    Dplate_norm = Dplate_curve ./ Dplate_curve(1);

    % Use the bending relation Kbend ~ E*h^3/[12(1-nu^2)].
    % With constant E and nu, h/h0 can be inferred from Kbend/Kbend0.
    h_norm = Dplate_norm.^(1/3);

    % Effective transverse stiffness is defined from the incremental center-
    % deflection compliance under a unit pressure load. Unlike Kbend, this
    % response quantity includes both flexural rigidity and in-plane tension.
    q_unit = 1.0; % Pa
    Ktrans_curve = NaN(size(eps_plot)); % Pa/m = N/m^3
    for i = 1:numel(eps_plot)
        [~, ~, w_unit] = loadedCavityState_fromMembraneLoad(q_unit, p.V00, ...
            a_load_curve(i), Dplate_curve(i), Tpre_eff_curve(i), p);
        if isfinite(w_unit) && abs(w_unit) > 0
            Ktrans_curve(i) = q_unit / abs(w_unit);
        end
    end
    Ktrans_norm = Ktrans_curve ./ Ktrans_curve(1);

    figure; hold on;

    yyaxis left;

    plot(eps_plot, h_norm, '-', ...
        'LineWidth', 2.5, ...
        'DisplayName', 'h(\epsilon) / h_0');

    plot(eps_plot, Dplate_norm, '--', ...
        'LineWidth', 2.5, ...
        'DisplayName', 'K_{bend}(\epsilon) / K_{bend}(0) (flexural rigidity)');

    xlabel('Effective Engineering Strain, \epsilon_{eff} (-)');
    ylabel('Normalized Thickness or Flexural Rigidity (-)');

    yyaxis right;
    plot(eps_plot, Ktrans_norm, '-.', ...
        'LineWidth', 2.8, ...
        'DisplayName', 'K_{trans}(\epsilon) / K_{trans}(0)');
    ylabel('Normalized Effective Transverse Stiffness (-)');

    title('Thickness, Flexural Rigidity, and Effective Transverse Stiffness');
    subtitle(sprintf(['K_{trans}=q/w_{max} includes bending and in-plane tension; ', ...
        'T_{pre} reaches %.1f N/m'], max(Tpre_eff_curve)));

    xlim([eps_plot(1) eps_plot(end)]);
    yyaxis left; ylim([0 1.05]);
    yyaxis right; ylim([0 1.05 * max(Ktrans_norm, [], 'omitnan')]);

    legend('Location', 'southwest');
    grid on; box on; formatAxes(gca);

    fprintf('\n--- combined transverse-stiffness diagnostic ---\n');
    fprintf('Ktrans = q/wmax from the unit-pressure center-deflection compliance.\n');
    fprintf('eps       Kbend/Kbend0    Tpre(N/m)    Ktrans/Ktrans0\n');
    reportStrains = [0, p.experimentalEffectiveStrains];
    for i = 1:numel(reportStrains)
        fprintf('%6.3f       %9.4f     %10.4f       %10.4f\n', ...
            reportStrains(i), ...
            interp1(eps_plot, Dplate_norm, reportStrains(i), 'linear'), ...
            interp1(eps_plot, Tpre_eff_curve, reportStrains(i), 'linear'), ...
            interp1(eps_plot, Ktrans_norm, reportStrains(i), 'linear'));
    end

    saveCurrentFigure(resultsDir, ['03_normalized_mechanics_', versionTag]);
end

function [kT0_rmse, kT0_wrmse] = runFocusedKTSweep(eps_plot, eps_data, y_data, y_err, p, kT0_sweep)
    kT0_rmse = zeros(numel(kT0_sweep), 1);
    kT0_wrmse = zeros(numel(kT0_sweep), 1);
    for j = 1:numel(kT0_sweep)
        pSweep = p;
        pSweep.preTensionMode = 'ideal';
        pSweep.kT0 = kT0_sweep(j);
        [S_sweep, ~] = runBareReferencedModelResponse(eps_plot, pSweep);
        [kT0_rmse(j), kT0_wrmse(j)] = computeModelErrors(eps_plot, S_sweep, eps_data, y_data, y_err);
    end
end

function plotKT0IntermediateDiagnostics(eps_plot, p, kT0_sweep, resultsDir, versionTag)
    dV_mat = NaN(numel(kT0_sweep), numel(eps_plot));
    w_mat = NaN(size(dV_mat));
    strain_avg_mat = NaN(size(dV_mat));
    strain_max_mat = NaN(size(dV_mat));
    conv_mat = false(size(dV_mat));
    vol_frac_mat = NaN(size(dV_mat));
    Psurf_diag = p.P_static + p.dP0/2;

    for ik = 1:numel(kT0_sweep)
        pSweep = p;
        pSweep.kT0 = kT0_sweep(ik);
        pSweep.preTensionMode = 'ideal';
        for ie = 1:numel(eps_plot)
            [Dplate, ~, Tpre_eff, ~, a_load, ~, ~, ~] = evaluateMembraneState(eps_plot(ie), pSweep);
            solverState.Pi = p.P_gas0;
            [~, ~, dV, wmax, ~, converged, ~] = solveCavityPressure(Psurf_diag, a_load, Dplate, Tpre_eff, pSweep, solverState);
            dV_mat(ik, ie) = dV;
            w_mat(ik, ie) = wmax;
            conv_mat(ik, ie) = converged;
            vol_frac_mat(ik, ie) = dV / pSweep.V00;
            strain_avg_mat(ik, ie) = (64/105) * (wmax / a_load)^2;
            strain_max_mat(ik, ie) = (32/27) * (wmax / a_load)^2;
        end
    end

    plotRawDebugQuantity(eps_plot, kT0_sweep, 1e9*dV_mat, conv_mat, '\Delta V (mm^3)', '\Delta V vs. Pre-Strain for k_{T0} Sweep', resultsDir, ['21_kT0_deltaV_diagnostic_', versionTag]);
    plotRawDebugQuantity(eps_plot, kT0_sweep, 1e3*w_mat, conv_mat, 'w_{max} (mm)', 'Center Deflection vs. Pre-Strain for k_{T0} Sweep', resultsDir, ['22_kT0_wmax_diagnostic_', versionTag]);
    plotRawDebugQuantity(eps_plot, kT0_sweep, vol_frac_mat, conv_mat, '\Delta V / V_0 (-)', 'Fractional Cavity Volume Change Diagnostic', resultsDir, ['24_kT0_volume_fraction_diagnostic_', versionTag]);

    idx0 = find(abs(kT0_sweep) < 1e-12, 1);
    if ~isempty(idx0)
        yAvg = strain_avg_mat(idx0,:);
        yMax = strain_max_mat(idx0,:);
        bad = ~conv_mat(idx0,:) | ~isfinite(yAvg) | ~isfinite(yMax) | abs(yAvg) > 1 | abs(yMax) > 1;
        yAvg(bad) = NaN;
        yMax(bad) = NaN;
        figure; hold on;
        plot(eps_plot, yAvg, '-', 'LineWidth', 2.5, 'DisplayName', 'Avg. deformation-induced strain');
        plot(eps_plot, yMax, '--', 'LineWidth', 2.5, 'DisplayName', 'Max local slope strain');
        plot(eps_plot, eps_plot, ':k', 'LineWidth', 2.2, 'DisplayName', 'Installed pre-strain');
        xlabel('Effective Engineering Strain, \epsilon_{eff} (-)');
        ylabel('Strain Estimate (-)');
        title('Pressure-Induced Strain Estimate for k_{T0}=0');
        legend('Location','northwest');
        grid on; box on; formatAxes(gca);
        set(gca,'YScale','log');
        ylim([1e-10 1]);
        saveCurrentFigure(resultsDir, ['23_kT0_zero_deformation_strain_', versionTag]);
    end
end

function plotKT0RawImplementationDiagnostics(eps_plot, p, kT0_debug_sweep, resultsDir, versionTag)
    [diagData, eps_report] = computeKT0Diagnostics(eps_plot, p, kT0_debug_sweep);

    plotRawDebugQuantity(eps_plot, kT0_debug_sweep, 1e3*diagData.wmax, diagData.conv, 'w_{max} (mm)', 'Raw Center Deflection, No Normalization', resultsDir, ['25_raw_kT0_wmax_', versionTag]);
    plotRawDebugQuantity(eps_plot, kT0_debug_sweep, diagData.dVfrac, diagData.conv, '\Delta V / V_0 (-)', 'Raw Fractional Volume Change, No Normalization', resultsDir, ['26_raw_kT0_dVfrac_', versionTag]);
    plotRawDebugQuantity(eps_plot, kT0_debug_sweep, diagData.Kbend, diagData.conv, 'K_{bend} (N m)', 'Raw Flexural Rigidity, No Normalization', resultsDir, ['27_raw_kT0_Kbend_', versionTag]);
    plotRawDebugQuantity(eps_plot, kT0_debug_sweep, diagData.Tpre, diagData.conv, 'T_{pre,eff} (N/m)', 'Raw Effective Pre-Tension, No Normalization', resultsDir, ['28_raw_kT0_Tpre_', versionTag]);
    plotRawDebugQuantity(eps_plot, kT0_debug_sweep, diagData.alpha, diagData.conv, '\alpha (-)', 'Raw Dimensionless Tension Parameter, No Normalization', resultsDir, ['29_raw_kT0_alpha_', versionTag]);
    plotRawDebugQuantity(eps_plot, kT0_debug_sweep, diagData.iter, diagData.conv, 'Solver Iterations (-)', 'Raw Solver Iteration Count', resultsDir, ['30_raw_kT0_iterations_', versionTag]);
    plotRawDebugQuantity(eps_plot, kT0_debug_sweep, diagData.dVraw_frac, diagData.conv, '\Delta V_{raw} / V_0 (-)', 'Raw Uncapped Fractional Volume Change', resultsDir, ['32_raw_kT0_dVraw_frac_', versionTag]);
    plotRawDebugQuantity(eps_plot, kT0_debug_sweep, diagData.Vi_frac, diagData.conv, 'V_i / V_0 (-)', 'Final Cavity Volume Fraction', resultsDir, ['33_raw_kT0_Vi_frac_', versionTag]);
    plotRawDebugQuantity(eps_plot, kT0_debug_sweep, diagData.q, diagData.conv, 'q = P_{surf} - P_i (Pa)', 'Final Membrane Pressure Load', resultsDir, ['34_raw_kT0_q_membrane_', versionTag]);

    figure; hold on;
    for ik = 1:numel(kT0_debug_sweep)
        plot(eps_plot, double(diagData.conv(ik, :)), '-', 'LineWidth', 2.2, 'DisplayName', sprintf('k_{T0}=%.2f', kT0_debug_sweep(ik)));
    end
    xlabel('Effective Engineering Strain, \epsilon_{eff} (-)');
    ylabel('Converged Flag (-)');
    title('Raw Solver Convergence Flag');
    ylim([-0.1 1.1]);
    legend('Location', 'best');
    grid on; box on; formatAxes(gca);
    saveCurrentFigure(resultsDir, ['31_raw_kT0_convergence_', versionTag]);

    printKT0DiagnosticsTable(eps_plot, kT0_debug_sweep, diagData, eps_report);
end

function [diagData, eps_report] = computeKT0Diagnostics(eps_plot, p, kT0_debug_sweep)
    nK = numel(kT0_debug_sweep);
    nE = numel(eps_plot);
    diagData.wmax       = NaN(nK, nE);
    diagData.dVfrac     = NaN(nK, nE);
    diagData.Kbend      = NaN(nK, nE);
    diagData.Tpre       = NaN(nK, nE);
    diagData.alpha      = NaN(nK, nE);
    diagData.iter       = NaN(nK, nE);
    diagData.conv       = false(nK, nE);
    diagData.dVraw_frac = NaN(nK, nE);
    diagData.Vi_frac    = NaN(nK, nE);
    diagData.q          = NaN(nK, nE);
    Psurf_diag = p.P_static + p.dP0/2;

    for ik = 1:nK
        pSweep = p;
        pSweep.kT0 = kT0_debug_sweep(ik);
        pSweep.preTensionMode = 'ideal';

        for ie = 1:nE
            eps_i = eps_plot(ie);
            [Dplate, ~, Tpre_eff, ~, a_load, ~, ~, ~] = evaluateMembraneState(eps_i, pSweep);
            solverState.Pi = pSweep.P_gas0;
            [~, Vi, dV, wmax, iter, converged, ~, dV_raw, q_last] = solveCavityPressure(Psurf_diag, a_load, Dplate, Tpre_eff, pSweep, solverState);

            diagData.Kbend(ik, ie)      = Dplate;
            diagData.Tpre(ik, ie)       = Tpre_eff;
            diagData.wmax(ik, ie)       = wmax;
            diagData.dVfrac(ik, ie)     = dV / pSweep.V00;
            diagData.dVraw_frac(ik, ie) = dV_raw / pSweep.V00;
            diagData.Vi_frac(ik, ie)    = Vi / pSweep.V00;
            diagData.q(ik, ie)          = q_last;
            diagData.iter(ik, ie)       = iter;
            diagData.conv(ik, ie)       = converged;

            if Tpre_eff > 0
                diagData.alpha(ik, ie) = a_load * sqrt(Tpre_eff / max(Dplate, 1e-30));
            else
                diagData.alpha(ik, ie) = 0;
            end
        end
    end

    eps_report = [0.000, 0.020, 0.111, 0.176, 0.250, 0.333];
end

function printKT0DiagnosticsTable(eps_plot, kT0_debug_sweep, diagData, eps_report)
    fprintf('\n--- raw kT0 implementation diagnostics ---\n');
    fprintf('kT0      eps      conv   iter      wmax(mm)      dV/V0        dVraw/V0     Vi/V0        q(Pa)       Kbend(Nm)      Tpre(N/m)      alpha\n');

    for ik = 1:numel(kT0_debug_sweep)
        for eps_i = eps_report
            w_i     = interp1(eps_plot, diagData.wmax(ik,:), eps_i, 'linear');
            vf_i    = interp1(eps_plot, diagData.dVfrac(ik,:), eps_i, 'linear');
            vfr_i   = interp1(eps_plot, diagData.dVraw_frac(ik,:), eps_i, 'linear');
            Vi_i    = interp1(eps_plot, diagData.Vi_frac(ik,:), eps_i, 'linear');
            q_i     = interp1(eps_plot, diagData.q(ik,:), eps_i, 'linear');
            K_i     = interp1(eps_plot, diagData.Kbend(ik,:), eps_i, 'linear');
            T_i     = interp1(eps_plot, diagData.Tpre(ik,:), eps_i, 'linear');
            a_i     = interp1(eps_plot, diagData.alpha(ik,:), eps_i, 'linear');
            iter_i  = interp1(eps_plot, diagData.iter(ik,:), eps_i, 'nearest');
            conv_i  = interp1(eps_plot, double(diagData.conv(ik,:)), eps_i, 'nearest');
            fprintf('%5.2f   %6.3f     %1.0f   %5.0f   %11.4e   %11.4e   %11.4e   %10.6f   %10.4e   %11.4e   %11.4e   %11.4e\n', ...
                kT0_debug_sweep(ik), eps_i, conv_i, iter_i, 1e3*w_i, vf_i, vfr_i, Vi_i, q_i, K_i, T_i, a_i);
        end
    end
end

function plotThicknessUpdateDiagnostics(eps_plot, p, resultsDir, versionTag)

    modes = {'constant', 'incompressible', 'poisson'};
    t_mat = NaN(numel(modes), numel(eps_plot));
    D_mat = NaN(size(t_mat));
    T_mat = NaN(size(t_mat));
    S_mat = NaN(size(t_mat));

    for im = 1:numel(modes)
        pSweep = p;
        pSweep.thicknessMode = modes{im};

        [T_cav, D_curve, ~, Tpre_curve] = runModelOverStrainRange(eps_plot, pSweep);

        Tref = interp1(eps_plot, T_cav, 0.0, 'linear');
        S_mat(im, :) = T_cav ./ Tref;
        D_mat(im, :) = D_curve;
        T_mat(im, :) = Tpre_curve;

        for ie = 1:numel(eps_plot)
            [~, ~, ~, ~, ~, t_plate] = evaluatePlateGeometry(eps_plot(ie), pSweep);
            t_mat(im, ie) = t_plate;
        end
    end

    figure; hold on;
    for im = 1:numel(modes)
        plot(eps_plot, 1e3*t_mat(im,:), '-', 'LineWidth', 2.2, ...
            'DisplayName', modes{im});
    end
    xlabel('Effective Engineering Strain, \epsilon_{eff} (-)');
    ylabel('Membrane Thickness, t (mm)');
    title('Thickness Update Diagnostic');
    legend('Location', 'best');
    grid on; box on; formatAxes(gca);
    saveCurrentFigure(resultsDir, ['50_thickness_update_', versionTag]);

    figure; hold on;
    for im = 1:numel(modes)
        plot(eps_plot, D_mat(im,:) ./ D_mat(im,1), '-', 'LineWidth', 2.2, ...
            'DisplayName', modes{im});
    end
    xlabel('Effective Engineering Strain, \epsilon_{eff} (-)');
    ylabel('D(\epsilon) / D(0) (-)');
    title('Bending Stiffness Sensitivity to Thickness Law');
    legend('Location', 'best');
    grid on; box on; formatAxes(gca);
    saveCurrentFigure(resultsDir, ['51_thickness_Dnorm_', versionTag]);

    figure; hold on;
    for im = 1:numel(modes)
        plot(eps_plot, S_mat(im,:), '-', 'LineWidth', 2.2, ...
            'DisplayName', modes{im});
    end
    xlabel('Effective Engineering Strain, \epsilon_{eff} (-)');
    ylabel('Normalized Transmission (-)');
    title('Model Sensitivity to Thickness Law');
    legend('Location', 'best');
    grid on; box on; formatAxes(gca);
    saveCurrentFigure(resultsDir, ['52_thickness_transmission_', versionTag]);

    fprintf('\n--- thickness update diagnostics ---\n');
    fprintf('mode             eps       t(mm)    Kbend/Kbend0    Tpre(N/m)      S/S0\n');

    for im = 1:numel(modes)
        for eps_i = [0.000 0.111 0.176 0.250 0.333]
            t_i = interp1(eps_plot, t_mat(im,:), eps_i, 'linear');
            D_i = interp1(eps_plot, D_mat(im,:), eps_i, 'linear');
            T_i = interp1(eps_plot, T_mat(im,:), eps_i, 'linear');
            S_i = interp1(eps_plot, S_mat(im,:), eps_i, 'linear');

            fprintf('%-15s  %6.3f   %10.4e   %10.4e   %11.4e   %8.4f\n', ...
                modes{im}, eps_i, 1e3*t_i, D_i/D_mat(im,1), T_i, S_i);
        end
    end
end

function plotKT0ThicknessReferenceFitSweep(eps_plot, eps_data, y_data, y_err, p, ...
    kT0_grid, c75, c80, c85, c90, markerSize, errorLineWidth, ...
    bareColor, resultsDir, versionTag)

    refModes = {'current', 'reference'};
    rmse_mat = NaN(numel(refModes), numel(kT0_grid));
    wrmse_mat = NaN(size(rmse_mat));

    bestS = cell(numel(refModes), 1);
    bestk = NaN(numel(refModes), 1);

    for im = 1:numel(refModes)
        for ik = 1:numel(kT0_grid)
            pSweep = p;
            pSweep.preTensionThicknessMode = refModes{im};
            pSweep.kT0 = kT0_grid(ik);

            [S_sweep, ~] = runBareReferencedModelResponse(eps_plot, pSweep);
            [rmse_mat(im, ik), wrmse_mat(im, ik)] = ...
                computeModelErrors(eps_plot, S_sweep, eps_data, y_data, y_err);
        end

        [~, bestIdx] = min(wrmse_mat(im, :));
        bestk(im) = kT0_grid(bestIdx);

        pBest = p;
        pBest.preTensionThicknessMode = refModes{im};
        pBest.kT0 = bestk(im);
        [bestS{im}, ~] = runBareReferencedModelResponse(eps_plot, pBest);
    end

    figure; hold on;
    for im = 1:numel(refModes)
        plot(kT0_grid, rmse_mat(im, :), '-o', 'LineWidth', 2.2, ...
            'MarkerSize', 7, 'DisplayName', refModes{im});
    end
    xlabel('k_{T0} (-)');
    ylabel('RMSE (-)');
    title('k_{T0} Fit Sweep: Thickness Reference');
    legend('Location', 'best');
    grid on; box on; formatAxes(gca);
    saveCurrentFigure(resultsDir, ['56_kT0_thickness_reference_RMSE_', versionTag]);

    figure; hold on;
    for im = 1:numel(refModes)
        plot(kT0_grid, wrmse_mat(im, :), '-o', 'LineWidth', 2.2, ...
            'MarkerSize', 7, 'DisplayName', refModes{im});
    end
    xlabel('k_{T0} (-)');
    ylabel('WRMSE (-)');
    title('Weighted k_{T0} Fit Sweep: Thickness Reference');
    legend('Location', 'best');
    grid on; box on; formatAxes(gca);
    saveCurrentFigure(resultsDir, ['57_kT0_thickness_reference_WRMSE_', versionTag]);

    figure; hold on;
    hBareMean = plotBareReference(bareColor);
    [h75, h80, h85, h90] = plotExperimentalData(eps_data, y_data, y_err, ...
        c75, c80, c85, c90, markerSize, errorLineWidth);

    hModel = gobjects(numel(refModes), 1);
    for im = 1:numel(refModes)
        hModel(im) = plot(eps_plot, bestS{im}, '-', 'LineWidth', 2.5, ...
            'DisplayName', sprintf('%s, best k_{T0}=%.4g', refModes{im}, bestk(im)));
    end

    xlabel('Effective Engineering Strain, \epsilon_{eff} (-)');
    ylabel('Normalized Sensitivity, S/S_0 (-)');
    title('Best k_{T0} Fit by Pre-Tension Thickness Reference');
    xlim([eps_plot(1) eps_plot(end)]);
    ylim([0.90 1.45]);
    grid on; box on; formatAxes(gca);
    legend([h75, h80, h85, h90, hBareMean, hModel'], ...
        {'75%', '80%', '85%', '90%', 'Bare-Port Reference', ...
         hModel(1).DisplayName, hModel(2).DisplayName}, ...
        'Location', 'northwest');
    saveCurrentFigure(resultsDir, ['58_kT0_thickness_reference_best_fit_', versionTag]);

    fprintf('\n--- kT0 thickness-reference fit sweep ---\n');
    fprintf('mode          best_kT0      min_RMSE      min_WRMSE\n');

    for im = 1:numel(refModes)
        [minWRMSE, bestIdx] = min(wrmse_mat(im, :));
        minRMSE = rmse_mat(im, bestIdx);
        fprintf('%-10s    %8.4f     %9.4f     %9.4f\n', ...
            refModes{im}, kT0_grid(bestIdx), minRMSE, minWRMSE);
    end
end

function plotSolverAndVolumeCapComparison(eps_plot, p, kT0_debug_sweep, resultsDir, versionTag)
    cases = {
        'root capped',       'root',       true;
        'root uncapped',     'root',       false;
        'fixedPoint capped', 'fixedPoint', true;
        'fixedPoint uncapped','fixedPoint',false
    };

    for ic = 1:size(cases, 1)
        pCase = p;
        pCase.solverMode = cases{ic, 2};
        pCase.useVolumeCap = cases{ic, 3};
        [diagData, ~] = computeKT0Diagnostics(eps_plot, pCase, kT0_debug_sweep);
        safeName = regexprep(cases{ic, 1}, '[^A-Za-z0-9]+', '_');

        plotRawDebugQuantity(eps_plot, kT0_debug_sweep, 1e3*diagData.wmax, diagData.conv, ...
            'w_{max} (mm)', ['Solver/Cap Comparison: ', cases{ic, 1}, ', w_{max}'], ...
            resultsDir, ['40_compare_', safeName, '_wmax_', versionTag]);

        plotRawDebugQuantity(eps_plot, kT0_debug_sweep, diagData.dVfrac, diagData.conv, ...
            '\Delta V / V_0 (-)', ['Solver/Cap Comparison: ', cases{ic, 1}, ', \Delta V/V_0'], ...
            resultsDir, ['41_compare_', safeName, '_dVfrac_', versionTag]);

        plotRawDebugQuantity(eps_plot, kT0_debug_sweep, diagData.iter, diagData.conv, ...
            'Solver Iterations (-)', ['Solver/Cap Comparison: ', cases{ic, 1}, ', iterations'], ...
            resultsDir, ['42_compare_', safeName, '_iterations_', versionTag]);
    end
end

function plotRawDebugQuantity(eps_plot, kT0_debug_sweep, Y, conv_mat, yLabelText, titleText, resultsDir, outName)
    figure; hold on;

    for ik = 1:numel(kT0_debug_sweep)
        y = Y(ik, :);
        bad = ~conv_mat(ik, :) | ~isfinite(y);
        yPlot = y;
        yPlot(bad) = NaN;

        plot(eps_plot, yPlot, '-', 'LineWidth', 2.2, 'DisplayName', sprintf('k_{T0}=%.4g', kT0_debug_sweep(ik)));

        if any(bad & isfinite(y))
            plot(eps_plot(bad), y(bad), 'x', 'LineWidth', 1.8, 'MarkerSize', 7, 'HandleVisibility', 'off');
        end
    end

    xlabel('Effective Engineering Strain, \epsilon_{eff} (-)');
    ylabel(yLabelText);
    title(titleText);
    legend('Location', 'best');
    grid on; box on; formatAxes(gca);
    saveCurrentFigure(resultsDir, outName);
end

function printSensitivitySweepSummary( ...
    kT0_sweep, kT0_rmse, kT0_wrmse)

    fprintf('\n--- local tensile-restraint sensitivity summary ---\n');
    fprintf('Bare-port normalization uses the fixed assumption R_bare.\n');
    fprintf('R_bare is not fitted to the experimental response.\n\n');

    fprintf('Ideal local tensile-restraint sweep:\n');
    for j = 1:numel(kT0_sweep)
        fprintf(['  kT0 = %.4f  | RMSE = %.4f', ...
                 ' | WRMSE = %.4f\n'], ...
            kT0_sweep(j), kT0_rmse(j), kT0_wrmse(j));
    end
end

function saveCurrentFigure(resultsDir, baseName)
    fig = gcf;
    formatManuscriptFigure(fig);
    % MATLAB can export the floating axes-toolbar background even when its
    % icons are not visible, masking vertical strips of plotted data. Hide
    % every axes toolbar explicitly before writing any output format.
    axList = findall(fig, 'Type', 'Axes');
    for ax = axList'
        try
            if ~isempty(ax.Toolbar)
                ax.Toolbar.Visible = 'off';
            end
        catch
            % Some specialized axes do not expose a Toolbar property.
        end
    end
    drawnow;
    exportgraphics(fig, fullfile(resultsDir, [baseName, '.png']), 'Resolution', 300);
    exportgraphics(fig, fullfile(resultsDir, [baseName, '.pdf']), 'ContentType', 'vector', 'BackgroundColor', 'white');
    savefig(fig, fullfile(resultsDir, [baseName, '.fig']));
    print(fig, fullfile(resultsDir, [baseName, '.eps']), '-depsc', '-vector');
end

function formatAxes(ax)
    if nargin < 1 || isempty(ax)
        ax = gca;
    end
    grid(ax, 'on');
    box(ax, 'on');
end

function formatManuscriptFigure(fig)
    if nargin < 1 || isempty(fig)
        fig = gcf;
    end

    set(fig, 'Color', 'white');
    axList = findall(fig, 'Type', 'Axes');

    isTransmissionAudit = strcmp(fig.Tag, 'transmissionAudit');
    isBucklingWeightedDiagnostic = strcmp(fig.Tag, 'bucklingWeightedDiagnostic');

    for ax = axList'
        ax.FontSize = 16;          % Tick labels (leave unchanged)
        ax.LineWidth = 1.4;
        ax.TickDir = 'in';
        ax.Box = 'on';
        ax.XGrid = 'on';
        ax.YGrid = 'on';

        % Slightly larger labels/titles
        ax.XLabel.FontSize = 20;   % was 18
        ax.YLabel.FontSize = 20;   % was 18
        ax.Title.FontSize  = 20;   % was 18
        ax.Title.FontWeight = 'bold';

        lines = findall(ax, 'Type', 'Line');
        for h = lines'
            if h.LineWidth < 2.2
                h.LineWidth = 2.2;
            end
        end

        errs = findall(ax, 'Type', 'ErrorBar');
        for h = errs'
            h.LineWidth = 2.0;
            h.MarkerSize = 9;
            h.CapSize = 8;
        end

        leg = ax.Legend;
        if ~isempty(leg)
            leg.FontSize = 14;
            leg.Box = 'on';
        end

        if isTransmissionAudit
            ax.FontSize = 11;
            ax.XLabel.FontSize = 13;
            ax.YLabel.FontSize = 13;
            ax.Title.FontSize = 13;
            if ~isempty(leg)
                leg.FontSize = 9;
            end
        end

        if isBucklingWeightedDiagnostic
            ax.FontSize = 11;
            ax.XLabel.FontSize = 14;
            ax.YLabel.FontSize = 14;
            ax.Title.FontSize = 15;
            if ~isempty(leg)
                leg.FontSize = 10;
            end
        end
    end
end

function plotPreTensionThicknessDiagnostics(eps_plot, p, resultsDir, versionTag)

    modes = {'current', 'reference'};
    S_mat = NaN(numel(modes), numel(eps_plot));
    T_mat = NaN(size(S_mat));
    alpha_mat = NaN(size(S_mat));

    for im = 1:numel(modes)
        pSweep = p;
        pSweep.preTensionThicknessMode = modes{im};

        [T_cav, D_curve, ~, Tpre_curve, ~, a_load_curve] = runModelOverStrainRange(eps_plot, pSweep);

        Tref = interp1(eps_plot, T_cav, 0.0, 'linear');
        S_mat(im, :) = T_cav ./ Tref;
        T_mat(im, :) = Tpre_curve;
        alpha_mat(im, :) = a_load_curve .* sqrt(max(Tpre_curve, 0) ./ max(D_curve, 1e-30));
    end

    figure; hold on;
    for im = 1:numel(modes)
        plot(eps_plot, T_mat(im,:), '-', 'LineWidth', 2.2, ...
            'DisplayName', modes{im});
    end
    xlabel('Effective Engineering Strain, \epsilon_{eff} (-)');
    ylabel('T_{pre,eff} (N/m)');
    title('Pre-Tension Sensitivity to Thickness Reference');
    legend('Location', 'best');
    grid on; box on; formatAxes(gca);
    saveCurrentFigure(resultsDir, ['53_Tpre_thickness_reference_', versionTag]);

    figure; hold on;
    for im = 1:numel(modes)
        plot(eps_plot, alpha_mat(im,:), '-', 'LineWidth', 2.2, ...
            'DisplayName', modes{im});
    end
    xlabel('Effective Engineering Strain, \epsilon_{eff} (-)');
    ylabel('\alpha = a(T/K_{bend})^{1/2} (-)');
    title('Dimensionless Tension Sensitivity to Thickness Reference');
    legend('Location', 'best');
    grid on; box on; formatAxes(gca);
    saveCurrentFigure(resultsDir, ['54_alpha_thickness_reference_', versionTag]);

    figure; hold on;
    for im = 1:numel(modes)
        plot(eps_plot, S_mat(im,:), '-', 'LineWidth', 2.2, ...
            'DisplayName', modes{im});
    end
    xlabel('Effective Engineering Strain, \epsilon_{eff} (-)');
    ylabel('Normalized Transmission (-)');
    title('Model Sensitivity to Pre-Tension Thickness Reference');
    legend('Location', 'best');
    grid on; box on; formatAxes(gca);
    saveCurrentFigure(resultsDir, ['55_transmission_Tpre_thickness_reference_', versionTag]);

    fprintf('\n--- pre-tension thickness-reference diagnostics ---\n');
    fprintf('mode          eps       Tpre(N/m)      alpha        S/S0\n');

    for im = 1:numel(modes)
        for eps_i = [0.000 0.111 0.176 0.250 0.333]
            T_i = interp1(eps_plot, T_mat(im,:), eps_i, 'linear');
            a_i = interp1(eps_plot, alpha_mat(im,:), eps_i, 'linear');
            S_i = interp1(eps_plot, S_mat(im,:), eps_i, 'linear');

            fprintf('%-10s  %6.3f   %11.4e   %10.4e   %8.4f\n', ...
                modes{im}, eps_i, T_i, a_i, S_i);
        end
    end
end

function plotStaticDepthDiagnostics(eps_plot, p, resultsDir, versionTag)

    rho_water = p.rho_water;
    g_water   = p.g_water;
    depth_sweep_m = p.depth_sweep_m;
    Pstatic_sweep = p.P_gas0 + rho_water * g_water * depth_sweep_m;

    nD = numel(depth_sweep_m);
    nE = numel(eps_plot);

    gain_mat = NaN(nD, nE);
    gain_rel_mat = NaN(nD, nE);
    Pi_mat = NaN(nD, nE);
    q_mat = NaN(nD, nE);
    wmax_mat = NaN(nD, nE);
    dVfrac_mat = NaN(nD, nE);
    Uplate_mat = NaN(nD, nE);
    Wgas_mat = NaN(nD, nE);
    conv_mat = false(nD, nE);

    for id = 1:nD
        pDepth = p;
        pDepth.P_static = Pstatic_sweep(id);
        pDepth.baselineDepth_m = depth_sweep_m(id);

        for ie = 1:nE
            eps_i = eps_plot(ie);

            [Dplate, ~, Tpre_eff, ~, a_load, ~, ~, ~] = evaluateMembraneState(eps_i, pDepth);

            solverState.Pi = pDepth.P_gas0;

            % First solve the hydrostatically preloaded equilibrium using the
            % original sealed-gas inventory P_gas0*V00.
            [Pi_i, Vi_i, dV_i, w_i, ~, conv_eq, staticState, ~, q_i] = ...
                solveCavityPressure(pDepth.P_static, a_load, Dplate, Tpre_eff, pDepth, solverState);

            % Then evaluate the tangent small-signal gain about that state.
            [gain_i, conv_i] = evaluateLocalP3GainFromState( ...
                pDepth.P_static, a_load, Dplate, Tpre_eff, pDepth, staticState);

            gain_mat(id, ie) = gain_i;
            Pi_mat(id, ie) = Pi_i;
            q_mat(id, ie) = q_i;
            wmax_mat(id, ie) = w_i;
            dVfrac_mat(id, ie) = dV_i / pDepth.V00;
            Uplate_mat(id, ie) = 0.5 * q_i * dV_i;
            Wgas_mat(id, ie) = pDepth.P_gas0 * pDepth.V00 * ...
                log(pDepth.V00 / Vi_i);
            conv_mat(id, ie) = conv_i && conv_eq;
        end

        refGain = interp1(eps_plot, gain_mat(id,:), 0.0, 'linear');
        gain_rel_mat(id,:) = gain_mat(id,:) ./ refGain;
    end

    figure; hold on;
    for id = 1:nD
        y = gain_rel_mat(id,:);
        y(~conv_mat(id,:)) = NaN;
        plot(eps_plot, y, '-', 'LineWidth', 2.2, ...
            'DisplayName', sprintf('depth = %.3g m', depth_sweep_m(id)));
    end
    xlabel('Effective Engineering Strain, \epsilon_{eff} (-)');
    ylabel('Normalized Small-Signal Gain (-)');
    title('Small-Signal Transmission vs. Static Depth');
    if p.includeKnoDeploymentDepths
        subtitle('12-15 m KNO cases are mechanistic extrapolations beyond small-deflection validity');
    end
    legend('Location', 'best');
    grid on; box on; formatAxes(gca);
    saveCurrentFigure(resultsDir, ['60_static_depth_gain_', versionTag]);

    figure; hold on;
    for id = 1:nD
        y = 1e3 * wmax_mat(id,:);
        y(~conv_mat(id,:)) = NaN;
        plot(eps_plot, y, '-', 'LineWidth', 2.2, ...
            'DisplayName', sprintf('depth = %.3g m', depth_sweep_m(id)));
    end
    xlabel('Effective Engineering Strain, \epsilon_{eff} (-)');
    ylabel('Static w_{max} (mm)');
    title('Hydrostatic Preload Deflection vs. Depth');
    if p.includeKnoDeploymentDepths
        subtitle('12-15 m KNO cases are mechanistic extrapolations beyond small-deflection validity');
    end
    legend('Location', 'best');
    grid on; box on; formatAxes(gca);
    saveCurrentFigure(resultsDir, ['61_static_depth_wmax_', versionTag]);

    figure; hold on;
    for id = 1:nD
        y = dVfrac_mat(id,:);
        y(~conv_mat(id,:)) = NaN;
        plot(eps_plot, y, '-', 'LineWidth', 2.2, ...
            'DisplayName', sprintf('depth = %.3g m', depth_sweep_m(id)));
    end
    xlabel('Effective Engineering Strain, \epsilon_{eff} (-)');
    ylabel('\Delta V_{static} / V_0 (-)');
    title('Hydrostatic Cavity Compression vs. Depth');
    if p.includeKnoDeploymentDepths
        subtitle('12-15 m KNO cases are mechanistic extrapolations beyond small-deflection validity');
    end
    legend('Location', 'best');
    grid on; box on; formatAxes(gca);
    saveCurrentFigure(resultsDir, ['62_static_depth_dVfrac_', versionTag]);

    figure; hold on;
    for id = 1:nD
        y = 1e3 * (Uplate_mat(id,:) + Wgas_mat(id,:));
        y(~conv_mat(id,:)) = NaN;
        plot(eps_plot, y, '-', 'LineWidth', 2.2, ...
            'DisplayName', sprintf('depth = %.3g m', depth_sweep_m(id)));
    end
    xlabel('Effective Engineering Strain, \epsilon_{eff} (-)');
    ylabel('Stored-energy diagnostics (mJ)');
    title('Static Plate Energy + Isothermal Gas Compression Work');
    if p.includeKnoDeploymentDepths
        subtitle('12-15 m KNO cases are mechanistic extrapolations beyond small-deflection validity');
    end
    legend('Location', 'best');
    grid on; box on; formatAxes(gca);
    saveCurrentFigure(resultsDir, ['63_static_depth_energy_', versionTag]);

    fprintf('\n--- two-stage hydrostatic preload diagnostics ---\n');
    fprintf('Gas inventory remains P_gas0*V00 at every depth.\n');
    fprintf('depth(m)   eps   conv   gain_rel   wstatic(mm)   dVstatic/V0   Pi-Pstatic(Pa)   qstatic(Pa)   Uplate(mJ)   Wgas(mJ)\n');

    eps_report = [0.000, 0.111, 0.176, 0.250, 0.333];

    for id = 1:nD
        for eps_i = eps_report
            g_i    = interp1(eps_plot, gain_rel_mat(id,:), eps_i, 'linear');
            w_i    = interp1(eps_plot, wmax_mat(id,:), eps_i, 'linear');
            vf_i   = interp1(eps_plot, dVfrac_mat(id,:), eps_i, 'linear');
            Pi_i   = interp1(eps_plot, Pi_mat(id,:), eps_i, 'linear');
            q_i    = interp1(eps_plot, q_mat(id,:), eps_i, 'linear');
            Up_i   = interp1(eps_plot, Uplate_mat(id,:), eps_i, 'linear');
            Wg_i   = interp1(eps_plot, Wgas_mat(id,:), eps_i, 'linear');
            conv_i = interp1(eps_plot, double(conv_mat(id,:)), eps_i, 'nearest');

            fprintf('%8.4f  %6.3f    %1.0f    %9.4f   %11.4e   %11.4e   %14.4e   %11.4e   %10.4e   %10.4e\n', ...
                depth_sweep_m(id), eps_i, conv_i, g_i, 1e3*w_i, vf_i, ...
                Pi_i - Pstatic_sweep(id), q_i, 1e3*Up_i, 1e3*Wg_i);
        end
    end

    % Export a compact deployment-focused table only when KNO cases are on.
    % Convergence indicates numerical equilibrium only; it does not override
    % the small-deflection limits represented by w/a and w/h.
    if p.includeKnoDeploymentDepths
    deploymentDepths = p.deploymentDepthRange_m;
    nRows = numel(deploymentDepths) * numel(eps_report);
    depth_m = NaN(nRows,1); strain = NaN(nRows,1);
    P_external_abs_Pa = NaN(nRows,1); P_cavity_abs_Pa = NaN(nRows,1);
    raw_gain = NaN(nRows,1); gain_rel_zero_strain = NaN(nRows,1);
    wstatic_mm = NaN(nRows,1); dVstatic_over_V0 = NaN(nRows,1);
    w_over_a = NaN(nRows,1); w_over_h = NaN(nRows,1);
    within_small_deflection = false(nRows,1); converged = false(nRows,1);
    row = 0;
    for depth_i = deploymentDepths
        id = find(abs(depth_sweep_m - depth_i) < 1e-12, 1);
        for eps_i = eps_report
            row = row + 1;
            [a_i, ~, ~, ~, ~, h_i] = evaluatePlateGeometry(eps_i, p);
            w_i = interp1(eps_plot, wmax_mat(id,:), eps_i, 'linear');
            depth_m(row) = depth_i;
            strain(row) = eps_i;
            P_external_abs_Pa(row) = Pstatic_sweep(id);
            P_cavity_abs_Pa(row) = interp1(eps_plot, Pi_mat(id,:), eps_i, 'linear');
            raw_gain(row) = interp1(eps_plot, gain_mat(id,:), eps_i, 'linear');
            gain_rel_zero_strain(row) = interp1(eps_plot, gain_rel_mat(id,:), eps_i, 'linear');
            wstatic_mm(row) = 1e3 * w_i;
            dVstatic_over_V0(row) = interp1(eps_plot, dVfrac_mat(id,:), eps_i, 'linear');
            w_over_a(row) = abs(w_i) / a_i;
            w_over_h(row) = abs(w_i) / h_i;
            converged(row) = interp1(eps_plot, double(conv_mat(id,:)), eps_i, 'nearest') > 0.5;
            within_small_deflection(row) = w_over_a(row) <= 0.1 && w_over_h(row) <= 1.0;
        end
    end
    knoTable = table(depth_m, strain, P_external_abs_Pa, P_cavity_abs_Pa, ...
        raw_gain, gain_rel_zero_strain, wstatic_mm, dVstatic_over_V0, ...
        w_over_a, w_over_h, within_small_deflection, converged);
    csvPath = fullfile(resultsDir, ['kno_depth_states_', versionTag, '.csv']);
    writetable(knoTable, csvPath);
    fprintf('Saved KNO deployment-depth table: %s\n', csvPath);
    end
end

function [Glocal, convergedBoth] = evaluateLocalP3GainFromState(P3_center, a_plate, Dplate, Tpre_eff, p, staticState)
    dP = p.dP_local;
    P3_plus  = P3_center + dP / 2;
    P3_minus = P3_center - dP / 2;

    [P1_plus, ~, ~, ~, ~, convPlus] = ...
        solveCavityPressure(P3_plus, a_plate, Dplate, Tpre_eff, p, staticState);
    [P1_minus, ~, ~, ~, ~, convMinus] = ...
        solveCavityPressure(P3_minus, a_plate, Dplate, Tpre_eff, p, staticState);

    Glocal = (P1_plus - P1_minus) / dP;
    convergedBoth = convPlus && convMinus && isfinite(Glocal);
end

function plotSpatialPreTensionDiagnostic(eps_plot, eps_data, y_data, y_err, p, ...
    c75, c80, c85, c90, markerSize, errorLineWidth, bareColor, resultsDir, versionTag)

    exponent_sweep = [0 2 4 6 8 10 12];
    k_global = 0.75;

    figure; hold on;
    hBareMean = plotBareReference(bareColor);
    [h75, h80, h85, h90] = plotExperimentalData(eps_data, y_data, y_err, ...
        c75, c80, c85, c90, markerSize, errorLineWidth);

    fprintf('\n--- spatial pre-tension coupling diagnostic ---\n');
    fprintf('kT0_global = %.3f\n', k_global);
    fprintf('n_exp      RMSE        WRMSE       S(0.111)   S(0.176)   S(0.250)   S(0.333)\n');

    hModel = gobjects(numel(exponent_sweep), 1);

    for i = 1:numel(exponent_sweep)
        pSweep = p;
        pSweep.preTensionCouplingMode = 'radiusPower';
        pSweep.kT0_global = k_global;
        pSweep.radiusPowerExponent = exponent_sweep(i);

        [S_sweep, ~] = runBareReferencedModelResponse(eps_plot, pSweep);
        [rmse_i, wrmse_i] = computeModelErrors(eps_plot, S_sweep, eps_data, y_data, y_err);

        S_pts = interp1(eps_plot, S_sweep, eps_data, 'linear');

        fprintf('%5.1f   %9.4f   %9.4f   %8.4f   %8.4f   %8.4f   %8.4f\n', ...
            exponent_sweep(i), rmse_i, wrmse_i, S_pts(1), S_pts(2), S_pts(3), S_pts(4));

        hModel(i) = plot(eps_plot, S_sweep, '-', 'LineWidth', 2.2, ...
            'DisplayName', sprintf('k_g=%.2f, n=%.0f', k_global, exponent_sweep(i)));
    end

    xlabel('Effective Engineering Strain, \epsilon_{eff} (-)');
    ylabel('Normalized Sensitivity, S/S_0 (-)');
    title('Spatial Coupling Test for Pre-Tension Formulation');
    xlim([eps_plot(1) eps_plot(end)]);
    ylim([0.90 1.45]);
    grid on; box on; formatAxes(gca);

    legend([h75, h80, h85, h90, hBareMean, hModel'], ...
        [{'75%', '80%', '85%', '90%', 'Bare-Port Reference'}, ...
        arrayfun(@(n) sprintf('k_g=%.2f, n=%.0f', k_global, n), exponent_sweep, 'UniformOutput', false)], ...
        'Location', 'northwest');

    saveCurrentFigure(resultsDir, ['70_spatial_pre_tension_coupling_', versionTag]);
end

function plotLocalTensionScaleDiagnostic(eps_plot, eps_data, y_data, y_err, p, ...
    c75, c80, c85, c90, markerSize, errorLineWidth, ...
    bareColor, resultsDir, versionTag)

    kT0_local_sweep = [0, 0.0075, 0.02, 0.05, 0.10, 0.25, 0.50, 0.75, 1.00];

    figure; hold on;
    hBareMean = plotBareReference(bareColor);
    [h75, h80, h85, h90] = plotExperimentalData(eps_data, y_data, y_err, ...
        c75, c80, c85, c90, markerSize, errorLineWidth);

    fprintf('\n--- local tension scale diagnostic ---\n');
    fprintf('This diagnostic separates global installed tension from the local tension entering the plate equation.\n');
    fprintf('kT0_local   RMSE      WRMSE     S(0.111)   S(0.176)   S(0.250)   S(0.333)\n');

    hModel = gobjects(numel(kT0_local_sweep), 1);

    for ik = 1:numel(kT0_local_sweep)
        pSweep = p;

        % Keep this diagnostic simple: local coupling is exactly kT0.
        pSweep.preTensionCouplingMode = 'direct';
        pSweep.kT0 = kT0_local_sweep(ik);

        [S_sweep, ~] = runBareReferencedModelResponse(eps_plot, pSweep);
        [rmse_i, wrmse_i] = computeModelErrors(eps_plot, S_sweep, eps_data, y_data, y_err);
        S_pts = interp1(eps_plot, S_sweep, eps_data, 'linear');

        fprintf('%9.4f   %7.4f   %8.4f   %8.4f   %8.4f   %8.4f   %8.4f\n', ...
            kT0_local_sweep(ik), rmse_i, wrmse_i, ...
            S_pts(1), S_pts(2), S_pts(3), S_pts(4));

        hModel(ik) = plot(eps_plot, S_sweep, '-', 'LineWidth', 2.2, ...
            'DisplayName', sprintf('k_{T,local}=%.4g', kT0_local_sweep(ik)));
    end

    xlabel('Effective Engineering Strain, \epsilon_{eff} (-)');
    ylabel('Normalized Sensitivity, S/S_0 (-)');
    title('Local Tension Scale Diagnostic');
    xlim([eps_plot(1) eps_plot(end)]);
    ylim([0.85 1.45]);
    grid on; box on; formatAxes(gca);

    legend([h75, h80, h85, h90, hBareMean, hModel'], ...
        [{'75%', '80%', '85%', '90%', 'Bare-Port Reference'}, ...
        arrayfun(@(k) sprintf('k_{T,local}=%.4g', k), kT0_local_sweep, 'UniformOutput', false)], ...
        'Location', 'northwest');

    saveCurrentFigure(resultsDir, ['71_local_tension_scale_', versionTag]);

    fprintf('\n--- local/global tension table ---\n');
    fprintf('kT0_local   eps      Tglobal(N/m)   Tlocal(N/m)    alpha_local    S/S0\n');

    eps_report = [0.000, 0.111, 0.176, 0.250, 0.333];

    for ik = 1:numel(kT0_local_sweep)
        pSweep = p;
        pSweep.preTensionCouplingMode = 'direct';
        pSweep.kT0 = kT0_local_sweep(ik);

        [S_sweep, ~] = runBareReferencedModelResponse(eps_plot, pSweep);

        for eps_i = eps_report
            [Dplate, E, Tlocal, ~, a_load, ~, ~, eps_pre] = evaluateMembraneState(eps_i, pSweep);

            [~, ~, ~, ~, ~, t_plate] = evaluatePlateGeometry(eps_i, pSweep);

            switch lower(pSweep.preTensionThicknessMode)
                case 'current'
                    t_for_Tpre = t_plate;
                case 'reference'
                    t_for_Tpre = pSweep.t_plate0;
                otherwise
                    error('Unknown preTensionThicknessMode: %s', pSweep.preTensionThicknessMode);
            end

            Tglobal = E * t_for_Tpre * eps_pre / (1 - pSweep.nu_plate);

            if Tlocal > 0
                alpha_local = a_load * sqrt(Tlocal / max(Dplate, 1e-30));
            else
                alpha_local = 0;
            end

            S_i = interp1(eps_plot, S_sweep, eps_i, 'linear');

            fprintf('%9.4f   %6.3f   %12.4e   %12.4e   %12.4e   %8.4f\n', ...
                kT0_local_sweep(ik), eps_i, Tglobal, Tlocal, alpha_local, S_i);
        end
    end
end

function plotTensionCapDiagnostic(eps_plot, eps_data, y_data, y_err, p, ...
    c75, c80, c85, c90, markerSize, errorLineWidth, ...
    bareColor, resultsDir, versionTag)

    Tcap_sweep = [0, 0.25, 0.50, 0.75, 1.00, 2.00, 5.00, 10.00];

    figure; hold on;
    hBareMean = plotBareReference(bareColor);
    [h75, h80, h85, h90] = plotExperimentalData(eps_data, y_data, y_err, ...
        c75, c80, c85, c90, markerSize, errorLineWidth);

    fprintf('\n--- tension cap diagnostic ---\n');
    fprintf('Tcap(N/m)   RMSE      WRMSE     S(0.111)   S(0.176)   S(0.250)   S(0.333)\n');

    hModel = gobjects(numel(Tcap_sweep), 1);

    for i = 1:numel(Tcap_sweep)
        pSweep = p;
        pSweep.preTensionCouplingMode = 'tensionCap';
        pSweep.Tpre_cap_Npm = Tcap_sweep(i);

        [S_sweep, ~] = runBareReferencedModelResponse(eps_plot, pSweep);
        [rmse_i, wrmse_i] = computeModelErrors(eps_plot, S_sweep, eps_data, y_data, y_err);
        S_pts = interp1(eps_plot, S_sweep, eps_data, 'linear');

        fprintf('%8.2f   %7.4f   %8.4f   %8.4f   %8.4f   %8.4f   %8.4f\n', ...
            Tcap_sweep(i), rmse_i, wrmse_i, ...
            S_pts(1), S_pts(2), S_pts(3), S_pts(4));

        hModel(i) = plot(eps_plot, S_sweep, '-', 'LineWidth', 2.2, ...
            'DisplayName', sprintf('T_{cap}=%.2g N/m', Tcap_sweep(i)));
    end

    xlabel('Effective Engineering Strain, \epsilon_{eff} (-)');
    ylabel('Normalized Sensitivity, S/S_0 (-)');
    title('Local Tension Cap Diagnostic');
    xlim([eps_plot(1) eps_plot(end)]);
    ylim([0.85 1.45]);
    grid on; box on; formatAxes(gca);

    legend([h75, h80, h85, h90, hBareMean, hModel'], ...
        [{'75%', '80%', '85%', '90%', 'Bare-Port Reference'}, ...
        arrayfun(@(T) sprintf('T_{cap}=%.2g N/m', T), Tcap_sweep, 'UniformOutput', false)], ...
        'Location', 'northwest');

    saveCurrentFigure(resultsDir, ['72_tension_cap_diagnostic_', versionTag]);
end

function plotModelComparisonDiagnostic(eps_plot, eps_data, y_data, y_err, p, ...
    c75, c80, c85, c90, markerSize, errorLineWidth, ...
    bareColor, resultsDir, versionTag)

    cases = {
        'Zero local tension',       'none',  'direct',     0.0075, 0.75,  6,  Inf;
        'Direct local tension',     'ideal', 'direct',     0.0075, 0.75,  6,  Inf;
        'Capped local tension',     'ideal', 'tensionCap', 0.0075, 0.75,  6,  0.25;
        'Spatial attenuation',      'ideal', 'radiusPower',0.0075, 0.75, 12,  Inf
    };

    figure; hold on;
    hBareMean = plotBareReference(bareColor);
    [h75, h80, h85, h90] = plotExperimentalData(eps_data, y_data, y_err, ...
        c75, c80, c85, c90, markerSize, errorLineWidth);

    fprintf('\n--- model comparison diagnostic ---\n');
    fprintf('case                                      RMSE      WRMSE     S(0.111)   S(0.176)   S(0.250)   S(0.333)\n');

    hModel = gobjects(size(cases,1), 1);

    for i = 1:size(cases,1)
        pSweep = p;
        pSweep.preTensionMode = cases{i,2};
        pSweep.preTensionCouplingMode = cases{i,3};
        pSweep.kT0 = cases{i,4};
        pSweep.kT0_global = cases{i,5};
        pSweep.radiusPowerExponent = cases{i,6};
        pSweep.Tpre_cap_Npm = cases{i,7};

        [S_sweep, ~] = runBareReferencedModelResponse(eps_plot, pSweep);
        [rmse_i, wrmse_i] = computeModelErrors(eps_plot, S_sweep, eps_data, y_data, y_err);
        S_pts = interp1(eps_plot, S_sweep, eps_data, 'linear');

        fprintf('%-40s  %7.4f   %8.4f   %8.4f   %8.4f   %8.4f   %8.4f\n', ...
            cases{i,1}, rmse_i, wrmse_i, S_pts(1), S_pts(2), S_pts(3), S_pts(4));

        hModel(i) = plot(eps_plot, S_sweep, '-', 'LineWidth', 2.5, ...
            'DisplayName', cases{i,1});
    end

    xlabel('Effective Engineering Strain, \epsilon_{eff} (-)');
    ylabel('Normalized Sensitivity, S/S_0 (-)');
    title('Local Tension Formulation Comparison');
    xlim([eps_plot(1) eps_plot(end)]);
    ylim([0.85 1.7]);
    grid on; box on; formatAxes(gca);

    legend([h75, h80, h85, h90, hModel', hBareMean], ...
        {'75%', ...
         '80%', ...
         '85%', ...
         '90%', ...
         'Zero local tension', ...
         'Direct local tension', ...
         'Capped local tension', ...
         'Spatial attenuation', ...
         'Bare-Port Reference'}, ...
        'Location', 'northwest');

    saveCurrentFigure(resultsDir, ['73_model_comparison_local_pretension_', versionTag]);
end

function plotModelClosureDiagnostic(eps_plot, eps_data, y_data, y_err, ...
    p, bareColor, resultsDir, versionTag)

    % Final physics-envelope audit. Every tensioned case retains kT0=1 and
    % the measured 5.5-mm unsupported radius. Parameter changes are explicit
    % bounds, not fitted values or candidate calibrated baselines.
    cases = {
        'Linear, isothermal baseline',   'linearEquibiaxial', 1.0, 1.00, 1.00;
        'Neo-Hookean, isothermal',       'neoHookean',        1.0, 1.00, 1.00;
        'Linear, adiabatic upper bound', 'linearEquibiaxial', 1.4, 1.00, 1.00;
        'Neo-Hookean, adiabatic bound',  'neoHookean',        1.4, 1.00, 1.00;
        'Linear, 50% modulus bound',     'linearEquibiaxial', 1.0, 0.50, 1.00;
        'Linear, 75% thickness bound',   'linearEquibiaxial', 1.0, 1.00, 0.75;
        'Combined optimistic envelope',  'neoHookean',        1.4, 0.50, 0.75
    };

    figure('Position', [100 100 1120 650]); hold on;
    errorbar(eps_data, y_data, y_err, 'ks', 'MarkerFaceColor', 'k', ...
        'LineWidth', 1.5, 'CapSize', 7, 'DisplayName', 'Experiment / bare port');
    yline(1.0, '--', 'Color', bareColor, 'LineWidth', 1.8, ...
        'DisplayName', sprintf('Assumed bare port, R_{bare}=%.2f', p.R_bare));

    fprintf('\n--- final model-closure physics envelope ---\n');
    fprintf(['All cases use a=5.5 mm and kT0=1. Parameter bounds are not fitted.\n', ...
        'case                                  n_gas  E/E0  h/h0  minRaw  ', ...
        'Rbare,max(all enhanced)  R/Rbare at eps [0.111 0.176 0.250 0.333]\n']);

    for i = 1:size(cases,1)
        pCase = p;
        pCase.preTensionMode = 'ideal';
        pCase.preTensionCouplingMode = 'direct';
        pCase.kT0 = 1.0;
        pCase.loadedRadiusMode = 'fixed';
        pCase.preTensionLaw = cases{i,2};
        pCase.gasExponent = cases{i,3};
        pCase.E_plate0 = p.E_plate0 * cases{i,4};
        pCase.t_plate0 = p.t_plate0 * cases{i,5};

        [rawCurve, ~, ~, ~, ~, ~, ~, ~, ~, convCurve] = ...
            runModelOverStrainRange(eps_plot, pCase);
        rawPts = interp1(eps_plot, rawCurve, eps_data, 'linear');
        relPts = rawPts ./ p.R_bare;
        minRaw = min(rawPts);
        allConverged = all(convCurve > 0.5);
        passive = all(rawCurve(convCurve > 0.5) >= 0 & rawCurve(convCurve > 0.5) <= 1);

        fprintf('%-37s %5.2f  %4.2f  %4.2f  %7.4f       %7.4f             %7.4f %7.4f %7.4f %7.4f  conv=%d passive=%d\n', ...
            cases{i,1}, cases{i,3}, cases{i,4}, cases{i,5}, minRaw, ...
            minRaw, relPts, allConverged, passive);

        plot(eps_plot, rawCurve ./ p.R_bare, '-', 'LineWidth', 2.2, ...
            'DisplayName', cases{i,1});
    end

    xlabel('Effective Engineering Strain, \epsilon_{eff} (-)');
    ylabel('R_{cav}/R_{bare} (-)');
    title('Final Non-Fitted Physics-Envelope Audit');
    xlim([eps_plot(1) eps_plot(end)]); grid on; box on; formatAxes(gca);
    legend('Location', 'best');
    saveCurrentFigure(resultsDir, ['75_model_closure_physics_envelope_', versionTag]);

    % Test whether finite differential amplitude changes the locked-tension
    % baseline. The sweep spans the original 5 Pa case, the 400 Pa nominal
    % total differential, and 800 Pa total (+400/-400 Pa at the ports).
    amplitudeSweep = [5, 10, 25, 50, 100, 200, 400, 800];
    fprintf('\n--- finite differential-pressure amplitude audit ---\n');
    fprintf('Full ideal linear tension, a=5.5 mm, isothermal gas.\n');
    fprintf('DeltaP(Pa)  R/Rbare at eps [0.111 0.176 0.250 0.333]\n');
    for dP = amplitudeSweep
        pAmp = p;
        pAmp.dP0 = dP;
        pAmp.preTensionMode = 'ideal';
        pAmp.preTensionCouplingMode = 'direct';
        pAmp.kT0 = 1.0;
        pAmp.loadedRadiusMode = 'fixed';
        pAmp.preTensionLaw = 'linearEquibiaxial';
        pAmp.gasExponent = 1.0;
        [relCurve, ~] = runBareReferencedModelResponse(eps_plot, pAmp);
        relPts = interp1(eps_plot, relCurve, eps_data, 'linear');
        fprintf('%10.1f      %7.4f %7.4f %7.4f %7.4f\n', dP, relPts);
    end
end
