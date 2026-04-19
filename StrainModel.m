%% sensor_interface_numerical_model_v18.m
clear; clc; close all;

% v18 
% v18 is an iterative reduced-order force-balance model for the
% cavity/interface system, developed to assess whether strain-dependent
% effective pre-tension, cavity compression, and modest geometry evolution
% are sufficient to explain the observed normalized sensitivity response.

% v18 assumptions:
% 1) axisymmetric reduced-order response
% 2) effective circular loaded region
% 3) small-deflection-style volume estimate regularized by saturation
% 4) pre-strain represented through a smoothed effective in-plane tension term
% 5) cavity response solved iteratively via ideal-gas compression
% 6) closure terms approximate higher-order mechanics not explicitly included
%    in the reduced-order formulation

%% -----------------------------
% Results folder
% ------------------------------
resultsDir = fullfile(pwd, 'results_v18');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

%% -----------------------------
% Plot / output toggles
% ------------------------------
makeMainFigure        = true;
makeTransmissionPlot  = true;
makeMechanicsPlot     = true;
makeGeometryPlot      = true;
makeJointRMSEPlot     = true;
makeJointWRMSEPlot    = true;
printSweepSummary     = true;

%% -----------------------------
% Experimental data / interface states
% ------------------------------
d0 = 21.978; % mm, nominal funnel diameter at zero engineering strain

scaleDiameters_raw = [16.484, 17.582, 18.681, 19.780, 21.978]; % mm
engStrains_raw     = (d0 - scaleDiameters_raw) ./ d0;          % (-)

[engStrains, sortIdx] = sort(engStrains_raw);
scaleDiameters = scaleDiameters_raw(sortIdx);

eps_data = [0.100, 0.150, 0.200, 0.250];
y_data   = [1.167, 1.328, 1.270, 1.224];
y_err    = [0.093, 0.082, 0.045, 0.042];

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
axisFontSize    = 16;
legendFontSize  = 12;

%% -----------------------------
% Model parameters
% ------------------------------
p = struct();

% Reference / perturbation pressures
p.P_ref = 101325;   % Pa
p.dP0   = 5.0;      % Pa

% Measured geometry states
p.d0_mm             = d0;
p.scaleDiameters_mm = scaleDiameters;
p.engStrains        = engStrains;

%% -----------------------------
% Cavity geometry (per sensing side)
% ------------------------------
h1 = 2.506;   % mm
h2 = 6.617;   % mm
r  = 5.5;     % mm

p.r_forced_m = r * 1e-3;
p.A_forced   = pi * p.r_forced_m^2;

p.V11 = pi * r^2 * h1;
p.V22 = pi * h2^2 * r - pi * h2^3 / 3;
V_total_mm3 = p.V11 + p.V22;
p.V00 = V_total_mm3 * 1e-9; % m^3

%% -----------------------------
% Plate / membrane model (v18)
% Reduced-order model with:
% 1) smoothed pretension engagement
% 2) smooth tanh saturation of dV
% 3) baseline bare-port scaling applied after T/T(0)
% 4) weak, capped strain-dependent loaded-radius correction
% ------------------------------
p.t_plate0 = 0.020 * 0.0254;  % m
p.nu_plate = 0.49;
p.E_plate0 = 6.0e5;           % Pa

p.kT0               = 0.250;
p.eps_char          = 0.160;

% Baseline multiplicative calibration that maps interface-relative
% transmission to the experimentally referenced bare-port sensitivity scale
p.S0_interface_bare = 1.60;

% Weak geometry-based loaded-radius correction
p.alpha_load         = 0.100;
p.a_load_cap_fraction = 0.10;

% Optional strain-dependent modulus
p.useStrainDependentE = false;
p.c1 = 0.0;
p.c2 = 0.0;

%% -----------------------------
% Iterative solver settings
% ------------------------------
p.tolP    = 1e-5;
p.maxIter = 1000;
p.relax   = 0.08;

%% -----------------------------
% Stability settings
% ------------------------------
p.dV_cap_fraction = 0.10;

%% -----------------------------
% Parameter sweeps
% Focused around current working region
% ------------------------------
eps_plot = linspace(0, 0.36, 300);

kT0_sweep       = [0.250, 0.275, 0.300];
S0_sweep        = [1.55, 1.60, 1.65, 1.70];
alpha_load_sweep = [0.050, 0.075, 0.100, 0.125];

kT0_joint_grid = 0.250:0.025:0.325;
S0_joint_grid   = 1.55:0.025:1.70;

%% -----------------------------
% Print configuration summary
% ------------------------------
fprintf('Nominal cavity volume V00: %.3e m^3\n', p.V00);
fprintf('Plate thickness: %.3e m\n', p.t_plate0);
fprintf('Baseline modulus: %.3e Pa\n', p.E_plate0);
fprintf('Poisson ratio: %.3f\n', p.nu_plate);
fprintf('Forced loading radius: %.3e m\n', p.r_forced_m);
fprintf('Pressure tolerance: %.3e Pa\n', p.tolP);
fprintf('Max iterations: %d\n', p.maxIter);
fprintf('Pressure relaxation factor: %.3f\n', p.relax);
fprintf('Pre-tension scale factor kT0: %.3f\n', p.kT0);
fprintf('Pretension smoothing strain eps_char: %.3f\n', p.eps_char);
fprintf('Baseline interface sensitivity S0_interface_bare: %.3f\n', p.S0_interface_bare);
fprintf('Loaded-radius weighting alpha_load: %.3f\n', p.alpha_load);
fprintf('Loaded-radius cap fraction: %.3f\n', p.a_load_cap_fraction);
fprintf('Volume-change saturation fraction: %.3f\n', p.dV_cap_fraction);
fprintf('Use strain-dependent E: %d\n\n', p.useStrainDependentE);

fprintf('Sampling interface diameters (mm): ');
fprintf('%.3f ', scaleDiameters);
fprintf('\n');

fprintf('Engineering strains (-):          ');
fprintf('%.3f ', engStrains);
fprintf('\n\n');

%% -----------------------------
% Evaluate baseline model over strain range
% ------------------------------
[T_cav, Dplate_curve, Eeff_curve, Tpre_eff_curve, phi_pre_curve, ...
 a_load_curve, a_geom_curve, D_installed_curve_mm, eps_pre_curve, ...
 converged_curve] = runModelOverStrainRange_v18(eps_plot, p);

Tref = interp1(eps_plot, T_cav, 0.00, 'linear');

if abs(Tref) < 1e-12
    error('Reference transmission for normalization is too small. Inspect raw T_cav first.');
end

T_cav_rel0   = T_cav ./ Tref;
S_bare_model = p.S0_interface_bare * T_cav_rel0;

%% -----------------------------
% Print summary at experimental points
% ------------------------------
fprintf('\n--- v18 model summary at experimental strain values ---\n');
fprintf('   eps      data      model    T_rel0      a_load(mm)   a_geom(mm)   D(Nm)   Tpre_eff(N/m)   phi_pre   resid    w_resid\n');

resid = zeros(size(eps_data));
wres  = zeros(size(eps_data));

for i = 1:numel(eps_data)
    Si_model = interp1(eps_plot, S_bare_model,       eps_data(i), 'linear');
    Ti_rel0  = interp1(eps_plot, T_cav_rel0,         eps_data(i), 'linear');
    ai_load  = interp1(eps_plot, 1e3 * a_load_curve, eps_data(i), 'linear');
    ai_geom  = interp1(eps_plot, 1e3 * a_geom_curve, eps_data(i), 'linear');
    Di       = interp1(eps_plot, Dplate_curve,       eps_data(i), 'linear');
    Tpre_effi = interp1(eps_plot, Tpre_eff_curve, eps_data(i), 'linear');
    phii     = interp1(eps_plot, phi_pre_curve,      eps_data(i), 'linear');

    resid(i) = Si_model - y_data(i);
    wres(i)  = resid(i) / y_err(i);

    fprintf('%7.3f   %7.3f   %7.3f   %7.3f   %10.3f   %10.3f   %8.3e   %8.3e   %7.3f   %7.3f   %7.3f\n', ...
        eps_data(i), y_data(i), Si_model, Ti_rel0, ai_load, ai_geom, Di, Tpre_effi, phii, resid(i), wres(i));
end

rmse  = sqrt(mean(resid .^ 2));
wrmse = sqrt(mean(wres  .^ 2));

fprintf('\nv18 RMSE  = %.4f\n', rmse);
fprintf('v18 WRMSE = %.4f\n', wrmse);
fprintf('v18 converged points: %d / %d\n', sum(converged_curve > 0.5), numel(converged_curve));

%% -----------------------------
% Main figure
% ------------------------------
if makeMainFigure
    figure; hold on;

    hBareMean = yline(1.0, '-', 'Color', bareColor, 'LineWidth', 1.8, 'DisplayName', 'Bare Port Reference');
    yline(1.04, '--', 'Color', bareColor, 'LineWidth', 1.2, 'HandleVisibility', 'off');
    yline(0.96, '--', 'Color', bareColor, 'LineWidth', 1.2, 'HandleVisibility', 'off');

    h75 = errorbar(eps_data(4), y_data(4), y_err(4), 's', 'Color', c75, 'MarkerFaceColor', c75, 'MarkerEdgeColor', c75, 'MarkerSize', markerSize, 'LineWidth', errorLineWidth, 'CapSize', 8, 'DisplayName', '75%');
    h80 = errorbar(eps_data(3), y_data(3), y_err(3), 's', 'Color', c80, 'MarkerFaceColor', c80, 'MarkerEdgeColor', c80, 'MarkerSize', markerSize, 'LineWidth', errorLineWidth, 'CapSize', 8, 'DisplayName', '80%');
    h85 = errorbar(eps_data(2), y_data(2), y_err(2), 's', 'Color', c85, 'MarkerFaceColor', c85, 'MarkerEdgeColor', c85, 'MarkerSize', markerSize, 'LineWidth', errorLineWidth, 'CapSize', 8, 'DisplayName', '85%');
    h90 = errorbar(eps_data(1), y_data(1), y_err(1), 's', 'Color', c90, 'MarkerFaceColor', c90, 'MarkerEdgeColor', c90, 'MarkerSize', markerSize, 'LineWidth', errorLineWidth, 'CapSize', 8, 'DisplayName', '90%');

    hModel = plot(eps_plot, S_bare_model, '-', 'Color', modelColor, 'LineWidth', modelLineWidth, ...
        'DisplayName', 'Iterative v18 Interface Transmission Model');

    xlabel('Engineering Strain');
    ylabel('Normalized Sensitivity (-)');
    xlim([0 0.36]);
    ylim([0.85 1.60]);
    grid on;
    box on;
    formatAxes(gca, 16);

    legend([h75, h80, h85, h90, hBareMean, hModel], ...
           {'75%', '80%', '85%', '90%', 'Bare Port Reference', 'Iterative v18 Interface Transmission Model'}, ...
           'Location', 'northeast', 'FontSize', legendFontSize);

    saveCurrentFigure(resultsDir, '01_main_normalized_sensitivity_v18');
end

%% -----------------------------
% Interface-relative transmission
% ------------------------------
if makeTransmissionPlot
    figure;
    plot(eps_plot, T_cav_rel0, '-', 'Color', physColor, 'LineWidth', 2.5); hold on;
    yline(1.0, ':k');

    xlabel('Engineering Strain');
    ylabel('Normalized Transmission (-)');
    title('v18 Cavity / Interface Transmission Relative to Zero-Strain Interface');
    legend('T_{cav,rel0}(\epsilon)', 'Reference = 1', 'Location', 'best');
    grid on;
    formatAxes(gca, 16);

    saveCurrentFigure(resultsDir, '02_transmission_rel0_v18');
end

%% -----------------------------
% Mechanics diagnostic
% ------------------------------
if makeMechanicsPlot
    figure;
    yyaxis left
    plot(eps_plot, Eeff_curve, 'k-', 'LineWidth', 2); hold on;
    ylabel('E_{eff}(\epsilon) [Pa]');

    yyaxis right
    plot(eps_plot, Dplate_curve, 'c--', 'LineWidth', 2); hold on;
    plot(eps_plot, Tpre_eff_curve, 'r-.', 'LineWidth', 2);
    plot(eps_plot, phi_pre_curve, 'm:', 'LineWidth', 2.2);
    ylabel('D_{plate}(\epsilon) [N m], T_{pre,eff}(\epsilon) [N/m], \phi_{pre}(\epsilon) [-]');
    
    xlabel('Engineering Strain');
    title('v18 Mechanics Diagnostics');
    legend('E_{eff}(\epsilon)', 'D_{plate}(\epsilon)', 'T_{pre,eff}(\epsilon)', '\phi_{pre}(\epsilon)', 'Location', 'best');
    grid on;
    formatAxes(gca, 16);

    saveCurrentFigure(resultsDir, '03_mechanics_v18');
end

%% -----------------------------
% Geometry / loaded-radius diagnostic
% ------------------------------
if makeGeometryPlot
    figure; hold on;
    plot(eps_plot, 1e3 * a_geom_curve, '-', 'Color', geomColor, 'LineWidth', 2.5, 'DisplayName', 'a_{geom}(\epsilon)');
    plot(eps_plot, 1e3 * a_load_curve, '--', 'Color', modelColor, 'LineWidth', 2.5, 'DisplayName', 'a_{load}(\epsilon)');
    yline(1e3 * p.r_forced_m, ':', 'Color', bareColor, 'LineWidth', 1.8, 'DisplayName', 'r_{forced}');

    xlabel('Engineering Strain');
    ylabel('Radius [mm]');
    title('v18 Geometry Diagnostics');
    legend('Location', 'best');
    grid on;
    formatAxes(gca, 16);

    saveCurrentFigure(resultsDir, '04_geometry_loaded_radius_v18');
end

%% -----------------------------
% Joint sweep: kT0 and S0_interface_bare
% ------------------------------
joint_rmse  = zeros(numel(S0_joint_grid), numel(kT0_joint_grid));
joint_wrmse = zeros(numel(S0_joint_grid), numel(kT0_joint_grid));

for ig = 1:numel(S0_joint_grid)
    for ik = 1:numel(kT0_joint_grid)
        pSweep = p;
        pSweep.S0_interface_bare = S0_joint_grid(ig);
        pSweep.kT0               = kT0_joint_grid(ik);

        [S_sweep, ~] = runBareReferencedModelResponse_v18(eps_plot, pSweep);
        [joint_rmse(ig, ik), joint_wrmse(ig, ik)] = ...
            computeModelErrors(eps_plot, S_sweep, eps_data, y_data, y_err);
    end
end

[minJointRMSE, idxJointRMSE] = min(joint_rmse(:));
[rowRMSE, colRMSE] = ind2sub(size(joint_rmse), idxJointRMSE);
bestS0_RMSE   = S0_joint_grid(rowRMSE);
bestkT0_RMSE  = kT0_joint_grid(colRMSE);

[minJointWRMSE, idxJointWRMSE] = min(joint_wrmse(:));
[rowWRMSE, colWRMSE] = ind2sub(size(joint_wrmse), idxJointWRMSE);
bestS0_WRMSE  = S0_joint_grid(rowWRMSE);
bestkT0_WRMSE = kT0_joint_grid(colWRMSE);

if makeJointRMSEPlot
    figure;
    contourf(kT0_joint_grid, S0_joint_grid, joint_rmse, 24, 'LineColor', 'none');
    hold on;
    colormap(parula);

    plot(p.kT0, p.S0_interface_bare, 'wo', ...
        'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', 'k', ...
        'MarkerSize', 8, ...
        'DisplayName', 'Baseline');

    plot(bestkT0_RMSE, bestS0_RMSE, 'rp', ...
        'MarkerFaceColor', 'r', ...
        'MarkerEdgeColor', 'k', ...
        'MarkerSize', 12, ...
        'DisplayName', 'Min RMSE');

    xlabel('k_{T0}');
    ylabel('S_{0,interface/bare}');
    title('Joint Calibration Sweep: RMSE');
    cb = colorbar;
    cb.Label.String = 'RMSE (-)';
    grid on;
    box on;
    formatAxes(gca, 16);
    legend('Location', 'best');

    saveCurrentFigure(resultsDir, '05_joint_sweep_rmse_v18');
end

if makeJointWRMSEPlot
    figure;
    contourf(kT0_joint_grid, S0_joint_grid, joint_wrmse, 24, 'LineColor', 'none');
    hold on;
    colormap(parula);

    plot(p.kT0, p.S0_interface_bare, 'wo', ...
        'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', 'k', ...
        'MarkerSize', 8, ...
        'DisplayName', 'Baseline');

    plot(bestkT0_WRMSE, bestS0_WRMSE, 'rp', ...
        'MarkerFaceColor', 'r', ...
        'MarkerEdgeColor', 'k', ...
        'MarkerSize', 12, ...
        'DisplayName', 'Min WRMSE');

    xlabel('k_{T0}');
    ylabel('S_{0,interface/bare}');
    title('Joint Calibration Sweep: WRMSE');
    cb = colorbar;
    cb.Label.String = 'WRMSE (-)';
    grid on;
    box on;
    formatAxes(gca, 16);
    legend('Location', 'best');

    saveCurrentFigure(resultsDir, '06_joint_sweep_wrmse_v18');
end

%% -----------------------------
% Focused 1D sweeps for console output
% ------------------------------
kT0_rmse  = zeros(numel(kT0_sweep), 1);
kT0_wrmse = zeros(numel(kT0_sweep), 1);

for j = 1:numel(kT0_sweep)
    pSweep = p;
    pSweep.kT0 = kT0_sweep(j);

    [S_sweep, ~] = runBareReferencedModelResponse_v18(eps_plot, pSweep);
    [kT0_rmse(j), kT0_wrmse(j)] = computeModelErrors(eps_plot, S_sweep, eps_data, y_data, y_err);
end

S0_rmse  = zeros(numel(S0_sweep), 1);
S0_wrmse = zeros(numel(S0_sweep), 1);

for j = 1:numel(S0_sweep)
    pSweep = p;
    pSweep.S0_interface_bare = S0_sweep(j);

    [S_sweep, ~] = runBareReferencedModelResponse_v18(eps_plot, pSweep);
    [S0_rmse(j), S0_wrmse(j)] = computeModelErrors(eps_plot, S_sweep, eps_data, y_data, y_err);
end

alpha_rmse  = zeros(numel(alpha_load_sweep), 1);
alpha_wrmse = zeros(numel(alpha_load_sweep), 1);

for j = 1:numel(alpha_load_sweep)
    pSweep = p;
    pSweep.alpha_load = alpha_load_sweep(j);

    [S_sweep, ~] = runBareReferencedModelResponse_v18(eps_plot, pSweep);
    [alpha_rmse(j), alpha_wrmse(j)] = computeModelErrors(eps_plot, S_sweep, eps_data, y_data, y_err);
end

%% -----------------------------
% Print sensitivity summary
% ------------------------------
if printSweepSummary
    fprintf('\n--- v18 sensitivity sweep summary ---\n');

    fprintf('\nkT0 sweep:\n');
    for j = 1:numel(kT0_sweep)
        fprintf('  kT0 = %.3f  | RMSE = %.4f | WRMSE = %.4f\n', ...
            kT0_sweep(j), kT0_rmse(j), kT0_wrmse(j));
    end

    fprintf('\nBaseline interface sensitivity sweep:\n');
    for j = 1:numel(S0_sweep)
        fprintf('  S0_interface_bare = %.3f  | RMSE = %.4f | WRMSE = %.4f\n', ...
            S0_sweep(j), S0_rmse(j), S0_wrmse(j));
    end

    fprintf('\nLoaded-radius weighting sweep:\n');
    for j = 1:numel(alpha_load_sweep)
        fprintf('  alpha_load = %.3f  | RMSE = %.4f | WRMSE = %.4f\n', ...
            alpha_load_sweep(j), alpha_rmse(j), alpha_wrmse(j));
    end

    fprintf('\nJoint kT0 / S0_interface_bare sweep:\n');
    fprintf('  kT0 range: %.3f to %.3f (%d cases)\n', ...
        kT0_joint_grid(1), kT0_joint_grid(end), numel(kT0_joint_grid));
    fprintf('  S0_interface_bare range: %.3f to %.3f (%d cases)\n', ...
        S0_joint_grid(1), S0_joint_grid(end), numel(S0_joint_grid));

    fprintf('  Minimum joint RMSE  = %.4f at kT0 = %.3f, S0_interface_bare = %.3f\n', ...
        minJointRMSE, bestkT0_RMSE, bestS0_RMSE);
    fprintf('  Minimum joint WRMSE = %.4f at kT0 = %.3f, S0_interface_bare = %.3f\n', ...
        minJointWRMSE, bestkT0_WRMSE, bestS0_WRMSE);
end

%% -----------------------------
% Local functions
% ------------------------------
function [T_cav, Dplate_curve, Eeff_curve, Tpre_eff_curve, phi_pre_curve, ...
          a_load_curve, a_geom_curve, D_installed_curve_mm, eps_pre_curve, ...
          converged_curve] = runModelOverStrainRange_v18(eps_plot, p)

    T_cav                = zeros(size(eps_plot));
    Dplate_curve         = zeros(size(eps_plot));
    Eeff_curve           = zeros(size(eps_plot));
    Tpre_eff_curve       = zeros(size(eps_plot));
    phi_pre_curve        = zeros(size(eps_plot));
    a_load_curve         = zeros(size(eps_plot));
    a_geom_curve         = zeros(size(eps_plot));
    D_installed_curve_mm = zeros(size(eps_plot));
    eps_pre_curve        = zeros(size(eps_plot));
    converged_curve      = zeros(size(eps_plot));

    solverState.Pi = p.P_ref;

    for k = 1:numel(eps_plot)
        eps_query = eps_plot(k);

        [T_cav(k), Dplate_curve(k), Eeff_curve(k), Tpre_eff_curve(k), ...
         phi_pre_curve(k), a_load_curve(k), a_geom_curve(k), ...
         D_installed_curve_mm(k), eps_pre_curve(k), converged_curve(k), ...
         solverState] = evaluateCavityTransmission_v18(eps_query, p, solverState);
    end
end

function [S_bare_model, T_cav_rel0] = runBareReferencedModelResponse_v18(eps_plot, p)
    [T_cav, ~, ~, ~, ~, ~, ~, ~, ~, ~] = runModelOverStrainRange_v18(eps_plot, p);

    Tref = interp1(eps_plot, T_cav, 0.00, 'linear');

    if abs(Tref) < 1e-12
        error('Reference transmission for normalization is too small. Inspect raw T_cav first.');
    end

    T_cav_rel0   = T_cav ./ Tref;
    S_bare_model = p.S0_interface_bare * T_cav_rel0;
end

function [rmse, wrmse] = computeModelErrors(eps_plot, S_model, eps_data, y_data, y_err)
    resid = zeros(size(eps_data));
    wres  = zeros(size(eps_data));

    for i = 1:numel(eps_data)
        Si_model = interp1(eps_plot, S_model, eps_data(i), 'linear');
        resid(i) = Si_model - y_data(i);
        wres(i)  = resid(i) / y_err(i);
    end

    rmse  = sqrt(mean(resid .^ 2));
    wrmse = sqrt(mean(wres  .^ 2));
end

function [Tcav, Dplate, Eeff, Tpre_eff, phi_pre, a_load, a_geom, ...
          D_installed_mm, eps_pre, converged, solverStateOut] = ...
          evaluateCavityTransmission_v18(eps_query, p, solverStateIn)

    [Dplate, Eeff, Tpre_eff, phi_pre, a_load, a_geom, D_installed_mm, eps_pre] = ...
        evaluateMembraneState_v18(eps_query, p);

    symmetricState = solverStateIn;
    if isempty(symmetricState)
        symmetricState.Pi = p.P_ref;
    end

    dPsens_plus  = sensorDifferentialResponse_v18(+p.dP0, a_load, Dplate, Tpre_eff, p, symmetricState);
    dPsens_minus = sensorDifferentialResponse_v18(-p.dP0, a_load, Dplate, Tpre_eff, p, symmetricState);

    Tcav = (dPsens_plus - dPsens_minus) / (2 * p.dP0);

    [~, ~, ~, ~, ~, converged, solverStateOut] = ...
        solveCavityPressure_iterative_v18(p.P_ref + p.dP0, a_load, Dplate, Tpre_eff, p, solverStateIn);
end

function [eps_pre, lambda_pre, D_installed_mm] = evaluateInstalledPrestretchState_v18(eps_query, p)
    D_installed_mm = interp1(p.engStrains, p.scaleDiameters_mm, eps_query, 'linear', 'extrap');
    eps_pre = max(eps_query, 0);
    lambda_pre = 1 + eps_pre;
end

function [a_load, a_geom, D_installed_mm, eps_pre, lambda_pre, t_plate] = ...
    evaluatePlateGeometry_v18(eps_query, p)

    [eps_pre, lambda_pre, D_installed_mm] = evaluateInstalledPrestretchState_v18(eps_query, p);

    a_geom = 0.5 * D_installed_mm * 1e-3;

    da_load = p.alpha_load * (a_geom - p.r_forced_m);
    da_load_max = p.a_load_cap_fraction * p.r_forced_m;
    da_load = max(min(da_load, da_load_max), -da_load_max);

    a_load = p.r_forced_m + da_load;

    t_plate = p.t_plate0 / (lambda_pre ^ 2);
end

function [Dplate, Eeff, Tpre_eff, phi_pre, a_load, a_geom, D_installed_mm, eps_pre] = ...
    evaluateMembraneState_v18(eps_query, p)

    [a_load, a_geom, D_installed_mm, eps_pre, ~, t_plate] = ...
        evaluatePlateGeometry_v18(eps_query, p);

    if p.useStrainDependentE
        Eeff = p.E_plate0 * (1 + p.c1 * eps_pre + p.c2 * eps_pre ^ 2);
    else
        Eeff = p.E_plate0;
    end

    Dplate = Eeff * t_plate ^ 3 / (12 * (1 - p.nu_plate ^ 2));

    phi_pre = 1 - exp(-(eps_pre / max(p.eps_char, 1e-12)) ^ 2);
    Tpre_eff = p.kT0 * phi_pre * Eeff * t_plate * eps_pre / (1 - p.nu_plate);
end

function dPsens = sensorDifferentialResponse_v18(dPsurf, a_plate, Dplate, Tpre_eff, p, solverStateIn)
    P3 = p.P_ref + dPsurf / 2;
    P4 = p.P_ref - dPsurf / 2;

    [P1, ~, ~, ~, ~, ~, ~] = solveCavityPressure_iterative_v18(P3, a_plate, Dplate, Tpre_eff, p, solverStateIn);
    [P2, ~, ~, ~, ~, ~, ~] = solveCavityPressure_iterative_v18(P4, a_plate, Dplate, Tpre_eff, p, solverStateIn);

    dPsens = P1 - P2;
end

function [Pi, Vi, dV, wmax, iter, converged, solverStateOut] = ...
    solveCavityPressure_iterative_v18(Psurf, a_plate, Dplate, Tpre_eff, p, solverStateIn)

    if nargin < 7 || isempty(solverStateIn)
        Pi = p.P_ref;
    else
        Pi = solverStateIn.Pi;
        if ~isfinite(Pi)
            Pi = p.P_ref;
        end
    end

    converged = false;
    Vi = p.V00;
    dV = 0;
    wmax = 0;

    for iter = 1:p.maxIter
        q_membrane = Psurf - Pi;

        [Vi_new, dV_new, wmax_new] = ...
            loadedCavityState_fromMembraneLoad_v18(q_membrane, p.V00, a_plate, Dplate, Tpre_eff, p);

        Pi_new_raw = p.P_ref * p.V00 / Vi_new;
        Pi_new = (1 - p.relax) * Pi + p.relax * Pi_new_raw;

        if abs(Pi_new - Pi) < p.tolP
            Pi = Pi_new;
            Vi = Vi_new;
            dV = dV_new;
            wmax = wmax_new;
            converged = true;
            solverStateOut.Pi = Pi;
            return;
        end

        Pi = Pi_new;
        Vi = Vi_new;
        dV = dV_new;
        wmax = wmax_new;
    end

    solverStateOut.Pi = Pi;
end

function [Vi, dV, wmax] = loadedCavityState_fromMembraneLoad_v18(q_membrane, V00, a, Dplate, Tpre_eff, p)
    denom = 64 * Dplate + 4 * Tpre_eff * a ^ 2;
    denom = max(denom, 1e-18);

    wmax_raw = q_membrane * a ^ 4 / denom;
    dV_raw   = pi * a ^ 2 * wmax_raw / 3;

    dV_cap = p.dV_cap_fraction * V00;
    dV     = dV_cap * tanh(dV_raw / max(dV_cap, 1e-18));

    Vi   = V00 - dV;
    Vi   = max(Vi, 1e-12);
    wmax = 3 * dV / (pi * a ^ 2);
end

function saveCurrentFigure(resultsDir, baseName)
    fig = gcf;
    ax = gca;

    if isprop(ax, 'Toolbar') && ~isempty(ax.Toolbar)
        ax.Toolbar.Visible = 'off';
    end

    set(fig, 'PaperPositionMode', 'auto');

    % Save raster version for quick viewing
    exportgraphics(fig, fullfile(resultsDir, [baseName, '.png']), 'Resolution', 300);

    % Save MATLAB figure for later editing
    savefig(fig, fullfile(resultsDir, [baseName, '.fig']));

    % Save vector version for manuscripts
    print(fig, fullfile(resultsDir, [baseName, '.eps']), '-depsc', '-painters');
end

function formatAxes(ax, fontSize)
    if nargin < 1 || isempty(ax)
        ax = gca;
    end
    if nargin < 2 || isempty(fontSize)
        fontSize = 16;
    end

    set(ax, ...
        'FontSize', fontSize, ...
        'LineWidth', 1.0, ...
        'Box', 'on', ...
        'Layer', 'top', ...
        'TickDir', 'in');

    grid(ax, 'on');
end