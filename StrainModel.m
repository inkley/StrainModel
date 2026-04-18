%% sensor_interface_numerical_model_v12b.m
clear; clc; close all;

%% -----------------------------
% Results folder
% ------------------------------
resultsDir = fullfile(pwd, 'results_v12b');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

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
p.dP0   = 5.0;      % Pa, applied differential perturbation for transmission estimate

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
p.V22 = pi * h2^2 * r - pi * h2^3/3;
V_total_mm3 = p.V11 + p.V22;
p.V00 = V_total_mm3 * 1e-9; % m^3

%% -----------------------------
% Plate / membrane model (v12B)
% ------------------------------
p.t_plate0 = 0.020 * 0.0254;  % m
p.nu_plate = 0.49;
p.E_plate0 = 6.0e5;           % Pa

% Optional scale factor on pre-tension term:
% T0 = kT0 * E*h*eps/(1-nu)
p.kT0 = 1.0;

% Optional strain-dependent modulus toggle
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
p.dV_cap_fraction = 0.10;   % max allowed |dV| / V00

%% -----------------------------
% Optional sensitivity sweep on kT0
% ------------------------------
kT0_sweep = [0.5, 1.0, 1.5];
kT0_grid  = linspace(0.25, 2.0, 36);

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
fprintf('Volume-change cap fraction: %.3f\n', p.dV_cap_fraction);
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
eps_plot = linspace(0, 0.36, 300);

[T_cav, Pi_test_curve, Vi_test_curve, dV_test_curve, Dplate_curve, ...
 Eeff_curve, T0_curve, wmax_test_curve, iter_test_curve, a_plate_curve, ...
 Dplate_geom_curve_mm, converged_curve] = ...
    runModelOverStrainRange_v12b(eps_plot, p);

eps_ref_norm = 0.10;
Tref = interp1(eps_plot, T_cav, eps_ref_norm, 'linear');

if abs(Tref) < 1e-12
    error('Reference transmission for normalization is too small. Inspect raw T_cav first.');
end

T_cav_norm = T_cav ./ Tref;
S_norm = T_cav_norm;

%% -----------------------------
% Print summary at experimental points
% ------------------------------
fprintf('\n--- v12B model summary at experimental strain values ---\n');
fprintf('   eps      data      model    T_cav,norm   a_plate(mm)   D(Nm)      T0(N/m)   resid    w_resid\n');

resid = zeros(size(eps_data));
wres  = zeros(size(eps_data));

for i = 1:numel(eps_data)
    Si_model = interp1(eps_plot, S_norm,            eps_data(i), 'linear');
    Ti_cav   = interp1(eps_plot, T_cav_norm,        eps_data(i), 'linear');
    ai_mm    = interp1(eps_plot, 1e3*a_plate_curve, eps_data(i), 'linear');
    Di       = interp1(eps_plot, Dplate_curve,      eps_data(i), 'linear');
    T0i      = interp1(eps_plot, T0_curve,          eps_data(i), 'linear');

    resid(i) = Si_model - y_data(i);
    wres(i)  = resid(i) / y_err(i);

    fprintf('%7.3f   %7.3f   %7.3f   %7.3f   %10.3f   %8.3e   %8.3e   %7.3f   %7.3f\n', ...
        eps_data(i), y_data(i), Si_model, Ti_cav, ai_mm, Di, T0i, resid(i), wres(i));
end

rmse  = sqrt(mean(resid.^2));
wrmse = sqrt(mean(wres.^2));

fprintf('\nv12B RMSE  = %.4f\n', rmse);
fprintf('v12B WRMSE = %.4f\n', wrmse);
fprintf('v12B converged points: %d / %d\n', sum(converged_curve > 0.5), numel(converged_curve));

%% -----------------------------
% Main figure
% ------------------------------
figure; hold on;

hBareMean = yline(1.0, '-', 'Color', bareColor, 'LineWidth', 1.8, 'DisplayName', 'Bare Port Reference');
yline(1.04, '--', 'Color', bareColor, 'LineWidth', 1.2, 'HandleVisibility', 'off');
yline(0.96, '--', 'Color', bareColor, 'LineWidth', 1.2, 'HandleVisibility', 'off');

h75 = errorbar(eps_data(4), y_data(4), y_err(4), 's', 'Color', c75, 'MarkerFaceColor', c75, 'MarkerEdgeColor', c75, 'MarkerSize', markerSize, 'LineWidth', errorLineWidth, 'CapSize', 8, 'DisplayName', '75%');
h80 = errorbar(eps_data(3), y_data(3), y_err(3), 's', 'Color', c80, 'MarkerFaceColor', c80, 'MarkerEdgeColor', c80, 'MarkerSize', markerSize, 'LineWidth', errorLineWidth, 'CapSize', 8, 'DisplayName', '80%');
h85 = errorbar(eps_data(2), y_data(2), y_err(2), 's', 'Color', c85, 'MarkerFaceColor', c85, 'MarkerEdgeColor', c85, 'MarkerSize', markerSize, 'LineWidth', errorLineWidth, 'CapSize', 8, 'DisplayName', '85%');
h90 = errorbar(eps_data(1), y_data(1), y_err(1), 's', 'Color', c90, 'MarkerFaceColor', c90, 'MarkerEdgeColor', c90, 'MarkerSize', markerSize, 'LineWidth', errorLineWidth, 'CapSize', 8, 'DisplayName', '90%');

hModel = plot(eps_plot, S_norm, '-', 'Color', modelColor, 'LineWidth', modelLineWidth, 'DisplayName', 'Iterative v12B Pre-Tension Model');

xlabel('Engineering Strain');
ylabel('Normalized Sensitivity (-)');
xlim([0 0.36]);
ylim([0.85 1.60]);
grid on;
box on;
set(gca, 'FontSize', axisFontSize, 'LineWidth', 1.0);

legend([h75, h80, h85, h90, hBareMean, hModel], ...
       {'75%', '80%', '85%', '90%', 'Bare Port Reference', 'Iterative v12B Pre-Tension Model'}, ...
       'Location', 'northeast', 'FontSize', legendFontSize);

saveCurrentFigure(resultsDir, '01_main_normalized_sensitivity_v12b');

%% -----------------------------
% Diagnostic: transmission only
% ------------------------------
figure;
plot(eps_plot, T_cav_norm, '-', 'Color', physColor, 'LineWidth', 2.5); hold on;
yline(1.0, ':k');

xlabel('Engineering Strain');
ylabel('Normalized Transmission (-)');
title('v12B Cavity / Interface Transmission');
legend('T_{cav,norm}(\epsilon)', 'Reference = 1', 'Location', 'best');
grid on;
set(gca, 'FontSize', 14);

saveCurrentFigure(resultsDir, '02_transmission_v12b');

%% -----------------------------
% Diagnostic: iterative test-load response
% ------------------------------
figure;
yyaxis left
plot(eps_plot, dV_test_curve, 'b-', 'LineWidth', 2); hold on;
plot(eps_plot, Vi_test_curve, 'g--', 'LineWidth', 2);
ylabel('\Delta V_{test}(\epsilon), V_{i,test}(\epsilon) [m^3]');

yyaxis right
plot(eps_plot, abs(wmax_test_curve), 'm-.', 'LineWidth', 2);
ylabel('|w_{max,test}(\epsilon)| [m]');

xlabel('Engineering Strain');
title('v12B Test-Load Volume Change and Deflection');
legend('\Delta V_{test}(\epsilon)', 'V_{i,test}(\epsilon)', '|w_{max,test}(\epsilon)|', 'Location', 'best');
grid on;
set(gca, 'FontSize', 14);

saveCurrentFigure(resultsDir, '03_test_load_volume_and_deflection_v12b');

%% -----------------------------
% Diagnostic: mechanics
% ------------------------------
figure;
yyaxis left
plot(eps_plot, Eeff_curve, 'k-', 'LineWidth', 2); hold on;
ylabel('E_{eff}(\epsilon) [Pa]');

yyaxis right
plot(eps_plot, Dplate_curve, 'c--', 'LineWidth', 2); hold on;
plot(eps_plot, T0_curve, 'r-.', 'LineWidth', 2);
plot(eps_plot, 1e3*a_plate_curve, 'm:', 'LineWidth', 2.2);
ylabel('D_{plate}(\epsilon) [N m], T_0(\epsilon) [N/m], a_{plate}(\epsilon) [mm]');

xlabel('Engineering Strain');
title('v12B Mechanics Diagnostics');
legend('E_{eff}(\epsilon)', 'D_{plate}(\epsilon)', 'T_0(\epsilon)', 'a_{plate}(\epsilon)', 'Location', 'best');
grid on;
set(gca, 'FontSize', 14);

saveCurrentFigure(resultsDir, '04_mechanics_v12b');

%% -----------------------------
% Diagnostic: installed interface diameter
% ------------------------------
figure;
plot(eps_plot, Dplate_geom_curve_mm, 'b-', 'LineWidth', 2); hold on;
plot(engStrains, scaleDiameters, 'ro', 'MarkerFaceColor', 'r');
xlabel('Engineering Strain');
ylabel('Installed Interface Diameter [mm]');
title('Sampling Interface Geometry vs Engineering Strain');
legend('Interpolated D_{plate}(\epsilon)', 'Measured interface states', 'Location', 'best');
grid on;
set(gca, 'FontSize', 14);

saveCurrentFigure(resultsDir, '05_installed_interface_diameter_v12b');

%% -----------------------------
% Diagnostic: iteration count
% ------------------------------
figure;
plot(eps_plot, iter_test_curve, 'k-', 'LineWidth', 2); hold on;
yline(p.maxIter, '--r', 'MaxIter');
xlabel('Engineering Strain');
ylabel('Iterations to Converge');
title('v12B Iterative Solver Diagnostics');
grid on;
set(gca, 'FontSize', 14);

saveCurrentFigure(resultsDir, '06_iteration_count_v12b');

%% -----------------------------
% Sensitivity sweep: kT0
% ------------------------------
kT0_results = zeros(numel(kT0_sweep), numel(eps_plot));
kT0_rmse    = zeros(numel(kT0_sweep), 1);
kT0_wrmse   = zeros(numel(kT0_sweep), 1);

for j = 1:numel(kT0_sweep)
    pSweep = p;
    pSweep.kT0 = kT0_sweep(j);

    [S_sweep, ~] = runNormalizedModelResponse_v12b(eps_plot, pSweep);

    kT0_results(j, :) = S_sweep;
    [kT0_rmse(j), kT0_wrmse(j)] = computeModelErrors(eps_plot, S_sweep, eps_data, y_data, y_err);
end

figure; hold on;
plot(eps_plot, kT0_results(1, :), '-', 'LineWidth', 2.0, 'DisplayName', sprintf('k_{T0} = %.2f', kT0_sweep(1)));
plot(eps_plot, kT0_results(2, :), '-', 'LineWidth', 2.5, 'DisplayName', sprintf('k_{T0} = %.2f', kT0_sweep(2)));
plot(eps_plot, kT0_results(3, :), '-', 'LineWidth', 2.0, 'DisplayName', sprintf('k_{T0} = %.2f', kT0_sweep(3)));

errorbar(eps_data(4), y_data(4), y_err(4), 's', 'Color', c75, 'MarkerFaceColor', c75, 'MarkerSize', markerSize, 'LineWidth', errorLineWidth, 'CapSize', 8, 'HandleVisibility', 'off');
errorbar(eps_data(3), y_data(3), y_err(3), 's', 'Color', c80, 'MarkerFaceColor', c80, 'MarkerSize', markerSize, 'LineWidth', errorLineWidth, 'CapSize', 8, 'HandleVisibility', 'off');
errorbar(eps_data(2), y_data(2), y_err(2), 's', 'Color', c85, 'MarkerFaceColor', c85, 'MarkerSize', markerSize, 'LineWidth', errorLineWidth, 'CapSize', 8, 'HandleVisibility', 'off');
errorbar(eps_data(1), y_data(1), y_err(1), 's', 'Color', c90, 'MarkerFaceColor', c90, 'MarkerSize', markerSize, 'LineWidth', errorLineWidth, 'CapSize', 8, 'HandleVisibility', 'off');
yline(1.0, '-', 'Color', bareColor, 'LineWidth', 1.6, 'HandleVisibility', 'off');

xlabel('Engineering Strain');
ylabel('Normalized Sensitivity (-)');
title('Sensitivity Sweep: Pre-Tension Scale Factor k_{T0}');
xlim([0 0.36]);
ylim([0.85 1.60]);
grid on;
box on;
set(gca, 'FontSize', 14);
legend('Location', 'northwest');

saveCurrentFigure(resultsDir, '07_sensitivity_kT0_v12b');

%% -----------------------------
% 2D sweep: kT0 vs strain contour
% ------------------------------
kT0_contour_results = zeros(numel(kT0_grid), numel(eps_plot));
kT0_contour_rmse    = zeros(numel(kT0_grid), 1);
kT0_contour_wrmse   = zeros(numel(kT0_grid), 1);

for j = 1:numel(kT0_grid)
    pSweep = p;
    pSweep.kT0 = kT0_grid(j);

    [S_sweep, ~] = runNormalizedModelResponse_v12b(eps_plot, pSweep);

    kT0_contour_results(j, :) = S_sweep;
    [kT0_contour_rmse(j), kT0_contour_wrmse(j)] = ...
        computeModelErrors(eps_plot, S_sweep, eps_data, y_data, y_err);
end

figure;
contourf(eps_plot, kT0_grid, kT0_contour_results, 24, 'LineColor', 'none');
hold on;
colormap(parula);

cb = colorbar;
cb.Label.String = 'Normalized Sensitivity, S_{norm} (-)';
cb.Label.FontSize = 14;

hBase = plot([eps_plot(1), eps_plot(end)], [p.kT0, p.kT0], 'k--', ...
    'LineWidth', 2.0, ...
    'DisplayName', sprintf('Baseline k_{T0} = %.2f', p.kT0));

for i = 1:numel(eps_data)
    plot(eps_data(i), p.kT0, 'wo', ...
        'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', 'k', ...
        'MarkerSize', 7, ...
        'HandleVisibility', 'off');
end

xlabel('Engineering Strain');
ylabel('k_{T0}');
title('2D Sensitivity Sweep: k_{T0} vs Engineering Strain');
xlim([0 0.36]);
ylim([kT0_grid(1), kT0_grid(end)]);
grid on;
box on;
set(gca, 'FontSize', 14);

legend(hBase, sprintf('Baseline k_{T0} = %.2f', p.kT0), ...
    'Location', 'northwest');

saveCurrentFigure(resultsDir, '08_sensitivity_kT0_contour_v12b');

%% -----------------------------
% Diagnostic: kT0 error metrics
% ------------------------------
figure;
yyaxis left
plot(kT0_grid, kT0_contour_rmse, 'b-', 'LineWidth', 2.5); hold on;
ylabel('RMSE (-)');

yyaxis right
plot(kT0_grid, kT0_contour_wrmse, 'r--', 'LineWidth', 2.5);
ylabel('WRMSE (-)');

xline(p.kT0, 'k--', 'LineWidth', 1.5, 'Label', 'Baseline', ...
    'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom');

xlabel('k_{T0}');
title('Fit Error vs Pre-Tension Scale Factor');
grid on;
box on;
set(gca, 'FontSize', 14);
legend('RMSE', 'WRMSE', 'Location', 'best');

saveCurrentFigure(resultsDir, '09_kT0_error_metrics_v12b');

%% -----------------------------
% Print sensitivity summary
% ------------------------------
fprintf('\n--- v12B sensitivity sweep summary ---\n');

fprintf('\nkT0 sweep:\n');
for j = 1:numel(kT0_sweep)
    fprintf('  kT0 = %.3f  | RMSE = %.4f | WRMSE = %.4f\n', ...
        kT0_sweep(j), kT0_rmse(j), kT0_wrmse(j));
end

fprintf('\n2D kT0 contour sweep:\n');
fprintf('  kT0 range: %.3f to %.3f (%d cases)\n', ...
    kT0_grid(1), kT0_grid(end), numel(kT0_grid));

[minRMSE, idxMinRMSE]   = min(kT0_contour_rmse);
[minWRMSE, idxMinWRMSE] = min(kT0_contour_wrmse);

fprintf('  Minimum RMSE  = %.4f at kT0 = %.3f\n', ...
    minRMSE, kT0_grid(idxMinRMSE));
fprintf('  Minimum WRMSE = %.4f at kT0 = %.3f\n', ...
    minWRMSE, kT0_grid(idxMinWRMSE));

%% -----------------------------
% Local functions
% ------------------------------
function [T_cav, Pi_test_curve, Vi_test_curve, dV_test_curve, Dplate_curve, ...
          Eeff_curve, T0_curve, wmax_test_curve, iter_test_curve, a_plate_curve, ...
          Dplate_geom_curve_mm, converged_curve] = ...
          runModelOverStrainRange_v12b(eps_plot, p)

    T_cav                = zeros(size(eps_plot));
    Pi_test_curve        = zeros(size(eps_plot));
    Vi_test_curve        = zeros(size(eps_plot));
    dV_test_curve        = zeros(size(eps_plot));
    Dplate_curve         = zeros(size(eps_plot));
    Eeff_curve           = zeros(size(eps_plot));
    T0_curve             = zeros(size(eps_plot));
    wmax_test_curve      = zeros(size(eps_plot));
    iter_test_curve      = zeros(size(eps_plot));
    a_plate_curve        = zeros(size(eps_plot));
    Dplate_geom_curve_mm = zeros(size(eps_plot));
    converged_curve      = zeros(size(eps_plot));

    solverState.Pi = p.P_ref;

    for k = 1:numel(eps_plot)
        epsn = eps_plot(k);

        [T_cav(k), Pi_test_curve(k), Vi_test_curve(k), dV_test_curve(k), ...
            Dplate_curve(k), Eeff_curve(k), T0_curve(k), wmax_test_curve(k), ...
            iter_test_curve(k), a_plate_curve(k), Dplate_geom_curve_mm(k), ...
            converged_curve(k), solverState] = ...
            evaluateCavityTransmission_v12b(epsn, p, solverState);
    end
end

function [S_norm, T_cav_norm] = runNormalizedModelResponse_v12b(eps_plot, p)
    [T_cav, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = runModelOverStrainRange_v12b(eps_plot, p);

    eps_ref_norm = 0.10;
    Tref = interp1(eps_plot, T_cav, eps_ref_norm, 'linear');

    if abs(Tref) < 1e-12
        error('Reference transmission for normalization is too small. Inspect raw T_cav first.');
    end

    T_cav_norm = T_cav ./ Tref;
    S_norm     = T_cav_norm;
end

function [rmse, wrmse] = computeModelErrors(eps_plot, S_norm, eps_data, y_data, y_err)
    resid = zeros(size(eps_data));
    wres  = zeros(size(eps_data));

    for i = 1:numel(eps_data)
        Si_model = interp1(eps_plot, S_norm, eps_data(i), 'linear');
        resid(i) = Si_model - y_data(i);
        wres(i)  = resid(i) / y_err(i);
    end

    rmse  = sqrt(mean(resid.^2));
    wrmse = sqrt(mean(wres.^2));
end

function [Tcav, Pi_test, Vi_test, dV_test, Dplate, Eeff, T0mem, wmax_test, ...
          iter_test, a_plate, Dplate_geom_mm, converged, solverStateOut] = ...
          evaluateCavityTransmission_v12b(epsn, p, solverStateIn)

    [Dplate, Eeff, T0mem, a_plate, Dplate_geom_mm] = evaluateMembraneState_v12b(epsn, p);

    % Forward test solve used for diagnostics and continuation
    Psurf_test = p.P_ref + p.dP0;
    [Pi_test, Vi_test, dV_test, wmax_test, iter_test, converged, solverStateOut] = ...
        solveCavityPressure_iterative_v12b(Psurf_test, a_plate, Dplate, T0mem, p, solverStateIn);

    % Symmetric differential evaluation
    symmetricState = solverStateIn;
    if isempty(symmetricState)
        symmetricState.Pi = p.P_ref;
    end

    dPsens_plus  = sensorDifferentialResponse_v12b(+p.dP0, a_plate, Dplate, T0mem, p, symmetricState);
    dPsens_minus = sensorDifferentialResponse_v12b(-p.dP0, a_plate, Dplate, T0mem, p, symmetricState);

    Tcav = (dPsens_plus - dPsens_minus) / (2 * p.dP0);
end

function [a_plate, D_plate_mm, t_plate] = evaluatePlateGeometry_v12b(epsn, p)
    D_plate_mm = interp1(p.engStrains, p.scaleDiameters_mm, epsn, 'linear', 'extrap');
    a_plate = 0.5 * D_plate_mm * 1e-3;
    t_plate = p.t_plate0;
end

function [Dplate, Eeff, T0mem, a_plate, D_plate_mm] = evaluateMembraneState_v12b(epsn, p)
    [a_plate, D_plate_mm, t_plate] = evaluatePlateGeometry_v12b(epsn, p);

    if p.useStrainDependentE
        Eeff = p.E_plate0 * (1 + p.c1 * epsn + p.c2 * epsn^2);
    else
        Eeff = p.E_plate0;
    end

    % Classical plate bending stiffness
    Dplate = Eeff * t_plate^3 / (12 * (1 - p.nu_plate^2));

    % Equibiaxial pre-strain -> in-plane tensile resultant
    % T0 = kT0 * E*h*eps/(1-nu)
    T0mem = p.kT0 * Eeff * t_plate * max(epsn, 0) / (1 - p.nu_plate);
end

function dPsens = sensorDifferentialResponse_v12b(dPsurf, a_plate, Dplate, T0mem, p, solverStateIn)
    % Differential sensor wrapper:
    % this function evaluates the two-sided IPS response by calling the
    % single-cavity solver once for each side.
    %
    % Side 1:
    %   q1 = P3 - P1
    %
    % Side 2:
    %   q2 = P4 - P2
    %
    % The final differential sensor output is:
    %   dPsens = P1 - P2
    P3 = p.P_ref + dPsurf / 2;
    P4 = p.P_ref - dPsurf / 2;

    stateA = solverStateIn;
    stateB = solverStateIn;

    [P1, ~, ~, ~, ~, ~, ~] = solveCavityPressure_iterative_v12b(P3, a_plate, Dplate, T0mem, p, stateA);
    [P2, ~, ~, ~, ~, ~, ~] = solveCavityPressure_iterative_v12b(P4, a_plate, Dplate, T0mem, p, stateB);

    dPsens = P1 - P2;
end

function [Pi, Vi, dV, wmax, iter, converged, solverStateOut] = ...
    solveCavityPressure_iterative_v12b(Psurf, a_plate, Dplate, T0mem, p, solverStateIn)

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
        % Single-cavity membrane loading:
        % q = P3 - P1
        % where:
        %   P3 = external pressure applied to this membrane/cavity
        %   P1 = internal cavity pressure
        P3 = Psurf;
        P1 = Pi;
    
        q_membrane = P3 - P1;
    
        [Vi_new, dV_new, wmax_new] = loadedCavityState_fromMembraneLoad_v12b(q_membrane, p.V00, a_plate, Dplate, T0mem, p);
    
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

function [Vi, dV, wmax] = loadedCavityState_fromMembraneLoad_v12b(q_membrane, V00, a, Dplate, T0mem, p)
    % v12B reduced-order hybrid plate + pre-tension model:
    %
    % Single-cavity membrane load:
    %   q_membrane = P3 - P1
    %
    % Reduced-order response:
    %   w(r) = W0 * (1 - r^2/a^2)^2
    %   W0   = q_membrane * a^4 / (64*D + 4*T0*a^2)
    %   dV   = pi*a^2/3 * W0
    %
    % Positive q_membrane means external pressure exceeds cavity pressure,
    % producing positive inward membrane deflection and a positive cavity
    % volume reduction term dV, so that:
    %   Vi = V00 - dV

    denom = 64 * Dplate + 4 * T0mem * a^2;
    denom = max(denom, 1e-18);
    
    wmax = q_membrane * a^4 / denom;
    dV   = pi * a^2 * wmax / 3;

    dV_max = p.dV_cap_fraction * V00;
    dV = max(min(dV, dV_max), -dV_max);

    Vi = V00 - dV;
    Vi = max(Vi, 1e-12);
end

function saveCurrentFigure(resultsDir, baseName) 
    fig = gcf;
    set(fig, 'PaperPositionMode', 'auto');
    exportgraphics(fig, fullfile(resultsDir, [baseName, '.png']), 'Resolution', 300);
    %savefig(fig, fullfile(resultsDir, [baseName, '.fig']));
end