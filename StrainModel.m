%% sensor_interface_numerical_model_v11.m
clear; clc; close all;

%% -----------------------------
% Results folder
% ------------------------------
resultsDir = fullfile(pwd, 'results_v11');
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
geomOnlyColor   = [0.10, 0.55, 0.85];
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

a_plate_ref    = 0.5 * d0 * 1e-3;
A_plate_ref    = pi * a_plate_ref^2;
p.geom_raw_ref = p.A_forced / A_plate_ref;

%% -----------------------------
% Plate / membrane model
% ------------------------------
p.t_plate0 = 0.020 * 0.0254;  % m
p.nu_plate = 0.49;
p.E_plate0 = 6.0e5;           % Pa

% Strain-dependent effective modulus:
% E_eff(eps) = E0 * (1 + c1*eps + c2*eps^2)
% c1 = linear stiffening coefficient
% c2 = quadratic stiffening coefficient
%
% These are retained for completeness, but v10 showed low model sensitivity
% to them over the current strain range. They are therefore held fixed in v11.
p.c1 = 0.9;
p.c2 = 1.4;

%% -----------------------------
% Iterative solver settings
% ------------------------------
p.tolP    = 1e-5;
p.maxIter = 1000;
p.relax   = 0.05;
p.relax_q = 0.05;

%% -----------------------------
% Stability / geometry coupling settings
% ------------------------------
p.beta_geom       = 0.25;
p.q_cap           = 50;     % Pa
p.dV_cap_fraction = 0.10;   % max allowed |dV| / V00

%% -----------------------------
% Sensitivity sweep settings
% ------------------------------
beta_geom_sweep = [0.15, 0.25, 0.35];
beta_geom_grid  = linspace(0.10, 0.40, 41);

%% -----------------------------
% Print configuration summary
% ------------------------------
fprintf('Nominal cavity volume V00: %.3e m^3\n', p.V00);
fprintf('Plate thickness: %.3e m\n', p.t_plate0);
fprintf('Baseline effective modulus: %.3e Pa\n', p.E_plate0);
fprintf('Forced loading radius: %.3e m\n', p.r_forced_m);
fprintf('Pressure tolerance: %.3e Pa\n', p.tolP);
fprintf('Max iterations: %d\n', p.maxIter);
fprintf('Pressure relaxation factor: %.3f\n', p.relax);
fprintf('Net-load relaxation factor: %.3f\n', p.relax_q);
fprintf('Geometry coupling parameter beta_geom: %.3f\n', p.beta_geom);
fprintf('Net-load cap q_cap: %.3f Pa\n', p.q_cap);
fprintf('Volume-change cap fraction: %.3f\n', p.dV_cap_fraction);
fprintf('Linear modulus coefficient c1: %.3f\n', p.c1);
fprintf('Quadratic modulus coefficient c2: %.3f\n\n', p.c2);

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
 Eeff_curve, wmax_test_curve, iter_test_curve, a_plate_curve, ...
 Dplate_geom_curve_mm, geom_gain_curve, G_eff_curve, converged_curve] = ...
    runModelOverStrainRange(eps_plot, p);

T0 = T_cav(1);
if abs(T0) < eps
    warning('Zero-strain transmission T0 is near zero. Using eps safeguard for normalization.');
    T0 = eps;
end

T_cav_norm = T_cav ./ T0;
S_norm     = T_cav_norm;

% Geometry-only comparison curve
S_geom_only = G_eff_curve;

%% -----------------------------
% Print summary at experimental points
% ------------------------------
fprintf('\n--- Baseline model summary at experimental strain values ---\n');
fprintf('   eps      data      model    geomOnly  T_cav,norm   a_plate(mm)  G_geom   G_eff    resid    w_resid\n');

resid = zeros(size(eps_data));
wres  = zeros(size(eps_data));

for i = 1:numel(eps_data)
    Si_model = interp1(eps_plot, S_norm,             eps_data(i), 'linear');
    Si_geom  = interp1(eps_plot, S_geom_only,        eps_data(i), 'linear');
    Ti_cav   = interp1(eps_plot, T_cav_norm,         eps_data(i), 'linear');
    ai_mm    = interp1(eps_plot, 1e3*a_plate_curve,  eps_data(i), 'linear');
    Gi       = interp1(eps_plot, geom_gain_curve,    eps_data(i), 'linear');
    Geffi    = interp1(eps_plot, G_eff_curve,        eps_data(i), 'linear');

    resid(i) = Si_model - y_data(i);
    wres(i)  = resid(i) / y_err(i);

    fprintf('%7.3f   %7.3f   %7.3f   %7.3f   %7.3f   %10.3f  %7.3f  %7.3f   %7.3f   %7.3f\n', ...
        eps_data(i), y_data(i), Si_model, Si_geom, Ti_cav, ai_mm, Gi, Geffi, resid(i), wres(i));
end

rmse  = sqrt(mean(resid.^2));
wrmse = sqrt(mean(wres.^2));

fprintf('\nBaseline RMSE  = %.4f\n', rmse);
fprintf('Baseline WRMSE = %.4f\n', wrmse);
fprintf('Baseline converged points: %d / %d\n', sum(converged_curve > 0.5), numel(converged_curve));

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

hGeomOnly = plot(eps_plot, S_geom_only, '--', 'Color', geomOnlyColor, 'LineWidth', 2.0, 'DisplayName', 'Geometry-Only Gain');
hModel    = plot(eps_plot, S_norm, '-', 'Color', modelColor, 'LineWidth', modelLineWidth, 'DisplayName', 'Iterative v11 Force-Balance Model');

xlabel('Engineering Strain');
ylabel('Normalized Sensitivity (-)');
xlim([0 0.36]);
ylim([0.85 1.60]);
grid on;
box on;
set(gca, 'FontSize', axisFontSize, 'LineWidth', 1.0);

legend([h75, h80, h85, h90, hBareMean, hGeomOnly, hModel], ...
       {'75%', '80%', '85%', '90%', 'Bare Port Reference', 'Geometry-Only Gain', 'Iterative v11 Force-Balance Model'}, ...
       'Location', 'northeast', 'FontSize', legendFontSize);

saveCurrentFigure(resultsDir, '01_main_normalized_sensitivity_v11');

%% -----------------------------
% Diagnostic: transmission and embedded geometry
% ------------------------------
figure;
plot(eps_plot, T_cav_norm, '-', 'Color', physColor, 'LineWidth', 2.0); hold on;
plot(eps_plot, geom_gain_curve, '--', 'Color', [0.0 0.2 0.8], 'LineWidth', 2.0);
plot(eps_plot, G_eff_curve, ':', 'Color', [0.6 0.1 0.8], 'LineWidth', 2.5);
plot(eps_plot, S_norm, 'k-', 'LineWidth', 2.5);
yline(1.0, ':k');

xlabel('Engineering Strain');
ylabel('Magnitude (-)');
title('Transmission and Embedded Geometry Effect');
legend('T_{cav,norm}(\epsilon)', 'G_{geom}(\epsilon)', 'G_{eff}(\epsilon)', 'S_{norm}(\epsilon)', 'Reference = 1', 'Location', 'best');
grid on;
set(gca, 'FontSize', 14);

saveCurrentFigure(resultsDir, '02_transmission_and_embedded_geometry_v11');

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
title('Iterative Test-Load Volume Change and Plate Deflection');
legend('\Delta V_{test}(\epsilon)', 'V_{i,test}(\epsilon)', '|w_{max,test}(\epsilon)|', 'Location', 'best');
grid on;
set(gca, 'FontSize', 14);

saveCurrentFigure(resultsDir, '03_test_load_volume_and_deflection_v11');

%% -----------------------------
% Diagnostic: plate mechanics and geometry
% ------------------------------
figure;
yyaxis left
plot(eps_plot, Eeff_curve, 'k-', 'LineWidth', 2); hold on;
ylabel('E_{eff}(\epsilon) [Pa]');

yyaxis right
plot(eps_plot, Dplate_curve, 'c--', 'LineWidth', 2); hold on;
plot(eps_plot, 1e3*a_plate_curve, 'm-.', 'LineWidth', 2);
ylabel('D_{plate}(\epsilon), a_{plate}(\epsilon) [mm]');

xlabel('Engineering Strain');
title('Plate Mechanics Diagnostics');
legend('E_{eff}(\epsilon)', 'D_{plate}(\epsilon)', 'a_{plate}(\epsilon)', 'Location', 'best');
grid on;
set(gca, 'FontSize', 14);

saveCurrentFigure(resultsDir, '04_plate_mechanics_and_geometry_v11');

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

saveCurrentFigure(resultsDir, '05_installed_interface_diameter_v11');

%% -----------------------------
% Diagnostic: iteration count
% ------------------------------
figure;
plot(eps_plot, iter_test_curve, 'k-', 'LineWidth', 2); hold on;
yline(p.maxIter, '--r', 'MaxIter');
xlabel('Engineering Strain');
ylabel('Iterations to Converge');
title('Iterative Solver Diagnostics');
grid on;
set(gca, 'FontSize', 14);

saveCurrentFigure(resultsDir, '06_iteration_count_v11');

%% -----------------------------
% Diagnostic: convergence flag
% ------------------------------
figure;
plot(eps_plot, converged_curve, 'b-', 'LineWidth', 2);
xlabel('Engineering Strain');
ylabel('Converged (1=yes, 0=no)');
title('Convergence Flag');
ylim([-0.1, 1.1]);
grid on;
set(gca, 'FontSize', 14);

saveCurrentFigure(resultsDir, '07_convergence_flag_v11');

%% -----------------------------
% Sensitivity sweep: beta_geom
% ------------------------------
beta_results = zeros(numel(beta_geom_sweep), numel(eps_plot));
beta_rmse    = zeros(numel(beta_geom_sweep), 1);
beta_wrmse   = zeros(numel(beta_geom_sweep), 1);

for j = 1:numel(beta_geom_sweep)
    pSweep = p;
    pSweep.beta_geom = beta_geom_sweep(j);

    [S_sweep, ~] = runNormalizedModelResponse(eps_plot, pSweep);

    beta_results(j, :) = S_sweep;
    [beta_rmse(j), beta_wrmse(j)] = computeModelErrors(eps_plot, S_sweep, eps_data, y_data, y_err);
end

figure; hold on;
plot(eps_plot, beta_results(1, :), '-', 'LineWidth', 2.0, 'DisplayName', sprintf('\\beta_{geom} = %.2f', beta_geom_sweep(1)));
plot(eps_plot, beta_results(2, :), '-', 'LineWidth', 2.5, 'DisplayName', sprintf('\\beta_{geom} = %.2f', beta_geom_sweep(2)));
plot(eps_plot, beta_results(3, :), '-', 'LineWidth', 2.0, 'DisplayName', sprintf('\\beta_{geom} = %.2f', beta_geom_sweep(3)));

errorbar(eps_data(4), y_data(4), y_err(4), 's', 'Color', c75, 'MarkerFaceColor', c75, 'MarkerSize', markerSize, 'LineWidth', errorLineWidth, 'CapSize', 8, 'HandleVisibility', 'off');
errorbar(eps_data(3), y_data(3), y_err(3), 's', 'Color', c80, 'MarkerFaceColor', c80, 'MarkerSize', markerSize, 'LineWidth', errorLineWidth, 'CapSize', 8, 'HandleVisibility', 'off');
errorbar(eps_data(2), y_data(2), y_err(2), 's', 'Color', c85, 'MarkerFaceColor', c85, 'MarkerSize', markerSize, 'LineWidth', errorLineWidth, 'CapSize', 8, 'HandleVisibility', 'off');
errorbar(eps_data(1), y_data(1), y_err(1), 's', 'Color', c90, 'MarkerFaceColor', c90, 'MarkerSize', markerSize, 'LineWidth', errorLineWidth, 'CapSize', 8, 'HandleVisibility', 'off');
yline(1.0, '-', 'Color', bareColor, 'LineWidth', 1.6, 'HandleVisibility', 'off');

xlabel('Engineering Strain');
ylabel('Normalized Sensitivity (-)');
title('Sensitivity Sweep: \beta_{geom}');
xlim([0 0.36]);
ylim([0.85 1.60]);
grid on;
box on;
set(gca, 'FontSize', 14);
legend('Location', 'northwest');

saveCurrentFigure(resultsDir, '08_sensitivity_beta_geom_v11');

%% -----------------------------
% 2D sweep: beta_geom vs strain contour
% ------------------------------
beta_contour_results = zeros(numel(beta_geom_grid), numel(eps_plot));
beta_contour_rmse    = zeros(numel(beta_geom_grid), 1);
beta_contour_wrmse   = zeros(numel(beta_geom_grid), 1);

for j = 1:numel(beta_geom_grid)
    pSweep = p;
    pSweep.beta_geom = beta_geom_grid(j);

    [S_sweep, ~] = runNormalizedModelResponse(eps_plot, pSweep);

    beta_contour_results(j, :) = S_sweep;
    [beta_contour_rmse(j), beta_contour_wrmse(j)] = ...
        computeModelErrors(eps_plot, S_sweep, eps_data, y_data, y_err);
end

figure;
contourf(eps_plot, beta_geom_grid, beta_contour_results, 24, 'LineColor', 'none');
hold on;
colormap(parula);

cb = colorbar;
cb.Label.String = 'Normalized Sensitivity, S_{norm} (-)';
cb.Label.FontSize = 14;

hBase = plot([eps_plot(1), eps_plot(end)], [p.beta_geom, p.beta_geom], 'k--', ...
    'LineWidth', 2.0, ...
    'DisplayName', sprintf('Baseline \\beta_{geom} = %.2f', p.beta_geom));

for i = 1:numel(eps_data)
    plot(eps_data(i), p.beta_geom, 'wo', ...
        'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', 'k', ...
        'MarkerSize', 7, ...
        'HandleVisibility', 'off');
end

xlabel('Engineering Strain');
ylabel('\beta_{geom}');
title('2D Sensitivity Sweep: \beta_{geom} vs Engineering Strain');
xlim([0 0.36]);
ylim([beta_geom_grid(1), beta_geom_grid(end)]);
grid on;
box on;
set(gca, 'FontSize', 14);

legend(hBase, sprintf('Baseline \\beta_{geom} = %.2f', p.beta_geom), ...
    'Location', 'northwest');

saveCurrentFigure(resultsDir, '09_sensitivity_beta_geom_contour_v11');

%% -----------------------------
% Diagnostic: beta_geom error metrics
% ------------------------------
figure;
yyaxis left
plot(beta_geom_grid, beta_contour_rmse, 'b-', 'LineWidth', 2.5); hold on;
ylabel('RMSE (-)');

yyaxis right
plot(beta_geom_grid, beta_contour_wrmse, 'r--', 'LineWidth', 2.5);
ylabel('WRMSE (-)');

xline(p.beta_geom, 'k--', 'LineWidth', 1.5, 'Label', 'Baseline', ...
    'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom');

xlabel('\beta_{geom}');
title('Fit Error vs Geometry Coupling Parameter');
grid on;
box on;
set(gca, 'FontSize', 14);
legend('RMSE', 'WRMSE', 'Location', 'best');

saveCurrentFigure(resultsDir, '10_beta_geom_error_metrics_v11');

%% -----------------------------
% Print sensitivity summary
% ------------------------------
fprintf('\n--- Sensitivity sweep summary ---\n');

fprintf('\nBeta_geom sweep:\n');
for j = 1:numel(beta_geom_sweep)
    fprintf('  beta_geom = %.3f  | RMSE = %.4f | WRMSE = %.4f\n', ...
        beta_geom_sweep(j), beta_rmse(j), beta_wrmse(j));
end

fprintf('\n2D beta_geom contour sweep:\n');
fprintf('  beta range: %.3f to %.3f (%d cases)\n', ...
    beta_geom_grid(1), beta_geom_grid(end), numel(beta_geom_grid));

[minRMSE, idxMinRMSE]   = min(beta_contour_rmse);
[minWRMSE, idxMinWRMSE] = min(beta_contour_wrmse);

fprintf('  Minimum RMSE  = %.4f at beta_geom = %.3f\n', ...
    minRMSE, beta_geom_grid(idxMinRMSE));
fprintf('  Minimum WRMSE = %.4f at beta_geom = %.3f\n', ...
    minWRMSE, beta_geom_grid(idxMinWRMSE));

%% -----------------------------
% Local functions
% ------------------------------
function [T_cav, Pi_test_curve, Vi_test_curve, dV_test_curve, Dplate_curve, ...
          Eeff_curve, wmax_test_curve, iter_test_curve, a_plate_curve, ...
          Dplate_geom_curve_mm, geom_gain_curve, G_eff_curve, converged_curve] = ...
          runModelOverStrainRange(eps_plot, p)

    T_cav                = zeros(size(eps_plot));
    Pi_test_curve        = zeros(size(eps_plot));
    Vi_test_curve        = zeros(size(eps_plot));
    dV_test_curve        = zeros(size(eps_plot));
    Dplate_curve         = zeros(size(eps_plot));
    Eeff_curve           = zeros(size(eps_plot));
    wmax_test_curve      = zeros(size(eps_plot));
    iter_test_curve      = zeros(size(eps_plot));
    a_plate_curve        = zeros(size(eps_plot));
    Dplate_geom_curve_mm = zeros(size(eps_plot));
    geom_gain_curve      = zeros(size(eps_plot));
    G_eff_curve          = zeros(size(eps_plot));
    converged_curve      = zeros(size(eps_plot));

    solverState.Pi = p.P_ref;
    solverState.q_net_old = 0.0;

    for k = 1:numel(eps_plot)
        epsn = eps_plot(k);

        [T_cav(k), Pi_test_curve(k), Vi_test_curve(k), dV_test_curve(k), ...
            Dplate_curve(k), Eeff_curve(k), wmax_test_curve(k), iter_test_curve(k), ...
            a_plate_curve(k), Dplate_geom_curve_mm(k), geom_gain_curve(k), ...
            G_eff_curve(k), converged_curve(k), solverState] = ...
            evaluateCavityTransmission_v11(epsn, p, solverState);
    end
end

function [S_norm, T_cav_norm] = runNormalizedModelResponse(eps_plot, p)
    [T_cav, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = runModelOverStrainRange(eps_plot, p);

    T0 = T_cav(1);
    if abs(T0) < eps
        T0 = eps;
    end

    T_cav_norm = T_cav ./ T0;
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

function [Tcav, Pi_test, Vi_test, dV_test, Dplate, Eeff, wmax_test, iter_test, ...
          a_plate, Dplate_geom_mm, geom_gain, G_eff, converged, solverStateOut] = ...
          evaluateCavityTransmission_v11(epsn, p, solverStateIn)

    [Dplate, Eeff, a_plate, ~, Dplate_geom_mm] = evaluatePlateStiffness(epsn, p);

    A_plate   = pi * a_plate^2;
    geom_raw  = p.A_forced / A_plate;
    geom_gain = geom_raw / p.geom_raw_ref;
    G_eff     = 1 + p.beta_geom * (geom_gain - 1);

    % Forward test solve used for diagnostics and continuation
    Psurf_test = p.P_ref + p.dP0;
    [Pi_test, Vi_test, dV_test, wmax_test, iter_test, converged, solverStateOut] = ...
        solveCavityPressure_iterative_v11(Psurf_test, epsn, a_plate, geom_gain, p, solverStateIn);

    % Symmetric differential evaluation:
    % use matched initial states for +dP0 and -dP0 to reduce path-dependent artifacts
    symmetricState = solverStateIn;
    if isempty(symmetricState)
        symmetricState.Pi = p.P_ref;
        symmetricState.q_net_old = 0.0;
    end

    dPsens_plus  = sensorDifferentialResponse_v11(+p.dP0, epsn, a_plate, geom_gain, p, symmetricState);
    dPsens_minus = sensorDifferentialResponse_v11(-p.dP0, epsn, a_plate, geom_gain, p, symmetricState);

    Tcav = (dPsens_plus - dPsens_minus) / (2 * p.dP0);
end

function [a_plate, D_plate_mm, t_plate] = evaluatePlateGeometry(epsn, p)
    D_plate_mm = interp1(p.engStrains, p.scaleDiameters_mm, epsn, 'linear', 'extrap');
    a_plate = 0.5 * D_plate_mm * 1e-3;
    t_plate = p.t_plate0;
end

function [Dplate, Eeff, a_plate, t_plate, D_plate_mm] = evaluatePlateStiffness(epsn, p)
    [a_plate, D_plate_mm, t_plate] = evaluatePlateGeometry(epsn, p);
    Eeff = p.E_plate0 * (1 + p.c1 * epsn + p.c2 * epsn^2);
    Dplate = Eeff * t_plate^3 / (12 * (1 - p.nu_plate^2));
end

function dPsens = sensorDifferentialResponse_v11(dPsurf, epsn, a_plate, geom_gain, p, solverStateIn)
    P3 = p.P_ref + dPsurf / 2;
    P4 = p.P_ref - dPsurf / 2;

    % Use identical initial states for both branches
    stateA = solverStateIn;
    stateB = solverStateIn;

    [P1, ~, ~, ~, ~, ~, ~] = solveCavityPressure_iterative_v11(P3, epsn, a_plate, geom_gain, p, stateA);
    [P2, ~, ~, ~, ~, ~, ~] = solveCavityPressure_iterative_v11(P4, epsn, a_plate, geom_gain, p, stateB);

    dPsens = P1 - P2;
end

function [Pi, Vi, dV, wmax, iter, converged, solverStateOut] = ...
    solveCavityPressure_iterative_v11(Psurf, epsn, a_plate, geom_gain, p, solverStateIn)

    [Dplate, ~, ~, ~, ~] = evaluatePlateStiffness(epsn, p);

    if nargin < 6 || isempty(solverStateIn)
        Pi = p.P_ref;
        q_net_old = 0;
    else
        Pi = solverStateIn.Pi;
        q_net_old = solverStateIn.q_net_old;
        if ~isfinite(Pi), Pi = p.P_ref; end
        if ~isfinite(q_net_old), q_net_old = 0; end
    end

    converged = false;
    Vi = p.V00;
    dV = 0;
    wmax = 0;

    for iter = 1:p.maxIter
        dPsurf = Psurf - p.P_ref;
        dPi    = Pi - p.P_ref;

        G_eff = 1 + p.beta_geom * (geom_gain - 1);

        q_net_raw = dPsurf * G_eff - dPi;
        q_net_raw = max(min(q_net_raw, p.q_cap), -p.q_cap);

        q_net = (1 - p.relax_q) * q_net_old + p.relax_q * q_net_raw;

        [Vi_new, dV_new, wmax_new] = loadedCavityState_fromDeltaP_v11(q_net, p.V00, a_plate, Dplate, p);

        Pi_new_raw = p.P_ref * p.V00 / Vi_new;
        Pi_new = (1 - p.relax) * Pi + p.relax * Pi_new_raw;

        if abs(Pi_new - Pi) < p.tolP
            Pi = Pi_new;
            Vi = Vi_new;
            dV = dV_new;
            wmax = wmax_new;
            converged = true;

            solverStateOut.Pi = Pi;
            solverStateOut.q_net_old = q_net;
            return;
        end

        Pi = Pi_new;
        q_net_old = q_net;
        Vi = Vi_new;
        dV = dV_new;
        wmax = wmax_new;
    end

    solverStateOut.Pi = Pi;
    solverStateOut.q_net_old = q_net_old;
end

function [Vi, dV, wmax] = loadedCavityState_fromDeltaP_v11(DeltaP, V00, a, Dplate, p)
    % Small-deflection circular plate response under uniform load
    wmax = DeltaP * a^4 / (64 * Dplate);
    dV   = DeltaP * pi * a^6 / (192 * Dplate);

    dV_max = p.dV_cap_fraction * V00;
    dV = max(min(dV, dV_max), -dV_max);

    Vi = V00 - dV;
    Vi = max(Vi, 1e-12);
end

function saveCurrentFigure(resultsDir, baseName)
    fig = gcf;
    set(fig, 'PaperPositionMode', 'auto');
    exportgraphics(fig, fullfile(resultsDir, [baseName, '.png']), 'Resolution', 300);
    savefig(fig, fullfile(resultsDir, [baseName, '.fig']));
end