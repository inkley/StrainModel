%% sensor_interface_numerical_model_v20.m
clear; clc; close all;

% v20
% v20 is a reduced-order force-balance model for the membrane-confined
% cavity/interface system. This version emphasizes physically defensible
% interpretation by using only the interface-relative model response,
% normalized by the zero-strain membrane condition.
%
% New in v20:
% 1) Removes geometry-amplification terms from the active model response.
% 2) Uses a fixed effective loaded radius independent of installed membrane
%    diameter.
% 3) Uses no empirical bare-port scale factor; model output is interpreted
%    as interface-relative transmission normalized by T(0).
% 4) Retains material, thickness, pressure-level, and calibration sweeps as
%    diagnostics rather than hidden fit terms.
%
% v20 assumptions:
% 1) axisymmetric reduced-order response
% 2) effective circular loaded region
% 3) small-deflection-style volume estimate regularized by saturation
% 4) fixed effective loaded radius independent of installed membrane diameter
% 5) strain-dependent thickness update based on incompressibility
% 6) pre-strain represented through a smoothed effective pre-tension term
% 7) cavity response solved iteratively using ideal-gas compression
% 8) closure terms approximate higher-order mechanics not explicitly resolved
%    in the reduced-order formulation

%% -----------------------------
% Results folder
% ------------------------------
resultsDir = fullfile(pwd, 'results_v20');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

%% -----------------------------
% Plot / output toggles
% ------------------------------
makeMainFigure                  = true;   % manuscript
makeTransmissionPlot            = true;   % manuscript/model interpretation
makeMechanicsPlot               = true;   % optional manuscript diagnostic
makeGeometryPlot                = false;

makeRMSEPlot                    = false;  % appendix/advisor only
makeWRMSEPlot                   = false;
makePressureSweepPlot           = false;
makePressureSurfacePlot         = false;
makePressureLevelMaterialPlots  = false;

makeMaterialSweepPlot           = true;  % appendix/advisor only
makeMaterialFitOverlayPlots     = true;  % appendix/advisor only

printSweepSummary               = true;
printPressureSweepSummary       = true;

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
% Plate / membrane model (v20)
% Reduced-order model with:
% 1) fixed effective loaded radius (no geometry amplification)
% 2) smoothed pre-tension engagement from installed strain
% 3) strain-dependent thickness update (incompressibility assumption)
% 4) smooth tanh regularization of cavity volume change (dV)
% 5) interface-relative transmission normalized by T(0)
% 6) no empirical bare-port scale factor
% 7) pressure-level and material sweeps used as diagnostics
% ------------------------------
p.t_plate0 = 0.020 * 0.0254;  % m
p.nu_plate = 0.49;
p.E_plate0 = 6.0e5;           % Pa

% Effective pre-tension scale factor.
% kT0 represents the fraction of the ideal equibiaxial membrane tension
% that is realized in the installed interface. This accounts for
% nonuniform strain during installation, clamping compliance,
% local seating effects, and departures from ideal linear-elastic behavior.
p.kT0 = 0.250;

% Characteristic strain for smooth activation of pre-tension.
p.eps_char = 0.160;

% No empirical bare-port scaling is applied.
% Model output is interpreted as interface-relative transmission,
% normalized by the zero-strain membrane condition.
p.S0_interface_bare = 1.0;

% Optional strain-dependent modulus (disabled in baseline model).
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

kT0_sweep = [0.025, 0.050, 0.100, 0.150, 0.200, 0.250, 0.300];
S0_sweep  = 1.0;

kT0_grid  = 0.025:0.025:0.400;
S0_grid   = 1.0;

%% -----------------------------
% pressure-level sweep settings
% ------------------------------
% These are the pre-strain cases where pressure-level dependence is checked.
% By default, this uses the same pre-strain values as the experimental data.
eps_pressure_cases = eps_data;

% Absolute external pressure on port 3 is swept around P_ref.
% Edit this range if you want to evaluate a wider pressure interval.
P3_offset_sweep = [-400, -200, 0, 200, 400]; % Pa relative to p.P_ref
P3_sweep        = p.P_ref + P3_offset_sweep;

% Local perturbation size used to estimate dP1/dP3 at each P3 level.
% This is intentionally small so the result behaves like a local sensitivity.
p.dP_local = p.dP0;

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
fprintf('Volume-change saturation fraction: %.3f\n', p.dV_cap_fraction);
fprintf('Use strain-dependent E: %d\n', p.useStrainDependentE);
fprintf('Pressure sweep local perturbation dP_local: %.3f Pa\n\n', p.dP_local);

fprintf('Sampling interface diameters (mm): ');
fprintf('%.3f ', scaleDiameters);
fprintf('\n');

fprintf('Engineering Strains, epsilon (-): ');
fprintf('%.3f ', engStrains);
fprintf('\n\n');

%% -----------------------------
% Evaluate baseline model over strain range
% ------------------------------
[T_cav, Dplate_curve, E_curve, Tpre_eff_curve, phi_pre_curve, ...
 a_load_curve, a_geom_curve, D_installed_curve_mm, eps_pre_curve, ...
 converged_curve] = runModelOverStrainRange(eps_plot, p);

Tref = interp1(eps_plot, T_cav, 0.00, 'linear');

if abs(Tref) < 1e-12
    error('Reference transmission for normalization is too small. Inspect raw T_cav first.');
end

T_cav_rel0 = T_cav ./ Tref;

S_model = T_cav_rel0;

% Keep this alias for downstream calculations
S_interface_model = S_model;

%% -----------------------------
% Print summary at experimental points
% ------------------------------
fprintf('\n--- model summary at experimental strain values ---\n');
fprintf('   eps      data      model    T_rel0      a_load(mm)   a_geom(mm)   D(Nm)   Tpre_eff(N/m)   phi_pre   resid    w_resid\n');

resid = zeros(size(eps_data));
wres  = zeros(size(eps_data));

for i = 1:numel(eps_data)
    Si_model  = interp1(eps_plot, S_interface_model,       eps_data(i), 'linear');
    Ti_rel0   = interp1(eps_plot, T_cav_rel0,         eps_data(i), 'linear');
    ai_load   = interp1(eps_plot, 1e3 * a_load_curve, eps_data(i), 'linear');
    ai_geom   = interp1(eps_plot, 1e3 * a_geom_curve, eps_data(i), 'linear');
    Di        = interp1(eps_plot, Dplate_curve,       eps_data(i), 'linear');
    Tpre_effi = interp1(eps_plot, Tpre_eff_curve,     eps_data(i), 'linear');
    phii      = interp1(eps_plot, phi_pre_curve,      eps_data(i), 'linear');

    resid(i) = Si_model - y_data(i);
    wres(i)  = resid(i) / y_err(i);

    fprintf('%7.3f   %7.3f   %7.3f   %7.3f   %10.3f   %10.3f   %8.3e   %8.3e   %7.3f   %7.3f   %7.3f\n', ...
        eps_data(i), y_data(i), Si_model, Ti_rel0, ai_load, ai_geom, Di, Tpre_effi, phii, resid(i), wres(i));
end

rmse  = sqrt(mean(resid .^ 2));
wrmse = sqrt(mean(wres  .^ 2));

fprintf('\n RMSE  = %.4f\n', rmse);
fprintf('WRMSE = %.4f\n', wrmse);
fprintf('converged points: %d / %d\n', sum(converged_curve > 0.5), numel(converged_curve));

%% -----------------------------
% pressure-level dependence check
% ------------------------------
[pressureGain, pressureGain_relLocal, pressureGain_relZeroStrain, ...
 S_pressure, pressureConverged] = runPressureLevelSweep(eps_pressure_cases, P3_sweep, P3_offset_sweep, p);

if printPressureSweepSummary
    fprintf('\n--- pressure-level dependence check ---\n');
    fprintf('Rows are pre-strain cases. Columns are P3 offsets from P_ref.\n\n');

    fprintf('P3 offsets (Pa): ');
    fprintf('%10.1f ', P3_offset_sweep);
    fprintf('\n');

    fprintf('\nRaw local gain dP1/dP3:\n');
    for ie = 1:numel(eps_pressure_cases)
        fprintf('eps = %.3f: ', eps_pressure_cases(ie));
        fprintf('%10.5f ', pressureGain(ie, :));
        fprintf('\n');
    end

    fprintf('\nGain normalized within each pre-strain case to P3 offset = 0:\n');
    for ie = 1:numel(eps_pressure_cases)
        fprintf('eps = %.3f: ', eps_pressure_cases(ie));
        fprintf('%10.5f ', pressureGain_relLocal(ie, :));
        fprintf('\n');
    end

    fprintf('\nGain normalized by zero-strain, zero-offset gain:\n');
    for ie = 1:numel(eps_pressure_cases)
        fprintf('eps = %.3f: ', eps_pressure_cases(ie));
        fprintf('%10.5f ', pressureGain_relZeroStrain(ie, :));
        fprintf('\n');
    end

    fprintf('\nInterface-relative sensitivity prediction across pressure levels:\n');
    for ie = 1:numel(eps_pressure_cases)
        fprintf('eps = %.3f: ', eps_pressure_cases(ie));
        fprintf('%10.5f ', S_pressure(ie, :));
        fprintf('\n');
    end

    fprintf('\nPressure sweep converged points: %d / %d\n', ...
        nnz(pressureConverged(:) > 0.5), numel(pressureConverged));
end

if makePressureSweepPlot
    figure; hold on; %#ok<UNRCH>

    for ie = 1:numel(eps_pressure_cases)
        plot(P3_offset_sweep, pressureGain_relLocal(ie, :), '-o', ...
            'LineWidth', 2.2, ...
            'MarkerSize', 7, ...
            'DisplayName', sprintf('\\epsilon = %.3f', eps_pressure_cases(ie)));
    end

    yline(1.0, ':k', 'Reference at P_3 - P_{ref} = 0', 'HandleVisibility', 'off');

    xlabel('P_3 - P_{ref} (Pa)');
    ylabel('Relative Local Gain (-)');
    title('Pressure-Level Dependence at Fixed Pre-Strain');
    legend('Location', 'best', 'FontSize', legendFontSize);
    grid on;
    box on;
    formatAxes(gca, axisFontSize);

    saveCurrentFigure(resultsDir, '07_pressure_level_dependence_v20');
end

if makePressureSurfacePlot
    figure; %#ok<UNRCH>
    contourf(P3_offset_sweep, eps_pressure_cases, S_pressure, 24, 'LineColor', 'none');
    hold on;
    colormap(parula);

    xlabel('P_3 - P_{ref} (Pa)');
    ylabel('Engineering Strains, \epsilon (-)');
    title('Interface-Relative Sensitivity vs. Pre-Strain and Pressure Level');
    cb = colorbar;
    cb.Label.String = 'Interface-Relative Sensitivity (-)';
    grid on;
    box on;
    formatAxes(gca, axisFontSize);

    saveCurrentFigure(resultsDir, '08_pressure_strain_sensitivity_map_v20');
end

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

    % Baseline model curve
    hModel = plot(eps_plot, S_model, '-', ...
        'Color', modelColor, ...
        'LineWidth', modelLineWidth, ...
        'DisplayName', 'Baseline Model, h_0 = 0.508 mm');
    
    % Best physically controlled thickness case
    pBest = p;
    pBest.t_plate0 = 1.50e-3; % m
    
    [S_bestFit, ~] = runBareReferencedModelResponse(eps_plot, pBest);
    
    hBestFit = plot(eps_plot, S_bestFit, '--', ...
        'Color', physColor, ...
        'LineWidth', modelLineWidth, ...
        'DisplayName', 'Physically Tuned Case, h_0 = 1.50 mm');

    xlabel('Engineering Strains, \epsilon (-)');
    ylabel('Normalized Sensitivity, S/S_0 (-)');
    title('Measured and Modeled Normalized Sensitivity');
    xlim([0 0.36]);
    ylim([0.85 1.5]);
    grid on;
    box on;
    formatAxes(gca, axisFontSize);

    legend([h75, h80, h85, h90, hBareMean, hModel, hBestFit], ...
           {'75%', '80%', '85%', '90%', 'Bare Port Reference', ...
            'Baseline Model, h_0 = 0.508 mm', ...
            'Best Physical Fit, h_0 = 1.50 mm'}, ...
           'Location', 'northeast', 'FontSize', legendFontSize);

    saveCurrentFigure(resultsDir, '01_main_normalized_sensitivity_v20');
end

%% -----------------------------
% Interface-relative transmission
% ------------------------------
if makeTransmissionPlot
    figure; hold on;

    % Model curve
    plot(eps_plot, T_cav_rel0, '-', ...
        'Color', physColor, ...
        'LineWidth', 2.5);

    % Bare port reference (consistent style)
    yline(1.0, '-', ...
        'Color', bareColor, ...
        'LineWidth', 1.8);

    xlabel('Engineering Strains, \epsilon (-)');
    ylabel('Normalized Transmission (-)');
    title('Cavity / Interface Transmission Relative to Zero-Strain Interface');

    legend('T_{cav,rel0}(\epsilon)', 'Bare Port Reference', ...
    'Location', 'east', ...
    'FontSize', legendFontSize);

    grid on;
    box on;
    formatAxes(gca, axisFontSize);

    saveCurrentFigure(resultsDir, '02_transmission_rel0_v20');
end

%% -----------------------------
% Material sensitivity sweep: E, nu, and thickness
% ------------------------------
E_sweep      = [3e5, 6e5, 9e5, 1.2e6];                          % Pa
nu_sweep     = [0.35, 0.45, 0.49];                              % (-)
t_sweep_mm   = [0.25, 0.508, 0.75, 1.00, 1.25, 1.50, 2.00];     % mm
t_sweep_m    = t_sweep_mm * 1e-3;                               % m

if makeMaterialSweepPlot
    figure; hold on;

    for iE = 1:numel(E_sweep)
        pSweep = p;
        pSweep.E_plate0 = E_sweep(iE);

        [~, Trel_E] = runBareReferencedModelResponse(eps_plot, pSweep);

        plot(eps_plot, Trel_E, '-', ...
            'LineWidth', 2.2, ...
            'DisplayName', sprintf('E = %.1e Pa', E_sweep(iE)));
    end

    yline(1.0, '-', 'Color', bareColor, 'LineWidth', 1.8, ...
        'DisplayName', 'Bare Port Reference');

    xlabel('Engineering Strains, \epsilon (-)');
    ylabel('Normalized Transmission (-)');
    title('Effect of Young''s Modulus on Transmission');
    legend('Location', 'west', 'FontSize', legendFontSize);
    grid on; box on;
    formatAxes(gca, axisFontSize);

    saveCurrentFigure(resultsDir, '09_E_sweep_transmission_v20_fixedLoadRadius');

    figure; hold on;

    for inu = 1:numel(nu_sweep)
        pSweep = p;
        pSweep.nu_plate = nu_sweep(inu);

        [~, Trel_nu] = runBareReferencedModelResponse(eps_plot, pSweep);

        plot(eps_plot, Trel_nu, '-', ...
            'LineWidth', 2.2, ...
            'DisplayName', sprintf('\\nu = %.2f', nu_sweep(inu)));
    end

    yline(1.0, '-', 'Color', bareColor, 'LineWidth', 1.8, ...
        'DisplayName', 'Bare Port Reference');

    xlabel('Engineering Strains, \epsilon (-)');
    ylabel('Normalized Transmission (-)');
    title('Effect of Poisson''s Ratio on Transmission');
    legend('Location', 'west', 'FontSize', legendFontSize);
    grid on; box on;
    formatAxes(gca, axisFontSize);

    saveCurrentFigure(resultsDir, '10_nu_sweep_transmission_v20_noGeomAmp');

    figure; hold on;

    for it = 1:numel(t_sweep_m)
        pSweep = p;
        pSweep.t_plate0 = t_sweep_m(it);

        [~, Trel_t] = runBareReferencedModelResponse(eps_plot, pSweep);

        plot(eps_plot, Trel_t, '-', ...
            'LineWidth', 2.2, ...
            'DisplayName', sprintf('h_0 = %.3f mm', t_sweep_mm(it)));
    end

    yline(1.0, '-', 'Color', bareColor, 'LineWidth', 1.8, ...
        'DisplayName', 'Bare Port Reference');

    xlabel('Engineering Strains, \epsilon (-)');
    ylabel('Normalized Transmission (-)');
    title('Effect of Membrane Thickness on Transmission');
    legend('Location', 'southeast', 'FontSize', legendFontSize);
    grid on; box on;
    formatAxes(gca, axisFontSize);

    saveCurrentFigure(resultsDir, '11_thickness_sweep_transmission_v20_noGeomAmp');

    if makePressureLevelMaterialPlots
        % Pressure-level sweep over strain
        P3_offset_curve_sweep = [-400, -200, 0, 200, 400]; %#ok<UNRCH>
        P3_curve_sweep        = p.P_ref + P3_offset_curve_sweep;
    
        figure; hold on;
    
        for ip = 1:numel(P3_curve_sweep)
            P3_center = P3_curve_sweep(ip);
            Trel_P3 = runLocalP3GainOverStrain(eps_plot, P3_center, p);
    
            plot(eps_plot, Trel_P3, '-', ...
                'LineWidth', 2.2, ...
                'DisplayName', sprintf('P_3 - P_{ref} = %+d Pa', P3_offset_curve_sweep(ip)));
        end
    
        yline(1.0, ':k', 'Reference = 1', 'HandleVisibility', 'off');
        xlabel('Engineering Strains, \epsilon (-)');
        ylabel('Normalized Local Transmission (-)');
        title('Effect of P_3 Pressure Level on Transmission');
        legend('Location', 'best', 'FontSize', legendFontSize);
        grid on; box on;
        formatAxes(gca, axisFontSize);
    
        saveCurrentFigure(resultsDir, '12_P3_sweep_transmission_v20_fixedLoadRadius');
    
        Trel_ref = runLocalP3GainOverStrain(eps_plot, p.P_ref, p);
    
        figure; hold on;
        for ip = 1:numel(P3_curve_sweep)
            Trel_P3 = runLocalP3GainOverStrain(eps_plot, P3_curve_sweep(ip), p);
            plot(eps_plot, Trel_P3 - Trel_ref, '-', ...
                'LineWidth', 2.2, ...
                'DisplayName', sprintf('P_3 - P_{ref} = %+d Pa', P3_offset_curve_sweep(ip)));
        end
    
        yline(0.0, ':k', 'Reference = 0', 'HandleVisibility', 'off');
        xlabel('Engineering Strains, \epsilon (-)');
        ylabel('\Delta Normalized Local Transmission (-)');
        title('Pressure-Level Effect Relative to P_3 = P_{ref}');
        legend('Location', 'best', 'FontSize', legendFontSize);
        grid on; box on;
        formatAxes(gca, axisFontSize);
    
        saveCurrentFigure(resultsDir, '13_P3_sweep_delta_transmission_v20_fixedLoadRadius');
    end
end

%% -----------------------------
% Material sweep mapped onto main sensitivity fit
% ------------------------------
if makeMaterialFitOverlayPlots

    % ---- E overlay ----
    figure; hold on;
    plotExperimentalData(eps_data, y_data, y_err, c75, c80, c85, c90, ...
        markerSize, errorLineWidth);

    for iE = 1:numel(E_sweep)
        pSweep = p;
        pSweep.E_plate0 = E_sweep(iE);

        [S_E, ~] = runBareReferencedModelResponse(eps_plot, pSweep);

        plot(eps_plot, S_E, '-', 'LineWidth', 2.2, ...
            'DisplayName', sprintf('E = %.1e Pa', E_sweep(iE)));
    end

    yline(1.0, '-', 'Color', bareColor, 'LineWidth', 1.8, ...
        'DisplayName', 'Bare Port Reference');
    xlabel('Engineering Strains, \epsilon (-)');
    ylabel('Normalized Sensitivity, S/S_0 (-)');
    title('Effect of Young''s Modulus on Model Fit');
    xlim([0 0.36]); ylim([0.85 1.70]);
    legend('Location', 'northeast', 'FontSize', legendFontSize);
    grid on; box on; formatAxes(gca, axisFontSize);
    saveCurrentFigure(resultsDir, '14_E_sweep_fit_overlay_v20_fixedLoadRadius');


    % ---- nu overlay ----
    figure; hold on;
    plotExperimentalData(eps_data, y_data, y_err, c75, c80, c85, c90, ...
        markerSize, errorLineWidth);

    for inu = 1:numel(nu_sweep)
        pSweep = p;
        pSweep.nu_plate = nu_sweep(inu);

        [S_nu, ~] = runBareReferencedModelResponse(eps_plot, pSweep);

        plot(eps_plot, S_nu, '-', 'LineWidth', 2.2, ...
            'DisplayName', sprintf('\\nu = %.2f', nu_sweep(inu)));
    end

    yline(1.0, '-', 'Color', bareColor, 'LineWidth', 1.8, ...
        'DisplayName', 'Bare Port Reference');
    xlabel('Engineering Strains, \epsilon (-)');
    ylabel('Normalized Sensitivity, S/S_0 (-)');
    title('Effect of Poisson''s Ratio on Model Fit');
    xlim([0 0.36]); ylim([0.85 1.70]);
    legend('Location', 'northeast', 'FontSize', legendFontSize);
    grid on; box on; formatAxes(gca, axisFontSize);
    saveCurrentFigure(resultsDir, '15_nu_sweep_fit_overlay_v20_fixedLoadRadius');


    % ---- thickness overlay ----
    figure; hold on;
    plotExperimentalData(eps_data, y_data, y_err, c75, c80, c85, c90, ...
        markerSize, errorLineWidth);

    for it = 1:numel(t_sweep_m)
        pSweep = p;
        pSweep.t_plate0 = t_sweep_m(it);

        [S_t, ~] = runBareReferencedModelResponse(eps_plot, pSweep);

        plot(eps_plot, S_t, '-', 'LineWidth', 2.2, ...
            'DisplayName', sprintf('h_0 = %.3f mm', t_sweep_mm(it)));
    end

    yline(1.0, '-', 'Color', bareColor, 'LineWidth', 1.8, ...
        'DisplayName', 'Bare Port Reference');
    xlabel('Engineering Strains, \epsilon (-)');
    ylabel('Normalized Sensitivity, S/S_0 (-)');
    title('Effect of Membrane Thickness on Model Fit');
    xlim([0 0.36]); ylim([0.85 1.90]);
    legend('Location', 'northeast', 'FontSize', legendFontSize);
    grid on; box on; formatAxes(gca, axisFontSize);
    saveCurrentFigure(resultsDir, '16_thickness_sweep_fit_overlay_v20_fixedLoadRadius');
end

%% -----------------------------
% Mechanics diagnostic
% ------------------------------
if makeMechanicsPlot
    Dplate_norm = Dplate_curve ./ Dplate_curve(1);

    if max(abs(Tpre_eff_curve)) > 1e-12
        Tpre_eff_norm = Tpre_eff_curve ./ max(abs(Tpre_eff_curve));
    else
        Tpre_eff_norm = Tpre_eff_curve;
    end

    phi_pre_norm = phi_pre_curve;

    figure; hold on;

    plot(eps_plot, Dplate_norm, '--', ...
        'Color', physColor, ...
        'LineWidth', 2.5, ...
        'DisplayName', 'D_{plate}(\epsilon)/D_{plate}(0)');

    plot(eps_plot, Tpre_eff_norm, '-.', ...
        'Color', modelColor, ...
        'LineWidth', 2.5, ...
        'DisplayName', 'T_{pre,eff}(\epsilon)/max(T_{pre,eff})');

    plot(eps_plot, phi_pre_norm, ':', ...
        'Color', geomColor, ...
        'LineWidth', 3.0, ...
        'DisplayName', '\phi_{pre}(\epsilon)');

    xlabel('Engineering Strains, \epsilon (-)');
    ylabel('Normalized Mechanical Quantity (-)');
    title('Normalized Mechanics Diagnostics');
    xlim([0 0.36]);
    ylim([0 1.10]);
    legend('Location', 'best', 'FontSize', legendFontSize);
    grid on;
    box on;
    formatAxes(gca, axisFontSize);

    saveCurrentFigure(resultsDir, '03_normalized_mechanics_v20');
end

%% -----------------------------
% Geometry / loaded-radius diagnostic
% ------------------------------
if makeGeometryPlot
    figure; hold on; %#ok<UNRCH>

    plot(eps_plot, 1e3 * a_geom_curve, '-', 'Color', geomColor, ...
        'LineWidth', 2.5, 'DisplayName', 'a_{installed}(\epsilon)');

    plot(eps_plot, 1e3 * a_load_curve, '--', 'Color', modelColor, ...
        'LineWidth', 2.5, 'DisplayName', 'a_{load}(\epsilon)');

    yline(1e3 * p.r_forced_m, ':', 'Color', bareColor, ...
        'LineWidth', 1.8, 'DisplayName', 'r_{forced}');

    xlabel('Engineering Strains, \epsilon (-)');
    ylabel('Radius (mm)');
    title('Geometry Diagnostics: Amplification Off');
    legend('Location', 'best');
    grid on;
    formatAxes(gca, axisFontSize);

    saveCurrentFigure(resultsDir, '04_geometry_loaded_radius_v20_noGeomAmp');
end

%% -----------------------------
% Sweep: kT0 only, with S0 fixed at 1.0
% ------------------------------
rmse  = zeros(numel(S0_grid), numel(kT0_grid));
wrmse = zeros(numel(S0_grid), numel(kT0_grid));

for ig = 1:numel(S0_grid)
    for ik = 1:numel(kT0_grid)
        pSweep = p;
        pSweep.S0_interface_bare = S0_grid(ig);
        pSweep.kT0               = kT0_grid(ik);

        [S_sweep, ~] = runBareReferencedModelResponse(eps_plot, pSweep);
        [rmse(ig, ik), wrmse(ig, ik)] = ...
            computeModelErrors(eps_plot, S_sweep, eps_data, y_data, y_err);
    end
end

[minRMSE, idxRMSE] = min(rmse(:));
[rowRMSE, colRMSE] = ind2sub(size(rmse), idxRMSE);
bestS0_RMSE   = S0_grid(rowRMSE);
bestkT0_RMSE  = kT0_grid(colRMSE);

[minWRMSE, idxWRMSE] = min(wrmse(:));
[rowWRMSE, colWRMSE] = ind2sub(size(wrmse), idxWRMSE);
bestS0_WRMSE  = S0_grid(rowWRMSE);
bestkT0_WRMSE = kT0_grid(colWRMSE);

if makeRMSEPlot
    figure; hold on; %#ok<UNRCH>

    plot(kT0_grid, rmse(1, :), '-o', ...
        'LineWidth', 2.2, ...
        'MarkerSize', 7, ...
        'DisplayName', 'S_0 = 1.0');

    plot(bestkT0_RMSE, minRMSE, 'rp', ...
        'MarkerFaceColor', 'r', ...
        'MarkerEdgeColor', 'k', ...
        'MarkerSize', 12, ...
        'DisplayName', 'Min RMSE');

    xline(p.kT0, '--k', ...
        'DisplayName', sprintf('Baseline k_{T0} = %.3f', p.kT0));

    xlabel('k_{T0}');
    ylabel('RMSE (-)');
    title('Pre-Tension Sensitivity Sweep: RMSE');
    grid on;
    box on;
    formatAxes(gca, axisFontSize);
    legend('Location', 'best', 'FontSize', legendFontSize);

    saveCurrentFigure(resultsDir, '05_kT0_sweep_rmse_v20');
end

if makeWRMSEPlot
    figure; hold on; %#ok<UNRCH>

    plot(kT0_grid, wrmse(1, :), '-o', ...
        'LineWidth', 2.2, ...
        'MarkerSize', 7, ...
        'DisplayName', 'S_0 = 1.0');

    plot(bestkT0_WRMSE, minWRMSE, 'rp', ...
        'MarkerFaceColor', 'r', ...
        'MarkerEdgeColor', 'k', ...
        'MarkerSize', 12, ...
        'DisplayName', 'Min WRMSE');

    xline(p.kT0, '--k', ...
        'DisplayName', sprintf('Baseline k_{T0} = %.3f', p.kT0));

    xlabel('k_{T0}');
    ylabel('WRMSE (-)');
    title('Pre-Tension Sensitivity Sweep: WRMSE');
    grid on;
    box on;
    formatAxes(gca, axisFontSize);
    legend('Location', 'best', 'FontSize', legendFontSize);

    saveCurrentFigure(resultsDir, '06_kT0_sweep_wrmse_v20');
end

%% -----------------------------
% Focused 1D sweeps for console output
% ------------------------------
kT0_rmse  = zeros(numel(kT0_sweep), 1);
kT0_wrmse = zeros(numel(kT0_sweep), 1);

for j = 1:numel(kT0_sweep)
    pSweep = p;
    pSweep.kT0 = kT0_sweep(j);

    [S_sweep, ~] = runBareReferencedModelResponse(eps_plot, pSweep);
    [kT0_rmse(j), kT0_wrmse(j)] = computeModelErrors(eps_plot, S_sweep, eps_data, y_data, y_err);
end

S0_rmse  = zeros(numel(S0_sweep), 1);
S0_wrmse = zeros(numel(S0_sweep), 1);

for j = 1:numel(S0_sweep)
    pSweep = p;
    pSweep.S0_interface_bare = S0_sweep(j);

    [S_sweep, ~] = runBareReferencedModelResponse(eps_plot, pSweep);
    [S0_rmse(j), S0_wrmse(j)] = computeModelErrors(eps_plot, S_sweep, eps_data, y_data, y_err);
end

%% -----------------------------
% Print sensitivity summary
% ------------------------------
if printSweepSummary
    fprintf('\n--- sensitivity sweep summary ---\n');

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

    fprintf('\nkT0 sweep with S0 fixed at 1.0:\n');
    fprintf('  kT0 range: %.3f to %.3f (%d cases)\n', ...
        kT0_grid(1), kT0_grid(end), numel(kT0_grid));
    fprintf('  S0_interface_bare range: %.3f to %.3f (%d cases)\n', ...
        S0_grid(1), S0_grid(end), numel(S0_grid));

    fprintf('  Minimum RMSE  = %.4f at kT0 = %.3f, S0_interface_bare = %.3f\n', ...
        minRMSE, bestkT0_RMSE, bestS0_RMSE);
    fprintf('  Minimum WRMSE = %.4f at kT0 = %.3f, S0_interface_bare = %.3f\n', ...
        minWRMSE, bestkT0_WRMSE, bestS0_WRMSE);
end

%% -----------------------------
% Local functions
% ------------------------------
function [T_cav, Dplate_curve, E_curve, Tpre_eff_curve, phi_pre_curve, ...
          a_load_curve, a_geom_curve, D_installed_curve_mm, eps_pre_curve, ...
          converged_curve] = runModelOverStrainRange(eps_plot, p)

    T_cav                = zeros(size(eps_plot));
    Dplate_curve         = zeros(size(eps_plot));
    E_curve           = zeros(size(eps_plot));
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

        [T_cav(k), Dplate_curve(k), E_curve(k), Tpre_eff_curve(k), ...
         phi_pre_curve(k), a_load_curve(k), a_geom_curve(k), ...
         D_installed_curve_mm(k), eps_pre_curve(k), converged_curve(k), ...
         solverState] = evaluateCavityTransmission(eps_query, p, solverState);
    end
end

function [S_bare_model, T_cav_rel0] = runBareReferencedModelResponse(eps_plot, p)
    [T_cav, ~, ~, ~, ~, ~, ~, ~, ~, ~] = runModelOverStrainRange(eps_plot, p);

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

function [Tcav, Dplate, E, Tpre_eff, phi_pre, a_load, a_geom, ...
          D_installed_mm, eps_pre, converged, solverStateOut] = ...
          evaluateCavityTransmission(eps_query, p, solverStateIn)

    [Dplate, E, Tpre_eff, phi_pre, a_load, a_geom, D_installed_mm, eps_pre] = ...
        evaluateMembraneState(eps_query, p);

    symmetricState = solverStateIn;
    if isempty(symmetricState)
        symmetricState.Pi = p.P_ref;
    end

    dPsens_plus  = sensorDifferentialResponse(+p.dP0, a_load, Dplate, Tpre_eff, p, symmetricState);
    dPsens_minus = sensorDifferentialResponse(-p.dP0, a_load, Dplate, Tpre_eff, p, symmetricState);

    Tcav = (dPsens_plus - dPsens_minus) / (2 * p.dP0);

    [~, ~, ~, ~, ~, converged, solverStateOut] = ...
        solveCavityPressure_iterative(p.P_ref + p.dP0, a_load, Dplate, Tpre_eff, p, solverStateIn);
end

function [eps_pre, lambda_pre, D_installed_mm] = evaluateInstalledPrestretchState(eps_query, p)
    D_installed_mm = interp1(p.engStrains, p.scaleDiameters_mm, eps_query, 'linear', 'extrap');
    eps_pre = max(eps_query, 0);
    lambda_pre = 1 + eps_pre;
end

function [a_load, a_geom, D_installed_mm, eps_pre, lambda_pre, t_plate] = ...
    evaluatePlateGeometry(eps_query, p)

    [eps_pre, lambda_pre, D_installed_mm] = evaluateInstalledPrestretchState(eps_query, p);

    a_geom = 0.5 * D_installed_mm * 1e-3;
    
    % Effective loaded radius used in the membrane deflection calculation.
    % Held fixed to isolate material, pre-tension, and cavity effects.
    a_load = p.r_forced_m;

    t_plate = p.t_plate0 / (lambda_pre ^ 2);
end

function [Dplate, E, Tpre_eff, phi_pre, a_load, a_geom, D_installed_mm, eps_pre] = ...
    evaluateMembraneState(eps_query, p)

    [a_load, a_geom, D_installed_mm, eps_pre, ~, t_plate] = ...
        evaluatePlateGeometry(eps_query, p);

    if p.useStrainDependentE
        E = p.E_plate0 * (1 + p.c1 * eps_pre + p.c2 * eps_pre ^ 2);
    else
        E = p.E_plate0;
    end

    Dplate = E * t_plate ^ 3 / (12 * (1 - p.nu_plate ^ 2));

    phi_pre = 1 - exp(-(eps_pre / max(p.eps_char, 1e-12)) ^ 2);
    Tpre_eff = p.kT0 * phi_pre * E * t_plate * eps_pre / (1 - p.nu_plate);
end

function dPsens = sensorDifferentialResponse(dPsurf, a_plate, Dplate, Tpre_eff, p, solverStateIn)
    P3 = p.P_ref + dPsurf / 2;
    P4 = p.P_ref - dPsurf / 2;

    [P1, ~, ~, ~, ~, ~, ~] = solveCavityPressure_iterative(P3, a_plate, Dplate, Tpre_eff, p, solverStateIn);
    [P2, ~, ~, ~, ~, ~, ~] = solveCavityPressure_iterative(P4, a_plate, Dplate, Tpre_eff, p, solverStateIn);

    dPsens = P1 - P2;
end

function [pressureGain, pressureGain_relLocal, pressureGain_relZeroStrain, ...
          S_pressure, pressureConverged] = ...
          runPressureLevelSweep(eps_pressure_cases, P3_sweep, P3_offset_sweep, p)

    pressureGain           = zeros(numel(eps_pressure_cases), numel(P3_sweep));
    pressureGain_relLocal  = zeros(size(pressureGain));
    pressureConverged      = zeros(size(pressureGain));

    for ie = 1:numel(eps_pressure_cases)
        eps_i = eps_pressure_cases(ie);

        [Dplate_i, ~, Tpre_eff_i, ~, a_load_i, ~, ~, ~] = ...
            evaluateMembraneState(eps_i, p);

        for ip = 1:numel(P3_sweep)
            P3_i = P3_sweep(ip);

            [pressureGain(ie, ip), pressureConverged(ie, ip)] = ...
                evaluateLocalP3Gain(P3_i, a_load_i, Dplate_i, Tpre_eff_i, p);
        end

        zeroOffsetIdx = find(abs(P3_offset_sweep) < 1e-12, 1);
        if isempty(zeroOffsetIdx)
            [~, zeroOffsetIdx] = min(abs(P3_offset_sweep));
        end

        localRef = pressureGain(ie, zeroOffsetIdx);
        if abs(localRef) < 1e-12
            warning('Local pressure-gain reference is near zero at eps = %.3f.', eps_i);
            localRef = NaN;
        end

        pressureGain_relLocal(ie, :) = pressureGain(ie, :) ./ localRef;
    end

    [Dplate_0, ~, Tpre_eff_0, ~, a_load_0, ~, ~, ~] = evaluateMembraneState(0.0, p);
    [G0_pressure, ~] = evaluateLocalP3Gain(p.P_ref, a_load_0, Dplate_0, Tpre_eff_0, p);

    if abs(G0_pressure) < 1e-12
        error('Zero-strain pressure-gain reference is too small. Inspect pressureGain before normalization.');
    end

    pressureGain_relZeroStrain = pressureGain ./ G0_pressure;
    S_pressure = p.S0_interface_bare * pressureGain_relZeroStrain;
end

function [Glocal, convergedBoth] = evaluateLocalP3Gain(P3_center, a_plate, Dplate, Tpre_eff, p)
    dP = p.dP_local;

    solverState.Pi = p.P_ref;

    P3_plus  = P3_center + dP / 2;
    P3_minus = P3_center - dP / 2;

    [P1_plus, ~, ~, ~, ~, convPlus, ~] = ...
        solveCavityPressure_iterative(P3_plus, a_plate, Dplate, Tpre_eff, p, solverState);

    [P1_minus, ~, ~, ~, ~, convMinus, ~] = ...
        solveCavityPressure_iterative(P3_minus, a_plate, Dplate, Tpre_eff, p, solverState);

    Glocal = (P1_plus - P1_minus) / dP;
    convergedBoth = convPlus && convMinus;
end

function [Pi, Vi, dV, wmax, iter, converged, solverStateOut] = ...
    solveCavityPressure_iterative(Psurf, a_plate, Dplate, Tpre_eff, p, solverStateIn)

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
            loadedCavityState_fromMembraneLoad(q_membrane, p.V00, a_plate, Dplate, Tpre_eff, p);

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

function [Vi, dV, wmax] = loadedCavityState_fromMembraneLoad(q_membrane, V00, a, Dplate, Tpre_eff, p)
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

function Trel_P3 = runLocalP3GainOverStrain(eps_plot, P3_center, p) %#ok<DEFNU>

    Glocal_curve = zeros(size(eps_plot));

    for k = 1:numel(eps_plot)
        eps_query = eps_plot(k);

        [Dplate, ~, Tpre_eff, ~, a_load, ~, ~, ~] = ...
            evaluateMembraneState(eps_query, p);

        [Glocal_curve(k), ~] = ...
            evaluateLocalP3Gain(P3_center, a_load, Dplate, Tpre_eff, p);
    end

    Gref = interp1(eps_plot, Glocal_curve, 0.00, 'linear');

    if abs(Gref) < 1e-12
        error('Reference local P3 gain is too small for normalization.');
    end

    Trel_P3 = Glocal_curve ./ Gref;
end

function plotExperimentalData(eps_data, y_data, y_err, c75, c80, c85, c90, markerSize, errorLineWidth)

    errorbar(eps_data(4), y_data(4), y_err(4), 's', ...
        'Color', c75, 'MarkerFaceColor', c75, 'MarkerEdgeColor', c75, ...
        'MarkerSize', markerSize, 'LineWidth', errorLineWidth, ...
        'CapSize', 8, 'DisplayName', '75%');

    errorbar(eps_data(3), y_data(3), y_err(3), 's', ...
        'Color', c80, 'MarkerFaceColor', c80, 'MarkerEdgeColor', c80, ...
        'MarkerSize', markerSize, 'LineWidth', errorLineWidth, ...
        'CapSize', 8, 'DisplayName', '80%');

    errorbar(eps_data(2), y_data(2), y_err(2), 's', ...
        'Color', c85, 'MarkerFaceColor', c85, 'MarkerEdgeColor', c85, ...
        'MarkerSize', markerSize, 'LineWidth', errorLineWidth, ...
        'CapSize', 8, 'DisplayName', '85%');

    errorbar(eps_data(1), y_data(1), y_err(1), 's', ...
        'Color', c90, 'MarkerFaceColor', c90, 'MarkerEdgeColor', c90, ...
        'MarkerSize', markerSize, 'LineWidth', errorLineWidth, ...
        'CapSize', 8, 'DisplayName', '90%');
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
    print(fig, fullfile(resultsDir, [baseName, '.eps']), '-depsc', '-vector');
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