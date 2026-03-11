%% sensor_interface_numerical_model_v0.m
% Reduced-order model for strain-dependent pressure transmission
% Tyler Inkley - updated first-pass journal model scaffold
%
% Notes:
% - S_norm(eps) = G(eps) * T(eps)
% - T(eps) comes from cavity compressibility + interface compliance
% - G(eps) is an effective strain-dependent gain termSho

clear; clc; close all;

%% -----------------------------
% Measured strain values from calibration script
% ------------------------------
d0              = 21.978;                                   % mm, unstretched reference diameter
scaleDiameters  = [16.484, 17.582, 18.681, 19.780, 21.978]; % mm
engStrains      = (d0 - scaleDiameters) ./ d0;
engStrains      = sort(engStrains);                         % [0, 0.111..., 0.176..., 0.25, 0.333...]

% Experimental summary points for direct comparison
eps_data        = [0.100, 0.150, 0.200, 0.250];
y_data          = [1.167, 1.328, 1.270, 1.224];

% Approximate error bars from champion figure / summary plot
y_err           = [0.093, 0.082, 0.045, 0.042];

%% -----------------------------
% Plot styling
% ------------------------------
modelColor      = [0.00, 0.20, 0.60];   % dark blue
bareColor       = [0.35, 0.35, 0.35];   % neutral gray
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
P_ref   = 101325;    % Pa, ambient pressure

% Geometry / cavity lumped parameter
V00     = 2.0e-7;    % m^3, nominal cavity volume per side

% Effective interface compliance model
Cm0     = 1.5e-11;   % m^3/Pa
a1      = 2.0;       % linear strain stiffening term
a2      = 4.0;       % quadratic strain stiffening term

% Optional cavity volume reduction with strain
alphaV  = 0.5;       % V0 = V00*(1 - alphaV*eps)

% Effective interface gain term
A_gain      = 0.58;  % peak gain above baseline of 1
eps_peak_G  = 0.15;  % gain peaks near 85% condition

% Small-signal surface pressure sweep
dPsurf_vec = linspace(-500, 500, 101);   % Pa

%% -----------------------------
% Run model
% ------------------------------
eps_plot        = linspace(0, 0.36, 200);
S_norm          = zeros(size(eps_plot));
T_only          = zeros(size(eps_plot));
G_curve         = zeros(size(eps_plot));
lambda_curve    = zeros(size(eps_plot));
V0_curve        = zeros(size(eps_plot));
Cm_curve        = zeros(size(eps_plot));

for k = 1:length(eps_plot)
    epsn = eps_plot(k);

    % Strain-dependent cavity volume
    V0 = V00 * (1 - alphaV * epsn);
    V0 = max(V0, 1e-12);   % guard against negative/singular volume

    % Strain-dependent interface compliance
    Cm = Cm0 / (1 + a1*epsn + a2*epsn^2);

    % Effective gain term with explicit peak strain
    if eps_peak_G > 0
        G = 1 + A_gain * (epsn / eps_peak_G) * exp(1 - epsn / eps_peak_G);
    else
        G = 1;
    end

    % Useful diagnostic lumped parameter
    lambda = P_ref * Cm / V0;

    dPsens_raw = zeros(size(dPsurf_vec));
    dPsens_eff = zeros(size(dPsurf_vec));

    for j = 1:length(dPsurf_vec)
        dPsurf = dPsurf_vec(j);

        % Symmetric applied differential loading
        P3 = P_ref + dPsurf/2;
        P4 = P_ref - dPsurf/2;

        % Internal cavity pressures
        P1 = solve_cavity_pressure_air(P3, P_ref, V0, Cm);
        P2 = solve_cavity_pressure_air(P4, P_ref, V0, Cm);

        % Raw sensor differential pressure
        dPsens = P1 - P2;

        % Apply effective gain term
        dPsens_raw(j) = dPsens;
        dPsens_eff(j) = G * dPsens;
    end

    % Small-signal slopes
    pfit_raw = polyfit(dPsurf_vec, dPsens_raw, 1);
    pfit_eff = polyfit(dPsurf_vec, dPsens_eff, 1);

    T_eff_raw   = pfit_raw(1);   % cavity/compliance transmission only
    T_eff_total = pfit_eff(1);   % full normalized model output

    % Store results
    T_only(k)       = T_eff_raw;
    S_norm(k)       = T_eff_total;
    G_curve(k)      = G;
    lambda_curve(k) = lambda;
    V0_curve(k)     = V0;
    Cm_curve(k)     = Cm;
end

%% -----------------------------
% Print quick summary at data points
% ------------------------------
fprintf('\n--- Model summary at experimental strain values ---\n');
fprintf('   eps        S_model      T_only       G(eps)\n');

for i = 1:length(eps_data)
    Si = interp1(eps_plot, S_norm, eps_data(i), 'linear');
    Ti = interp1(eps_plot, T_only, eps_data(i), 'linear');
    Gi = interp1(eps_plot, G_curve, eps_data(i), 'linear');
    fprintf('%8.4f    %8.4f    %8.4f    %8.4f\n', eps_data(i), Si, Ti, Gi);
end

%% -----------------------------
% Publication-style model figure
% ------------------------------
figure; hold on;

% Bare-port reference and bounds
hBareMean = yline(1.0, '-', 'Color', bareColor, 'LineWidth', 1.8, ...
    'DisplayName', 'Bare Port Reference');
yline(1.04, '--', 'Color', bareColor, 'LineWidth', 1.2, 'HandleVisibility', 'off');
yline(0.96, '--', 'Color', bareColor, 'LineWidth', 1.2, 'HandleVisibility', 'off');

% Experimental points with error bars, plotted individually for legend control
h75 = errorbar(eps_data(4), y_data(4), y_err(4), 's', ...
    'Color', c75, 'MarkerFaceColor', c75, 'MarkerEdgeColor', c75, ...
    'MarkerSize', markerSize, 'LineWidth', errorLineWidth, ...
    'CapSize', 8, 'DisplayName', '75%');

h80 = errorbar(eps_data(3), y_data(3), y_err(3), 's', ...
    'Color', c80, 'MarkerFaceColor', c80, 'MarkerEdgeColor', c80, ...
    'MarkerSize', markerSize, 'LineWidth', errorLineWidth, ...
    'CapSize', 8, 'DisplayName', '80%');

h85 = errorbar(eps_data(2), y_data(2), y_err(2), 's', ...
    'Color', c85, 'MarkerFaceColor', c85, 'MarkerEdgeColor', c85, ...
    'MarkerSize', markerSize, 'LineWidth', errorLineWidth, ...
    'CapSize', 8, 'DisplayName', '85%');

h90 = errorbar(eps_data(1), y_data(1), y_err(1), 's', ...
    'Color', c90, 'MarkerFaceColor', c90, 'MarkerEdgeColor', c90, ...
    'MarkerSize', markerSize, 'LineWidth', errorLineWidth, ...
    'CapSize', 8, 'DisplayName', '90%');

% Reduced-order model curve
hModel = plot(eps_plot, S_norm, '-', ...
    'Color', modelColor, 'LineWidth', modelLineWidth, ...
    'DisplayName', 'Reduced-Order Model');

xlabel('Engineering Strain');
ylabel('Normalized Sensitivity (-)');
title('Normalized Sensitivity vs. Engineering Strain');

% No title for journal figure
xlim([0 0.36]);
ylim([0.85 1.45]);
grid on;
box on;
set(gca, 'FontSize', axisFontSize, 'LineWidth', 1.0);

legend([h75, h80, h85, h90, hBareMean, hModel], ...
       {'75%', '80%', '85%', '90%', 'Bare Port Reference', 'Reduced-Order Model'}, ...
       'Location', 'northeast', 'FontSize', legendFontSize);

%% -----------------------------
% Diagnostic plot: decomposition
% ------------------------------
figure;
plot(eps_plot, T_only, 'b-', 'LineWidth', 2); hold on;
plot(eps_plot, G_curve, 'r--', 'LineWidth', 2);
plot(eps_plot, S_norm, '-', 'Color', modelColor, 'LineWidth', 2.5);
yline(1.0, ':k');
xlabel('Engineering Strain');
ylabel('Magnitude (-)');
title('Model Decomposition');
legend('T(\epsilon): cavity/compliance only', ...
       'G(\epsilon): effective gain term', ...
       'S_{norm}(\epsilon) = G(\epsilon)T(\epsilon)', ...
       'Reference = 1', ...
       'Location', 'best');
grid on;
set(gca, 'FontSize', 14);

%% -----------------------------
% Diagnostic plot: lambda, V0, Cm
% ------------------------------
figure;
yyaxis left
plot(eps_plot, lambda_curve, 'k-', 'LineWidth', 2); hold on;
ylabel('\lambda(\epsilon) = P_{ref} C_m / V_0');

yyaxis right
plot(eps_plot, V0_curve, 'g--', 'LineWidth', 2);
plot(eps_plot, Cm_curve, 'c-.', 'LineWidth', 2);
ylabel('V_0(\epsilon), C_m(\epsilon)');

xlabel('Engineering Strain');
title('Diagnostic Parameters vs Engineering Strain');
legend('\lambda(\epsilon)', 'V_0(\epsilon)', 'C_m(\epsilon)', 'Location', 'best');
grid on;
set(gca, 'FontSize', 14);

%% -----------------------------
% Optional: overlay on saved .fig
% ------------------------------
figFile = 'NSES_long.fig';
if isfile(figFile)
    openfig(figFile, 'new', 'visible');
    ax = gca;
    hold(ax, 'on');

    % Add reduced-order model to existing experimental figure
    plot(ax, eps_plot, S_norm, '-', ...
        'Color', modelColor, 'LineWidth', modelLineWidth, ...
        'DisplayName', 'Reduced-Order Model');

    % Optional cleanup for the imported figure
    xlim(ax, [0 0.36]);
    ylim(ax, [0.85 1.50]);
    set(ax, 'FontSize', axisFontSize, 'LineWidth', 1.0);

    lgd = legend(ax, 'show');
    set(lgd, 'FontSize', legendFontSize, 'Location', 'northeast');
end

%% -----------------------------
% Local function
% ------------------------------
function Pi = solve_cavity_pressure_air(Pext, Pref, V0, Cm)
    % Isothermal cavity model:
    %
    % Pi * [V0 - Cm*(Pext - Pi)] = Pref * V0
    %
    % Rearranged:
    % Cm*Pi^2 + (V0 - Cm*Pext)*Pi - Pref*V0 = 0

    A = Cm;
    B = V0 - Cm*Pext;
    C = -Pref * V0;

    r = roots([A B C]);

    % Keep positive real roots only
    r = r(abs(imag(r)) < 1e-10);
    r = real(r);
    r = r(r > 0);

    if isempty(r)
        Pi = Pref;
    else
        [~, idx] = min(abs(r - Pref));
        Pi = r(idx);
    end
end