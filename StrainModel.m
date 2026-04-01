%% sensor_interface_numerical_model_v1.m
% Reduced-order model for strain-dependent pressure transmission
% Improved structure with explicit normalization and local sensitivity
%
% S_norm(eps) = T_cav,norm(eps) * eta_cpl(eps)
% where:
%   T_cav,norm(eps) = cavity/compressibility transmission normalized by its
%                     zero-strain value
%   eta_cpl(eps)    = strain-dependent interface coupling efficiency
%
% Interpretation:
% - T_cav,norm captures trapped-air compression and effective cavity loading
% - eta_cpl captures improved pressure transmission with moderate pre-strain,
%   followed by tapering as interface stiffening becomes dominant
%
% Tyler Inkley - revised scaffold

clear; clc; close all;

%% -----------------------------
% Experimental data
% ------------------------------
d0              = 21.978;                                   % mm
scaleDiameters  = [16.484, 17.582, 18.681, 19.780, 21.978]; % mm
engStrains      = sort((d0 - scaleDiameters) ./ d0);

eps_data        = [0.100, 0.150, 0.200, 0.250];
y_data          = [1.167, 1.328, 1.270, 1.224];
y_err           = [0.093, 0.082, 0.045, 0.042];

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

% Reference conditions
p.P_ref = 101325;   % Pa, ambient/reference pressure

% Cavity geometry (per sensing side)
% All geometry defined in mm, volumes computed in mm^3 and converted to m^3

h1 = 2.506;   % mm, top funnel cylindrical section height
h2 = 6.617;   % mm, bottom funnel spherical-cap height
r  = 5.5;     % mm, sphere radius (used in spherical-cap formulation)

% Volume components (mm^3)
p.V11 = pi * r^2 * h1;                 % mm^3, cylindrical volume
p.V22 = pi * h2^2 * r - pi * h2^3/3;   % mm^3, spherical-cap volume

% Total cavity volume (per side)
V_total_mm3 = p.V11 + p.V22;           % mm^3

% Convert to SI units for model use
p.V00 = V_total_mm3 * 1e-9;            % m^3 (1 mm^3 = 1e-9 m^3)

%% -----------------------------
% Compliance model selection
% ------------------------------
% Options:
%   'lumped'   -> purely phenomenological effective compliance model
%   'material' -> material-informed compliance surrogate using E and t
p.complianceModel = 'material';

%% -----------------------------
% Lumped compliance model parameters
% ------------------------------
% Cm represents an effective lumped interface compliance, intended to
% capture the combined effects of membrane material stiffness, thickness,
% geometry, and strain-dependent tensioning.
p.Cm0 = 1.5e-11;   % m^3/Pa, nominal interface compliance at zero strain
p.a1  = 1.2;       % (-), linear strain-dependent stiffening coefficient
p.a2  = 2.2;       % (-), quadratic strain-dependent stiffening coefficient

%% -----------------------------
% Explicit membrane material parameters
% ------------------------------
p.interfaceMaterial = 'Hygienic latex';
p.t_mem             = 0.020 * 0.0254;   % m, membrane thickness (0.020 in)
p.E_mem             = 1.2e6;            % Pa, nominal effective elastic modulus
p.nu_mem            = 0.49;             % (-), approximate Poisson ratio for latex-like material
p.a_mem             = r * 1e-3;         % m, effective membrane radius/span from funnel radius

% Material-model strain stiffening modifiers
p.b1 = 0.8;        % (-), linear strain-dependent stiffening coefficient
p.b2 = 1.3;        % (-), quadratic strain-dependent stiffening coefficient

% Optional scale factor for material-based compliance
% This allows the material-informed model to be tuned without changing E or t
p.k_cm = 0.0075;   % (-), dimensionless scaling factor

fprintf('Compliance model: %s\n', p.complianceModel);
fprintf('Coupling term enabled: %d\n', p.useCouplingTerm);
fprintf('Nominal cavity volume per side: %.3e m^3\n', p.V00);
fprintf('Membrane thickness: %.3e m\n', p.t_mem);
fprintf('Membrane modulus: %.3e Pa\n\n', p.E_mem);

%% -----------------------------
% Optional strain-dependent cavity volume reduction
% ------------------------------
p.alphaV = 0.06;   % (-), fractional reduction coefficient: V0 = V00*(1 - alphaV*eps)

%% -----------------------------
% Interface coupling model
% ------------------------------
% eta_cpl represents strain-dependent pressure-coupling efficiency of the
% stretched interface. The intent is to capture the idea that moderate
% pre-strain can improve pressure transmission, while excessive pre-strain
% leads to tapering of that benefit.
p.useCouplingTerm = true;   % enable/disable strain-dependent coupling term
p.A_cpl           = 0.39;   % (-), peak coupling enhancement above unity
p.eps_peak_cpl    = 0.18;   % (-), strain at which coupling term peaks

%% -----------------------------
% Numerical settings
% ------------------------------
p.dP0 = 1.0;   % Pa, small-signal perturbation for centered derivative

%% -----------------------------
% Evaluate model over strain range
% ------------------------------
eps_plot        = linspace(0, 0.36, 300);

T_cav          = zeros(size(eps_plot));   % cavity transmission sensitivity
eta_cpl        = zeros(size(eps_plot));   % interface coupling efficiency

lambda_curve   = zeros(size(eps_plot));
V0_curve       = zeros(size(eps_plot));
Cm_curve       = zeros(size(eps_plot));

for k = 1:numel(eps_plot)
    epsn = eps_plot(k);

    [T_cav(k), V0_curve(k), Cm_curve(k), lambda_curve(k)] = evaluateCavityTransmission(epsn, p);
    eta_cpl(k) = evaluateCouplingEfficiency(epsn, p);
end

% Normalize cavity-only transmission relative to zero-strain value
T0 = T_cav(1);
T_cav_norm = T_cav ./ T0;

% Full physics-informed reduced-order model
S_norm = T_cav_norm .* eta_cpl;

%% -----------------------------
% Print summary at experimental points
% ------------------------------
fprintf('\n--- Model summary at experimental strain values ---\n');
fprintf('   eps      data      model    T_cav,norm   eta_cpl   resid    w_resid\n');

resid = zeros(size(eps_data));
wres  = zeros(size(eps_data));

for i = 1:numel(eps_data)
    Si_model = interp1(eps_plot, S_norm,      eps_data(i), 'linear');
    Ti_cav   = interp1(eps_plot, T_cav_norm,  eps_data(i), 'linear');
    Ei_cpl   = interp1(eps_plot, eta_cpl,     eps_data(i), 'linear');

    resid(i) = Si_model - y_data(i);
    wres(i)  = resid(i) / y_err(i);

    fprintf('%7.3f   %7.3f   %7.3f   %7.3f   %7.3f   %7.3f   %7.3f\n', ...
        eps_data(i), y_data(i), Si_model, Ti_cav, Ei_cpl, resid(i), wres(i));
end

rmse = sqrt(mean(resid.^2));
wrmse = sqrt(mean(wres.^2));

fprintf('\nRMSE  = %.4f\n', rmse);
fprintf('WRMSE = %.4f\n', wrmse);

%% -----------------------------
% Publication-style comparison figure
% ------------------------------
figure; hold on;

hBareMean = yline(1.0, '-', 'Color', bareColor, 'LineWidth', 1.8, ...
    'DisplayName', 'Bare Port Reference');
yline(1.04, '--', 'Color', bareColor, 'LineWidth', 1.2, 'HandleVisibility', 'off');
yline(0.96, '--', 'Color', bareColor, 'LineWidth', 1.2, 'HandleVisibility', 'off');

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

hPhys = plot(eps_plot, T_cav_norm, '--', ...
    'Color', physColor, 'LineWidth', 2.2, ...
    'DisplayName', 'Cavity Transmission Only');

hModel = plot(eps_plot, S_norm, '-', ...
    'Color', modelColor, 'LineWidth', modelLineWidth, ...
    'DisplayName', 'Physics-Informed Reduced-Order Model');

xlabel('Engineering Strain');
ylabel('Normalized Sensitivity (-)');
xlim([0 0.36]);
ylim([0.85 1.45]);
grid on;
box on;
set(gca, 'FontSize', axisFontSize, 'LineWidth', 1.0);

legend([h75, h80, h85, h90, hBareMean, hPhys, hModel], ...
       {'75%', '80%', '85%', '90%', 'Bare Port Reference', ...
        'Cavity Transmission Only', 'Physics-Informed Reduced-Order Model'}, ...
       'Location', 'northeast', 'FontSize', legendFontSize);

%% -----------------------------
% Diagnostic decomposition
% ------------------------------
figure;
plot(eps_plot, T_cav_norm, 'LineStyle', '--', 'Color', physColor, 'LineWidth', 2.2); hold on;
plot(eps_plot, eta_cpl, 'r--', 'LineWidth', 2.0);
plot(eps_plot, S_norm, '-', 'Color', modelColor, 'LineWidth', 2.8);
yline(1.0, ':k');

xlabel('Engineering Strain');
ylabel('Magnitude (-)');
title('Model Decomposition');
legend('T_{cav,norm}(\epsilon)', ...
       '\eta_{cpl}(\epsilon)', ...
       'S_{norm}(\epsilon)', ...
       'Reference = 1', ...
       'Location', 'best');
grid on;
set(gca, 'FontSize', 14);

%% -----------------------------
% Diagnostic parameter plot
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
% Local functions
% ------------------------------
function [Tcav, V0, Cm, lambda] = evaluateCavityTransmission(epsn, p)
    % Strain-dependent cavity volume
    V0 = p.V00 * (1 - p.alphaV * epsn);
    V0 = max(V0, 1e-12);

    % Strain-dependent effective interface compliance
    Cm = evaluateInterfaceCompliance(epsn, p);

    % Lumped parameter
    lambda = p.P_ref * Cm / V0;

    % Small-signal centered derivative
    dP0 = p.dP0;

    dPsens_plus  = sensorDifferentialResponse(+dP0, p.P_ref, V0, Cm);
    dPsens_minus = sensorDifferentialResponse(-dP0, p.P_ref, V0, Cm);

    Tcav = (dPsens_plus - dPsens_minus) / (2 * dP0);
end

function Cm = evaluateInterfaceCompliance(epsn, p)
    % Returns effective interface compliance according to the selected model.
    %
    % complianceModel options:
    %   'lumped'   -> phenomenological compliance model
    %   'material' -> material-informed compliance surrogate using explicit
    %                 membrane thickness, modulus, and span

    switch lower(p.complianceModel)
        case 'lumped'
            % Purely phenomenological effective compliance model
            Cm = p.Cm0 / (1 + p.a1 * epsn + p.a2 * epsn^2);

        case 'material'
            % Material-informed effective compliance surrogate
            %
            % The present material model captures the cavity-side compliance trend only.
            % Additional strain-dependent transmission improvement is represented
            % separately through eta_cpl(eps).
            %
            % Compliance trend:
            %   larger span   -> more compliant
            %   thicker film  -> less compliant
            %   larger modulus-> less compliant

            A_eff = pi * p.a_mem^2;   % m^2, effective loaded area

            % Material-based zero-strain compliance scale
            Cm_mat0 = p.k_cm * (A_eff * p.a_mem) / (p.E_mem * p.t_mem);

            % Apply strain-dependent stiffening
            Cm = Cm_mat0 / (1 + p.b1 * epsn + p.b2 * epsn^2);

        otherwise
            error('Unknown compliance model: %s. Use ''lumped'' or ''material''.', p.complianceModel);
    end
end

function dPsens = sensorDifferentialResponse(dPsurf, Pref, V0, Cm)
    P3 = Pref + dPsurf / 2;
    P4 = Pref - dPsurf / 2;

    P1 = solveCavityPressureAir(P3, Pref, V0, Cm);
    P2 = solveCavityPressureAir(P4, Pref, V0, Cm);

    dPsens = P1 - P2;
end

function eta = evaluateCouplingEfficiency(epsn, p)
    % Strain-dependent interface coupling efficiency.
    %
    % Intended interpretation:
    % - At low pre-strain, pressure transmission is inefficient
    % - At moderate pre-strain, coupling improves
    % - At high pre-strain, the benefit tapers as stiffness dominates

    if ~p.useCouplingTerm
        eta = 1.0;
        return;
    end

    if p.eps_peak_cpl > 0
        eta = 1 + p.A_cpl * (epsn / p.eps_peak_cpl) * exp(1 - epsn / p.eps_peak_cpl);
    else
        eta = 1.0;
    end
end

function Pi = solveCavityPressureAir(Pext, Pref, V0, Cm)
    % Isothermal trapped-air cavity model based on ideal gas compression
    %
    % Assumes a sealed air cavity where pressure-volume behavior follows:
    %     P * V = constant  (isothermal ideal gas)
    %
    % External pressure (Pext) loads the compliant interface, modifying the
    % effective cavity volume via Cm, and resulting in a new internal pressure Pi.
    %
    % Governing relation:
    %     Pi * [V0 - Cm*(Pext - Pi)] = Pref * V0
    %
    % Rearranged to quadratic form:
    %     Cm*Pi^2 + (V0 - Cm*Pext)*Pi - Pref*V0 = 0

    A = Cm;
    B = V0 - Cm * Pext;
    C = -Pref * V0;

    rts = roots([A B C]);
    rts = rts(abs(imag(rts)) < 1e-10);
    rts = real(rts);
    rts = rts(rts > 0);
    
    if isempty(rts)
        Pi = Pref;
    else
        [~, idx] = min(abs(rts - Pref));
        Pi = rts(idx);
    end
end