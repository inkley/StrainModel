%% strainModel_v40.m -- standalone manuscript model
% Run this file in MATLAB R2025b (the release used for verification).
% No CSV, MAT, formatter, or helper files are required. Outputs are written
% beside this script in results_v40; Figures 7/8 are also exported formatted.
% This is the frozen v40 presentation model, not its archived exploratory sweeps.
% Experimental summaries are embedded below; raw-trial reduction is a separate
% calibration workflow. Editing a normalization does not change equilibrium.
% The 50% thickness-restoration and all-inward redistribution rules are
% illustrative assumptions, not measured or fitted mechanical properties.
% Each R_bare is conditional on its model, not a separate measured reference.
% The nominal R_bare is selected; the other two are least-squares level fits.
clear; clc;
root = fileparts(mfilename('fullpath'));
out = fullfile(root, 'results_v40');
if ~isfolder(out), mkdir(out); end

%% Physical parameters (SI units unless noted)
p.P_gas0 = 101325;                 % sealed absolute gas pressure, Pa
p.rho_water = 1000;                % assumed freshwater density, kg/m^3
p.g_water = 9.81;                  % m/s^2
p.baselineDepth_m = 0.3048;        % membrane-center depth, m
p.P_static = p.P_gas0 + p.rho_water*p.g_water*p.baselineDepth_m;
p.dP0 = 800;                      % total differential, +/-400 Pa
p.gasExponent = 1;                % isothermal air
p.r_forced_m = 5.5e-3;            % washer inner/pressure-loaded radius
p.strainWeightingOuterRadius_m = 21.978e-3/2;
p.strainWeightingCenterRadius_m = p.r_forced_m;
% Preserve the v39 funnel-CAD volume calculation, in mm^3 before conversion.
r = p.r_forced_m*1e3;
h1 = 2.506; h2 = 6.617;
p.V00 = (pi*r^2*h1 + pi*h2^2*r - pi*h2^3/3)*1e-9;
% Tubing, sensor, and fitting volumes were not measured or added.
p.t_plate0 = 0.020*0.0254;         % unstretched latex thickness, m
p.nu_plate = 0.49;
p.E_plate0 = 6e5;                 % Pa, constant effective modulus
p.tolP = 1e-5;                    % root residual or pressure bracket tolerance
p.maxIter = 1000;
p.relax = 0.08;                   % fixed-point numerical check only
p.solverMode = 'root';
p.minVolumeFraction = 1e-6;        % numerical guard, inactive in reported states
cfg.nominalRbare = 0.290;
cfg.restoredThicknessFraction = 0.5;
cfg.inwardFraction = 1.0;

%% Embedded experimental summary (frozen calibration rebuild)
% Source: calibration_rebuild/both_good/experimental_summary.csv from v39.
% Means are individual membrane-trial slopes / mean same-day bare slopes.
% SEM = sample SD of these ratios / sqrt(n); shared-reference uncertainty
% is not separately propagated. Exact archived precision is retained.
% Four additional bare trials appear only in response plots (17 total;
% 13 supply same-day references). No raw trial data are embedded here.
obs = table([90;85;80;75], [11;6;6;10], [10;6;3;7], ...
    [1.1524000577985434;1.3274895411289143;1.2704967098145357;1.2245493081181904], ...
    [0.08493972983326892;0.08193755157780681;0.044504994837231936;0.042193622437757594], ...
    'VariableNames', {'scale','raw_n','normalized_n','mean','sem'});
writetable(obs, fullfile(out,'embedded_experimental_summary.csv'));

%% Build current manuscript outputs and run numerical checks
buildPresentation(out,p,cfg,obs);
runSolverAudit(out,p,cfg);
fprintf('\nv40 complete. Results: %s\nFormatted Figures 7/8: %s\n',out,fullfile(out,'formatted'));

function buildPresentation(out,p,cfg,obs)
% Single-panel, three-case comparison. The supplied state solver sets the
% explicitly assumed clamp geometry; only normalization factors are fitted.
scales=[90;85;80;75];[ok,j]=ismember(scales,obs.scale);assert(all(ok));obs=obs(j,:);
en=[.111;.176;.250;.333];
% Include exact experimental states in a dense solve; do not interpolate fits.
strainGrid=unique([linspace(0,.36,241)';en]);
raw=zeros(numel(strainGrid),3); eff=nan(size(strainGrid)); map=cell(numel(strainGrid),1); energy=zeros(numel(strainGrid),3);
for i=1:numel(strainGrid)
 st=solvePresentationState(strainGrid(i),p,cfg);raw(i,:)=st.raw;eff(i)=st.effectiveStrain;
 map{i}=st.mapping;energy(i,:)=st.energy;
 assert(st.converged && st.volumeError<1e-12);
end
assert(all(isfinite(raw(:))) && all(raw(:)>=0 & raw(:)<=1));
assert(abs(eff(1))<1e-12 && all(diff(eff)>0),'Effective-strain axis must start at zero and increase.');
[ok,j]=ismember(en,strainGrid);assert(all(ok));md=raw(j,:);y=obs.mean;sem=obs.sem;
R=zeros(1,3);Rw=R;names=["Full nominal-strain tension","Effective-strain tension (volume-conserving)","Locally relaxed + deformation-induced tension"];
fitRows=cell(3,7);
for k=1:3
 if k<3,idx=2:4;objective="Postpeak response levels";else,idx=1:4;objective="All-four response levels";end
 a=md(idx,k);target=y(idx);w=1./sem(idx).^2;
 scale=max(1,(a'*target)/(a'*a));scaleW=max(1,sum(w.*a.*target)/sum(w.*a.^2));
 R(k)=1/scale;Rw(k)=1/scaleW;
 % Both optima checked against a fine numerical sweep.
 rg=linspace(.01,1,9901);L=sum((a./rg-target).^2,1);
 assert(sum((a*scale-target).^2)<=min(L)+1e-11);
 leastSquaresR=R(k);
 if k==1 && isfield(cfg,'nominalRbare')
  R(k)=cfg.nominalRbare;assert(R(k)>0 && R(k)<=1);
  objective="Illustrative placement below 85% mean";
  assert(md(2,k)/R(k)<y(2));
 end
 fitRows(k,:)={names(k),objective,R(k),leastSquaresR,Rw(k),sqrt(mean((a/R(k)-target).^2)),md(2,k)/R(k)-md(4,k)/R(k)};
end
fits=cell2table(fitRows,'VariableNames',{'case_name','objective','R_bare_used','R_bare_unweighted_fit','R_bare_SEM_weighted','fit_subset_RMSE','predicted_85_to_75_drop'});
writetable(fits,fullfile(out,'champion_fit_summary.csv'));
M=struct2table(vertcat(map{:}));writetable(M,fullfile(out,'champion_geometry.csv'));
curves=table(strainGrid,eff,raw(:,1),raw(:,2),raw(:,3),raw(:,1)/R(1),raw(:,2)/R(2),raw(:,3)/R(3),'VariableNames',{'nominal_strain','effective_strain','nominal_raw','volume_raw','relaxed_raw','nominal_normalized','volume_normalized','relaxed_normalized'});
writetable(curves,fullfile(out,'champion_curves.csv'));
points=table(scales,en,eff(j),y,sem,md(:,1)/R(1),md(:,2)/R(2),md(:,3)/R(3),'VariableNames',{'scale','nominal_strain','effective_strain','mean','SEM','nominal_prediction','volume_prediction','relaxed_prediction'});
writetable(points,fullfile(out,'champion_experimental_comparison.csv'));
colors=[0 .447 .741;.85 .325 .098;.929 .694 .125];
f=figure('Visible','off','Theme','light','Color','w','Position',[100 100 1250 760]);ax=axes(f,'Position',[.10 .26 .87 .60]);hold(ax,'on');h=gobjects(1,5);
for k=1:3
 h(k)=plot(ax,eff,raw(:,k)/R(k),'Color',colors(k,:),'LineWidth',2.2);
end
h(4)=errorbar(ax,eff(j),y,sem,'ko','MarkerFaceColor','k','MarkerSize',8,'LineStyle','none','LineWidth',1.6,'CapSize',9);
h(5)=yline(ax,1,':','Color',[.4 .4 .4],'LineWidth',1.4);
xlabel(ax,'Effective engineering strain, \epsilon_{eff} (-)');ylabel(ax,'Bare-port-normalized response (-)');
xlim(ax,[0 eff(end)]);
% Include full fitted curves at zero rather than silently clipping them.
normalized=raw./R; vals=[normalized(:);y+sem]; upper=ceil(max(vals)*10)/10;
% Explicit lower limit from normalized response matrix (column order).
ylim(ax,[max(0,floor(min([normalized(:);y-sem])*10)/10-.1),upper+.05]);
set(ax,'FontSize',13,'LineWidth',1.2,'Box','on','TickDir','in','XGrid','on','YGrid','on','GridAlpha',.15);
title(ax,'Tension-Model Cases and Experimental Response','FontSize',18,'FontWeight','bold');
labels={sprintf('Full nominal-strain tension; R_{bare} = %.3f',R(1)),sprintf('Effective-strain tension (volume-conserving); R_{bare} = %.3f',R(2)),sprintf('Locally relaxed + deformation-induced tension; R_{bare} = %.3f',R(3)),'Experiment (bare-port normalized, mean \pm SEM)','Bare-port reference'};
lg=legend(ax,h,labels,'NumColumns',2,'Location','southoutside','FontSize',10,'Interpreter','tex');lg.Position=[.10 .065 .87 .12];
savefig(f,fullfile(out,'fig08_model_experiment_comparison.fig'));
exportgraphics(f,fullfile(out,'fig08_model_experiment_comparison.png'),'Resolution',180);
exportgraphics(f,fullfile(out,'fig08_model_experiment_comparison.pdf'),'ContentType','vector');
exportgraphics(f,fullfile(out,'fig08_model_experiment_comparison.eps'),'ContentType','vector');
close(f);
formatModelFigure(fullfile(out,'fig08_model_experiment_comparison.fig'),fullfile(out,'formatted','fig08_model_experiment_comparison'),true);
% Matched energy output for the same center state, independent of R_bare.
writetable(table(strainGrid,eff,energy(:,1),energy(:,2),energy(:,3),'VariableNames',{'nominal_strain','effective_strain','U_low_uJ','U_hydrostatic_uJ','U_high_uJ'}),fullfile(out,'champion_energy.csv'));
f=figure('Visible','off','Theme','light','Color','w');hold on;
plot(eff,energy,'LineWidth',1.8);xlabel('Effective engineering strain, \epsilon_{eff} (-)');ylabel('Structural energy (microjoules)');grid on;
legend('Low side','Hydrostatic','High side','Location','best');title('Energy diagnostic: same volume-conserving installed state');
savefig(f,fullfile(out,'champion_energy.fig'));exportgraphics(f,fullfile(out,'champion_energy.png'),'Resolution',180);close(f);
exportPotentialEnergyVsTension(out,M,energy);
save(fullfile(out,'champion_workspace.mat'),'p','cfg','fits','curves','points','M','energy');
fid=fopen(fullfile(out,'champion_validation.txt'),'w');fprintf(fid,'All dense-strainGrid states converged, passive raw response, volume error <1e-12.\nEffective coordinate starts at zero and is strictly increasing.\nExact experimental states included in solve. Least-squares optima checked against grid; nominal display normalization is a separate illustrative setting.\nGeometry assumption is not calibrated. SEM-weighted alternatives reported without confidence intervals.\n');fclose(fid);
disp(fits);disp(points);
end

function exportPotentialEnergyVsTension(out,M,energy)
% Figure 7 companion: same installed state as the volume-conserving case.
assert(size(energy,1)==height(M));
[ok,markers]=ismember([0;.111;.176;.250;.333],M.nominal_strain);assert(all(ok));
last=markers(end);idx=1:last;
x=M.center_tension_Npm(idx);assert(all(diff(x)>0));
f=figure('Visible','off','Theme','light','Color','w','Position',[100 100 880 650]);
ax=axes(f,'Position',[.13 .14 .83 .75]);hold(ax,'on');
colors=[0 .447 .741;.85 .325 .098;.929 .694 .125];
cols=[2 3 1];symbols={'o','s','d'};h=gobjects(1,3);
for k=1:3
 h(k)=plot(ax,x,energy(idx,cols(k)),'Color',colors(k,:),'LineWidth',2, ...
  'Marker',symbols{k},'MarkerIndices',markers,'MarkerSize',7,'MarkerFaceColor','none');
end
xlabel(ax,'Effective installed tension, T_0 (N/m)');
ylabel(ax,'Structural potential energy, U_{struct} (\muJ)');
title(ax,'Potential Energy vs. Effective Installed Tension','FontWeight','bold','FontSize',16);
legend(ax,h,{'Hydrostatic equilibrium','High side: P_{static} + \DeltaP/2','Low side: P_{static} - \DeltaP/2'},'Location','southeast','FontSize',12);
set(ax,'FontSize',13,'LineWidth',1.2,'Box','on','TickDir','in','XGrid','on','YGrid','on','GridAlpha',.18);
xlim(ax,[0 ceil(max(x)/5)*5]);ylim(ax,[1 10]);
base=fullfile(out,'fig07_potential_energy_vs_tension');
savefig(f,[base '.fig']);exportgraphics(f,[base '.png'],'Resolution',180);
exportgraphics(f,[base '.pdf'],'ContentType','vector');exportgraphics(f,[base '.eps'],'ContentType','vector');
close(f);
formatModelFigure([base '.fig'],fullfile(out,'formatted','fig07_potential_energy_vs_tension'),false);
end

function [D,E,T,a] = installedNominalState(strain,p)
    a=p.r_forced_m;
    h=p.t_plate0*(1+strain)^(-2*p.nu_plate/(1-p.nu_plate));
    E=p.E_plate0;
    D=E*h^3/(12*(1-p.nu_plate^2));
    T=E*h*strain/(1-p.nu_plate);
end

function st=solvePresentationState(enom,p,cfg)
% Original nominal/relaxed material closure retained for those two cases.
% Volume case uses incompressible thickness and a separate equivalent tension.
r=(p.strainWeightingOuterRadius_m/p.strainWeightingCenterRadius_m)^2-1;
hn=p.t_plate0/(1+enom)^2;
hc=hn+cfg.restoredThicknessFraction*(p.t_plate0-hn);
c=(hc/hn-1)/(cfg.inwardFraction*r);
hw=hn*(1-c);
outward=(1-cfg.inwardFraction)*r*c*hn;
volumeError=abs((1+r)*hn-r*hw-hc-outward)/((1+r)*hn);
eff=sqrt(p.t_plate0/hc)-1;
assert(hw>0 && c>=-1e-12 && c<1 && eff>=-1e-12,'Invalid tensile geometry.');
eff=max(0,eff); % only floating-point roundoff at the exact zero limit
[D,E,T,a]=installedNominalState(enom,p);
hOriginal=p.t_plate0*(1+enom)^(-2*p.nu_plate/(1-p.nu_plate));
Dc=E*hc^3/(12*(1-p.nu_plate^2));Tc=E*hc*eff/(1-p.nu_plate);
external=p.P_static+[-p.dP0/2 p.dP0/2];
piNom=zeros(1,2);piVol=piNom;piRelax=piNom;ok=true;
for k=1:2
 [piNom(k),~,~,~,~,convN]=solveCavityPressure(external(k),a,D,T,p,struct('Pi',p.P_gas0));
 [piVol(k),~,~,~,~,convV]=solveCavityPressure(external(k),a,Dc,Tc,p,struct('Pi',p.P_gas0));
 nl=solveConservativeNonlinearState(external(k),a,D,0,E,hOriginal,.5,1,p);
 piRelax(k)=nl.Pi;ok=ok && convN && convV && nl.converged;
end
energy=zeros(1,3);
for k=1:3
 nl=solveConservativeNonlinearState(p.P_static+(k-2)*p.dP0/2,a,Dc,Tc,E,hc,.5,1,p);
 energy(k)=1e6*nl.Ustruct;ok=ok && nl.converged;
end
st.raw=[diff(piNom) diff(piVol) diff(piRelax)]/p.dP0;
st.effectiveStrain=eff;st.energy=energy;st.converged=ok;st.volumeError=volumeError;
st.mapping=struct('nominal_strain',enom,'effective_strain',eff,'restored_thinning_fraction',cfg.restoredThicknessFraction, ...
 'inward_fraction',cfg.inwardFraction,'washer_compression',c,'h_nom_mm',hn*1e3,'h_washer_mm',hw*1e3, ...
 'h_center_mm',hc*1e3,'center_tension_Npm',Tc,'volume_relative_error',volumeError);
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

function F = nonlinearEnergyResidual(w, Psurf, A, B, areaFactor, p)
    V = p.V00 - areaFactor*w;
    if V <= 0
        F = -Inf;
        return;
    end
    Pi = p.P_gas0 * (p.V00/V)^p.gasExponent;
    F = Psurf - Pi - A*w - B*w^3;
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

    dV = dV_raw;

    minVolume = max(p.minVolumeFraction * V00, realmin);
    Vi = max(V00 - dV, minVolume);
end

function runSolverAudit(out,p,cfg)
% Two linear retained-tension cases; nonlinear case has no fixed-point solver.
% Root ignores starting-state input: route equality is a consistency check,
% not independent validation of material path independence.
    rows=cell(30,13); row=0;
    for en=[0 .111 .176 .250 .333]
        [D,E,T,a]=installedNominalState(en,p);
        hn=p.t_plate0/(1+en)^2;
        hc=hn+cfg.restoredThicknessFraction*(p.t_plate0-hn);
        ec=sqrt(p.t_plate0/hc)-1;
        Dc=E*hc^3/(12*(1-p.nu_plate^2)); Tc=E*hc*ec/(1-p.nu_plate);
        for cas=1:2
            if cas==1, dk=D;tk=T;name='nominal';else,dk=Dc;tk=Tc;name='effective';end
            [~,~,~,~,~,ch,hs]=solveCavityPressure(p.P_static,a,dk,tk,p,struct('Pi',p.P_gas0));
            for off=[-p.dP0/2 0 p.dP0/2]
                pe=p.P_static+off;
                [pr,~,vr,wr,~,cr]=solveCavityPressure(pe,a,dk,tk,p,struct('Pi',p.P_gas0));
                [ps,~,vs,ws,~,cs]=solveCavityPressure(pe,a,dk,tk,p,hs);
                pf=p;pf.solverMode='fixedPoint';
                [fp,~,fv,fw,it,cf]=solveCavityPressure(pe,a,dk,tk,pf,struct('Pi',p.P_gas0));
                residual=cavityPressureResidual(fp,pe,a,dk,tk,p);
                row=row+1;
                rows(row,:)={name,en,off,ch&&cr&&cs,abs(pr-ps),abs(vr-vs),abs(wr-ws),cf,it,abs(pr-fp),residual,abs(vr-fv),abs(wr-fw)};
            end
        end
    end
    audit=cell2table(rows,'VariableNames',{'case_name','nominal_strain','port_offset_Pa','root_converged','path_dPi_Pa','path_dV_m3','path_dw_m','fixed_converged','fixed_iterations','root_fixed_dPi_Pa','fixed_residual_Pa','root_fixed_dV_m3','root_fixed_dw_m'});
    writetable(audit,fullfile(out,'solver_comparison.csv'));
    assert(all(audit.root_converged & audit.fixed_converged),'Solver convergence check failed.');
    assert(all(audit.path_dPi_Pa==0 & audit.path_dV_m3==0 & audit.path_dw_m==0));
    assert(max(audit.root_fixed_dPi_Pa)<8e-5,'Solver agreement exceeds the manuscript audit threshold.');
    fid=fopen(fullfile(out,'solver_validation.txt'),'w');
    assert(fid~=-1,'Cannot write validation report.');
    fprintf(fid,'Both solvers converged in %d comparisons.\nMaximum pressure difference %.12g Pa.\nMaximum fixed-point residual %.12g Pa.\n',height(audit),max(audit.root_fixed_dPi_Pa),max(abs(audit.fixed_residual_Pa)));
    fprintf(fid,'Bisection ignores initial pressure; direct/staged equality is by construction.\nFixed-point termination tests relaxed iterate change, not equilibrium residual.\nThis audit does not include a fixed-point solver for the nonlinear relaxed case.\n');
    fclose(fid);
    fprintf('Solver cross-check: %d/30 converged; max difference %.6g Pa.\n',sum(audit.fixed_converged),max(audit.root_fixed_dPi_Pa));
end

function formatModelFigure(inputFig,outputBase,isComparison)
% Embedded Figure 7/8 subset of formatCalibrationFigure.m (v39 styling).
    fig=openfig(inputFig,'invisible');
    cleanup=onCleanup(@() close(fig));
    folder=fileparts(outputBase);
    if ~isfolder(folder), mkdir(folder); end
    fig.Color='white'; fig.Units='inches';
    if isComparison, fig.Position=[1 1 7.16 4.3];else,fig.Position=[1 1 3.5 2.55];end
    fig.PaperPositionMode='auto';
    ax=findall(fig,'Type','Axes');
    for a=reshape(ax,1,[])
        set(a,'FontName','Helvetica','FontSize',8,'TickLabelInterpreter','tex', ...
            'LineWidth',.8,'TickDir','in','TickLength',[.018 .018],'Box','on', ...
            'Layer','bottom','XGrid','on','YGrid','on','XMinorGrid','off', ...
            'YMinorGrid','off','GridColor',[.82 .82 .82],'GridAlpha',.55,'MinorGridAlpha',.30);
        set(a.XLabel,'FontName','Helvetica','FontSize',9,'FontWeight','normal','Interpreter','tex');
        set(a.YLabel,'FontName','Helvetica','FontSize',9,'FontWeight','normal','Interpreter','tex');
        if isComparison, titleSize=11;else,titleSize=9;end
        set(a.Title,'FontName','Helvetica','FontSize',titleSize,'FontWeight','bold','Interpreter','tex');
        for h=reshape(findall(a,'Type','Line'),1,[])
            h.LineWidth=1.1;
            if ~strcmp(h.Marker,'none'),h.MarkerSize=3.5;end
        end
        for h=reshape(findall(a,'Type','ErrorBar'),1,[])
            set(h,'LineStyle','none','LineWidth',1,'MarkerSize',5,'CapSize',5);
            uistack(h,'top');
        end
        leg=a.Legend;
        set(leg,'FontName','Helvetica','FontSize',7.5,'Box','on','LineWidth',.6,'Interpreter','tex');
        if isComparison
            leg.Location='southoutside';leg.NumColumns=2;
            a.Units='normalized';a.Position=[.105 .25 .875 .65];
        else
            leg.Location='southeast';
        end
    end
    drawnow;
    exportgraphics(fig,[outputBase '.png'],'Resolution',600,'BackgroundColor','white');
    exportgraphics(fig,[outputBase '.pdf'],'ContentType','vector','BackgroundColor','white');
    exportgraphics(fig,[outputBase '.eps'],'ContentType','vector','BackgroundColor','white');
    savefig(fig,[outputBase '.fig']);
    clear cleanup
end
