function renderFigure06(pointsFile,out)
% Measured voltage, equal-weight trial means, and between-trial SD.
p=readtable(pointsFile);
if ~isfolder(out), mkdir(out); end
    scales=[75 80 85 90 100];
    labels={'75%','80%','85%','90%','Bare-Port Reference'};
    markers={'s','d','o','h','+'};
    colors=[0 .447 .741; .85 .325 .098; .929 .694 .125; .494 .184 .556; .466 .674 .188];
    beta=unique(p.beta);
    pressure=1000*9.81*.040004*sin(beta);
    fig=figure('Visible','off','Color','w');
    cleanup=onCleanup(@()close(fig));
if isprop(fig,'Theme'), fig.Theme='light'; end
    fig.Position=[100 100 900 600];
    ax=axes(fig); hold(ax,'on'); handles=gobjects(1,5);
for g=1:5
    mu=zeros(size(beta)); sd=mu;
    for k=1:numel(beta)
        values=p.voltage(p.scale==scales(g) & p.beta==beta(k));
        mu(k)=mean(values); sd(k)=std(values);
    end
    coeff=polyfit(pressure,mu,1);
    xx=linspace(min(pressure),max(pressure),100);
    plot(ax,xx,polyval(coeff,xx),'Color',colors(g,:),'LineWidth',1.1,'HandleVisibility','off');
    handles(g)=scatter(ax,pressure,mu,24,markers{g},'MarkerEdgeColor',colors(g,:), ...
        'MarkerFaceColor',colors(g,:),'LineWidth',.8,'DisplayName',labels{g});

    % Repeated pressures occur at different angles. Preserve every mean
    % marker, but use the outer pointwise SD envelope at each pressure.
    [boundPressure,~,idx]=unique(round(pressure,8));
    lower=accumarray(idx,mu-sd,[],@min);
    upper=accumarray(idx,mu+sd,[],@max);
    plot(ax,boundPressure,lower,'--','Color',colors(g,:),'Tag','SDbound', ...
        'LineWidth',.7,'HandleVisibility','off');
    plot(ax,boundPressure,upper,'--','Color',colors(g,:),'Tag','SDbound', ...
        'LineWidth',.7,'HandleVisibility','off');
end

xlabel(ax,'\Delta P (Pa)'); ylabel(ax,'Mean Output Voltage (V)');
title(ax,'Longitudinal Voltage–Pressure Response');
xlim(ax,[-400 400]); xticks(ax,-400:200:400);
ylim(ax,[.4 1.8]); yticks(ax,.4:.2:1.8);
leg=legend(ax,handles,labels,'Location','southwest'); leg.AutoUpdate='off';
fig.UserData=struct('selection','both_good','offsetCorrect',false, ...
    'trialIDs',unique(p.trial),'bounds','Outer envelope of pointwise mean +/- between-trial SD at repeated pressures', ...
    'fitDefinition','Linear fit to configuration mean recording voltages; not used for Figure 9');
set(ax,'FontName','Times New Roman','FontSize',14,'LineWidth',1);
box(ax,'on');
savefig(fig,fullfile(out,'fig06_voltage_pressure_response.fig'));
exportgraphics(fig,fullfile(out,'fig06_voltage_pressure_response.png'),'Resolution',300);
exportgraphics(fig,fullfile(out,'fig06_voltage_pressure_response.pdf'),'ContentType','vector');
end
