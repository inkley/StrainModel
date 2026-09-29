function renderFigure05(pointsFile,out)
% Figure 5 uses measured recording means, with no offset subtraction.
p=readtable(pointsFile);
if ~isfolder(out), mkdir(out); end
    scales=[75 80 85 90 100];
    labels={'75%','80%','85%','90%','Bare-Port Reference'};
    markers={'s','d','o','h','+'};
    colors=[0 .447 .741; .85 .325 .098; .929 .694 .125; .494 .184 .556; .466 .674 .188];
    beta=unique(p.beta);
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
    handles(g)=scatter(ax,beta,mu,24,markers{g},'MarkerEdgeColor',colors(g,:), ...
        'MarkerFaceColor',colors(g,:),'LineWidth',.8,'DisplayName',labels{g});
    plot(ax,beta,mu-sd,'--','Color',colors(g,:),'Tag','SDbound', ...
        'LineWidth',.7,'HandleVisibility','off');
    plot(ax,beta,mu+sd,'--','Color',colors(g,:),'Tag','SDbound', ...
        'LineWidth',.7,'HandleVisibility','off');
end

xlabel(ax,'Orientation \beta (radians)'); ylabel(ax,'Mean Output Voltage (V)');
title(ax,'Longitudinal Voltage–Orientation Response');
xticks(ax,0:pi/4:2*pi);
xticklabels(ax,{'0','\pi/4','\pi/2','3\pi/4','\pi','5\pi/4','3\pi/2','7\pi/4','2\pi'});
xlim(ax,[0 2*pi]); ylim(ax,[.4 1.8]); yticks(ax,.4:.2:1.8);
leg=legend(ax,handles,labels,'Location','southeast'); leg.AutoUpdate='off';
fig.UserData=struct('selection','both_good','offsetCorrect',false, ...
    'trialIDs',unique(p.trial),'bounds','Pointwise mean +/- between-trial sample SD of recording means');
set(ax,'FontName','Times New Roman','FontSize',14,'LineWidth',1);
box(ax,'on');
savefig(fig,fullfile(out,'fig05_rotation_response.fig'));
exportgraphics(fig,fullfile(out,'fig05_rotation_response.png'),'Resolution',300);
exportgraphics(fig,fullfile(out,'fig05_rotation_response.pdf'),'ContentType','vector');
end
