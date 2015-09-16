%PLOT_EXP4 script generates Fig. 10 of the manuscript ''Sensitivity 
%encoding for aligned multishot magnetic resonance reconstruction'', L. 
%Cordero-Grande, R. P. A. G. Teixeira, E. J. Hughes, J. Hutter, A. N. 
%Price, and J. V. Hajnal using the results of running the script 
%alignedSENSE_Exp4

%Study case considered
pathOu='./data';%Data path
nameExp=sprintf('%s/Exp4',pathOu);
load(nameExp)

%NS=[2 4 8 16 32];%Number of shots
S=length(NS);
V=length(theta);

close all
cont=1;
%pathFig=sprintf('%s/Fig10',pathOu);
%mkdir(pathFig)

FontSize=20;
FontSizeA=16;
FontSizeB=13;
ColorSet=varycolor(S);
markers = {'o','s','d','^','>','v','<','p','h','x'};
for v=1:V        
    figure(cont)
    for s=1:S
        NSP=length(errF{v}{1}{s});
        plot(log10(1:NSP),log10(errF{v}{1}{s}),'Color',ColorSet(s,:),'Marker',markers{s},'LineWidth',2);
        hold on
        NSVal{s}=sprintf('S=%d',NS(s));
    end   
    axis([0 3.5 -6 8])
    xlabel('$\log_{10}(i)$','interpreter','latex','FontSize',FontSizeA)
    ylabel('$\log_{10}(f)$','interpreter','latex','FontSize',FontSizeA)
    title(sprintf('$\\Delta\\theta=$%d$^{\\circ}$',theta(v)),'interpreter','latex','Fontsize',FontSize)
    AX=legend(NSVal,'Location','SouthEast');
    LEG = findobj(AX,'type','text');
    set(AX,'Position',get(AX,'Position')+[0.025 -0.000 0 0])%To the right and to the bottom
    set(LEG,'FontSize',FontSizeB)
    legend('boxoff')
    grid on
    set(gca,'fontsize',FontSizeB)
    set(gcf,'Color',[1 1 1])
    %%set(gcf, 'Position', get(0,'Screensize'))        
    %export_fig(sprintf('%s/SubFigNew%02d.pdf',pathFig,cont))
    cont=cont+1;
end
