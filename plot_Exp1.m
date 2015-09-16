%PLOT_EXP1 script generates Fig. 5 of the manuscript ''Sensitivity encoding
%for aligned multishot magnetic resonance reconstruction'', L. 
%Cordero-Grande, R. P. A. G. Teixeira, E. J. Hughes, J. Hutter, A. N. 
%Price, and J. V. Hajnal using the results of running the script 
%alignedSENSE_Exp1

%Study case considered
pathOu='./data';%Data path
nameExp=sprintf('%s/Exp1',pathOu);
load(nameExp)

EncMeth{1}='Cartesian Sequential';
EncMeth{2}='Radial Sequential';
EncMeth{3}='Cartesian Parallel 1D';
EncMeth{4}='Cartesian Parallel 2D';
EncMeth{5}='Random';

E=length(EncMeth);
V=length(theta);

close all
cont=1;
%pathFig=sprintf('%s/Fig5',pathOu);
%mkdir(pathFig)

FontSize=20;
FontSizeA=16;
FontSizeB=13;
ColorSet=varycolor(E);    
markers = {'o','s','d','^','>','v','<','p','h','x'};
for v=1:V
    figure(cont)        
    for e=1:E
        NS=length(errF{v}{e}{1});
        plot(log10(1:NS),log10(errF{v}{e}{1}),'Color',ColorSet(e,:),'Marker',markers{e},'LineWidth',2);
        hold on
    end   
    axis([0 3.5 -6 8])
    xlabel('$\log_{10}(i)$','interpreter','latex','FontSize',FontSizeA)
    ylabel('$\log_{10}(f)$','interpreter','latex','FontSize',FontSizeA)
    title(sprintf('$\\Delta\\theta=$%d$^{\\circ}$',theta(v)),'interpreter','latex','Fontsize',FontSize)
    AX=legend(EncMeth,'Location','SouthEast');    
    LEG = findobj(AX,'type','text');
    set(AX,'Position',get(AX,'Position')+[0.015 -0.000 0 0])%To the right and to the bottom
    set(LEG,'FontSize',FontSizeB)           
    legend('boxoff')
    grid on
    set(gca,'fontsize',FontSizeB)
    set(gcf,'Color',[1 1 1])
    %%set(gcf, 'Position', get(0,'Screensize'))
    %export_fig(sprintf('%s/SubFigNew%02d.pdf',pathFig,cont))
    cont=cont+1;
end
