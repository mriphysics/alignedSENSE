%PLOT_EXP2 script generates Fig. 8 of the manuscript ''Sensitivity encoding
%for aligned multishot magnetic resonance reconstruction'', L. 
%Cordero-Grande, R. P. A. G. Teixeira, E. J. Hughes, J. Hutter, A. N. 
%Price, and J. V. Hajnal using the results of running the script 
%alignedSENSE_Exp2

%Study case considered
pathOu='./data';%Data path
nameExp=sprintf('%s/Exp2',pathOu);
load(nameExp)

EncMeth{1}='Cartesian Sequential';
EncMeth{2}='Cartesian Parallel 1D';

E=length(EncMeth);
V=length(theta);

close all
cont=1;
%pathFig=sprintf('%s/Fig8',pathOu);
%mkdir(pathFig)

FontSize=20;
FontSizeA=16;
FontSizeB=13;
ColorSet=varycolor(V); 
markers = {'o','s','d','^','>','v','<','p','h','x'};
for e=1:E
    figure(cont)
    for v=1:V
        NS=length(errFT{v}{e}{1});
        plot(log10(1:NS),log10(errFT{v}{e}{1}),'Color',ColorSet(v,:),'Marker',markers{v},'LineWidth',2);
        hold on
        ThVal{v}=sprintf('\\Delta\\theta=%d%c',theta(v),char(176));
    end   
    xlabel('$\log_{10}(i)$','interpreter','latex','FontSize',FontSizeA)
    ylabel('$\log_{10}(f)$','interpreter','latex','FontSize',FontSizeA)
    title(sprintf('%s',EncMeth{e}),'interpreter','latex','Fontsize',FontSize)
    AX=legend(ThVal,'Location','SouthWest');
    LEG = findobj(AX,'type','text');
    set(LEG,'FontSize',FontSizeB)
    legend('boxoff')
    grid on
    set(gca,'fontsize',FontSizeB)
    set(gcf,'Color',[1 1 1])
    %%set(gcf, 'Position', get(0,'Screensize'))
    %export_fig(sprintf('%s/SubFigNew%02d.pdf',pathFig,cont))
    cont=cont+1;
end

