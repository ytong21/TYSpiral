figure(90)
set(gcf,'color','w','InvertHardcopy','off')
set(gcf,'units','centimeters','position',[4 4 40 20],'paperunits','centimeters','paperposition',[0 0 40 20])
clf
PlotFALim = [0 20]; FontSize = 12;   CLB = cell(10,1);   CLBFontSize = 12;
Figs = cell(4,3);

Figs{2,1} = subplot(3,2,3);
imagesc(ImgToPlot.b1CP(plotRolArray,plotCowArray,SliceIdx)',[0 5]); 
ylabel('DREAM B1 Map','FontSize',FontSize,'FontWeight','bold');axis off;
Figs{2,1}.YLabel.Visible = 'on';
CLB{1} = colorbar('FontSize',CLBFontSize);
CLB{1}.YLabel.String = 'Hz/V';    CLB{1}.YLabel.Rotation = 0;     
CLB{1}.YLabel.Position = [0.5 5.6 0];
hold on;colormap(Figs{2,1},'gray')
visboundaries(maskedMaps.mask(plotRolArray,plotCowArray,SliceIdx)','Color','r','LineWidth',0.4,...
    'EnhanceVisibility', true,'LineStyle','-');axis off;
nudge(Figs{2,1},[0 0.05 0 0])

Figs{3,1} = subplot(3,2,5); 
imagesc(ImgToPlot.b0(plotRolArray,plotCowArray,SliceIdx)',[-200 200]); 
ylabel('B0 Map','FontSize',FontSize,'FontWeight','bold');axis off;hold on;
Figs{3,1}.YLabel.Visible = 'on';
CLB{2} = colorbar('FontSize',CLBFontSize);
CLB{2}.YLabel.String = 'Hz';    CLB{2}.YLabel.Rotation = 0;     
%CLB{2}.YLabel.Position = [0 235 0];clc
CLB{2}.YLabel.Position = [0.5 245 0];
colormap(Figs{3,1},'gray')
visboundaries(maskedMaps.mask(plotRolArray,plotCowArray,SliceIdx)','Color','r','LineWidth',0.4,...
    'EnhanceVisibility', true,'LineStyle','-');axis off;
nudge(Figs{3,1},[0 0.1 0 0])

Figs{1,1} = subplot(3,2,1); 
imagesc(maskedMaps.localiser(plotRolArray,plotCowArray,SliceIdx)'); 
ylabel('TOF image','FontSize',FontSize,'FontWeight','bold');axis off;colormap(Figs{1,1},'gray')
Figs{1,1}.YLabel.Visible = 'on';
hold on
visboundaries(maskedMaps.mask(plotRolArray,plotCowArray,SliceIdx)','Color','r','LineWidth',0.4,...
    'EnhanceVisibility', true,'LineStyle','-');axis off;
hp3 = get(subplot(3,2,1),'Position');
set(Figs{1,1},'Position',[hp3(1) hp3(2) Figs{3,1}.Position(3) hp3(4)]);



HoriOffset = 0.038;
Figs{1,2} = subplot(3,2,2);
imagesc(FullImage.CP(plotRolArray,plotCowArray,2)',PlotFALim); axis off;
title('Flip angle maps','FontSize',FontSize);
ylabel('CP Mode','FontSize',FontSize,'FontWeight','bold');
Figs{1,2}.YLabel.Visible = 'on';
set(Figs{1,2},'Position',[Figs{1,2}.Position(1) Figs{1,2}.Position(2) Figs{1,2}.Position(3) Figs{1,2}.Position(4)]);
colormap(Figs{1,2},'hot')

Figs{2,2} = subplot(3,2,4);
imagesc(FullImage.PhaseOnly(plotRolArray,plotCowArray,2)',PlotFALim);axis off;
ylabel('Phase Only','FontSize',FontSize,'FontWeight','bold');
Figs{2,2}.YLabel.Visible = 'on';
set(Figs{2,2},'Position',[Figs{2,2}.Position(1) Figs{2,2}.Position(2) Figs{2,2}.Position(3) Figs{2,2}.Position(4)]);
colormap(Figs{2,2},'hot')
nudge(Figs{2,2},[0 0.05 0 0]);

Figs{3,2} = subplot(3,2,6);
imagesc(FullImage.AS(plotRolArray,plotCowArray,2)',PlotFALim);axis off;
set(Figs{3,2},'Position',[Figs{3,2}.Position(1) Figs{3,2}.Position(2) Figs{3,2}.Position(3) Figs{3,2}.Position(4)]);
ylabel('Full RF Shimming','FontSize',FontSize,'FontWeight','bold');
Figs{3,2}.YLabel.Visible = 'on';
CLB{6} = colorbar('southoutside','FontSize',CLBFontSize);
CLB{6}.YLabel.String = 'Degrees (°)';    CLB{6}.YLabel.Rotation = 0;   
colormap(Figs{3,2},'hot')
nudge(Figs{3,2},[0 0.1 0 0])
nudge(CLB{6},[0 0.02 0 0])
%


%

% Chopping the middle sections of the FA plots
%HoriOffset = 0.038;
 ToChop = 20;
      %plotCowArrayReduced = [col(1):((col(1)+col(end))/2-ToChop) (ceil((col(1)+col(end))/2)+ToChop):col(end)];
      plotRowArrayReduced = [row(1):(((row(1)+row(end))/2)-ToChop) (ceil((row(1)+row(end))/2)+ToChop):row(end)];
 FA_map_positions = cell(3,1);
for iDx = 1:3
    FA_map_positions{iDx} = get(Figs{iDx,2},'position');
    FA_map_positions{iDx}(3) = FA_map_positions{iDx}(3)*(numel(plotRowArrayReduced)/numel(plotRolArray));
    delete(Figs{iDx,2});
end
Figs{1,3} = subplot('Position',FA_map_positions{1});
imagesc(FullImage.CP(plotRowArrayReduced,plotCowArray,2)',PlotFALim); axis off;
title('Flip angle maps','FontSize',FontSize);
ylabel('CP Mode','FontSize',FontSize,'FontWeight','bold');
Figs{1,3}.YLabel.Visible = 'on';
%set(Figs{1,3},'Position',[Figs{1,3}.Position(1) Figs{1,3}.Position(2) Figs{1,3}.Position(3) Figs{1,3}.Position(4)]);
colormap(Figs{1,3},'hot')

Figs{2,3} = subplot('Position',FA_map_positions{2});
imagesc(FullImage.PhaseOnly(plotRowArrayReduced,plotCowArray,2)',PlotFALim);axis off;
ylabel('Phase Only','FontSize',FontSize,'FontWeight','bold');
Figs{2,3}.YLabel.Visible = 'on';
%set(Figs{2,3},'Position',[Figs{2,3}.Position(1) Figs{2,3}.Position(2) Figs{2,3}.Position(3) Figs{2,3}.Position(4)]);
colormap(Figs{2,3},'hot')
%nudge(Figs{2,3},[0 0.05 0 0]);

Figs{3,3} = subplot('Position',FA_map_positions{3});
imagesc(FullImage.AS(plotRowArrayReduced,plotCowArray,2)',PlotFALim);axis off;
%set(Figs{3,3},'Position',[Figs{3,3}.Position(1) Figs{3,3}.Position(2) Figs{3,3}.Position(3) Figs{3,3}.Position(4)]);
ylabel('Full RF Shimming','FontSize',FontSize,'FontWeight','bold');
Figs{3,3}.YLabel.Visible = 'on';
CLB{9} = colorbar('eastoutside','FontSize',CLBFontSize);
CLB{9}.YLabel.String = 'Degrees (°)';    %CLB{9}.YLabel.Rotation = 0;   
colormap(Figs{3,3},'hot')
%nudge(Figs{3,3},[0 0.1 0 0])

set(CLB{9},'Position',[ 0.878059964726631         0.209876543209877         0.017636684303351         0.715437545388525]);
nudge(CLB{9},[-0.04 0 0 0])
%%
% PlotDiffLim = [-20 20];
% Figs{1,3} = subplot(3,4,4);
% imagesc(DiffImage.CP(plotRolArray,plotCowArray,2)',PlotDiffLim); 
% title('Difference images','FontSize',FontSize);axis off;
% set(Figs{1,3},'Position',[Figs{1,3}.Position(1)-HoriOffset Figs{1,3}.Position(2) ...
%     Figs{1,3}.Position(3) Figs{1,3}.Position(4)]);
% colormap(Figs{1,3},'parula')
% 
% Figs{2,3} = subplot(3,4,8);
% imagesc(DiffImage.PhaseOnly(plotRolArray,plotCowArray,2)',PlotDiffLim); axis off;
% set(Figs{2,3},'Position',[Figs{2,3}.Position(1)-HoriOffset Figs{2,3}.Position(2) ...
%     Figs{2,3}.Position(3) Figs{2,3}.Position(4)]);
% colormap(Figs{2,3},'parula')
% nudge(Figs{2,3},[0 0.05 0 0])
% 
% 
% Figs{3,3} = subplot(3,4,12);
% imagesc(DiffImage.AS(plotRolArray,plotCowArray,2)',PlotDiffLim);axis off
% set(Figs{3,3},'Position',[Figs{3,3}.Position(1)-HoriOffset Figs{3,3}.Position(2) ...
%     Figs{3,3}.Position(3) Figs{3,3}.Position(4)]);
% CLB{7} = colorbar('southoutside');
% CLB{7} = colorbar('southoutside','FontSize',CLBFontSize);
% colormap(Figs{3,3},'parula')
% nudge(Figs{3,3},[0 0.1 0 0])
% nudge(CLB{7},[0 0.02 0 0])
% 
% CLB{7}.YLabel.String = 'Percentage error (%)';    CLB{7}.YLabel.Rotation = 0;     
% %CLB{7}.YLabel.Position = [0.5 700 0];
% 
% hp4 = get(Figs{3,3},'Position');
% %colorbar('Position', [hp4(1)+hp4(3)+0.015  hp4(2)  0.02  0.815],'FontSize',FontSize);
