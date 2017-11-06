% Plotting
  if ~exist('plotRolArray','var')
      figure(56)
      imagesc(maskedMaps.localiser(:,:,SliceIdx));
      title('Select Rect Regions to Plot')
      h2plot = imrect;  wait(h2plot);
      plotMask = logical(createMask(h2plot));
      [rol,cow] = find(plotMask);
      plotRolArray = rol(1):rol(end);
      plotCowArray = cow(1):cow(end);     
  end
 % Other plots
figure(90)
set(gcf,'color','w','InvertHardcopy','off')
set(gcf,'units','centimeters','position',[4 4 40 20],'paperunits','centimeters','paperposition',[0 0 40 20])
clf
PlotFALim = [0 30]; FontSize = 6;   CLB = cell(5,1); 
Figs = cell(3,3);
Figs{2,1} = subplot(3,2,3);
imagesc(ImgToPlot.b1CP(plotRolArray,plotCowArray,SliceIdx)'); title('B1 Map','FontSize',FontSize);axis off;
CLB{1} = colorbar('FontSize',FontSize);
CLB{1}.YLabel.String = 'Hz';    CLB{1}.YLabel.Rotation = 0;     CLB{1}.YLabel.Position = [0. 0.28 0];
hold on;colormap(Figs{2,1},'gray')
visboundaries(maskedMaps.mask(plotRolArray,plotCowArray,SliceIdx)','Color','r','LineWidth',0.4,...
    'EnhanceVisibility', true,'LineStyle','-');axis off;
nudge(Figs{2,1},[0 0.05 0 0])

Figs{3,1} = subplot(3,2,5); 
imagesc(ImgToPlot.b0(plotRolArray,plotCowArray,SliceIdx)'); title('B0 Map','FontSize',FontSize);axis off;hold on;
CLB{1} = colorbar('FontSize',FontSize);
colormap(Figs{3,1},'gray')
visboundaries(maskedMaps.mask(plotRolArray,plotCowArray,SliceIdx)','Color','r','LineWidth',0.4,...
    'EnhanceVisibility', true,'LineStyle','-');axis off;

Figs{1,1} = subplot(3,2,1); 
imagesc(maskedMaps.localiser(plotRolArray,plotCowArray,SliceIdx)'); title('TOF','FontSize',FontSize);axis off;colormap(Figs{1,1},'gray')
hold on
visboundaries(maskedMaps.mask(plotRolArray,plotCowArray,SliceIdx)','Color','r','LineWidth',0.4,...
    'EnhanceVisibility', true,'LineStyle','-');axis off;
hp3 = get(subplot(3,2,1),'Position');
set(Figs{1,1},'Position',[hp3(1) hp3(2) Figs{3,1}.Position(3) hp3(4)]);

HoriOffset = 0.05;
Figs{1,2} = subplot(3,4,3);
imagesc(FullImage.CP(plotRolArray,plotCowArray,2)',PlotFALim); title('CP Mode','FontSize',FontSize);axis off;
set(Figs{1,2},'Position',[Figs{1,2}.Position(1)-HoriOffset Figs{1,2}.Position(2) Figs{1,2}.Position(3) Figs{1,2}.Position(4)]);
colormap(Figs{1,2},'hot')

Figs{2,2} = subplot(3,4,7);
imagesc(FullImage.PhaseOnly(plotRolArray,plotCowArray,2)',PlotFALim); ylabel('Phase Only','FontSize',FontSize);axis off;
set(Figs{2,2},'Position',[Figs{2,2}.Position(1)-HoriOffset Figs{2,2}.Position(2) Figs{2,2}.Position(3) Figs{2,2}.Position(4)]);
colormap(Figs{2,2},'hot')
nudge(Figs{2,2},[0 0.05 0 0])
Figs{3,2} = subplot(3,4,11);
imagesc(FullImage.AS(plotRolArray,plotCowArray,2)',PlotFALim); ylabel('Full RF Shimming','FontSize',FontSize);axis off;
set(Figs{3,2},'Position',[Figs{3,2}.Position(1)-HoriOffset Figs{3,2}.Position(2) Figs{3,2}.Position(3) Figs{3,2}.Position(4)]);
colormap(Figs{3,2},'hot')

% = get(Figs{3,2},'Position');

%CLB{3} = colorbar('Position', [hp4(1)+hp4(3)+0.015  hp4(2)  0.02  0.815],'FontSize',FontSize);
%CLB{3}.YLabel.String = 'Hz';    CLB{3}.YLabel.Rotation = 0;     CLB{3}.YLabel.Position = [1 1 0];

PlotDiffLim = [-20 20];
Figs{1,3} = subplot(3,4,4);
imagesc(DiffImage.CP(plotRolArray,plotCowArray,2)',PlotDiffLim); 
title('CP Mode','FontSize',FontSize);axis off;
%set(Figs{2,2},'Position',[Figs{2,2}.Position(1) Figs{2,2}.Position(2) Figs{1,1}.Position(3) Figs{2,2}.Position(4)]);
colormap(Figs{1,3},'parula')

Figs{2,3} = subplot(3,4,8);
imagesc(DiffImage.PhaseOnly(plotRolArray,plotCowArray,2)',PlotDiffLim); 
title('Phase Only','FontSize',FontSize);axis off;
%set(Figs{3,2},'Position',[Figs{3,2}.Position(1) Figs{3,2}.Position(2) Figs{1,1}.Position(3) Figs{3,2}.Position(4)]);
colormap(Figs{2,3},'parula')

Figs{3,3} = subplot(3,4,12);
imagesc(DiffImage.AS(plotRolArray,plotCowArray,2)',PlotDiffLim); 
title('Full RF Shimming','FontSize',FontSize);axis off;
%set(Figs{1,2},'Position',[Figs{1,2}.Position(1) Figs{1,2}.Position(2) Figs{1,1}.Position(3) Figs{1,2}.Position(4)]);
colormap(Figs{3,3},'parula')

hp4 = get(Figs{3,3},'Position');
%colorbar('Position', [hp4(1)+hp4(3)+0.015  hp4(2)  0.02  0.815],'FontSize',FontSize);