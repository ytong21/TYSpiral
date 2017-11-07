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
PlotFALim = [0 30]; FontSize = 8;   CLB = cell(7,1); 
Figs = cell(3,3);

Figs{2,1} = subplot(3,2,3);
imagesc(ImgToPlot.b1CP(plotRolArray,plotCowArray,SliceIdx)'); title('B1 Map','FontSize',FontSize);axis off;
CLB{1} = colorbar('FontSize',FontSize);
CLB{1}.YLabel.String = 'Hz';    CLB{1}.YLabel.Rotation = 0;     
CLB{1}.YLabel.Position = [0.5 280 0];
hold on;colormap(Figs{2,1},'gray')
visboundaries(maskedMaps.mask(plotRolArray,plotCowArray,SliceIdx)','Color','r','LineWidth',0.4,...
    'EnhanceVisibility', true,'LineStyle','-');axis off;
nudge(Figs{2,1},[0 0.05 0 0])

Figs{3,1} = subplot(3,2,5); 
imagesc(ImgToPlot.b0(plotRolArray,plotCowArray,SliceIdx)',[-200 200]); title('B0 Map','FontSize',FontSize);axis off;hold on;
CLB{2} = colorbar('FontSize',FontSize);
CLB{2}.YLabel.String = 'Hz/V';    CLB{2}.YLabel.Rotation = 0;     
CLB{2}.YLabel.Position = [0.5 700 0];
colormap(Figs{3,1},'gray')
visboundaries(maskedMaps.mask(plotRolArray,plotCowArray,SliceIdx)','Color','r','LineWidth',0.4,...
    'EnhanceVisibility', true,'LineStyle','-');axis off;
nudge(Figs{3,1},[0 0.1 0 0])

Figs{1,1} = subplot(3,2,1); 
imagesc(maskedMaps.localiser(plotRolArray,plotCowArray,SliceIdx)'); title('TOF','FontSize',FontSize);axis off;colormap(Figs{1,1},'gray')
hold on
visboundaries(maskedMaps.mask(plotRolArray,plotCowArray,SliceIdx)','Color','r','LineWidth',0.4,...
    'EnhanceVisibility', true,'LineStyle','-');axis off;
hp3 = get(subplot(3,2,1),'Position');
set(Figs{1,1},'Position',[hp3(1) hp3(2) Figs{3,1}.Position(3) hp3(4)]);

HoriOffset = 0.038;
Figs{1,2} = subplot(3,4,3);
imagesc(FullImage.CP(plotRolArray,plotCowArray,2)',PlotFALim); axis off;
title('Flip angle maps','FontSize',FontSize);
ylabel('CP Mode','FontSize',FontSize,'FontWeight','bold');
Figs{1,2}.YLabel.Visible = 'on';
set(Figs{1,2},'Position',[Figs{1,2}.Position(1) Figs{1,2}.Position(2) Figs{1,2}.Position(3) Figs{1,2}.Position(4)]);
colormap(Figs{1,2},'hot')

Figs{2,2} = subplot(3,4,7);
imagesc(FullImage.PhaseOnly(plotRolArray,plotCowArray,2)',PlotFALim);axis off;
ylabel('Phase Only','FontSize',FontSize,'FontWeight','bold');
Figs{2,2}.YLabel.Visible = 'on';
set(Figs{2,2},'Position',[Figs{2,2}.Position(1) Figs{2,2}.Position(2) Figs{2,2}.Position(3) Figs{2,2}.Position(4)]);
colormap(Figs{2,2},'hot')
nudge(Figs{2,2},[0 0.05 0 0])

Figs{3,2} = subplot(3,4,11);
imagesc(FullImage.AS(plotRolArray,plotCowArray,2)',PlotFALim);axis off;
set(Figs{3,2},'Position',[Figs{3,2}.Position(1) Figs{3,2}.Position(2) Figs{3,2}.Position(3) Figs{3,2}.Position(4)]);
ylabel('Full RF Shimming','FontSize',FontSize,'FontWeight','bold');
Figs{3,2}.YLabel.Visible = 'on';
CLB{6} = colorbar('southoutside','FontSize',FontSize);
CLB{6}.YLabel.String = 'Degrees (°)';    CLB{6}.YLabel.Rotation = 0;   
colormap(Figs{3,2},'hot')
nudge(Figs{3,2},[0 0.1 0 0])
nudge(CLB{6},[0 0.02 0 0])

PlotDiffLim = [-20 20];
Figs{1,3} = subplot(3,4,4);
imagesc(DiffImage.CP(plotRolArray,plotCowArray,2)',PlotDiffLim); 
title('Difference images','FontSize',FontSize);axis off;
set(Figs{1,3},'Position',[Figs{1,3}.Position(1)-HoriOffset Figs{1,3}.Position(2) ...
    Figs{1,3}.Position(3) Figs{1,3}.Position(4)]);
colormap(Figs{1,3},'parula')

Figs{2,3} = subplot(3,4,8);
imagesc(DiffImage.PhaseOnly(plotRolArray,plotCowArray,2)',PlotDiffLim); axis off;
set(Figs{2,3},'Position',[Figs{2,3}.Position(1)-HoriOffset Figs{2,3}.Position(2) ...
    Figs{2,3}.Position(3) Figs{2,3}.Position(4)]);
colormap(Figs{2,3},'parula')
nudge(Figs{2,3},[0 0.05 0 0])


Figs{3,3} = subplot(3,4,12);
imagesc(DiffImage.AS(plotRolArray,plotCowArray,2)',PlotDiffLim);axis off
set(Figs{3,3},'Position',[Figs{3,3}.Position(1)-HoriOffset Figs{3,3}.Position(2) ...
    Figs{3,3}.Position(3) Figs{3,3}.Position(4)]);
CLB{7} = colorbar('southoutside');
CLB{7} = colorbar('southoutside','FontSize',FontSize);
colormap(Figs{3,3},'parula')
nudge(Figs{3,3},[0 0.1 0 0])
nudge(CLB{7},[0 0.02 0 0])

CLB{7}.YLabel.String = 'Percentage error (%)';    CLB{7}.YLabel.Rotation = 0;     
%CLB{7}.YLabel.Position = [0.5 700 0];

hp4 = get(Figs{3,3},'Position');
%colorbar('Position', [hp4(1)+hp4(3)+0.015  hp4(2)  0.02  0.815],'FontSize',FontSize);