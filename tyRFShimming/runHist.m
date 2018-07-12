figure(104)
set(gcf,'color','w','InvertHardcopy','off')
set(gcf,'units','centimeters','position',[4 4 40 20],'paperunits','centimeters','paperposition',[0 0 40 20])
%histogram(abs(FAFinal.AS),20,'Normalization','probability');hold on;
%[f1,xi] = ksdensity(abs(FAFinal.AS),linspace(17,22,100));
%plot(xi,f1,'b');

% Histo{2,1} = histfit(abs(FAFinal.CP),10,'kernel'); delete(Histo{2,1}(1));  Histo{2,1}(2).Color = 'c';  hold on
% Histo{3,1} = histfit(abs(FAFinal.PhaseOnly),10,'kernel');  delete(Histo{3,1}(1));  Histo{3,1}(2).Color = 'b';
% Histo{1,1} = histfit(abs(FAFinal.AS),10,'kernel');  delete(Histo{1,1}(1));  Histo{1,1}(2).Color = 'r'; 

edges = 0:0.1:20;
subplot(1,2,1)
Histo2{1,1} = histogram(abs(HistToPlot{2}.CP),edges,'FaceColor','b');    Histo2{1,1}.Normalization = 'probability'; hold on
Histo2{2,1} = histogram(abs(HistToPlot{2}.PhaseOnly),edges,'FaceColor','g');    Histo2{2,1}.Normalization = 'probability';
Histo2{3,1} = histogram(abs(HistToPlot{2}.AS),edges,'FaceColor','y');    Histo2{3,1}.Normalization = 'probability';
box off
xlabel('Flip angle (°)')
ylabel('Probability'); ylim([0 0.06])
title('Subject 1: flip angle achieved in all vessels');
lgd = legend('CP mode','Phase only shimming','Full RF shimming');legend('boxoff')   
lgd.FontSize = 13;set(gca, 'FontSize', 16)
subplot(1,2,2)
Histo2{1,2} = histogram(abs(HistToPlot{3}.CP),edges,'FaceColor','b');    Histo2{1,2}.Normalization = 'probability'; hold on
Histo2{2,2} = histogram(abs(HistToPlot{3}.PhaseOnly),edges,'FaceColor','g');    Histo2{2,2}.Normalization = 'probability';
Histo2{3,2} = histogram(abs(HistToPlot{3}.AS),edges,'FaceColor','y');    Histo2{3,2}.Normalization = 'probability';
box off
xlabel('Flip angle (°)')
ylabel('Probability'); ylim([0 0.06])
title('Subject 2: flip angle achieved in all vessels');
lgd = legend('CP mode','Phase only shimming','Full RF shimming');legend('boxoff')   
lgd.FontSize = 13;set(gca, 'FontSize', 16)

%%
figure(105)
set(gcf,'color','w','InvertHardcopy','off')
set(gcf,'units','centimeters','position',[4 4 40 20],'paperunits','centimeters','paperposition',[0 0 40 20])

subplot(1,2,1)
  Mode = {'CP Mode','Phase only','Full RF shimming'};   VesselName = {'RICA','RVA','LICA','LVA','Total'};
  BarChart = bar(1:5,BarMtxToPlot{1});  lgdTmp = legend(Mode); box off;  ylim([0 0.8]);  ylabel('Labelling efficiency');  
  lgdTmp.FontSize = 13; set(gca,'XTickLabel',VesselName);set(gca, 'FontSize', 16);legend('boxoff')  
  %pause(0.1);
  hold on
  for iDx = 1:numel(BarChart)
    xData = BarChart(iDx).XData+BarChart(iDx).XOffset;
    errorbar(xData,BarMtxToPlot{1}(:,iDx),StdMtxToPlot{1}(:,iDx),'r.')
  end
  BarChart(1).FaceColor = 'b'; BarChart(2).FaceColor = 'g'; BarChart(3).FaceColor = 'y';
  title('Subject 1: labeling efficiency')
  hold off
  subplot(1,2,2)
  BarChart = bar(1:5,BarMtxToPlot{2});  lgdTmp = legend(Mode); box off;  ylim([0 0.8]);  ylabel('Labelling efficiency');  
  lgdTmp.FontSize = 13; set(gca,'XTickLabel',VesselName);set(gca, 'FontSize', 16);legend('boxoff')  
  %pause(0.1);
  hold on
  for iDx = 1:numel(BarChart)
    xData = BarChart(iDx).XData+BarChart(iDx).XOffset;
    errorbar(xData,BarMtxToPlot{2}(:,iDx),StdMtxToPlot{2}(:,iDx),'r.')
  end
  BarChart(1).FaceColor = 'b'; BarChart(2).FaceColor = 'g'; BarChart(3).FaceColor = 'y';
  title('Subject 2: labeling efficiency')
  hold off