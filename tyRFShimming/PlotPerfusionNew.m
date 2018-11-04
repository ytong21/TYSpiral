function PlotPerfusionNew(dt,SerieNumberArray)

V_Shimmed_SeriesNumber = SerieNumberArray(1);
G_CP_SeriesNumber = SerieNumberArray(3);
G_Shimmed_SeriesNumber = SerieNumberArray(2);

matched_G_Shimmed = dt.search('target','instance','query',@(inst,ser,stu) ser.SeriesNumber==G_Shimmed_SeriesNumber);
matched_G_CP = dt.search('target','instance','query',@(inst,ser,stu) ser.SeriesNumber==G_CP_SeriesNumber);
matched_V_Shimmed = dt.search('target','instance','query',@(inst,ser,stu) ser.SeriesNumber==V_Shimmed_SeriesNumber);
if isempty(matched_G_Shimmed) || isempty(matched_G_CP) || isempty(matched_V_Shimmed)
    disp('Please check series number. Not images found');
else
    G_CP_Img = Spectro.dicomImage({matched_G_CP(1).Filename});
    G_Shimmed_Img = Spectro.dicomImage({matched_G_Shimmed(1).Filename});
    V_Shimmed_Img = Spectro.dicomImage({matched_V_Shimmed(1).Filename});
    AllPoints = nonzeros([G_Shimmed_Img.image(:);G_CP_Img.image(:);V_Shimmed_Img.image(:)]);
    G_CPtoPlot = (max(AllPoints) - G_CP_Img.image)';
    G_ShimmedtoPlot = (max(AllPoints) - G_Shimmed_Img.image)';
    V_ShimmedtoPlot = (max(AllPoints) - V_Shimmed_Img.image)';
    G_ShimmedtoPlot(G_ShimmedtoPlot == 4095) = 0;
    G_CPtoPlot(G_CPtoPlot == 4095) = 0;
    V_ShimmedtoPlot(V_ShimmedtoPlot == 4095) = 0;
    CutOff = prctile(max(AllPoints)-AllPoints,99);
    G_ShimmedtoPlot(331:440,441:550) =CutOff+1;
    G_CPtoPlot(331:440,441:550) =CutOff+1;
    V_ShimmedtoPlot(331:440,441:550) =CutOff+1;
    figure(99)
    set(gcf,'color','w','InvertHardcopy','off')
    Width = 40; Length = 13;LowerBound = 1800;FontSizeTmp = 16;
    set(gcf,'units','centimeters','position',[4 4 Width Length],'paperunits',...
        'centimeters','paperposition',[0 0 Width Length])
    clf
    SubPlotsCell = cell(3,1);
    SubPlotsCell{1} = subplot(1,3,1);
    imagesc(G_CPtoPlot(1:440,:),[LowerBound CutOff]);colormap('gray');axis equal;axis off
    title('Gaussian CP mode','FontSize',FontSizeTmp)
    SubPlotsCell{2} = subplot(1,3,2);
    imagesc(G_ShimmedtoPlot(1:440,:),[LowerBound CutOff]);colormap('gray');axis equal;axis off;
    title('Gaussian Shimmed','FontSize',FontSizeTmp)
    SubPlotsCell{3} = subplot(1,3,3);
    imagesc(V_ShimmedtoPlot(1:440,:),[LowerBound CutOff]);colormap('gray');axis equal;axis off;
    title('VERSE Shimmed','FontSize',FontSizeTmp)
    Amplify = 1.35;
    for iDx = 1:numel(SubPlotsCell)
        SubPlotsCell{iDx}.Position = [SubPlotsCell{iDx}.Position(1)-0.1+0.035*(iDx-1) SubPlotsCell{iDx}.Position(2)-0.2...
        SubPlotsCell{iDx}.Position(3)*Amplify SubPlotsCell{iDx}.Position(4)*Amplify];
        SubPlotsCell{iDx}.Title.Units = 'normalized';
        SubPlotsCell{iDx}.Title.Position(2) = SubPlotsCell{iDx}.Title.Position(2)-0.15;
        %nudge(SubPlotsCell{iDx}.Title,[0 -0.02 0 0])
    end
end