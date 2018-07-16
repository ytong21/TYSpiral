function PlotPerfusion(dt,CPSeriesNumber,pTxSeriesNumber)

matched_ptx = dt.search('target','instance','query',@(inst,ser,stu) ser.SeriesNumber==pTxSeriesNumber);
matched_CP = dt.search('target','instance','query',@(inst,ser,stu) ser.SeriesNumber==CPSeriesNumber);
if isempty(matched_ptx) || isempty(matched_CP)
    disp('Please check series number. Not images found');
else
    CP_Img = Spectro.dicomImage({matched_CP(1).Filename});
    pTx_Img = Spectro.dicomImage({matched_ptx(1).Filename});
    AllPoints = nonzeros([pTx_Img.image(:);CP_Img.image(:)]);
    CPtoPlot = (max(AllPoints) - CP_Img.image)';
    pTxtoPlot = (max(AllPoints) - pTx_Img.image)';
    pTxtoPlot(pTxtoPlot == 4095) = 0;
    CPtoPlot(CPtoPlot == 4095) = 0;
    CutOff = prctile(max(AllPoints)-AllPoints,99);
    figure(99)
    set(gcf,'color','w','InvertHardcopy','off')
    set(gcf,'units','centimeters','position',[4 4 40 20],'paperunits','centimeters','paperposition',[0 0 40 20])
    clf
    subplot(1,2,1)
    imagesc(CPtoPlot,[1900 CutOff]);colormap('gray');axis off;title('CP mode','FontSize',16)
    subplot(1,2,2)
    imagesc(pTxtoPlot,[1900 CutOff]);colormap('gray');axis off;title('RF shimming','FontSize',16)

end