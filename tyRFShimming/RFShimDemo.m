RFStruct_Example = tyMakeGaussian(600,5);

ShimArray = 2*rand(1,8);
%%
color_array = [255, 0, 0;255, 127, 0;164, 205, 0;0, 255, 0;0, 0, 255;139, 0, 255;0, 72, 29;123, 79, 201]/255;

figure(60)
for iDx = 1:8
    subplot(8,2,(iDx-1)*2+1)
    plot(1:numel(RFStruct_Example.RF_pulse),RFStruct_Example.RF_pulse,'k','LineWidth',3)
    ylim([0 2]);
    xlim([0 numel(RFStruct_Example.RF_pulse)])
    %axes('Color','none','YColor','none');
    set(gca,'ytick',[],'Ycolor','w','box','off','xtick',[])
    %set(gca,'xtick',[],'ytick',[])
    subplot(8,2,iDx*2)
    plot(1:numel(RFStruct_Example.RF_pulse),ShimArray(iDx)*RFStruct_Example.RF_pulse,...
        'color',color_array(iDx,:),'LineWidth',3);
    ylim([0 2]);
    xlim([0 numel(RFStruct_Example.RF_pulse)])
    set(gca,'ytick',[],'Ycolor','w','box','off','xtick',[])

end