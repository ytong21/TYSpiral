    addpath('/Users/ytong/Documents/MATLAB/For_James_Larkin');
    cd('/Users/ytong/Documents/MATLAB/tong-acptx/tySpiral/tyRFShimming/')
RFStruct_Example = tyMakeGaussian(600,5);
b = [0 0 RFStruct_Example.RF_pulse 0 0];  % arbitrary unit, peak=1
g = [0 0 6*ones(size(RFStruct_Example.RF_pulse)) 0 0];            % in mT/m
dt = 10e-3;                     % in ms
bmax = 0.5;
gmax = 22;                      % in mT/m
smax = 200;                     % in mT/(m*ms)
dtout = 10e-3;                  % in ms
emax = sum(b.^2)*dt*0.5;        % in units of b*b*dt
[bv,gv] = mintverse(b,g,dt,bmax,gmax,smax,dtout,emax);
ShimArray = 2*rand(1,8);


%%
color_array = [255, 0, 0;255, 127, 0;164, 205, 0;0, 255, 0;0, 0, 255;139, 0, 255;0, 72, 29;123, 79, 201]/255;
MyBox = cell(9,1);
SubPlotsCell = cell(9,3);
b2plot = [0 b 0];
g2plot = [0 2 4 6*ones(size(RFStruct_Example.RF_pulse)) 4 2 0];
FontSizeTmp = 13;
figure(60)
clf
set(gcf,'color','w','InvertHardcopy','off')
set(gcf,'units','centimeters','position',[4 4 13 40],'paperunits','centimeters','paperposition',[0 0 13 40])
set ( 0, 'DefaultFigureColor', [1 0 0] )
for iDx = 1:8
    SubPlotsCell{iDx,1} = subplot(9,3,(iDx-1)*3+1);
    plot(1:numel(b2plot),b2plot,'k','LineWidth',3)
    ylim([0 2]);
    xlim([0 numel(b2plot)])
    pos_temp = SubPlotsCell{iDx,1}.Position;
    xpos = pos_temp(1); ypos = pos_temp(2); 
    %axes('Color','none','YColor','none');
    set(gca,'ytick',[],'Ycolor','w','box','off','xtick',[])
    %ylabel(['RF Chan ',num2str(iDx)])
    %MyBox{iDx} = text(2,1,['RF Chan ',num2str(iDx)], 'clipping', 'off');
    %nudge(MyBox{iDx},[-0.01 0 0 0])
     MyBox{iDx} = axes;
     set(MyBox{iDx},'FontSize',12,'Ycolor','w','Xcolor','w',...
         'Unit','normalized','Position',[xpos-0.003,ypos-0.012,0.01,0.07]);
    MyBox{iDx}.YLabel.String =  ['RF Chan ',num2str(iDx)];  MyBox{iDx}.YLabel.FontSize = FontSizeTmp;
    MyBox{iDx}.YLabel.Color = [0,0,0]; MyBox{iDx}.YLabel.FontWeight = 'Bold';
    %set(MyBox,'String',['RF Chan ',num2str(iDx)],'FontSize',14);
   
    SubPlotsCell{iDx,2} = subplot(9,3,(iDx-1)*3+2);
    plot(1:numel(b2plot),ShimArray(iDx)*b2plot,...
        'color',color_array(iDx,:),'LineWidth',3);
    ylim([0 2]);
    xlim([0 numel(b2plot)])
    set(gca,'ytick',[],'Ycolor','w','box','off','xtick',[])
    
    
    SubPlotsCell{iDx,3} = subplot(9,3,iDx*3);
    plot(1:numel(bv),ShimArray(iDx)*bv,...
        'color',color_array(iDx,:),'LineWidth',3);
    ylim([0 2]);
    xlim([0 numel(bv)])
    set(gca,'ytick',[],'Ycolor','w','box','off','xtick',[])
    %PosTmp = [SubPlotsCell{iDx,3}.Position(1) SubPlotsCell{iDx,3}.Position(2)...
        %SubPlotsCell{iDx,3}.Position(3) SubPlotsCell{iDx,3}.Position(4)*950/600]
    SubPlotsCell{iDx,3}.Position = [SubPlotsCell{iDx,3}.Position(1) SubPlotsCell{iDx,3}.Position(2)...
        SubPlotsCell{iDx,3}.Position(3)*840/600 SubPlotsCell{iDx,3}.Position(4)];
    %pbaspect([1 0.835164835164835 0.835164835164835])
end
PulseType = cell(1,3);

for iDx = 1:3
    PulseType{iDx} = axes;
    set(PulseType{iDx},'FontSize',12,'Ycolor','w','Xcolor','w',...
         'Unit','normalized','Position',[SubPlotsCell{1,iDx}.Position(1)+0.07 SubPlotsCell{1,iDx}.Position(2)+0.1...
        0.07 0.005]);
    switch iDx 
        case 1 
            PulseType{iDx}.XLabel.String = 'CP Mode';
        case 2
            PulseType{iDx}.XLabel.String = 'Gaussian Shimmed';
        case 3
            PulseType{iDx}.XLabel.String = 'VERSE Shimmed';
    end
    PulseType{iDx}.XLabel.FontSize = FontSizeTmp;
    PulseType{iDx}.XLabel.Color = [0,0,0]; PulseType{iDx}.XLabel.FontWeight = 'Bold';
    
    SubPlotsCell{9,iDx} = subplot(9,3,24+iDx);
    switch iDx 
        case 1
            plot(1:numel(g2plot),g2plot,'color','k','LineWidth',3);xlim([0 numel(g2plot)])
        case 2
            plot(1:numel(g2plot),g2plot,'color','k','LineWidth',3);xlim([0 numel(g2plot)])
        case 3
            plot(1:numel(gv),gv,'color','k','LineWidth',3);xlim([0 numel(gv)])
    end
        ylim([0 12]);
        set(gca,'ytick',[],'Ycolor','w','box','off','xtick',[])
        
end
 SubPlotsCell{9,3}.Position = [SubPlotsCell{9,3}.Position(1) SubPlotsCell{9,3}.Position(2)...
        SubPlotsCell{9,3}.Position(3)*840/600 SubPlotsCell{9,3}.Position(4)];
nudge(PulseType{3},[0.045 0 0 0])
    

     MyBox{9} = axes;
     set(MyBox{9},'FontSize',12,'Ycolor','w','Xcolor','w','Unit','normalized',...
         'Position',[SubPlotsCell{9,1}.Position(1)-0.003,SubPlotsCell{9,1}.Position(1)-0.03,0.00001,0.07]);
    MyBox{9}.YLabel.String =  'Gradient';  MyBox{9}.YLabel.FontSize = FontSizeTmp;
    MyBox{9}.YLabel.Color = [0,0,0]; MyBox{9}.YLabel.FontWeight = 'Bold';