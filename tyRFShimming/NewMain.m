    %%   Loading the data and masking
    addpath('/Users/ytong/Documents/MATLAB/For_James_Larkin');
    cd('/Users/ytong/Documents/MATLAB/tong-acptx/tySpiral/tyRFShimming/')
    %pTxPath = '/Volumes/Data/DICOM/2018-10/20181019_F7T_2013_50_097';
    %pTxPath = '/Volumes/Data/DICOM/2018-10/20181022_F7T_2013_50_098';
    %pTxPath = '/Volumes/Data/DICOM/2018-10/20181026_F7T_2013_50_099';
    %pTxPath = '/Volumes/Data/DICOM/2018-10/20181030_F7T_2013_50_100';
    %pTxPath = '/Users/ytong/Documents/Data/ForISMRMabstract/20181022_F7T_2013_50_098';
    pTxPath = '/Volumes/Data/DICOM/2019-03/20190325_F7T_2013_50_112';
    dt = Spectro.dicomTree('dir',pTxPath,'recursive',false);
    %%
    ptxFMObj = DicomFM.WTCpTxFieldmaps(dt,'B1String','dt_dream_wIce_60deg_160VRef_tagging_32meas_higher__B1',...
        'B0String','fieldmap_ptx7t_iso4mm_trans_RL','LocaliserString','fl_tof','InterpTarget',...
        'Localiser');
%     ptxFMObj = DicomFM.WTCpTxFieldmaps(dt,'B1String','dt_dream_wIce_60deg_150VRef_tagging__B1',...
%         'B0String','fieldmap_ptx7t_iso4mm_trans_RL','LocaliserString','fl_tof','InterpTarget',...
%         'Localiser');  
%     ptxFMObj = DicomFM.WTCpTxFieldmaps(dt,'B1String','dt_dream_wIce_60deg_90VRef__B1',...
%         'B0String','fieldmap_ptx7t_iso4mm_trans_RL','LocaliserString','localizer_3D_moved','InterpTarget',...
%         'Localiser');
    %%
    SliceIdx = 39;  % Decided by visual inspection of the tof images
    % 18 for subject 98
    ptxFMObj.createMask(@(x) DicomFM.maskFunctions.VEellipseMask(x,SliceIdx,[1 1 1 1]),true);  
    
    %% Defineding a verse waveform
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
    bv = bv/max(bv);
     %%  Specifying parameters
    param.targetFlipAngle = 20;
    param.numCh = 8;
    param.TR = 6;% sec
    param.CGtikhonov = 1e-6;
    param.tol = 1e-5;
    param.MaxEvaluation = 40000;
    param.FOX = 25;%25cm
    param.RFsep = 1.2e-3;% sec
    param.tagging_duration = 1; %   sec
    param.num_EPI_slices = 17;
    param.PLD = 1;%sec
    param.RefVol = 80;  % reference voltage. Not applied to tagging pulses
    %param.RFdur = 600e-6;%sec
    %param.RFdur_in_us = param.RFdur*1e6;
    %%
    load Vessel.mat
    
    maskedMaps.b1SensMasked = ptxFMObj.getB1PerV('uT','delete');
    maskedMaps.b0MapMasked = ptxFMObj.getB0('Hz','delete');
    maskedMaps.b0MapMaskedRad = (2*pi)*maskedMaps.b0MapMasked;
    maskedMaps.posVox = ptxFMObj.getPositions('cm','delete').';
    maskedMaps.localiser = ptxFMObj.getLoc('none');
    maskedMaps.b1SensMaskedHz = ptxFMObj.getB1PerV('Hz','delete');
    maskedMaps.mask = ptxFMObj.getMask();
    maskedMaps.TargetMasked = Vessel.TargetMasked;
       
    ImgToPlot.b0 = ptxFMObj.getB0('Hz','none');
    ImgToPlot.b1 = ptxFMObj.getB1PerV('Hz','none');
    ImgToPlot.b1CP = abs(sum(ImgToPlot.b1,4));    

    %RFStruct = tyMakeGaussian(param.RFdur_in_us,4);
       
  figure(96)
  LocToPlot = maskedMaps.localiser(:,:,SliceIdx);
  LocMax = max(LocToPlot(:));
  %[dim1, dim2] = size(LocToPlot);
  %set(gcf,'units','centimeters','position',[4 4 30 30*(dim2/dim1)],'paperunits','centimeters','paperposition',[0 0 40 20])
  for iDx = 1:numel(VesselMask)
          LocTemp = LocToPlot;
          LocTemp(VesselMask{iDx}) = LocMax;
          imagesc(LocTemp,[min(LocToPlot(:)) max(LocToPlot(:))]); axis off;drawnow;
          pause(0.5)
  end
  close(96)
  clear LocToPlot LocMax LocTemp
  %% Running gaussian
     %RFStruct = struct('RF_pulse',bv);
    disp('Calculating system matrix...');
    SysMatMode = 'Full';
    if strcmp(SysMatMode,'Approximation')
        AFull_gaussian = getAMatSimp(RFStruct_Example,maskedMaps.b1SensMasked,maskedMaps.b0MapMasked); %Nv-(2*Nc) sys mtx, rad/V
    elseif strcmp(SysMatMode,'Full')
        gr = zeros(numel(RFStruct_Example.RF_pulse),3);
        %gr = zeros(numel(gv),3);
        %gr(:,3) = gv;
        rfOn = true(size(gr,1),1);
        tp = 10E-6 * ones(size(rfOn)); % in seconds
        DA = genAMatFull(tp,rfOn,gr,maskedMaps.b1SensMaskedHz,...
            maskedMaps.b0MapMasked,maskedMaps.posVox);
     %  To take the RF pulse shape into account. Construct a pseudo-diagonal matrix
        RFDiag = zeros(size(DA,2),8);
        %RFSize = numel(RFStruct.RF_pulse);
        RFSize = numel(RFStruct_Example.RF_pulse);
        for iDx = 1:8
            RFDiag(((1:RFSize)+(iDx-1)*RFSize),iDx) = RFStruct_Example.RF_pulse';
            %RFDiag(((1:RFSize)+(iDx-1)*RFSize),iDx) = bv;
        end
        AFullMag = DA*RFDiag;  %    magnetization /V
        AFull_gaussian = asin(AFullMag);   % rad/V  Dupas paper in rad/V
    end
    %   Loading PulsesFromDSV.mat
    %   A struct containing information of tagging, WET, inversion, fat sat
    %   and imaging pulses.
    if(exist('PulsesFromDSV','var') ~= 1)
        load PulsesFromDSV.mat
    end
    disp('Complete')
    %  Running RF shimming Step 1: Variable-exchange method
    disp('Running RF shimming variable-exchange optimization...')
    tic
    tikhonovArray = power(10,-7:0.5:-3);
    bVE = zeros(8,numel(tikhonovArray));
    for iDx = 1:numel(tikhonovArray)
        param.CGtikhonov = tikhonovArray(iDx);
        [bVE(:,iDx),NRMSETmp,~,~] = runVE(AFull_gaussian,param,maskedMaps);
        disp(NRMSETmp)  
    end
    toc
    disp('Complete')
    %  Running RF shimming Step 2: Active-set method
    disp('Running RF shimming active-set optimization...')
    tic
    bAS = zeros(16,numel(tikhonovArray));
    NRMSE = zeros(size(tikhonovArray));
    output = cell(size(NRMSE));
    exitflag = zeros(size(NRMSE));
    for iDx = 1:numel(tikhonovArray)
        param.CGtikhonov = tikhonovArray(iDx);
        [bAS(:,iDx),NRMSE(iDx),exitflag(iDx),output{iDx}] = runAS(bVE(:,iDx),RFStruct_Example,maskedMaps,param,AFull_gaussian);
    end
    [~, minIndex] = min(NRMSE);
    toc
        bmin_gaussain = bAS(1:8,minIndex).*exp(1i*bAS(9:16,minIndex));

    disp('Complete')
   %%  Running verse
   % Constructing system matrix
   RFStruct = struct('RF_pulse',bv);
    disp('Calculating system matrix...');
    SysMatMode = 'Full';
    if strcmp(SysMatMode,'Approximation')
        AFull_verse = getAMatSimp(RFStruct,maskedMaps.b1SensMasked,maskedMaps.b0MapMasked); %Nv-(2*Nc) sys mtx, rad/V
    elseif strcmp(SysMatMode,'Full')
        %gr = zeros(numel(RFStruct.RF_pulse),3);
        gr = zeros(numel(gv),3);
        gr(:,3) = gv;
        rfOn = true(size(gr,1),1);
        tp = 10E-6 * ones(size(rfOn)); % in seconds
        DA = genAMatFull(tp,rfOn,gr,maskedMaps.b1SensMaskedHz,...
            maskedMaps.b0MapMasked,maskedMaps.posVox);
     %  To take the RF pulse shape into account. Construct a pseudo-diagonal matrix
        RFDiag = zeros(size(DA,2),8);
        %RFSize = numel(RFStruct.RF_pulse);
        RFSize = numel(bv);
        for iDx = 1:8
            %RFDiag(((1:RFSize)+(iDx-1)*RFSize),iDx) = RFStruct.RF_pulse';
            RFDiag(((1:RFSize)+(iDx-1)*RFSize),iDx) = bv;
        end
        AFullMag = DA*RFDiag;  %    magnetization /V
        AFull_verse = asin(AFullMag);   % rad/V  Dupas paper in rad/V
    end
    %   Loading PulsesFromDSV.mat
    %   A struct containing information of tagging, WET, inversion, fat sat
    %   and imaging pulses.
    if(exist('PulsesFromDSV','var') ~= 1)
        load PulsesFromDSV.mat
    end
    disp('Complete')
    %  Running RF shimming Step 1: Variable-exchange method
    disp('Running RF shimming variable-exchange optimization...')
    tic
    tikhonovArray = power(10,-7:0.5:-3);
    bVE = zeros(8,numel(tikhonovArray));
    for iDx = 1:numel(tikhonovArray)
        param.CGtikhonov = tikhonovArray(iDx);
        [bVE(:,iDx),NRMSETmp,~,~] = runVE(AFull_verse,param,maskedMaps);
        disp(NRMSETmp)  
    end
    toc
    disp('Complete')
    %  Running RF shimming Step 2: Active-set method
    disp('Running RF shimming active-set optimization...')
    tic
    bAS = zeros(16,numel(tikhonovArray));
    NRMSE = zeros(size(tikhonovArray));
    output = cell(size(NRMSE));
    exitflag = zeros(size(NRMSE));
    for iDx = 1:numel(tikhonovArray)
        param.CGtikhonov = tikhonovArray(iDx);
        [bAS(:,iDx),NRMSE(iDx),exitflag(iDx),output{iDx}] = runAS(bVE(:,iDx),RFStruct,maskedMaps,param,AFull_verse);
    end
    [~, minIndex] = min(NRMSE);
    toc
        bmin_verse = bAS(1:8,minIndex).*exp(1i*bAS(9:16,minIndex));

    disp('Complete')
     %%  Bloch Simulation
    disp('Running Bloch simulation...')
    %bmin = bVE(:,minIndex);
    %RFToSim = bsxfun(@times,repmat(RFStruct.RF_pulse.',1,8),bmin.');    % Note the diff between .' and '
    RFToSim_gaussian = bsxfun(@times,repmat(RFStruct_Example.RF_pulse.',1,8),bmin_gaussain.');
    GToSim = zeros(numel(RFToSim_gaussian)/8,3);
    %GToSim(:,3) = gv;
    %GToSim(:,3) = 425.8*RFStruct.G_amp*ones(numel(RFToSim)/8,1);     %from mT/m to Hz/cm
    B1ToSim = ptxFMObj.getB1PerV('Hz','delete');
    magnetization_gaussian = make_blochSim(RFToSim_gaussian,B1ToSim,...
        maskedMaps.b0MapMasked,GToSim,10E-6,maskedMaps.posVox,maskedMaps.mask);
    disp('Bloch simulation complete')
    
    RFToSim_verse = bsxfun(@times,repmat(bv,1,8),bmin_verse.');
    GToSim = zeros(numel(RFToSim_verse)/8,3);
    GToSim(:,3) = gv;
    %GToSim(:,3) = 425.8*RFStruct.G_amp*ones(numel(RFToSim)/8,1);     %from mT/m to Hz/cm
    B1ToSim = ptxFMObj.getB1PerV('Hz','delete');
    magnetization_verse = make_blochSim(RFToSim_verse,B1ToSim,...
        maskedMaps.b0MapMasked,GToSim,10E-6,maskedMaps.posVox,maskedMaps.mask);
    disp('Bloch simulation complete')
    
    
    %  Finding out what CP mode can do
    %bCP = runCP(AFull,param,RFStruct);
    bCP = ones(8,1)*sqrt(sum(abs(bmin_gaussain).^2)/8);
    %  Running phase only optimization
    tic
        [bPhaseTmp, ErrorOut] = runPhaseOnly(AFull_gaussian,param,RFStruct);
    toc
    [~, minPhaseIndex] = min(ErrorOut); 
    bPhase = bPhaseTmp(1,minPhaseIndex)*exp(1i*bPhaseTmp(2:9,minPhaseIndex));

    %  Calculate magnetization and NRMSE
  rmse = @(x,xref) sqrt(immse(x,xref));
  nrmse = @(x,xref) rmse(x,xref)/mean(x);
  % FA results in degrees in a masked vector form.
  FAFinal = struct;
  FAFinal.CP = rad2deg(AFull_gaussian*bCP);
  FAFinal.AS = rad2deg(AFull_gaussian*bmin_gaussain);
  FAFinal.PhaseOnly = rad2deg(AFull_gaussian*bPhase);

  FullImage = struct;
  % JawSlicesMask was originally designed to contain 3 slices (selected slice plus 2 adjacent slices)
  % Attempting to reduce the design to one slice only. 
  JawSlicesMask = maskedMaps.mask(:,:,SliceIdx);
  FullImage.CP = zeros(size(JawSlicesMask));
  FullImage.CP(JawSlicesMask) = abs(FAFinal.CP);
  FullImage.AS = zeros(size(JawSlicesMask));
  FullImage.AS(JawSlicesMask) = abs(FAFinal.AS);  
  FullImage.ASBloch = asind(magnetization_gaussian.mxy);
  FullImage.PhaseOnly = zeros(size(JawSlicesMask));
  FullImage.PhaseOnly((JawSlicesMask)) = abs(FAFinal.PhaseOnly);
  FullImage.verse = zeros(size(JawSlicesMask));
  FullImage.verse(JawSlicesMask) = rad2deg(abs(AFull_verse*bmin_verse));
  
  DiffImage.CP = NaN(size(JawSlicesMask));
  DiffImage.CP(JawSlicesMask) = 100*(abs(FAFinal.CP) - ones(size(FAFinal.CP))*param.targetFlipAngle)/param.targetFlipAngle;
  DiffImage.AS = NaN(size(JawSlicesMask));
  DiffImage.AS(JawSlicesMask) = 100*(abs(FAFinal.AS) - ones(size(FAFinal.AS))*param.targetFlipAngle)/param.targetFlipAngle;  
  DiffImage.ASBloch = asind(magnetization_gaussian.mxy);
  DiffImage.PhaseOnly = NaN(size(JawSlicesMask));
  DiffImage.PhaseOnly((JawSlicesMask)) = 100*(abs(FAFinal.PhaseOnly) - ones(size(FAFinal.PhaseOnly))*param.targetFlipAngle)/param.targetFlipAngle;
  
  Error = struct;
  Error.CP = nrmse(abs(FAFinal.CP),ones(size(FAFinal.CP))*param.targetFlipAngle);
  Error.AS = nrmse(abs(FAFinal.AS),ones(size(FAFinal.CP))*param.targetFlipAngle);
  Error.PhaseOnly = nrmse(abs(FAFinal.PhaseOnly),ones(size(FAFinal.CP))*param.targetFlipAngle); 

  OneSlice = FullImage;
  %OneSlice.CP = FullImage.CP(:,:,2);      OneSlice.PhaseOnly = FullImage.PhaseOnly(:,:,2);  OneSlice.Full = FullImage.AS(:,:,2);
  RICA{1} = OneSlice.CP(VesselMask{1});   RICA{2} = OneSlice.PhaseOnly(VesselMask{1});  RICA{3} = OneSlice.AS(VesselMask{1});
  RVA{1} = OneSlice.CP(VesselMask{2});   RVA{2} = OneSlice.PhaseOnly(VesselMask{2});  RVA{3} = OneSlice.AS(VesselMask{2});
  LICA{1} = OneSlice.CP(VesselMask{3});   LICA{2} = OneSlice.PhaseOnly(VesselMask{3});  LICA{3} = OneSlice.AS(VesselMask{3});
  LVA{1} = OneSlice.CP(VesselMask{4});   LVA{2} = OneSlice.PhaseOnly(VesselMask{4});  LVA{3} = OneSlice.AS(VesselMask{4});
  Full{1} = abs(FAFinal.CP);       Full{2} = abs(FAFinal.PhaseOnly);  Full{3} = abs(FAFinal.AS);
    %% Plotting
    % Plotting
  if ~exist('plotRolArray','var')
      figure(56)
      imagesc(maskedMaps.localiser(:,:,SliceIdx));
      title('Select Rect Regions to Plot')
      h2plot = imrect;  wait(h2plot);
      plotMask = logical(createMask(h2plot));
      [row,col] = find(plotMask);
      plotRolArray = row(1):row(end);
      plotCowArray = col(1):col(end);
      close(56)
  end
    figure(92)
set(gcf,'color','w','InvertHardcopy','off')
set(gcf,'units','centimeters','position',[4 4 40 20],'paperunits','centimeters','paperposition',[0 0 40 20])
clf
PlotFALim = [0 20]; FontSize = 16;   CLB = cell(10,1);   CLBFontSize = 14;
Figs = cell(4,3);

Figs{2,1} = subplot(3,2,3);
imagesc(ImgToPlot.b1CP(plotRolArray,plotCowArray,SliceIdx)',[0 5]); 
ylabel('DREAM B1 Map','FontSize',FontSize,'FontWeight','bold');axis off;
Figs{2,1}.YLabel.Visible = 'on';
CLB{1} = colorbar('FontSize',CLBFontSize);
CLB{1}.YLabel.String = 'Hz/V';    CLB{1}.YLabel.Rotation = 0;     
CLB{1}.YLabel.Position = [0.5 5.7 0];
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
CLB{2}.YLabel.Position = [0.5 252 0];
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
imagesc(FullImage.CP(plotRolArray,plotCowArray)',PlotFALim); axis off;
title('Flip angle maps','FontSize',FontSize);
ylabel('Gaussian CP Mode','FontSize',FontSize,'FontWeight','bold');
Figs{1,2}.YLabel.Visible = 'on';
set(Figs{1,2},'Position',[Figs{1,2}.Position(1) Figs{1,2}.Position(2) Figs{1,2}.Position(3) Figs{1,2}.Position(4)]);
colormap(Figs{1,2},'hot')

Figs{2,2} = subplot(3,2,4);
imagesc(FullImage.verse(plotRolArray,plotCowArray)',PlotFALim);axis off;
ylabel('VERSE Shimmed','FontSize',FontSize,'FontWeight','bold');
Figs{2,2}.YLabel.Visible = 'on';
set(Figs{2,2},'Position',[Figs{2,2}.Position(1) Figs{2,2}.Position(2) Figs{2,2}.Position(3) Figs{2,2}.Position(4)]);
colormap(Figs{2,2},'hot')
nudge(Figs{2,2},[0 0.05 0 0]);

Figs{3,2} = subplot(3,2,6);
imagesc(FullImage.AS(plotRolArray,plotCowArray)',PlotFALim);axis off;
set(Figs{3,2},'Position',[Figs{3,2}.Position(1) Figs{3,2}.Position(2) Figs{3,2}.Position(3) Figs{3,2}.Position(4)]);
ylabel('Gaussian Shimmed','FontSize',FontSize,'FontWeight','bold');
Figs{3,2}.YLabel.Visible = 'on';
%CLB{6} = colorbar('southoutside','FontSize',CLBFontSize);
%CLB{6}.YLabel.String = 'Degrees (°)';    CLB{6}.YLabel.Rotation = 0;   
%colormap(Figs{3,2},'hot')
nudge(Figs{3,2},[0 0.1 0 0])
nudge(CLB{6},[0 0.02 0 0])
%


% Chopping the middle sections of the FA plots
%HoriOffset = 0.038;
 ToChop = 12;
      %plotCowArrayReduced = [col(1):((col(1)+col(end))/2-ToChop) (ceil((col(1)+col(end))/2)+ToChop):col(end)];
      plotRowArrayReduced = [row(1):(((row(1)+row(end))/2)-ToChop) (ceil((row(1)+row(end))/2)+ToChop):row(end)];
 FA_map_positions = cell(3,1);
for iDx = 1:3
    FA_map_positions{iDx} = get(Figs{iDx,2},'position');
    FA_map_positions{iDx}(3) = FA_map_positions{iDx}(3)*(numel(plotRowArrayReduced)/numel(plotRolArray));
    delete(Figs{iDx,2});
end
Figs{1,3} = subplot('Position',FA_map_positions{1});
imagesc(FullImage.CP(plotRowArrayReduced,plotCowArray)',PlotFALim); axis off;
title('Flip angle maps','FontSize',FontSize);
ylabel({'Gaussian', 'CP Mode'},'FontSize',FontSize,'FontWeight','bold');
Figs{1,3}.YLabel.Visible = 'on';
%set(Figs{1,3},'Position',[Figs{1,3}.Position(1) Figs{1,3}.Position(2) Figs{1,3}.Position(3) Figs{1,3}.Position(4)]);
colormap(Figs{1,3},'hot')

%FA_map_positions{3}(4) = FA_map_positions{1}(4);
%FA_map_positions{3}
Figs{3,3} = subplot('Position',FA_map_positions{3});
imagesc(FullImage.verse(plotRowArrayReduced,plotCowArray)',PlotFALim);axis off;
ylabel({'VERSE','Shimmed'},'FontSize',FontSize,'FontWeight','bold');
Figs{3,3}.YLabel.Visible = 'on';
%set(Figs{2,3},'Position',[Figs{2,3}.Position(1) Figs{2,3}.Position(2) Figs{2,3}.Position(3) Figs{2,3}.Position(4)]);
colormap(Figs{3,3},'hot')
%nudge(Figs{2,3},[0 0.05 0 0]);

Figs{2,3} = subplot('Position',FA_map_positions{2});
imagesc(FullImage.AS(plotRowArrayReduced,plotCowArray)',PlotFALim);axis off;
%set(Figs{3,3},'Position',[Figs{3,3}.Position(1) Figs{3,3}.Position(2) Figs{3,3}.Position(3) Figs{3,3}.Position(4)]);
ylabel({'Gaussian','Shimmed'},'FontSize',FontSize,'FontWeight','bold');
Figs{2,3}.YLabel.Visible = 'on';
CLB{9} = colorbar('eastoutside','FontSize',CLBFontSize);
CLB{9}.YLabel.String = 'Degrees (°)';    %CLB{9}.YLabel.Rotation = 0;   
colormap(Figs{2,3},'hot')
%nudge(Figs{3,3},[0 0.1 0 0])

set(CLB{9},'Position',[ 0.878059964726631         0.209876543209877         0.017636684303351         0.715437545388525]);
%nudge(CLB{9},[-0.04 0 0 0])

  %% Writing the pulse
  makeRFToWrite = @(b) [abs(b)/max(abs(b)), angle(b)];
  vCP = sqrt(sum(abs(bmin_verse).^2)/8)*ones(size(bCP));
  
  
  ToWrite.Full = makeRFToWrite(bmin_gaussain);
  ToWrite.CP = makeRFToWrite(bCP);
  ToWrite.PhaseOnly = makeRFToWrite(bPhase);
  ToWrite.VERSE = makeRFToWrite(bmin_verse);
  ToWrite.vCP = makeRFToWrite(vCP);
  for iDx = 1:8
      if ToWrite.Full(iDx,2) < 0
          ToWrite.Full(iDx,2) = ToWrite.Full(iDx,2) + (2*pi);
      elseif ToWrite.PhaseOnly(iDx,2) < 0
          ToWrite.PhaseOnly(iDx,2) = ToWrite.PhaseOnly(iDx,2) + (2*pi);           
      elseif ToWrite.CP(iDx,2) < 0
          ToWrite.CP(iDx,2) = ToWrite.CP(iDx,2) + (2*pi);
      elseif ToWrite.VERSE(iDx,2) < 0
          ToWrite.VERSE(iDx,2) = ToWrite.VERSE(iDx,2) + (2*pi);
      elseif ToWrite.vCP(iDx,2) < 0
          ToWrite.vCP(iDx,2) = ToWrite.vCP(iDx,2) + (2*pi);          
      end
  end
  RFAmp.Full = max(abs(bmin_gaussain));   RFAmp.CP = max(abs(bCP));    RFAmp.PhaseOnly = max(abs(bPhase));
  RFAmp.VERSE = max(abs(bmin_verse)); RFAmp.vCP = max(abs(vCP));
  %%
  PulseToWrite = 'VERSE';
  switch PulseToWrite
      case 'Full'
        RFShimWrite(ToWrite.Full);
        fprintf('The max voltage is %f Volts.\n',RFAmp.Full);
      case 'CP'
        RFShimWrite(ToWrite.CP);
        fprintf('The max voltage is %f Volts.\n',RFAmp.CP);
      case 'PhaseOnly'
        RFShimWrite(ToWrite.PhaseOnly); 
        fprintf('The max voltage is %f Volts.\n',RFAmp.PhaseOnly);
      case 'VERSE'
        RFShimWrite(ToWrite.VERSE); 
        fprintf('The max voltage is %f Volts.\n',RFAmp.VERSE); 
      case 'vCP'
        RFShimWrite(ToWrite.vCP); 
        fprintf('The max voltage is %f Volts.\n',RFAmp.vCP);          
  end
if isfolder('/Volumes/Disk_C')
    copyfile('/Users/ytong/Documents/MATLAB/tong-acptx/tySpiral/tyRFShimming/pTXVEPCASLShim.ini','/Volumes/Disk_C/MedCom/MriCustomer/seq/RFPulses/pTXVEPCASLShim.ini')
end
if isfolder('/Volumes/Disk_C-1')
    copyfile('/Users/ytong/Documents/MATLAB/tong-acptx/tySpiral/tyRFShimming/pTXVEPCASLShim.ini','/Volumes/Disk_C-1/MedCom/MriCustomer/seq/RFPulses/pTXVEPCASLShim.ini')
end


%% Simulating verse labelling efficiency
v = 0.3; zmax = 0.1;
T = 2*zmax / v;
delta_t = 10e-6;

% Round up T to the nearest dt and calculate the number of time steps
T = T - mod(T,delta_t) + delta_t;
N = T/delta_t;
%disp(['Number of time steps required is: ' num2str(N) ] );

% Generate the time array
t = 0:delta_t:T;

% Generate the off resonance frequency array
ORFreq = 0;
OffRes = ones(1,size(t,2))*ORFreq;

% Set up the position array of the spin: starting at -zmax and ending at
% +zmax
Ps = [0 0];
P = generate_position('const_vel', [Ps(1) Ps(2) -zmax]', [0 0 v]', t);
bv_norm = bv/(max(bv(:)));
%gv_norm = gv/(max(gv(:)));
[mx,my,mz] = bloch_CTR_Hz(42.57*bv_norm,(gv/10)*42.57E3,10e-6,inf,inf,0,0,0); % Per uT
mxy = abs(mx+1i*my);
MagOnRes = max(mxy(:));
FArBloch = asin(MagOnRes);
FAinDeg =  12.4262954448681;    % Will's for plot
%FAinDeg = 15;
FAinRad = deg2rad(FAinDeg);
RFAmpInuT = FAinRad/FArBloch;
RFAmpInmT = RFAmpInuT/1000;
T1 = inf;
T2 = 0.3;
RF = zeros(2,numel(t));
G = zeros(3,numel(t));
mean_tag_amp = 0.8;
Rephase_Moment = mean_tag_amp*120-trapz(gv);
% Making a triangle gradient
PeakTriag = Rephase_Moment/((120-95));
%Rephase_Grad = linspace(0,PeakTriag,13);
%Rephase_Grad = [Rephase_Grad, flip(Rephase_Grad(2:end))];
for iDx = 1:floor(T/1.2E-3)
    RF(1,((iDx-1)*120+1):((iDx-1)*120+numel(bv))) = RFAmpInmT*bv_norm;
    G(3,((iDx-1)*120+1):((iDx-1)*120+numel(bv))) = gv;
    G(3,((iDx-1)*120+numel(bv)):(iDx*120)) = PeakTriag;
end
% This bloch sim nees mT
M = bloch_sim(P, G, RF, OffRes, t, T1, T2);
(1-M(3,end))/2

%% Calculate tSNR 

   
    DirStructArray = dir('/Users/ytong/Documents/Data/ForISMRMabstract/2018*');
    DicomTree = cell(numel(DirStructArray,1));
    SeriesArray = cell(numel(DirStructArray,1));
    % 1. VERSE Shimmed 2. Gaussian Shimmed 3. Gaussian CP
    % ALL EPI FILES! Not Perf-weighted
    SeriesArray{1} = [18,22,26];    % 097
    SeriesArray{2} = [18,22,26];    % 098
    SeriesArray{3} = [22,26,30];    % 099
    SeriesArray{4} = [20,24,28];    % 100
    SeriesArray{5} = [19,27,23];    % 101
    SeriesArray{6} = [18,22,32];    % 102
    tSNR = cell(numel(DirStructArray,1));
    PerfMask = cell(numel(DirStructArray,1));
for iDx = 1:numel(DirStructArray)
        PathTemp = fullfile('/Users/ytong/Documents/Data/ForISMRMabstract', DirStructArray(iDx).name) ;
        DicomTree{iDx} = Spectro.dicomTree('dir',PathTemp,'recursive',false);
        [tSNR{iDx},PerfMask{iDx}] = CalcTSNR(DicomTree{iDx},SeriesArray{iDx});
end    
    
    
%% 
AvgTSNR = zeros(6,3);
StdTSNR = zeros(6,3);
for iDx = 1:numel(tSNR)
    for jDx = 1:3
        AvgTSNR(iDx,jDx) = mean(tSNR{iDx}(:,jDx));
        StdTSNR(iDx,jDx) = std(tSNR{iDx}(:,jDx));
    end
end
%%
tSNR_all_subjects = [tSNR{1};tSNR{2};tSNR{3};tSNR{4};tSNR{5};tSNR{6}];

%%
tSNR_Img = cell(3,1);
for iDx = 1:3
    tSNR_Img{iDx} = zeros(110,110);
    tSNR_Img{iDx}(logical(PerfMask{6})) = tSNR{6}(:,iDx);
    tSNR_Img{iDx}(91:110,:) = [];
    tSNR_Img{iDx}(:,end-9:end) = [];
    tSNR_Img{iDx}(:,1:10) = [];
end
Fig298 = figure(298);
set(gcf,'color','w','InvertHardcopy','off')
set(gcf,'units','centimeters','position',[4 4 30 12],'paperunits','centimeters','paperposition',[0 0 30 12])
clf
SPTemp = cell(1,3);
TitleArray = {'VERSE Shimmed','Gaussian Shimmed','Gaussian CP Mode'};
for iDx = 1:3
    SPTemp{iDx} = subplot(1,3,iDx);
    imagesc(tSNR_Img{iDx},[0 1.5]);axis equal;axis off;colormap('hot')
    Title_Temp = title(TitleArray{iDx});
    Title_Temp.FontSize = 14;
    Title_Temp.Position(2) = Title_Temp.Position(2)+23;
    nudge(SPTemp{iDx},[-0.05 -0.15 0 0]);
    SPTemp{iDx}.Position(3) = SPTemp{iDx}.Position(3)*1.2;
    SPTemp{iDx}.Position(4) = SPTemp{iDx}.Position(4)*1.2;
    if iDx == 3
        clb_obj = colorbar('FontSize',11);
        %clb_obj.Position(4) = SPTemp{3}.Position(2);
        nudge(clb_obj,[0.08 0.055 0 0]);
        clb_obj.Position(4) = clb_obj.Position(4)*0.85;
        %clb_obj.YLabel.String = 'tSNR';
    end
end
saveas(Fig298,'/Users/ytong/Documents/ISMRM/Main2019/tSNR.png')
%%
F299= figure(299);
set(gcf,'color','w','InvertHardcopy','off')
set(gcf,'units','centimeters','position',[4 4 17 20],'paperunits','centimeters','paperposition',[0 0 17 20])
clf
hold on
Mean2Plot = mean(AvgTSNR,1);
Std2Plot = std(StdTSNR,1);
BarChart = cell(1,3);
for iDx = 1:3
    BarChart{iDx} = bar(iDx,Mean2Plot(iDx),0.6);
    if iDx == 1
        set(BarChart{iDx},'FaceColor','k');
    elseif iDx == 2
        set(BarChart{iDx},'FaceColor','k');
    elseif iDx == 3
        set(BarChart{iDx},'FaceColor','k');
    end
end
  for iDx = 1:numel(BarChart)
    xData = BarChart{iDx}.XData+BarChart{iDx}.XOffset;
    ErrB_obj = errorbar(xData,Mean2Plot(iDx),Std2Plot(iDx),'-r');
    ErrB_obj.LineWidth = 3;
  end
box off;xtick off;ylim([0 0.55])
yl = ylabel('tSNR'); yl.FontSize = 20;
%Title_tSNR = title({'Mean tSNR of the Central','Slice of All Subjects'}); Title_tSNR.FontSize = 18;
Title_tSNR = title('Mean tSNR'); Title_tSNR.FontSize = 18;
%lgdTmp = legend({'VERSE Shimmed','Gaussian Shimmed','Gaussian CP Mode'});
%lgdTmp.FontSize = 18;
xt = get(gca, 'YTick');
set(gca, 'FontSize', 18)
set(gca, 'XTick', [1 2 3])
labels = {'VERSE Shimmed','Gaussian Shimmed','Gaussian CP-Mode'};
labels = cellfun(@(x) strrep(x,' ','\newline'), labels,'UniformOutput',false);
set(gca, 'XTickLabel', labels)
%set(gca, 'Xtickangle', 45)
saveas(gcf,'/Users/ytong/Documents/ISMRM/Main2019/Bar.png')

%%
Mask_Full_Img = maskedMaps.localiser(:,:,SliceIdx)';
figure(71);imagesc(Mask_Full_Img);axis equal;axis off
h_vessel = drawellipse('Color','w');

ellipse_mask = h_vessel.createMask();
%Vessel_Rect_Img = zeros(size(rect_mask));
%Vessel_Rect_Img(:) = Mask_Full_Img(rect_mask);
%Vessel_Rect_Img(Vessel_Rect_Img==0) = [];

Vessel_Rect_Img = ellipse_mask.*Mask_Full_Img;
%figure(72);imagesc(Vessel_Rect_Img);axis equal;axis off
thresh = multithresh(Vessel_Rect_Img,2);
seg_img = imquantize(Vessel_Rect_Img,thresh);
figure(72);imagesc(seg_img);axis equal;axis off


%%
% Trying to get the mean tSNR on all subjects after perfusion analysis
StrArray = {'099','100','101','102','106'};
tSNR = cell(5,1);   dir_7T = cell(5,1);    dir_3T = cell(5,1);
load('/Users/ytong/Documents/MATLAB/tong-acptx/tySpiral/tyRFShimming/SubjectInfo.mat');
clear PCASL
for iDx = 1:numel(StrArray)
    SubjectInfoSingle = FindSubject(strcat('F7T_2013_50_',StrArray{iDx}), 7);
    dir_7T{iDx} = strcat('/Users/ytong/Documents/Data/ForMRM/7T/F7T_2013_50_',StrArray{iDx});
    dir_3T{iDx} = strcat('/Users/ytong/Documents/Data/ForMRM/3T/F3T_2013_50_',num2str(SubjectInfoSingle.num_3T));
    cd(dir_7T{iDx});
    PCASL = dir('*PCASL*');
    for iDy = 1:numel(PCASL)
        if isfolder(PCASL(iDy).name)
            cd(PCASL(iDy).name)
                if  contains(PCASL(iDy).name,'43')
                        tSNR_file = dir('*tSNR_masked.nii.gz*');
                        [tSNR{iDx}.mean(1),tSNR{iDx}.std(1)] = Find_Mean_Std(tSNR_file(1).name);                      
                elseif (contains(PCASL(iDy).name,'53') ||    contains(PCASL(iDy).name,'52'))
                        tSNR_file = dir('*tSNR_masked.nii.gz*');
                        [tSNR{iDx}.mean(2),tSNR{iDx}.std(2)] = Find_Mean_Std(tSNR_file(1).name);
                elseif contains(PCASL(iDy).name,'23')   
                        tSNR_file = dir('*tSNR_masked.nii.gz*');
                        [tSNR{iDx}.mean(3),tSNR{iDx}.std(3)] = Find_Mean_Std(tSNR_file(1).name);
                elseif contains(PCASL(iDy).name,'28')   
                        tSNR_file = dir('*tSNR_masked.nii.gz*');
                        [tSNR{iDx}.mean(4),tSNR{iDx}.std(4)] = Find_Mean_Std(tSNR_file(1).name); 
                end
            cd ..
        end
    end
    cd(dir_3T{iDx});
    PCASL = dir('*PCASL*');
    for iDy = 1:numel(PCASL)
        if isfolder(PCASL(iDy).name)
            cd(PCASL(iDy).name)
            filename_parts = strsplit(PCASL(iDy).name,'_');
            filename__wo_prefix = filename_parts{end};
            switch  filename__wo_prefix
                case 'toep2dPCASLmatch7TMatchFlip'
                        tSNR_file = dir('*tSNR_masked.nii.gz*');
                        [tSNR{iDx}.mean(5),tSNR{iDx}.std(5)] = Find_Mean_Std(tSNR_file(1).name);                      
                case 'toep2dPCASLmatch7T'
                        tSNR_file = dir('*tSNR_masked.nii.gz*');
                        [tSNR{iDx}.mean(6),tSNR{iDx}.std(6)] = Find_Mean_Std(tSNR_file(1).name);
                case 'toep2dPCASL3Ttagging' 
                        tSNR_file = dir('*tSNR_masked.nii.gz*');
                        [tSNR{iDx}.mean(7),tSNR{iDx}.std(7)] = Find_Mean_Std(tSNR_file(1).name);
                case 'toep2dPCASL3Toptimized'   
                        tSNR_file = dir('*tSNR_masked.nii.gz*');
                        [tSNR{iDx}.mean(8),tSNR{iDx}.std(8)] = Find_Mean_Std(tSNR_file(1).name); 
            end
            cd ..
        end
    end    
end
%%
tSNR_mean_std = zeros(10,8);
for ii = 1:numel(tSNR)
    tSNR_mean_std((ii-1)*2+1,:) = tSNR{ii}.mean;
    tSNR_mean_std(ii*2,:) = tSNR{ii}.std;
end
%%
[image_SNR,Background] = calculate_image_SNR;
%%
Sub2_3T = find_perf_per_measurement('/Users/ytong/Documents/Data/ForMRM/3T/F3T_2013_50_104/images_011_toep2dPCASLmatch7TMatchFlip/images_011_toep2dPCASLmatch7TMatchFlip_diff.nii.gz');
Sub2_7T = find_perf_per_measurement('/Users/ytong/Documents/Data/ForMRM/7T/F7T_2013_50_100/images_020_toVEPCASLVERSE28/images_020_toVEPCASLVERSE28_diff.nii.gz');
%%
figure(222)
PlotRange = [0 700];
for iDx = 2:2:8
    subplot(2,4,iDx/2)
    imagesc(Sub2_7T(:,:,9,iDx*2),PlotRange);axis off;PlotTitle = sprintf('%s meas',num2str(iDx*2));title(PlotTitle);
    subplot(2,4,4+iDx/2)
    imagesc(Sub2_3T(:,:,10,iDx*2),PlotRange);axis off;title(PlotTitle);
end

%%  Calculate B1 correction
correction_map_new = get_signal_correction('/Users/ytong/Documents/Data/7TDicom/20181122_F7T_2013_50_106');

cd('/Users/ytong/Documents/Data/7TNifti/F7T_2013_50_106/')
niftiwrite(correction_map_new,'correction.nii');
system('fslcpgeom images_010_dtdreamwIce60deg90VRef1001.nii.gz correction.nii');
%%
FA_map_CP = get_signal_correction('/Users/ytong/Documents/Data/7TDicom/20181122_F7T_2013_50_106');
cd('/Users/ytong/Documents/Data/7TNifti/F7T_2013_50_106/')
niftiwrite(FA_map_CP,'FA_map.nii');
system('fslmaths images_010_dtdreamwIce60deg90VRef1001.nii.gz -Tmean B1Map_Tmean.nii.gz')
system('fslcpgeom B1Map_Tmean.nii.gz FA_map.nii');
%%
close all
plot_tsnr_trend(tSNR_mean_std_new)
%%
rund_perf_plot
%%
plot_hist(rad2deg((AFull_gaussian*bCP)),rad2deg((AFull_gaussian*bmin_gaussain)),...
    rad2deg((AFull_verse*bmin_verse)));
%%
plot_stacked_hist_all(Sim_FA_array);
%%
function plot_stacked_hist_all(Sim_FA_array)
    Sim_FA_cell = cell(5,5);
    for ii = 1:numel(Sim_FA_array)
        Sim_FA_cell{ii,1} = Sim_FA_array(ii).Gaussian_CP;
        Sim_FA_cell{ii,2} = Sim_FA_array(ii).Gaussian_shimmed;
        Sim_FA_cell{ii,3} = Sim_FA_array(ii).VERSE_shimmed;
        Sim_FA_cell{ii,4} = Sim_FA_array(ii).VERSE_relax_2;
        Sim_FA_cell{ii,5} = Sim_FA_array(ii).VERSE_relax_10;
    end
    Nbins = [20 30 32 33 15];
    YLIM = [0 150;0 150;0 150;0 150;0 400];
    TitleStr = {'Gaussian CP','Gaussian Shimmed','VERSE Shimmed','VERSE relaxed by 2','VERSE relaxed by 10'};
    EDGES = 0:0.5:25;
    for ii = 1:size(Sim_FA_cell,1)
        temp = {Sim_FA_cell{1,ii},Sim_FA_cell{2,ii},Sim_FA_cell{3,ii},...
            Sim_FA_cell{4,ii},Sim_FA_cell{5,ii}};
        plot_stacked_single(temp,ii,YLIM,TitleStr,Nbins,EDGES);
    end
    function plot_stacked_single(FA_cell,Order,YLIM,TitleStr,Nbins,EDGES)
        colormat =    [0.2078    0.1843    0.5294;
                       0.1804    0.6824    0.6706;
                       1         0.549     0];
        FA_all = [FA_cell{1};FA_cell{2};FA_cell{3};FA_cell{4};FA_cell{5}];
        %[~,edges] = histcounts(FA_all,Nbins(Order));
        [n1] = histcounts(FA_cell{1},EDGES);
        [n2] = histcounts(FA_cell{2},EDGES);
        [n3] = histcounts(FA_cell{3},EDGES);
        [n4] = histcounts(FA_cell{4},EDGES);
        [n5] = histcounts(FA_cell{5},EDGES);
        centres = EDGES(1:(end-1)) + diff(EDGES);
        subplot(5,1,Order)
        H = bar(centres,[n1.' n2.' n3.' n4.' n5'],'stacked','barwidth',1,'edgecolor','k','linewidth',0.01);
        xlabel('Flip-angle (°)')
        xlim([0 25])
        %xlim(XLIM(Order,:))
        ylim(YLIM(Order,:))
        title(TitleStr{Order})
    end
end
%%
function plot_hist(GCP,Gshimmed,Vshimmed)
figure(104)
set(gcf,'color','w','InvertHardcopy','off')
set(gcf,'units','centimeters','position',[4 4 30 20],'paperunits','centimeters','paperposition',[0 0 40 20])
%histogram(abs(FAFinal.AS),20,'Normalization','probability');hold on;
%[f1,xi] = ksdensity(abs(FAFinal.AS),linspace(17,22,100));
%plot(xi,f1,'b');
edges = 5:0.5:25;
Histo2{1,1} = histogram(abs(GCP),edges,'FaceColor','b');    Histo2{1,1}.Normalization = 'probability'; hold on
Histo2{2,1} = histogram(abs(Gshimmed),edges,'FaceColor','g');    Histo2{2,1}.Normalization = 'probability';
Histo2{3,1} = histogram(abs(Vshimmed),edges,'FaceColor','y');    Histo2{3,1}.Normalization = 'probability';
box off
xlabel('Flip angle (°)');xlim([5 25])
ylabel('Probability'); ylim([0 0.15])
title('Flip angle achieved in all vessels');
lgd = legend('CP mode','Gaussian shimmed','VERSE shimmed');legend('boxoff')   
lgd.FontSize = 17;set(gca, 'FontSize', 20)


end
%%
function rund_perf_plot
perf_to_plot = cell(2,2);
perf_to_plot{2,1} = niftiread('/Users/ytong/Documents/Data/ForMRM/7T/F7T_2013_50_099/images_022_toVEPCASLVERSE28/images_022_toVEPCASLVERSE28_mean.nii.gz');
perf_to_plot{1,1} = niftiread('/Users/ytong/Documents/Data/ForMRM/7T/F7T_2013_50_099/images_030_toVEPCASLGCP43/images_030_toVEPCASLGCP43_mean.nii.gz');
perf_to_plot{1,2} = niftiread('/Users/ytong/Documents/Data/ForMRM/7T/F7T_2013_50_099/images_026_toVEPCASLGShimmed53/images_026_toVEPCASLGShimmed53_mean.nii.gz');
perf_to_plot{2,2} = niftiread('/Users/ytong/Documents/Data/ForMRM/3T/F3T_2013_50_105/images_011_toep2dPCASLmatch7TMatchFlip/images_011_toep2dPCASLmatch7TMatchFlip_mean.nii.gz');
slices_to_plot = 9:2:13;
plot_range = [-1.5 11];
%close(77)
figure(77)
%Title_Offset = [-120,60];
Subplot_Offset = [-0.13 -0.07 0.19 0.19];
%perf_title = {{'Gaussian';'CP Mode'},'Gaussian Shimmed','VERSE Shimmed','3T Matched Seq & FA'};
for iii = 1:size(perf_to_plot,1)
    for jjj = 1:size(perf_to_plot,2)
        PLT = subplot(3,2,(iii-1)*2+jjj);
        img_temp = perf_to_plot{iii,jjj};
        img_temp_three = horzcat(rot90(img_temp(:,:,slices_to_plot(1))),...
            rot90(img_temp(:,:,slices_to_plot(2))),rot90(img_temp(:,:,slices_to_plot(3))));
        imagesc(img_temp_three,plot_range);axis off;axis equal
        colormap(PLT,gray)
        nudge(PLT,Subplot_Offset);
%         Ttl = title(perf_title{(iii-1)*2+jjj});
%         Ttl.FontSize = 16;      Ttl.Color = 'k';
%         Ttl.Position(1) = Ttl.Position(1)+Title_Offset(1);
%         Ttl.Position(2) = Ttl.Position(2)+Title_Offset(2);
    end
end
cbf_to_plot = cell(2,1);
cbf_to_plot{2} = niftiread('/Users/ytong/Documents/Data/ForMRM/3T/F3T_2013_50_105/images_017_toep2dPCASL3Toptimized/native_space/perfusion_calib.nii.gz');
cbf_to_plot{1} = niftiread('/Users/ytong/Documents/Data/ForMRM/7T/F7T_2013_50_099/images_022_toVEPCASLVERSE28/native_space/perfusion_calib.nii.gz');
CBF_range = [19 89];
%CBF_title = {'e','f'};
for iii = 1:2
        PLT = subplot(3,2,4+iii);
        img_temp = cbf_to_plot{iii};
        img_temp_three = horzcat(rot90(img_temp(:,:,slices_to_plot(1))),...
            rot90(img_temp(:,:,slices_to_plot(2))),rot90(img_temp(:,:,slices_to_plot(3))));
        imagesc(img_temp_three,CBF_range);axis off;axis equal
        colormap(PLT,jet)
        nudge(PLT,Subplot_Offset);
        if iii == 2
            clb_obj = colorbar('FontSize',18);
            clb_obj.XLabel.String = 'ml/100g/min';
            nudge(clb_obj,[0.02 0.0015 0 0]);
            clb_obj.Position(4) = clb_obj.Position(4)*0.97;
            clb_obj.Position(2) = clb_obj.Position(2)+0.002;
        end
%         Ttl = title(CBF_title{iii});
%         Ttl.FontSize = 20;
%         Ttl.Position(1) = Ttl.Position(1)+Title_Offset(1);
%         Ttl.Position(2) = Ttl.Position(2)+Title_Offset(2);
end
Position_vec = [4 4 40 20];
set(gcf,'color','w','InvertHardcopy','off')
set(gcf,'units','centimeters','position',Position_vec,'paperunits','centimeters','paperposition',Position_vec)
end
%%
function [image_SNR,Background] = calculate_image_SNR
load('/Users/ytong/Documents/MATLAB/tong-acptx/tySpiral/tyRFShimming/SubjectInfo.mat');
StrArray = {'099','100','101','102','106'};
dir_7T = cell(5,1);    dir_3T = cell(5,1);
image_SNR = cell(5,1);
Background = cell(5,2);
for ii = 1:numel(StrArray)
    dir_7T{ii} = strcat('/Users/ytong/Documents/Data/ForMRM/7T/F7T_2013_50_',StrArray{ii});
    SubjectInfoSingle = FindSubject(strcat('F7T_2013_50_',StrArray{ii}), 7);
    dir_3T{ii} = strcat('/Users/ytong/Documents/Data/ForMRM/3T/F3T_2013_50_',num2str(SubjectInfoSingle.num_3T));
end
for ii = 1:numel(StrArray)
    cd(dir_7T{ii});%PCASL = dir('*PCASL*');
    mask_cmd_7 = sprintf('bet EPI_unwarp/M0.nii.gz M0_masked');
    system(mask_cmd_7);
    [image_SNR{ii}.M0_7T,Background{ii,1}] = find_img_SNR('EPI_unwarp/M0','M0_masked.nii.gz');
%     temp_snr = zeros(4,3);
%     for iDy = 1:numel(PCASL)
%         if isfolder(PCASL(iDy).name)
%             cd(PCASL(iDy).name)
%                 if  contains(PCASL(iDy).name,'43')
%                         temp_snr(1,:) = find_SNR;                      
%                 elseif (contains(PCASL(iDy).name,'53') ||    contains(PCASL(iDy).name,'52'))
%                         temp_snr(2,:) = find_SNR;
%                 elseif contains(PCASL(iDy).name,'23')   
%                         temp_snr(3,:) = find_SNR;
%                 elseif contains(PCASL(iDy).name,'28')   
%                         temp_snr(4,:) = find_SNR; 
%                 end
%             cd ..
%         end
%     end
%     image_SNR{ii}.perf_7T = temp_snr;  
    cd(dir_3T{ii});%PCASL = dir('*PCASL*');
    mask_cmd_3 = 'bet M0_0000.nii.gz M0_masked';
    system(mask_cmd_3);
    [image_SNR{ii}.M0_3T,Background{ii,2}] = find_img_SNR('M0_0000.nii.gz','M0_masked.nii.gz');
%     temp_snr = zeros(4,3);
%     for iDy = 1:numel(PCASL)
%         if isfolder(PCASL(iDy).name)
%             cd(PCASL(iDy).name)
%             filename_parts = strsplit(PCASL(iDy).name,'_');
%             filename__wo_prefix = filename_parts{end};
%             switch  filename__wo_prefix
%                 case 'toep2dPCASLmatch7TMatchFlip'
%                         temp_snr(1,:) = find_SNR;                    
%                 case 'toep2dPCASLmatch7T'
%                         temp_snr(2,:) = find_SNR;   
%                 case 'toep2dPCASL3Ttagging' 
%                         temp_snr(3,:) = find_SNR;   
%                 case 'toep2dPCASL3Toptimized'   
%                         temp_snr(4,:) = find_SNR;   
%             end
%             cd ..
%         end
%     end
%     image_SNR{ii}.perf_3T = temp_snr;  
end
end
%%
function SNR = find_SNR
         perfusion_mean = dir('*PCASL*_mean.nii.gz');
                cell_temp = strsplit(perfusion_mean(1).name,'.');
                filename_mean = cell_temp{1};
                filename_mean_masked = strcat(filename_mean,'_masked');
                filename_base = filename_mean(1:end-5);
%                 mask_cmd{1} = sprintf('fslmaths %s -mul ../GM_pve_ASL_space_thresh_mask.nii.gz %s',filename_mean,...
%                     filename_mean_masked);
%                 split_cmd = sprintf('asl_file --data=%s --ibf=tis --iaf=tc --spairs --out=%s --ntis=1',...
%                     strcat(filename_base,'_mcf'),filename_base);
%                 mask_cmd{2} = sprintf('fslmaths %s -mul ../GM_pve_ASL_space_thresh_mask.nii.gz %s',strcat(filename_base,'_odd'),...
%                     strcat(filename_base,'_odd_masked'));
%                 mask_cmd{3} = sprintf('fslmaths %s -mul ../GM_pve_ASL_space_thresh_mask.nii.gz %s',strcat(filename_base,'_even'),...
%                     strcat(filename_base,'_even_masked'));
%                 system(mask_cmd{1});system(split_cmd);system(mask_cmd{2});system(mask_cmd{3})
        SNR = zeros(1,3);
        SNR(1,1) = find_img_SNR(filename_mean,filename_mean_masked);
        SNR(1,2) = find_img_SNR(strcat(filename_base,'_odd'),strcat(filename_base,'_odd_masked'));
        SNR(1,3) = find_img_SNR(strcat(filename_base,'_even'),strcat(filename_base,'_even_masked'));
end
%%
function [SNR,BackgroundArray] = find_img_SNR(Unmasked,Masked)
        [MEAN, ~] = Find_Mean_Std(Masked);
        NiiImg = niftiread(Unmasked);
        ColIndex = [1:5 106:110];
        TwoColumns = NiiImg(ColIndex,:,:);
        BackgroundArray = single(TwoColumns(:));
        Background_std = std(BackgroundArray);
        SNR = MEAN/Background_std;    
end
%%
function SubjectInfoSingle = FindSubject(FileName, FieldStrength)
    load SubjectInfo.mat SubjectInfo
    PathNameSplitted = strsplit(FileName,'_');
    SubjectNum = str2double(PathNameSplitted{end});

    for iDx = 1:numel(SubjectInfo)
        if FieldStrength == 3
            NumFromStorage = SubjectInfo(iDx).num_3T;
        elseif FieldStrength == 7
            NumFromStorage = SubjectInfo(iDx).num_7T;
        end
        if SubjectNum == NumFromStorage
            SubjectInfoSingle = SubjectInfo(iDx);
            return
        end
    end
end
function [MEAN, STD] = Find_Mean_Std(filename)
                mean_cmd = sprintf('fslstats %s -M',filename);
                [~,MEAN] = system(mean_cmd);
                MEAN = str2double(MEAN);
                std_cmd = sprintf('fslstats %s -S',filename);
                [~,STD] = system(std_cmd); 
                STD = str2double(STD);
end
%%
function output = find_perf_per_measurement(perf_mcf)
    input = niftiread(perf_mcf);
    [dim1,dim2,dim3,dim4] = size(input);
    output = zeros(dim1,dim2,dim3,dim4/2);
    for iDx = 1:dim4/2
        output(:,:,:,iDx) = rot90(sum(input(:,:,:,1:iDx),4));
    end
end
function FA_map_CP_final = get_signal_correction(foldername)
    STR_CELL = strsplit(foldername,'/');
    subject_dir_name = STR_CELL{end};
    dt = Spectro.dicomTree('dir',foldername,'recursive',false);
    ptxFMObj = DicomFM.WTCpTxFieldmaps(dt,'B1String','dt_dream_wIce_60deg_90VRef__B1',...
        'B0String','fieldmap_ptx7t_iso4mm_trans_RL','LocaliserString','localizer_3D_moved');
    B1_perV_map = ptxFMObj.getB1PerV('Hz','none');
    B1_abs_map = B1_perV_map*40;                    % 90 deg pulse equivalent to 40V hard pulse (1ms)
    rotation_map = B1_abs_map*1E-3;                 % Multiply by 1ms to get number of rotations
    FA_map_8_chan = rotation_map*360;                      % Multiply by 360 to get FA in degs
    FA_map_CP = abs(sum(FA_map_8_chan,4));
    %FA_map_intermediate = flip(FA_map_CP,1);
    FA_map_CP_final = flip(FA_map_CP,2);
    %correction_map = 1/sind(FA_map_CP);         
    signal_map = sind(FA_map_CP_final);             % Take the sine of FA map (in degrees)
    correction_map = 1/signal_map;                  % Take the inverse
    cd(strcat('/Users/ytong/Documents/Data/ForMRM/7T/',subject_dir_name))
    niftiwrite(correction_map,'correction.nii');
    niftiwrite(FA_map_CP_final,'FA_map.nii');
    whole_brain_B1 = dir('*dtdreamwIce60deg90VRef*');
    TMEAN_cmd = sprintf('fslmaths %s -Tmean B1Map_Tmean.nii.gz',whole_brain_B1.name);
    system(TMEAN_cmd);
    system('fslcpgeom B1Map_Tmean.nii.gz correction.nii');
    system('fslcpgeom B1Map_Tmean.nii.gz FA_map.nii');
    flirt_cmd = sprintf('flirt -applyxfm -usesqform -in correction.nii -ref %s -out correction_in_EPi',...
        'EPI_unwarp/M0.nii.gz');
    system(flirt_cmd);
end
%%
function plot_tsnr_trend(tSNR_mean_std)

figure(55)
hold on
%TickCell = {'Gaussian CP','Gaussian shimmed','VERSE CP','VERSE shimmed','Matched 3T matched FA'};
    for iDx = 1:2:9
%         errorbar(1:8,tSNR_mean_std(iDx,:),tSNR_mean_std(iDx+1,:),...
%             '-d','MarkerSize',8)
        plot_temp = plot(1:8,tSNR_mean_std(iDx,:),...
            '-d','MarkerSize',11);
        plot_temp.LineWidth = 3.5;
    end
hold off
line([4.5 4.5],[0 0.6],'Color','red','LineStyle','--','LineWidth',2.5)
xlim([1 8])
xticks(1:8)
A1 = annotation('textbox',[0.278 0.6 0.2 0.3],'String','7 T Shimming Approaches','FitBoxToText','on');
A1.FontSize = 32;   A1.FontWeight = 'Bold';
A2 = annotation('textbox',[0.542 0.6 0.2 0.3],'String','3 T Comparions','FitBoxToText','on');
A2.FontSize = 32;   A2.FontWeight = 'Bold';
legend({'Subject 1','Subject 2','Subject 3','Subject 4','Subject 5'},'Location','northwest')
% xticklabels({'Gaussian CP','Gaussian shimmed','VERSE CP','VERSE shimmed',...
%     'Matched seq & FA',strcat('Matched seq (20',char(176),')'),'Optimised labelling','3T optimised'})
x_ticks = cell(8,1);
x_ticks{1} = sprintf('Gaussian\nCP Mode');
x_ticks{2} = sprintf('Gaussian\nShimmed');
x_ticks{3} = sprintf('VERSE\nCP Mode');
x_ticks{4} = sprintf('VERSE\nShimmed');
x_ticks{5} = sprintf('Matched\nSeq & FA');
x_ticks{6} = sprintf('Matched\nseq & 20%s',char(176));
x_ticks{7} = sprintf('Optmised\nLabelling');
x_ticks{8} = sprintf('3 T\nOptmised');
[hx,~] = format_ticks(gca,x_ticks);
for iDx = 1:numel(hx)
    hx(iDx).FontSize = 23;
    hx(iDx).FontWeight = 'Bold';
end
set(gca,'FontSize',23,'FontWeight','Bold')
ylabel('tSNR')
set(gcf,'color','w','InvertHardcopy','off')
set(gcf,'units','centimeters','position',[4 4 72 24],'paperunits','centimeters','paperposition',[4 4 72 24])
end