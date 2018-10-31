    %%   Loading the data and masking
    addpath('/Users/ytong/Documents/MATLAB/For_James_Larkin');
    cd('/Users/ytong/Documents/MATLAB/tong-acptx/tySpiral/tyRFShimming/')
    %pTxPath = '/Volumes/Data/DICOM/2018-10/20181019_F7T_2013_50_097';
    %pTxPath = '/Volumes/Data/DICOM/2018-10/20181022_F7T_2013_50_098';
    %pTxPath = '/Volumes/Data/DICOM/2018-10/20181026_F7T_2013_50_099';
    pTxPath = '/Volumes/Data/DICOM/2018-10/20181030_F7T_2013_50_100';
    dt = Spectro.dicomTree('dir',pTxPath,'recursive',false);
    %%
    ptxFMObj = DicomFM.WTCpTxFieldmaps(dt,'B1String','dt_dream_wIce_60deg_160VRef_tagging_32meas__B1',...
        'B0String','fieldmap_ptx7t_iso4mm_trans_RL','LocaliserString','fl_tof','InterpTarget',...
        'Localiser');
    
%     ptxFMObj = DicomFM.WTCpTxFieldmaps(dt,'B1String','dt_dream_wIce_60deg_90VRef__B1',...
%         'B0String','fieldmap_ptx7t_iso4mm_trans_RL','LocaliserString','localizer_3D_moved','InterpTarget',...
%         'Localiser');
    %%
    SliceIdx = 18;  % Decided by visual inspection of the tof images
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
    if strcmp(SysMatMode,'Approximation');
        AFull_gaussian = getAMatSimp(RFStruct_Example,maskedMaps.b1SensMasked,maskedMaps.b0MapMasked); %Nv-(2*Nc) sys mtx, rad/V
    elseif strcmp(SysMatMode,'Full');
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
    if strcmp(SysMatMode,'Approximation');
        AFull_verse = getAMatSimp(RFStruct,maskedMaps.b1SensMasked,maskedMaps.b0MapMasked); %Nv-(2*Nc) sys mtx, rad/V
    elseif strcmp(SysMatMode,'Full');
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
ylabel('CP Mode','FontSize',FontSize,'FontWeight','bold');
Figs{1,2}.YLabel.Visible = 'on';
set(Figs{1,2},'Position',[Figs{1,2}.Position(1) Figs{1,2}.Position(2) Figs{1,2}.Position(3) Figs{1,2}.Position(4)]);
colormap(Figs{1,2},'hot')

Figs{2,2} = subplot(3,2,4);
imagesc(FullImage.verse(plotRolArray,plotCowArray)',PlotFALim);axis off;
ylabel('VERSE','FontSize',FontSize,'FontWeight','bold');
Figs{2,2}.YLabel.Visible = 'on';
set(Figs{2,2},'Position',[Figs{2,2}.Position(1) Figs{2,2}.Position(2) Figs{2,2}.Position(3) Figs{2,2}.Position(4)]);
colormap(Figs{2,2},'hot')
nudge(Figs{2,2},[0 0.05 0 0]);

Figs{3,2} = subplot(3,2,6);
imagesc(FullImage.AS(plotRolArray,plotCowArray)',PlotFALim);axis off;
set(Figs{3,2},'Position',[Figs{3,2}.Position(1) Figs{3,2}.Position(2) Figs{3,2}.Position(3) Figs{3,2}.Position(4)]);
ylabel('Full RF Shimming','FontSize',FontSize,'FontWeight','bold');
Figs{3,2}.YLabel.Visible = 'on';
CLB{6} = colorbar('southoutside','FontSize',CLBFontSize);
CLB{6}.YLabel.String = 'Degrees (°)';    CLB{6}.YLabel.Rotation = 0;   
colormap(Figs{3,2},'hot')
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
ylabel('CP Mode','FontSize',FontSize,'FontWeight','bold');
Figs{1,3}.YLabel.Visible = 'on';
%set(Figs{1,3},'Position',[Figs{1,3}.Position(1) Figs{1,3}.Position(2) Figs{1,3}.Position(3) Figs{1,3}.Position(4)]);
colormap(Figs{1,3},'hot')

Figs{3,3} = subplot('Position',FA_map_positions{3});
imagesc(FullImage.verse(plotRowArrayReduced,plotCowArray)',PlotFALim);axis off;
ylabel('VERSE','FontSize',FontSize,'FontWeight','bold');
Figs{3,3}.YLabel.Visible = 'on';
%set(Figs{2,3},'Position',[Figs{2,3}.Position(1) Figs{2,3}.Position(2) Figs{2,3}.Position(3) Figs{2,3}.Position(4)]);
colormap(Figs{3,3},'hot')
%nudge(Figs{2,3},[0 0.05 0 0]);

Figs{2,3} = subplot('Position',FA_map_positions{2});
imagesc(FullImage.AS(plotRowArrayReduced,plotCowArray)',PlotFALim);axis off;
%set(Figs{3,3},'Position',[Figs{3,3}.Position(1) Figs{3,3}.Position(2) Figs{3,3}.Position(3) Figs{3,3}.Position(4)]);
ylabel('Full RF Shimming','FontSize',FontSize,'FontWeight','bold');
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
  PulseToWrite = 'vCP';
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
if isdir('/Volumes/Disk_C')
    copyfile('/Users/ytong/Documents/MATLAB/tong-acptx/tySpiral/tyRFShimming/pTXVEPCASLShim.ini','/Volumes/Disk_C/MedCom/MriCustomer/seq/RFPulses/pTXVEPCASLShim.ini')
end
if isdir('/Volumes/Disk_C-1')
    copyfile('/Users/ytong/Documents/MATLAB/tong-acptx/tySpiral/tyRFShimming/pTXVEPCASLShim.ini','/Volumes/Disk_C-1/MedCom/MriCustomer/seq/RFPulses/pTXVEPCASLShim.ini')
end