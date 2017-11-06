    %%   Loading the data and masking
    addpath('/Users/ytong/Documents/MATLAB/For_James_Larkin');
    pTxPath = '/Users/ytong/Documents/MATLAB/Temp/20171013_F7T_2013_40_387';
    dt = Spectro.dicomTree('dir',pTxPath,'recursive',false);
    ptxFMObj = DicomFM.WTCpTxFieldmaps(dt,'B1String','dt_dream_wIce_60deg_150VRef__B1',...
        'B0String','fieldmap_ptx7t_iso4mm_trans_RL','LocaliserString','fl_tof');
    %ptxFMObj = DicomFM.WTCpTxFieldmaps(pTxPath);
    SliceIdx = input('Please enter the slice number: \n');
    ptxFMObj.interpolateTo('Localiser');
    ptxFMObj.createMask(@(x) DicomFM.maskFunctions.ellipseMask(x,10),true);  
    %ptxFMObj.setSlice(10);
    %%  Specifying parameters
    param.targetFlipAngle = 20;
    param.numCh = 8;
    param.TR = 1e-3;% sec
    param.CGtikhonov = 1e-6;
    param.tol = 1e-5;
    param.MaxEvaluation = 25000;
    param.FOX = 25;%25cm

    maskedMaps.b1SensMasked = ptxFMObj.getB1PerV('uT','delete');
    maskedMaps.b0MapMasked = ptxFMObj.getB0('Hz','delete');
    maskedMaps.b0MapMaskedRad = (2*pi)*maskedMaps.b0MapMasked;
    maskedMaps.posVox = ptxFMObj.getPositions('cm','delete').';
    maskedMaps.localiser = ptxFMObj.getLoc('none');
    maskedMaps.b1SensMaskedHz = ptxFMObj.getB1PerV('Hz','delete');
    maskedMaps.mask = ptxFMObj.getMask();
    
    ImgToPlot.b0 = ptxFMObj.getB0('Hz','none');
    ImgToPlot.b1 = ptxFMObj.getB1('Hz','none');
    ImgToPlot.b1CP = abs(sum(ImgToPlot.b1,4));    

    RFStruct = tyMakeHanning(600,5);
    %%  Constructing system matrix
    disp('Calculating system matrix...');
    SysMatMode = 'Full';
    if strcmp(SysMatMode,'Approximation');
        AFull = getAMatSimp(RFStruct,maskedMaps.b1SensMasked,maskedMaps.b0MapMasked); %Nv-(2*Nc) sys mtx, rad/V
    elseif strcmp(SysMatMode,'Full');
        gr = zeros(numel(RFStruct.RF_pulse),3);
        rfOn = true(size(gr,1),1);
        tp = 10E-6 * ones(size(rfOn));
        DA = genAMatFull(tp,rfOn,gr,maskedMaps.b1SensMaskedHz,...
            maskedMaps.b0MapMasked,maskedMaps.posVox);
     %  To take the RF pulse shape into account. Construct a pseudo-diagonal matrix
        RFDiag = zeros(size(DA,2),8);
        RFSize = numel(RFStruct.RF_pulse);
        for iDx = 1:8
            RFDiag(((1:RFSize)+(iDx-1)*RFSize),iDx) = RFStruct.RF_pulse';
        end
        AFullMag = DA*RFDiag;  %    magnetization /V
        AFull = asin(AFullMag);   % rad/V  Dupas paper in rad/V
    end
    disp('Complete')
    %%  Running RF shimming Step 1: Variable-exchange method
    disp('Running RF shimming variable-exchange optimization...')
    tikhonovArray = power(10,-6:-1);
    bVE = zeros(8,numel(tikhonovArray));
    for iDx = 1:numel(tikhonovArray)
        param.CGtikhonov = tikhonovArray(iDx);
        [bVE(:,iDx),NRMSETmp,~,~] = runVE(AFull,param,maskedMaps);
        disp(NRMSETmp)
    end
    disp('Complete')
    %%  Running RF shimming Step 2: Active-set method
    disp('Running RF shimming active-set optimization...')
    bAS = zeros(16,numel(tikhonovArray));
    NRMSE = zeros(size(tikhonovArray));
    output = cell(size(NRMSE));
    exitflag = zeros(size(NRMSE));
    for iDx = 1:numel(tikhonovArray)
        param.CGtikhonov = tikhonovArray(iDx);
        [bAS(:,iDx),NRMSE(iDx),exitflag(iDx),output{iDx}] = runAS(bVE(:,iDx),RFStruct,maskedMaps,param,AFull);
    end
    [~, minIndex] = min(NRMSE);
    disp('Complete')
    
    %%  Bloch Simulation
    disp('Running Bloch simulation...')
    %bmin = bVE(:,minIndex);
    bmin = bAS(1:8,minIndex).*exp(1i*bAS(9:16,minIndex));
    RFToSim = bsxfun(@times,repmat(RFStruct.RF_pulse.',1,8),bmin.');    % Note the diff between .' and '
    GToSim = zeros(numel(RFToSim)/8,3);
    %GToSim(:,3) = 425.8*RFStruct.G_amp*ones(numel(RFToSim)/8,1);     %from mT/m to Hz/cm
    B1ToSim = ptxFMObj.getB1PerV('Hz','delete');
    magnetization = make_blochSim(RFToSim,B1ToSim,...
        maskedMaps.b0MapMasked,GToSim,10E-6,maskedMaps.posVox,maskedMaps.mask);
    disp('Bloch simulation complete')
    
    %%  Finding out what CP mode can do
    bCP = runCP(AFull,param,RFStruct);
    %%  Running phase only optimization
    tic
        [bPhaseTmp, ErrorOut] = runPhaseOnly(AFull,param,RFStruct);
    toc
    [~, minPhaseIndex] = min(ErrorOut);
    %[~, maxPhaseIndex] = max(ErrorOut);    
    
    bPhase = bPhaseTmp(1,minPhaseIndex)*exp(1i*bPhaseTmp(2:9,minPhaseIndex));
    
    %%  Calculate magnetization and NRMSE
  rmse = @(x,xref) sqrt(immse(x,xref));
  nrmse = @(x,xref) rmse(x,xref)/mean(x);
  % FA results in degrees in a masked vector form.
  FAFinal = struct;
  FAFinal.CP = rad2deg(AFull*bCP);
  FAFinal.AS = rad2deg(AFull*bmin);
  FAFinal.PhaseOnly = rad2deg(AFull*bPhase);
  
  FullImage = struct;
  JawSlicesMask = maskedMaps.mask(:,:,SliceIdx-1:SliceIdx+1);
  FullImage.CP = zeros(size(JawSlicesMask));
  FullImage.CP(JawSlicesMask) = abs(FAFinal.CP);
  FullImage.AS = zeros(size(JawSlicesMask));
  FullImage.AS(JawSlicesMask) = abs(FAFinal.AS);  
  FullImage.ASBloch = asind(magnetization.mxy);
  FullImage.PhaseOnly = zeros(size(JawSlicesMask));
  FullImage.PhaseOnly((JawSlicesMask)) = abs(FAFinal.PhaseOnly);  
  
  DiffImage.CP = NaN(size(JawSlicesMask));
  DiffImage.CP(JawSlicesMask) = 100*(abs(FAFinal.CP) - ones(size(FAFinal.CP))*param.targetFlipAngle)/param.targetFlipAngle;
  DiffImage.AS = NaN(size(JawSlicesMask));
  DiffImage.AS(JawSlicesMask) = 100*(abs(FAFinal.AS) - ones(size(FAFinal.AS))*param.targetFlipAngle)/param.targetFlipAngle;  
  DiffImage.ASBloch = asind(magnetization.mxy);
  DiffImage.PhaseOnly = NaN(size(JawSlicesMask));
  DiffImage.PhaseOnly((JawSlicesMask)) = 100*(abs(FAFinal.PhaseOnly) - ones(size(FAFinal.PhaseOnly))*param.targetFlipAngle)/param.targetFlipAngle;
  
  
  Error = struct;
  Error.CP = nrmse(abs(FAFinal.CP),ones(size(FAFinal.CP))*param.targetFlipAngle);
  Error.AS = nrmse(abs(FAFinal.AS),ones(size(FAFinal.CP))*param.targetFlipAngle);
  Error.PhaseOnly = nrmse(abs(FAFinal.PhaseOnly),ones(size(FAFinal.CP))*param.targetFlipAngle); 
  
  %%    Write RF pulse into file.
  makeRFToWrite = @(b) [abs(b)/max(abs(b)), angle(b)];

  ToWrite.Full = makeRFToWrite(bmin);
  ToWrite.CP = makeRFToWrite(bCP);
  ToWrite.PhaseOnly = makeRFToWrite(bPhase);
  for iDx = 1:8
      if ToWrite.Full(iDx,2) < 0
          ToWrite.Full(iDx,2) = ToWrite.Full(iDx,2) + (2*pi);
      elseif ToWrite.PhaseOnly(iDx,2) < 0
          ToWrite.PhaseOnly(iDx,2) = ToWrite.PhaseOnly(iDx,2) + (2*pi);           
      elseif ToWrite.CP(iDx,2) < 0
          ToWrite.CP(iDx,2) = ToWrite.CP(iDx,2) + (2*pi);         
      end
  end
  RFAmp.Full = max(abs(bmin));   RFAmp.CP = max(abs(bCP));    RFAmp.PhaseOnly = max(abs(bPhase));
  RFShimWrite(ToWrite.Full);
    %%  Plotting
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
 
%% Other plots
figure(90)
clf
 PlotFALim = [0 30];   
Figs = cell(3,3);
Figs{2,1} = subplot(3,2,3);
FontSize = 6;
imagesc(ImgToPlot.b1CP(plotRolArray,plotCowArray,SliceIdx)'); title('B1 Map','FontSize',FontSize);axis off;colorbar('FontSize',FontSize);hold on;
colormap(Figs{2,1},'gray')
visboundaries(maskedMaps.mask(plotRolArray,plotCowArray,SliceIdx)','Color','r','LineWidth',0.4,...
    'EnhanceVisibility', true,'LineStyle','-');axis off;

Figs{3,1} = subplot(3,2,5); 
imagesc(ImgToPlot.b0(plotRolArray,plotCowArray,SliceIdx)'); title('B0 Map','FontSize',FontSize);axis off;hold on;colorbar('FontSize',FontSize);
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
%imagesc(maskedMaps.mask(plotRolArray,plotCowArray,SliceIdx)');title('TOF');axis off;

%hh = figure(91);
HoriOffset = 0.05;
Figs{1,2} = subplot(3,4,3);
imagesc(FullImage.CP(plotRolArray,plotCowArray,2)',PlotFALim); title('CP Mode','FontSize',FontSize);axis off;
set(Figs{1,2},'Position',[Figs{1,2}.Position(1)-HoriOffset Figs{1,2}.Position(2) Figs{1,2}.Position(3) Figs{1,2}.Position(4)]);
colormap(Figs{1,2},'hot')

Figs{2,2} = subplot(3,4,7);
imagesc(FullImage.PhaseOnly(plotRolArray,plotCowArray,2)',PlotFALim); title('Phase Only','FontSize',FontSize);axis off;
set(Figs{2,2},'Position',[Figs{2,2}.Position(1)-HoriOffset Figs{2,2}.Position(2) Figs{2,2}.Position(3) Figs{2,2}.Position(4)]);
colormap(Figs{2,2},'hot')

Figs{3,2} = subplot(3,4,11);
imagesc(FullImage.AS(plotRolArray,plotCowArray,2)',PlotFALim); title('Full RF Shimming','FontSize',FontSize);axis off;
set(Figs{3,2},'Position',[Figs{3,2}.Position(1)-HoriOffset Figs{3,2}.Position(2) Figs{3,2}.Position(3) Figs{3,2}.Position(4)]);
colormap(Figs{3,2},'hot')

hp4 = get(Figs{3,2},'Position');
colorbar('Position', [hp4(1)+hp4(3)+0.015  hp4(2)  0.02  0.815],'FontSize',FontSize);

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
colorbar('Position', [hp4(1)+hp4(3)+0.015  hp4(2)  0.02  0.815],'FontSize',FontSize);
%%
 b1map = ptxFMObj.getB1PerV('Hz','NaN');
 b1mapCP = squeeze(sum(b1map,4));
 figure(54);imagesc(abs(b1map(:,:,11)))
 colorbar