    %%   Loading the data and masking
    addpath('/Users/ytong/Documents/MATLAB/For_James_Larkin');
    pTxPath = '/Volumes/DICOM/2017-11/20171107_F7T_2013_50_083';
    dt = Spectro.dicomTree('dir',pTxPath,'recursive',false);
    ptxFMObj = DicomFM.WTCpTxFieldmaps(dt,'B1String','dt_dream_wIce_60deg_100VRef__B1',...
        'B0String','fieldmap_ptx7t_iso4mm_trans_RL','LocaliserString','fl_tof');
    %ptxFMObj = DicomFM.WTCpTxFieldmaps(pTxPath);
    SliceIdx = input('Please enter the slice number: \n');
    ptxFMObj.interpolateTo('Localiser');
    ptxFMObj.createMask(@(x) DicomFM.maskFunctions.ellipseMask(x,SliceIdx),true);  
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
  PulseToWrite = 'Full';
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
  end
if isdir('/Volumes/Disk_C')
    copyfile('/Users/ytong/Documents/MATLAB/tong-acptx/tySpiral/tyRFShimming/pTXVEPCASLShim.ini?,?/Volumes/Disk_C/MedCom/MriCustomer/seq/RFPulses/pTXVEPCASLShim.ini')
end
if isdir('/Volumes/Disk_C-1')
    copyfile('/Users/ytong/Documents/MATLAB/tong-acptx/tySpiral/tyRFShimming/pTXVEPCASLShim.ini?,?/Volumes/Disk_C-1/MedCom/MriCustomer/seq/RFPulses/pTXVEPCASLShim.ini')
end

%%
runPlot;
%%
 b1map = ptxFMObj.getB1PerV('Hz','NaN');
 b1mapCP = squeeze(sum(b1map,4));
 figure(54);imagesc(abs(b1map(:,:,11)))
 colorbar