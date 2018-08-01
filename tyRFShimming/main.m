    %%   Loading the data and masking
    addpath('/Users/ytong/Documents/MATLAB/For_James_Larkin');
    cd('/Users/ytong/Documents/MATLAB/tong-acptx/tySpiral/tyRFShimming/')
    %pTxPath = '/Volumes/Data/DICOM/2018-06/20180601_F7T_2013_50_092';
    pTxPath = '/Users/ytong/Documents/Data/20180601_F7T_2013_50_092';
    %pTxPath = '/Users/ytong/Documents/Data/20171107_F7T_2013_50_083';
    %pTxPath = '/Users/ytong/Documents/Data/20171107_F7T_2013_50_083';
    %pTxPath = '/Users/ytong/Documents/MATLAB/Temp/20171013_F7T_2013_40_387';
    %pTxPath = '/Users/ytong/Documents/Data/20180525_F7T_2013_50_091';
    %%
    dt = Spectro.dicomTree('dir',pTxPath,'recursive',false);
    %
    ptxFMObj = DicomFM.WTCpTxFieldmaps(dt,'B1String','dt_dream_wIce_60deg_150VRef_tagging__B1',...
        'B0String','fieldmap_ptx7t_iso4mm_trans_RL','LocaliserString','fl_tof','InterpTarget',...
        'Localiser');
    %ptxFMObj = DicomFM.WTCpTxFieldmaps(pTxPath);
    
    SliceIdx = 15;%input('Please enter the slice number: \n');
    
    %ptxFMObj.interpolateTo('Localiser');
    
    ptxFMObj.createMask(@(x) DicomFM.maskFunctions.VEellipseMask(x,SliceIdx,[1 1 1 1]),true);  
    %ptxFMObj.setSlice(10);

    
    %RFFreqShift = @(RF,tVec,dfInRad) RF.*((1i)*dfInRad*tVec);
    %% 
    load Vessel.mat
  figure(96)
  for iDx = 1:numel(VesselMask)
    imagesc(VesselMask{iDx}); drawnow;
    pause(1)
  end
  close(96)
    %%  Specifying parameters
    param.targetFlipAngle = 20;
    param.numCh = 8;
    param.TR = 1.1e-3;% sec
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
    maskedMaps.TargetMasked = Vessel.TargetMasked;
    
    
    ImgToPlot.b0 = ptxFMObj.getB0('Hz','none');
    ImgToPlot.b1 = ptxFMObj.getB1PerV('Hz','none');
    ImgToPlot.b1CP = abs(sum(ImgToPlot.b1,4));    

    %RFStruct = tyMakeHanning(600,5);
    RFStruct = tyMakeGaussian(600,4);
    %%  Constructing system matrix
    disp('Calculating system matrix...');
    SysMatMode = 'Full';
    if strcmp(SysMatMode,'Approximation');
        AFull = getAMatSimp(RFStruct,maskedMaps.b1SensMasked,maskedMaps.b0MapMasked); %Nv-(2*Nc) sys mtx, rad/V
    elseif strcmp(SysMatMode,'Full');
        gr = zeros(numel(RFStruct.RF_pulse),3);
        rfOn = true(size(gr,1),1);
        tp = 10E-6 * ones(size(rfOn)); % in seconds
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
    %  Running RF shimming Step 1: Variable-exchange method
    disp('Running RF shimming variable-exchange optimization...')
    tic
    tikhonovArray = power(10,-6:-1);
    bVE = zeros(8,numel(tikhonovArray));
    for iDx = 1:numel(tikhonovArray)
        param.CGtikhonov = tikhonovArray(iDx);
        [bVE(:,iDx),NRMSETmp,~,~] = runVE(AFull,param,maskedMaps);
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
        [bAS(:,iDx),NRMSE(iDx),exitflag(iDx),output{iDx}] = runAS(bVE(:,iDx),RFStruct,maskedMaps,param,AFull);
    end
    [~, minIndex] = min(NRMSE);
    toc
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
    
    %  Finding out what CP mode can do
    bCP = runCP(AFull,param,RFStruct);
    %  Running phase only optimization
    tic
        [bPhaseTmp, ErrorOut] = runPhaseOnly(AFull,param,RFStruct);
    toc
    [~, minPhaseIndex] = min(ErrorOut); 
    bPhase = bPhaseTmp(1,minPhaseIndex)*exp(1i*bPhaseTmp(2:9,minPhaseIndex));


    %  Calculate magnetization and NRMSE
  rmse = @(x,xref) sqrt(immse(x,xref));
  nrmse = @(x,xref) rmse(x,xref)/mean(x);
  % FA results in degrees in a masked vector form.
  FAFinal = struct;
  FAFinal.CP = rad2deg(AFull*bCP);
  FAFinal.AS = rad2deg(AFull*bmin);
  FAFinal.PhaseOnly = rad2deg(AFull*bPhase);

  FullImage = struct;
  % JawSlicesMask was originally designed to contain 3 slices (selected slice plus 2 adjacent slices)
  % Attempting to reduce the design to one slice only. 
  JawSlicesMask = maskedMaps.mask(:,:,SliceIdx);
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

  OneSlice = FullImage;
  %OneSlice.CP = FullImage.CP(:,:,2);      OneSlice.PhaseOnly = FullImage.PhaseOnly(:,:,2);  OneSlice.Full = FullImage.AS(:,:,2);
  RICA{1} = OneSlice.CP(VesselMask{1});   RICA{2} = OneSlice.PhaseOnly(VesselMask{1});  RICA{3} = OneSlice.AS(VesselMask{1});
  RVA{1} = OneSlice.CP(VesselMask{2});   RVA{2} = OneSlice.PhaseOnly(VesselMask{2});  RVA{3} = OneSlice.AS(VesselMask{2});
  LICA{1} = OneSlice.CP(VesselMask{3});   LICA{2} = OneSlice.PhaseOnly(VesselMask{3});  LICA{3} = OneSlice.AS(VesselMask{3});
  LVA{1} = OneSlice.CP(VesselMask{4});   LVA{2} = OneSlice.PhaseOnly(VesselMask{4});  LVA{3} = OneSlice.AS(VesselMask{4});
  Full{1} = abs(FAFinal.CP);       Full{2} = abs(FAFinal.PhaseOnly);  Full{3} = abs(FAFinal.AS);
  %%    Calculate labelling efficiency 
  %EffArray = LabelEff(15:0.05:25,RFStruct);
  if ~exist('Efficiency','var')
    [Efficiency,FinalMag] = LabelEff(0:0.05:25,RFStruct);
  end
  %%
  BarMtx = zeros(3,5); StdMtx= zeros(3,5); FAToSim = 0:0.05:25;
  CalcEff = @(x) mean(interp1(FAToSim,Efficiency,x));
  getStd = @(x) std(interp1(FAToSim,Efficiency,x));
  for iDx = 1:numel(RICA)   
      BarMtx(iDx,1) = CalcEff(RICA{iDx});    BarMtx(iDx,2) = CalcEff(RVA{iDx});   BarMtx(iDx,3) = CalcEff(LICA{iDx});  
      BarMtx(iDx,4) = CalcEff(LVA{iDx});   BarMtx(iDx,5) = CalcEff(Full{iDx});
      StdMtx(iDx,1) = getStd(RICA{iDx});      StdMtx(iDx,2) = getStd(RVA{iDx});   StdMtx(iDx,3) = getStd(LICA{iDx});
      StdMtx(iDx,4) = getStd(LVA{iDx});      StdMtx(iDx,5) = getStd(Full{iDx});
      
  end
  BarMtx = BarMtx'; StdMtx = StdMtx';

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
  %%
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
    copyfile('/Users/ytong/Documents/MATLAB/tong-acptx/tySpiral/tyRFShimming/pTXVEPCASLShim.ini','/Volumes/Disk_C/MedCom/MriCustomer/seq/RFPulses/pTXVEPCASLShim.ini')
end
if isdir('/Volumes/Disk_C-1')
    copyfile('/Users/ytong/Documents/MATLAB/tong-acptx/tySpiral/tyRFShimming/pTXVEPCASLShim.ini','/Volumes/Disk_C-1/MedCom/MriCustomer/seq/RFPulses/pTXVEPCASLShim.ini')
end

%%
runPlot;
%%
  b1map = ptxFMObj.getB1PerV('Hz','NaN');
  b1mapCP = squeeze(sum(b1map,4));
  figure(54);imagesc(abs(b1map(:,:,16)))
  colorbar
 
%%
RF_duration_sec_vec = (600:25:750)*(1e-6);
Efficiency_cell = cell(numel(RF_duration_sec_vec),1);
FinalMag_cell = cell(numel(RF_duration_sec_vec),1);
for iDx = 1:numel(RF_duration_sec_vec)
    RFStruct = tyMakeGaussian(RF_duration_sec_vec(iDx)*1E6,4);
    [Efficiency_cell{iDx},FinalMag_cell{iDx}] = LabelEff(0:0.2:30,RFStruct,RF_duration_sec_vec(iDx), 1100e-6);
end

figure(572)
hold on 

for iDx = 1:numel(RF_duration_sec_vec)
    plot(0:0.2:30,Efficiency_cell{iDx})

end
%% Writing the workspace into a .mat file

FilePartsCell = strsplit(pTxPath,filesep);
FileNameStr = FilePartsCell{end};
StructToSave = struct('ptxFMObj',ptxFMObj,'param',param,'maskedMaps',maskedMaps,'RFStruct',RFStruct);
StructToSave.shims = struct('bAS',bAS,'bCP',bCP,'bmin',bmin,'bPhase',bPhase);


