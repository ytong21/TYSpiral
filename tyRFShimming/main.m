    %%   Loading the data and masking
    addpath('/Users/ytong/Documents/MATLAB/For_James_Larkin');
    pTxPath = '/Users/ytong/Documents/MATLAB/Temp/20171013_F7T_2013_40_387';
    dt = Spectro.dicomTree('dir',pTxPath,'recursive',false);
    ptxFMObj = DicomFM.WTCpTxFieldmaps(dt,'B1String','dt_dream_wIce_60deg_150VRef__B1',...
        'B0String','fieldmap_ptx7t_iso4mm_trans_RL','LocaliserString','fl_tof');
    %ptxFMObj = DicomFM.WTCpTxFieldmaps(pTxPath);
    SliceIdx = input('Please enter the slice number: \n');
    ptxFMObj.interpolateTo('Localiser');
    ptxFMObj.createMask(@(x) DicomFM.maskFunctions.ellipseMask(x),true);  
    %ptxFMObj.setSlice(10);
    %%  Specifying parameters
    param.targetFlipAngle = 15;
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
    RFStruct = tyMakeHanning(600,5);

    
    %%  Constructing system matrix
    disp('Calculating system matrix...');
    SysMatMode = 'Full';
    if strcmp(SysMatMode,'Approximation');
        AFull = getAMatSimp(RFStruct,maskedMaps.b1SensMasked,maskedMaps.b0MapMasked,...
        maskedMaps.posVox); %Nv-(2*Nc) sys mtx, rad/V
    elseif strcmp(SysMatMode,'Full');
        gr = zeros(numel(RFStruct.RF_pulse),3);
        rfOn = true(size(gr,1),1);
        tp = (0:(numel(RFStruct.RF_pulse)-1))*10E-6;
        DA = genAMatFull(tp',rfOn,gr,maskedMaps.b1SensMaskedHz,...
            maskedMaps.b0MapMasked,maskedMaps.b0MapMasked);
     %  To take the RF pulse shape into account. Construct a
     %  pseudo-diagonal matrix
        RFDiag = zeros(size(DA,2),8);
        RFSize = numel(RFStruct.RF_pulse);
        for iDx = 1:8
            RFDiag(((1:RFSize)+(iDx-1)*RFSize),iDx) = RFStruct.RF_pulse';
        end
        AFull = DA*RFDiag;
    end
    disp('complete')
    %%  Running RF shimming Step 1: Variable-exchange method
    disp('Running RF shimming variable-exchange optimization...')

    tikhonovArray = power(10,-8:-1);
    bVE = zeros(8,numel(tikhonovArray));
    for iDx = 1:numel(tikhonovArray)
        param.CGtikhonov = tikhonovArray(iDx);
        [bVE(:,iDx),~,~,~] = runVE(AFull,param,maskedMaps);
    end
    disp('Variable-exchange optimization finished')
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
    disp('Active-set optimization finished')    
    %%  Bloch Simulation
    disp('Running Bloch simulation...')
    bmin = bAS(1:8,minIndex).*exp(1i*bAS(9:16,minIndex));
    RFToSim = bsxfun(@times, bmin, repmat(RFStruct.RF_pulse,8,1));
    RFToSim = RFToSim';
    GToSim = zeros(numel(RFToSim)/8,3);
    %GToSim(:,3) = 425.8*RFStruct.G_amp*ones(numel(RFToSim)/8,1);     %from mT/m to Hz/cm
    B1ToSim = ptxFMObj.getB1PerV('Hz','delete');
    magnetization = make_blochSim(RFToSim,B1ToSim,...
        maskedMaps.b0MapMasked,GToSim,10E-6,maskedMaps.posVox,maskedMaps.mask);
    disp('Bloch simulation complete')
    %%  Plotting
% 
%  figure(56)
%   imagesc(maskedMaps.localiser(:,:,SliceIdx));
%   title('Select Rect Regions to Plot')
%   h2plot = imrect;
%   wait(h2plot);
%   plotMask = logical(createMask(h2plot));
%   [rol,cow] = find(plotMask);
%   plotRolArray = rol(1):rol(end);
%   plotCowArray = cow(1):cow(end);   
    %bTmp = bAS(1:8,1).*exp(1i*bAS(9:16,1));
  mVec = AFull*bmin;
  m = zeros(size(maskedMaps.mask));
  m(maskedMaps.mask) = mVec;
  
  
  figure(57)
  
  subplot(1,3,1);imagesc(abs(magnetization.mxy(plotRolArray,plotCowArray,SliceIdx-1)));title(sprintf('Slice %d',SliceIdx-1));colorbar
  subplot(1,3,2);imagesc(abs(magnetization.mxy(plotRolArray,plotCowArray,SliceIdx)));title(sprintf('Slice %d',SliceIdx));colorbar 
  subplot(1,3,3);imagesc(abs(magnetization.mxy(plotRolArray,plotCowArray,SliceIdx+1)));title(sprintf('Slice %d',SliceIdx+1));colorbar
   
  figure(58)
  
  subplot(1,3,1);imagesc(abs(m(plotRolArray,plotCowArray,SliceIdx-1)));title(sprintf('Slice %d',SliceIdx-1));colorbar
  subplot(1,3,2);imagesc(abs(m(plotRolArray,plotCowArray,SliceIdx)));title(sprintf('Slice %d',SliceIdx));colorbar 
  subplot(1,3,3);imagesc(abs(m(plotRolArray,plotCowArray,SliceIdx+1)));title(sprintf('Slice %d',SliceIdx+1));colorbar
  
%   figure(59)
%   subplot(1,2,1);imagesc(abs(m(plotRolArray,plotCowArray,SliceIdx)));title('Central Slice');colorbar
%   subplot(1,2,2);imagesc(abs(abs(m(plotRolArray,plotCowArray,SliceIdx)) - deg2rad(param.targetFlipAngle)*abs(maskedMaps.mask(:,:,10))))
%   title('Central Slice Error');colorbar
%     
    
%%
 b1map = ptxFMObj.getB1PerV('uT','NaN');
 b1map = squeeze(sum(b1map,4));
 figure(54);imagesc(abs(b1map(plotRolArray,plotCowArray,10)))
 colorbar