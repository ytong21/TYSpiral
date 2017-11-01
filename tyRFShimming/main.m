    %%   Loading the data and masking
    addpath('/Users/ytong/Documents/MATLAB/For_James_Larkin');
    pTxPath = '/Users/ytong/Documents/MATLAB/Temp/20171013_F7T_2013_40_387';
    dt = Spectro.dicomTree('dir',pTxPath,'recursive',false);
    ptxFMObj = DicomFM.WTCpTxFieldmaps(dt,'B1String','dt_dream_wIce_60deg_150VRef__B1',...
        'B0String','fieldmap_ptx7t_iso4mm_trans_RL','LocaliserString','fl_tof');
    %ptxFMObj = DicomFM.WTCpTxFieldmaps(pTxPath);
    SliceIdx = input('Please enter the slice number: \n');
    ptxFMObj.interpolateTo('Localiser');
    ptxFMObj.createMask(@(x) DicomFM.maskFunctions.ellipseMask(x,SliceIdx),true);  
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
    maskedMaps.posVox = ptxFMObj.getPositions('cm','delete').';
    maskedMaps.localiser = ptxFMObj.getLoc('none');
    
    maskedMaps.mask = ptxFMObj.getMask();
    RFStruct = tyMakeHanning(600,5);
    %%  Running RF shimming Step 1: Variable-exchange method
    disp('Running RF shimming variable-exchange optimization...')
    AFull = getAMatSimp(RFStruct,maskedMaps.b1SensMasked,maskedMaps.b0MapMasked,...
    maskedMaps.posVox); %Nv-(2*Nc) sys mtx, rad/V
    tikhonovArray = power(10,-5:-1);
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
    %%  Bloch Simulation and Plotting

  figure(56)
  title('Select Rect Regions to Plot')
  imagesc(maskedMaps.localiser(:,:,SliceIdx));
  h2plot = imrect;
  wait(h2plot);
  close(56)
  plotMask = logical(createMask(h2plot));
  [rol,cow] = find(plotMask);
  plotRolArray = rol(1):rol(end);
  plotCowArray = cow(1):cow(end);
    
  bmin = bAS(1:8,minIndex).*exp(1i*bAS(9:16,minIndex));
  mVec = AFull*bmin;
  m = zeros(size(maskedMaps.mask));
  m(maskedMaps.mask) = mVec;
  
  
  figure(57)
  NonZeroIdx = find(maskedMaps.mask);
  subplot(1,3,1);imagesc(abs(m(plotRolArray,plotCowArray,SliceIdx-1)));title(sprintf('Slice %d',SliceIdx-1));colorbar
  subplot(1,3,2);imagesc(abs(m(plotRolArray,plotCowArray,SliceIdx)));title(sprintf('Slice %d',SliceIdx));colorbar 
  subplot(1,3,3);imagesc(abs(m(plotRolArray,plotCowArray,SliceIdx+1)));title(sprintf('Slice %d',SliceIdx+1));colorbar
  
  figure(58)
  subplot(1,2,1);imagesc(abs(m(plotRolArray,plotCowArray,SliceIdx)));title('Central Slice');colorbar
  subplot(1,2,2);imagesc(abs(abs(m(plotRolArray,plotCowArray,SliceIdx)) - deg2rad(param.targetFlipAngle)*abs(maskedMaps.mask(:,:,10))))
  title('Central Slice Error');colorbar
    
    