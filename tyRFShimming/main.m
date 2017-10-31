    %%   Loading the data and masking
    addpath('/Users/ytong/Documents/MATLAB/For_James_Larkin');
    pTxPath = '/Users/ytong/Documents/MATLAB/Temp/20171013_F7T_2013_40_387';
    dt = Spectro.dicomTree('dir',pTxPath,'recursive',false);
    ptxFMObj = DicomFM.WTCpTxFieldmaps(dt,'B1String','dt_dream_wIce_60deg_150VRef__B1',...
        'B0String','fieldmap_ptx7t_iso4mm_trans_RL','LocaliserString','fl_tof');
    %ptxFMObj = DicomFM.WTCpTxFieldmaps(pTxPath);
    ptxFMObj.interpolateTo('Localiser');
    ptxFMObj.createMask(@(x) DicomFM.maskFunctions.ellipseMask(x),true);  
    %ptxFMObj.setSlice(10);
    %%  Specifying parameters
    param.targetFlipAngle = 15;
    param.numCh = 8;
    param.TR = 5e-3;% sec
    param.CGtikhonov = 1e-6;
    param.tol = 1e-5;
    param.MaxEvaluation = 25000;
    param.FOX = 25;%25cm

    maskedMaps.b1SensMasked = ptxFMObj.getB1PerV('uT','delete');
    maskedMaps.b0MapMasked = ptxFMObj.getB0('Hz','delete');
    maskedMaps.posVox = ptxFMObj.getPositions('cm','delete').';
    
    maskedMaps.mask = ptxFMObj.getMask();
    RFStruct = tyMakeHanning(1000,5);
    %%  Running RF shimming Step 1: Variable-exchange method
    AFull = getAMatSimp(RFStruct,maskedMaps.b1SensMasked,maskedMaps.b0MapMasked,...
    maskedMaps.posVox); %Nv-(2*Nc) sys mtx, rad/V
    tikhonovArray = power(10,-7:-2);
    bVE = zeros(8,numel(tikhonovArray));
    for iDx = 1:numel(tikhonovArray)
        param.CGtikhonov = tikhonovArray(iDx);
        [bVE(:,iDx),~,~,~] = runVE(AFull,param,maskedMaps);
    end
    %%  Running RF shimming Step 2: Active-set method
    bAS = zeros(16,numel(tikhonovArray));
    for iDx = 1:numel(tikhonovArray)
        param.CGtikhonov = tikhonovArray(iDx);
        [bAS(:,iDx),~,~,~] = runAS(bVE(:,iDx),RFStruct,maskedMaps,param,AFull);
    end    
    
    
    