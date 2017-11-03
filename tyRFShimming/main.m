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
    param.TR = 1.5e-3;% sec
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
    SysMatMode = 'Approximation';
    if strcmp(SysMatMode,'Approximation');
        AFull = getAMatSimp(RFStruct,maskedMaps.b1SensMasked,maskedMaps.b0MapMasked); %Nv-(2*Nc) sys mtx, rad/V
    elseif strcmp(SysMatMode,'Full');
        gr = zeros(numel(RFStruct.RF_pulse),3);
        rfOn = true(size(gr,1),1);
        tp = 10E-6 * ones(size(rfOn));
        DA = genAMatFull(tp,rfOn,gr,maskedMaps.b1SensMaskedHz,...
            maskedMaps.b0MapMasked,maskedMaps.posVox);
     %  To take the RF pulse shape into account. Construct a
     %  pseudo-diagonal matrix
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

    tikhonovArray = power(10,-8:-1);
    bVE = zeros(8,numel(tikhonovArray));
    for iDx = 1:numel(tikhonovArray)
        param.CGtikhonov = tikhonovArray(iDx);
        [bVE(:,iDx),NRMSETmp,~,~] = runVE(AFull,param,maskedMaps);
        disp(NRMSETmp)
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
    bmin = bVE(:,minIndex);
    %bmin = bAS(1:8,minIndex).*exp(1i*bAS(9:16,minIndex));
    RFToSim = bsxfun(@times, (bmin), repmat(RFStruct.RF_pulse,8,1));
    RFToSim = RFToSim';
    GToSim = zeros(numel(RFToSim)/8,3);
    %GToSim(:,3) = 425.8*RFStruct.G_amp*ones(numel(RFToSim)/8,1);     %from mT/m to Hz/cm
    B1ToSim = ptxFMObj.getB1PerV('Hz','delete');
    magnetization = make_blochSim(RFToSim,B1ToSim,...
        maskedMaps.b0MapMasked,GToSim,10E-6,maskedMaps.posVox,maskedMaps.mask);
    disp('Bloch simulation complete')
    %  Plotting
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
  mVec = AFull*(bmin);
  m = zeros(size(maskedMaps.mask));
  m(maskedMaps.mask) = mVec;
  % Scaling
  m = rad2deg(m);
  %mBloch = asind(abs(magnetization.mxy));
  mBloch = magnetization.mxy;
  
  figure(57)
  PlotFALim = [10 25];
  
  subplot(1,3,1);imagesc(abs(mBloch(plotRolArray,plotCowArray,SliceIdx-1)));title(sprintf('Slice %d',SliceIdx-1));colorbar
  subplot(1,3,2);imagesc(abs(mBloch(plotRolArray,plotCowArray,SliceIdx)));title(sprintf('Slice %d',SliceIdx));colorbar 
  subplot(1,3,3);imagesc(abs(mBloch(plotRolArray,plotCowArray,SliceIdx+1)));title(sprintf('Slice %d',SliceIdx+1));colorbar
   
  figure(58)
  
  subplot(1,3,1);imagesc(abs(m(plotRolArray,plotCowArray,SliceIdx-1)),PlotFALim);title(sprintf('Slice %d',SliceIdx-1));colorbar
  subplot(1,3,2);imagesc(abs(m(plotRolArray,plotCowArray,SliceIdx)),PlotFALim);title(sprintf('Slice %d',SliceIdx));colorbar 
  subplot(1,3,3);imagesc(abs(m(plotRolArray,plotCowArray,SliceIdx+1)),PlotFALim);title(sprintf('Slice %d',SliceIdx+1));colorbar
  
%   figure(59)
%   subplot(1,2,1);imagesc(abs(m(plotRolArray,plotCowArray,SliceIdx)));title('Central Slice');colorbar
%   subplot(1,2,2);imagesc(abs(abs(m(plotRolArray,plotCowArray,SliceIdx)) - deg2rad(param.targetFlipAngle)*abs(maskedMaps.mask(:,:,10))))
%   title('Central Slice Error');colorbar
%     
    
%%
 b1map = ptxFMObj.getB1PerV('Hz','NaN');
 b1mapCP = squeeze(sum(b1map,4));
 figure(54);imagesc(abs(b1map(:,:,11)))
 colorbar
 %%
 subplot(4,2,1);imagesc(abs(b1map(:,:,10)))
 
 %%
 B1ToSim = ptxFMObj.getB1PerV('Hz','delete');
 
 DA = genAMatFull(600e-6,true,false,maskedMaps.b1SensMaskedHz,...
            maskedMaps.b0MapMasked,maskedMaps.posVox);
        singleSLiceMask = maskedMaps.mask(:,:,9:11);
  fullImg = zeros(size(singleSLiceMask));      
 figure (1234)
 clf
 
 for iDx = 1
    %voltvec = 50*ones(1,8);
    %voltvec(iDx) = 10;
    voltvec = bmin.';
    RFToSim = bsxfun(@times,repmat(RFStruct.RF_pulse.',1,8),voltvec);
    magnetization1 = make_blochSim(RFToSim,B1ToSim,... 
        maskedMaps.b0MapMasked,zeros(size(RFToSim,1),3),10E-6,maskedMaps.posVox,maskedMaps.mask);
    
    subplot(1,3,1)
    imagesc(abs(magnetization1.mxy(:,:,10)))
    caxis([0 sind(30)])
  
    
    subplot(1,3,2)
    amatMag = DA*(voltvec.'*(sum(RFStruct.RF_pulse)/numel(RFStruct.RF_pulse)));
    fullImg(singleSLiceMask(:)) = amatMag;
       imagesc(abs(fullImg(:,:,2)))
    caxis([0 sind(30)])
    
    subplot(1,3,3)
  
       imagesc(abs(fullImg(:,:,2))-abs(magnetization1.mxy(:,:,10)))
       colorbar
    caxis([0 sind(30)/20])
    drawnow()
    pause(0.1)
 end
 
 
 %%
     RFToSim = bsxfun(@times,repmat(RFStruct.RF_pulse.',1,8),voltvec);
    magnetization = make_blochSim(RFToSim,B1ToSim,... 
        maskedMaps.b0MapMasked,zeros(size(RFToSim,1),3),10E-6,maskedMaps.posVox,maskedMaps.mask);
    
    
    
        magnetization = make_blochSim(RFToSim,B1ToSim,...
        maskedMaps.b0MapMasked,GToSim,10E-6,maskedMaps.posVox,maskedMaps.mask);
    
    
    B1ToSim = ptxFMObj.getB1PerV('Hz','delete');    %WTC
     B1ToSim = ptxFMObj.getB1PerV('Hz','delete');   %TY