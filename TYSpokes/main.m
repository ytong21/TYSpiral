    %%    
   
    pTxPath = '/Users/ytong/Documents/MATLAB/Temp/20170509_F7T_2013_40_346';
    ptxFMObj = DicomFM.WTCpTxFieldmaps(pTxPath,'B1String','dt_dream_wIce_60deg_80VRef__B1','B0String','fieldmap_ptx7t_iso4mm_trans_RL');
    ptxFMObj.interpolateTo('B1');
    ptxFMObj.createMask(@DicomFM.maskFunctions.multithreshMask,true);
    ptxFMObj.setSlice(16);
    %%
    param.targetFlipAngle = 15;
    param.numCh = 8;
    param.TR = 5e-3;% sec
    param.CGtikhonov = 1e-6;
    param.tol = 1e-2;
    param.MaxEvaluation = 3000;

    maskedMaps.b1SensMasked = ptxFMObj.getB1PerV('uT','delete');
    maskedMaps.b0MapMasked = ptxFMObj.getB0('Hz','delete');
    maskedMaps.posVox = ptxFMObj.getPositions('cm','delete').';
    maskedMaps.mask = ptxFMObj.getMask();
    %maskedMaps.targetAngle = ones(size(maskedMaps.b0MapMasked))*deg2rad(options.targetFlipAngle);
    
    %%
    deltaXAdj = -3:0.1:3;
    deltaYAdj = -3:0.1:3;
        StructToSave = struct('deltaX',deltaXAdj,'deltaY',deltaYAdj);
    FOX = 25; %25cm
        StructToSave.FOX = FOX;     StructToSave.CGtikhonov = param.CGtikhonov;
        StructToSave.tolerance = param.tol;     StructToSave.MaxEvaluation = param.MaxEvaluation;
    deltaX = deltaXAdj*(2*pi)/FOX;
    deltaY = deltaYAdj*(2*pi)/FOX;
    %
    Duration = 1000;        %us
    NumZeros = 6;           %zero crossings
    SliceThickness = 5;     %mm
    SINC = tyMakeSinc(Duration,NumZeros,SliceThickness);    %struct

     NumIt = 0;
     deltaKVec = zeros(2,numel(deltaX)*numel(deltaY));
     NRMSE = zeros(numel(deltaX)*numel(deltaY),1);
     bOut = zeros(numel(deltaX)*numel(deltaY),2*param.numCh);
for xDx = 1:numel(deltaX)
    for yDx = 1:numel(deltaY)
        NumIt = NumIt+1;
        deltaKVec(:,NumIt) = [deltaX(xDx);deltaY(yDx)];
    end
end
OptimType = 'b';
tic
parfor iDx = 1:size(deltaKVec,2)
    [bOut(iDx,:),NRMSE(iDx,:)] = VE_AS(OptimType,deltaKVec(:,iDx),SINC,maskedMaps,param);
end
toc
StructToSave.ComputationTime = toc;
    bOut = reshape(bOut,[numel(deltaX),numel(deltaY),2*param.numCh]);
    NRMSE = reshape(NRMSE,[numel(deltaX),numel(deltaY)]);
    StructToSave.bOut = bOut;   StructToSave.NRMSE = NRMSE;
    StructFileName = sprintf('Brain_Lamda_%.1d',param.CGtikhonov);
%% 
figure(102)
clf
imagesc(deltaXAdj,deltaYAdj,StructToSave.NRMSE)
colorbar

hold on
plot([0 1],[0 1],'ok', 'MarkerSize', 8);
%%
save(StructFileName,'StructToSave');
%%
         bComplex = (bOut(1:2*param.numCh)).*(cos(bOut(1+2*param.numCh:4*param.numCh))...
             +1i*sin(bOut(1+2*param.numCh:4*param.numCh)));
         finalMag = AFullSpokes*bComplex;   
         m = zeros(size(maskedMaps.mask));
         m(~maskedMaps.mask==0) = finalMag;
         figure(21);imagesc(abs(m))

%%
finalFA = AFullSpokes*100*[ones(4,1);zeros(12,1)];
%finalFA = AFullSpokes*100*[ones(8,1);ones(8,1)];
m(~maskedMaps.mask==0) = finalFA;
figure(21);imagesc(abs(m)); colorbar
%%
b1map = ptxFMObj.getB1PerV('uT','none');
b1map = squeeze(sum(b1map,3));
imagesc(abs(b1map))
colorbar