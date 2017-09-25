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
    param.FOX = 25;%25cm

    maskedMaps.b1SensMasked = ptxFMObj.getB1PerV('uT','delete');
    maskedMaps.b0MapMasked = ptxFMObj.getB0('Hz','delete');
    maskedMaps.posVox = ptxFMObj.getPositions('cm','delete').';
    maskedMaps.mask = ptxFMObj.getMask();
    %maskedMaps.targetAngle = ones(size(maskedMaps.b0MapMasked))*deg2rad(options.targetFlipAngle);
    SINC = tyMakeSinc(1000,6,5); %Duration = 1000us; number of zero crossings = 6; SliceThickness = 5 %mm

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