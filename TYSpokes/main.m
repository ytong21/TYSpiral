    %%    
   
    pTxPath = '/Users/ytong/Documents/Data/20170509_F7T_2013_40_345';
    ptxFMObj = DicomFM.WTCpTxFieldmaps(pTxPath,'B1String','dt_dream_wIce_60deg_80VRef__B1','B0String','fieldmap_ptx7t_iso4mm_trans_RL');
    %ptxFMObj = DicomFM.WTCpTxFieldmaps(pTxPath);
    ptxFMObj.interpolateTo('B1');
    ptxFMObj.createMask(@DicomFM.maskFunctions.multithreshMask,true);
    %%
    
    ptxFMObj.setSlice(16);
    %%
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
    %maskedMaps.targetAngle = ones(size(maskedMaps.b0MapMasked))*deg2rad(options.targetFlipAngle);
    SINC = tyMakeSinc(1000,6,5); %Duration = 1000us; number of zero crossings = 6; SliceThickness = 5 %mm

    %%
    tic
        [Landscape,OutputStruct] = landscape_run(0.1,param,maskedMaps,SINC);
    toc
    
    %%
    
    tic;Cart11 = cartesian_run(11,param,maskedMaps,SINC);toc
    %%
    tic;[Cart5,Cart5_Output] = cartesian_run(5,param,maskedMaps,SINC);toc
    %%
    figure(102)
    clf;imagesc(-3:0.05:3,-3:0.05:3,Landscape.NRMSE);colorbar
    hold on; plot(Cart11.kIn(1,:),Cart11.kIn(2,:),'.k', 'MarkerSize', 18);
    hold on; plot(Cart11.kOut(1,:),Cart11.kOut(2,:),'+r', 'MarkerSize', 10);
    title('11 by 11 Cartesian Grid Optimization')
    %%
    figure(107)
    clf;imagesc(-3:0.05:3,-3:0.05:3,Landscape.NRMSE);colorbar
    hold on; plot(Cart5.kIn(1,:),Cart5.kIn(2,:),'.k', 'MarkerSize', 18);
    cvals = zeros(6,3);
    cvals(2,:) = [1 0 0];
    cvals(6,:) = [1 1 0];
    for iDx = 1:numel(Cart5.ExitFlag)
    colorVecs(iDx,:) = cvals(Cart5.ExitFlag(iDx)+1,:);
    
    plot(Cart5.kOut(1,iDx).',Cart5.kOut(2,iDx).','+', 'MarkerSize', 10,'color',colorVecs(iDx,:));
    end
    for ii = 1:size(Cart5.kOut,2)
        plot([Cart5.kIn(1,ii); Cart5.kOut(1,ii)],[Cart5.kIn(2,ii); Cart5.kOut(2,ii)],'-g', 'MarkerSize', 18);
    end
    title('5 by 5 Cartesian Grid Optimization')
    %%
     Prob = zeros(20,1);
     global_min = min(Landscape.NRMSE(:));
     NRMSE_Vec = zeros(1000,30);
     %%
     tic
     for kk = 20
        
         NRMSE_Vec(:,kk) = rand_run(1000,kk,param,maskedMaps,SINC);
         Ref = 1.05*ones(size(NRMSE_Vec))*global_min;
         Prob(kk) = nnz(Ref>NRMSE_Vec)/numel(Ref);
     end
     toc
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
b1map = ptxFMObj.getB1PerV('uT','NaN');
b1map = squeeze(sum(b1map,3));
imagesc(abs(b1map))
colorbar

%%
loc = ptxFMObj.getLoc('none');
figure(97)
imagesc(loc)