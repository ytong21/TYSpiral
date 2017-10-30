    %%    
    addpath('/Users/ytong/Documents/MATLAB/For_James_Larkin');
    pTxPath = '/Users/ytong/Documents/MATLAB/Temp/20171013_F7T_2013_40_387';
    dt = Spectro.dicomTree('dir',pTxPath,'recursive',false);
    ptxFMObj = DicomFM.WTCpTxFieldmaps(dt,'B1String','dt_dream_wIce_60deg_150VRef__B1',...
        'B0String','fieldmap_ptx7t_iso4mm_trans_RL','LocaliserString','fl_tof');
    %ptxFMObj = DicomFM.WTCpTxFieldmaps(pTxPath);
    ptxFMObj.interpolateTo('Localiser');
    ptxFMObj.createMask(@(x) DicomFM.maskFunctions.ellipseMask(x),true);  
    %ptxFMObj.setSlice(10);
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
    

    %%
    tic
        SINC_10kHz = tyMakeSinc(1000,10,5);
        [Landscape_10kHz,OutputStruct_10kHz] = landscape_run(0.05,param,maskedMaps,SINC_10kHz);
    toc
    %%
    tic
        SINC_8kHz = tyMakeSinc(1000,8,5);
        [Landscape_8kHz,OutputStruct_8kHz] = landscape_run(0.05,param,maskedMaps,SINC_8kHz);
    toc
    %%
    tic
        SINC_6kHz = tyMakeSinc(1000,6,5);
        [Landscape_6kHz,OutputStruct_6kHz] = landscape_run(0.05,param,maskedMaps,SINC_6kHz);
    toc
    %%
    tic
        SINC_2kHz = tyMakeSinc(1000,2,5);
        [Landscape_2kHz,OutputStruct_2kHz] = landscape_run(0.05,param,maskedMaps,SINC_2kHz);
    toc
        
    
    %%
    
    [Cart5_10kHz,Cart5_10kHz_Output] = cartesian_run(5,param,maskedMaps,SINC_10kHz);
    [Cart5_8kHz,Cart5_8kHz_Output] = cartesian_run(5,param,maskedMaps,SINC_8kHz);
    [Cart5_6kHz,Cart5_6kHz_Output] = cartesian_run(5,param,maskedMaps,SINC_6kHz);
    %%
    tic;[Cart5_2kHz,Cart5_2kHz_Output] = cartesian_run(5,param,maskedMaps,SINC_2kHz);toc
    %%
    
    tic;Cart11 = cartesian_run(11,param,maskedMaps,SINC);toc
    
    tic;[Cart5,Cart5_Output] = cartesian_run(5,param,maskedMaps,SINC);toc
    %%
    figure(102)
    clf;imagesc(-3:0.05:3,-3:0.05:3,Landscape.NRMSE);colorbar
    hold on; plot(Cart11.kIn(1,:),Cart11.kIn(2,:),'.k', 'MarkerSize', 18);
    hold on; plot(Cart11.kOut(1,:),Cart11.kOut(2,:),'+r', 'MarkerSize', 10);
    title('11 by 11 Cartesian Grid Optimization')
    %%
    figure(107)
    clf;imagesc(-3:0.05:3,-3:0.05:3,Landscape_2kHz.NRMSE);colorbar
    hold on; plot(Cart5_2kHz.kIn(1,:),Cart5_2kHz.kIn(2,:),'.k', 'MarkerSize', 18);
    cvals = zeros(6,3);
    cvals(2,:) = [1 0 0];
    cvals(6,:) = [1 1 0];
    for iDx = 1:numel(Cart5_2kHz.ExitFlag)
    colorVecs(iDx,:) = cvals(Cart5_2kHz.ExitFlag(iDx)+1,:);
    
    plot(Cart5_2kHz.kOut(1,iDx).',Cart5_2kHz.kOut(2,iDx).','+', 'MarkerSize', 10,'color',colorVecs(iDx,:));
    end
    for ii = 1:size(Cart5_2kHz.kOut,2)
        plot([Cart5_2kHz.kIn(1,ii); Cart5_2kHz.kOut(1,ii)],[Cart5_2kHz.kIn(2,ii); Cart5_2kHz.kOut(2,ii)],'-g', 'MarkerSize', 18);
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
    ptxFMObj.setSlice(9);
 b1map = ptxFMObj.getB1PerV('uT','NaN');
 b1map = squeeze(sum(b1map,3));
 figure(54);imagesc(abs(b1map))
 colorbar


%%
% To draw 2 ellipsoid ROIs
loc = ptxFMObj.getLoc('none');
figure(97)
imagesc(loc)

h1 = imellipse;
wait(h1);
h2 = imellipse;
wait(h2);
MaskElip = logical(createMask(h1)+createMask(h2));
