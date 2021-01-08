    %% Load in the data and interporlate to B0
    pTxPath = '/Users/ytong/Documents/Data/TRUST/20200723/20200723_Phantom_23_07';
    dt = Spectro.dicomTree('dir',pTxPath,'recursive',false);
    ptxFMObj = DicomFM.WTCpTxFieldmaps(dt,'B1String','dt_dream_wIce_60deg_90VRef__B1',...
        'B0String','fieldmap_ptx7t_iso4mm_trans_RL','LocaliserString','localizer','InterpTarget',...
        'B0');
    %%
    nii_dir = '/Users/ytong/Documents/Data/TRUST/20200316_nii';
    dicm2nii(pTxPath,'/Users/ytong/Documents/Data/TRUST/20200316_nii')
    
    %% prep steps for fsl
    PATH = getenv('PATH');
    setenv('PATH', [PATH ':/usr/local/bin:/usr/local/fsl/bin']);
    setenv('FSLOUTPUTTYPE','NIFTI_GZ')
    setenv('FSLMULTIFILEQUIT','TRUE')
    setenv('FSLTCLSH','/usr/local/fsl/bin/fsltclsh')
    setenv('FSLWISH','/usr/local/fsl/bin/fslwish')
    setenv('FSLDIR','/usr/local/fsl')
    %% Creating a mask using betted nii file
    %B0_nii = niftiread('/Users/ytong/Documents/Data/20191105_B0_nii/OilPhantom.nii.gz');
    %B0_Hz_nii = niftiread('/Users/ytong/Documents/Data/20191105_B0_nii/dt_fieldmap_ptx7t_iso4mm_trans_RL_s007.nii.gz');
    % For some reason the second dimension is flipped
    % This is established by comparing B0_Hz_nii and ptxFMObj.getB0Raw()
    
    % fsl cannot handle .gz and it nees to be decompressed
    B0_filename = 'dt_fieldmap_ptx7t_iso4mm_trans_RL_s034.nii.gz';
    %gunzip(B0_filename);
    system(sprintf('bet %s%s%s %s -m','/Users/ytong/Documents/Data/TRUST/20200316_nii',...
        filesep,strcat(nii_dir, filesep,'OilPhan0316.nii'),...
        'OilPhan0316'));
    
    %%
    mask_nii = niftiread('/Users/ytong/Documents/Data/TRUST/20200723/0723_nii/oil_0723_mask.nii.gz');
    mask_nii = flip(mask_nii,2);
    
    %% Specify the slice number and make a mask
    Slice_Array = 31:40;
    ptxFMObj.createMask(@(x) DicomFM.maskFunctions.TRUST_mask(x,Slice_Array,3,mask_nii),true);
    
    %ptxFMObj.createMask(@(x) DicomFM.maskFunctions.TRUST_mask(x,Slice_Array,3,mask_nii),true);
   
    %% Load essential information
    xz_res = 4;
    y_res = 4;    
    maskedMaps.b1SensMasked = ptxFMObj.getB1PerV('uT','delete');
    maskedMaps.b0MapMasked = ptxFMObj.getB0('Hz','delete');
    maskedMaps.b0MapMaskedRad = (2*pi)*maskedMaps.b0MapMasked;
    maskedMaps.posVox = ptxFMObj.getPositions('cm','delete').';
    maskedMaps.localiser = ptxFMObj.getLoc('none');
    maskedMaps.b1SensMaskedHz = ptxFMObj.getB1PerV('Hz','delete');
    maskedMaps.b0 = ptxFMObj.getB0('Hz','NaN');
    maskedMaps.mask = ptxFMObj.getMask();
    maskedMaps.mask_one_slice = maskedMaps.mask(:,:,Slice_Array);
    
    %% Solve the issue of NANs in b1 map
    [Nan_Row, Nan_Col] = find(isnan(maskedMaps.b1SensMaskedHz));
    for iDx = 1:numel(Nan_Row)
        maskedMaps.b1SensMaskedHz(Nan_Row(iDx),Nan_Col(iDx)) = ...
            maskedMaps.b1SensMaskedHz(Nan_Row(iDx)-1,Nan_Col(iDx));
        maskedMaps.b1SensMasked(Nan_Row(iDx),Nan_Col(iDx)) = ...
            maskedMaps.b1SensMasked(Nan_Row(iDx)-1,Nan_Col(iDx));
    end
    
    %% Constructing an excitation target
    
    PlotSingleSlice(squeeze(maskedMaps.b0(:,:,Slice_Array(5))),xz_res,y_res,90,true);
    %%
    PlotSingleSlice(squeeze(maskedMaps.b0(:,:,Slice_Array(5))),xz_res,y_res,0,false);
    %%
    hTmp = drawellipse; %Do not close this figure until next section
    
    %%
    CircleMask = createMask(hTmp);
    CircleMaskFiltered = imgaussfilt(single(CircleMask),[1 1]);
    %%
    TargetTemp = zeros(size(maskedMaps.b0));
    for iDx = 1:numel(Slice_Array)
        TargetTemp(:,:,Slice_Array(iDx)) = squeeze(maskedMaps.mask(:,:,Slice_Array(iDx))) - CircleMaskFiltered;
    end
    %%
    %PlotSingleSlice(squeeze(TargetTemp(:,:,Slice_Array)),xz_res,y_res,90,true);
    maskedMaps.target_masked = TargetTemp(maskedMaps.mask);
    % maskedMaps.target_masked is the target mag in the slices in interest
    % most of them will be 1, few voxes will be close to 0;
    % YT change 21/07/2020. Select only a circle for debugging.
    %maskedMaps.target_masked = 1-TargetTemp(maskedMaps.mask);
    %% creating a mask for the ROI whose blood T2 we are interested in
    
    maskedMaps.ROImask = zeros(size(maskedMaps.b0));
    for iDx = 1:numel(Slice_Array)
        maskedMaps.ROImask(:,:,Slice_Array(iDx)) = CircleMaskFiltered;
    end
    maskedMaps.ROImask = logical(maskedMaps.ROImask);
    %%
    temp1 = ptxFMObj.getB0('Hz','none');        % Hz
    MLEVRoi.b0 = temp1(maskedMaps.ROImask);
    temp2 = ptxFMObj.getB1PerV('Hz','none');    % Hz/V
    for iDx = 1:size(temp2,4)
        temp3 = temp2(:,:,:,iDx);
        MLEVRoi.b1(:,iDx) = temp3(maskedMaps.ROImask);
    end

   %%
   T_init = makeTraj();
   tyPlotTrajectory( T_init ); set(gcf,'name','Initial k-space Trajectory')
   fprintf('Pulse duration = %.2f ms\n', T_init.t(end));
   %plot3(kx_vec,ky_vec,kz_vec);
   
   %% Calculate the RF pulse using 1us resolution with VE method
   Gradient = 425.77*downsample([T_init.Sam.Gx; T_init.Sam.Gy; T_init.Sam.Gz]',10);
   Time_Vec = downsample(T_init.Sam.t,10);
   %% AFull magnetisation per V
   AFullMag = genAMatFull(1E-5*ones(numel(Time_Vec),1),ones(numel(Time_Vec),1),Gradient,maskedMaps.b1SensMaskedHz,...
            maskedMaps.b0MapMasked,maskedMaps.posVox);
   AFull = asin(AFullMag);
   clear AFullMag
   %% Extending the RF pulse while keeping the number of points the same.
   Time_Vec_extended = Time_Vec*2;
   AFull_Extended = genAMatFull(1E-5*ones(numel(Time_Vec_extended),1),ones(numel(Time_Vec_extended),1),Gradient,maskedMaps.b1SensMaskedHz,...
            maskedMaps.b0MapMasked,maskedMaps.posVox);
   %% Building ALambda matrix. This is time-consuming. This is taken outside of the VE function
    param.targetFlipAngle = 90;
    param.numCh = 8;
    param.tol = 1e-5;
    %param.CGtikhonov = 1e-6;
    %tikhonovArray = power(10,-6:1:-4);
    tikhonovArray = power(10,-6);       % -6 for 0707
    ALambda_cell = cell(size(tikhonovArray));
    %% 
    % run time: "\" < lsqminnorm < pinv
    for iDx = 1:numel(tikhonovArray)
        ALambda_cell{iDx} = ((AFull'*AFull + tikhonovArray(iDx)*eye(size(AFull,2)))\sparse(eye(5736)))...
            *AFull';
    end
    %%
    OutCell = cell(size(tikhonovArray));
    for iDx = 1:numel(tikhonovArray)
        param.CGtikhonov = tikhonovArray(iDx);
        
        OutCell{iDx} = run_variable_exchange(AFull,param,maskedMaps,ALambda_cell{iDx});
        %OutCell{3} = run_variable_exchange(AFull,param,maskedMaps,ALambda_cell{3});
    end
%     %% Plotting the target and final mag
%     close all
%     PlotSingleSlice(TargetTemp(:,:,Slice_Array),xz_res,xz_res,90,true);
%     PlotSingleSlice(abs(OutCell{1}.finalMag),xz_res,xz_res,90,true);
%     PlotSingleSlice(TargetTemp(:,:,Slice_Array)-1-abs(finalMag),xz_res,xz_res,90,true);
%     %%
%     [U1,S1,V1] = svd(AFull);
%     %%
%     SysMat_SVD.U = U1(:,1:600); SysMat_SVD.S = sparse(S1(1:600,1:600)); SysMat_SVD.V = V1(:,1:600);
%     %%
%     out = run_active_set(OutCell{2}.bOut,maskedMaps,param,AFull);
%     
%     %%  
%     out_singval = run_SingVal(OutCell{2}.bOut,maskedMaps,param,AFull);
%     
%% Trying to smooth the RF pulse
[bSmooth,bSmooth_by_chan] = fileter_RF(OutCell{1}.bOut);
Smooth_finalmag = MakeFullImg(maskedMaps.mask_one_slice,AFull*bSmooth);
%%
run_comparison_plot(CircleMaskFiltered,Smooth_finalmag,maskedMaps);

%%
plotRF(bSmooth)
%%
writeIniFile_ty(bSmooth_by_chan,gradIn);
    %%
    run_comparison_plot(CircleMaskFiltered,OutCell{1}.finalMag,maskedMaps);
    %%
      %% Write the pulse into .ini file
    %rfIn = 1.2*reshape(OutCell{1}.bOut,[numel(OutCell{1}.bOut)/8 8]);
    rfIn = reshape(OutCell{1}.bOut,[numel(OutCell{1}.bOut)/8 8]);
    gradIn = Gradient;
    
    %%
    writeIniFile_ty(rfIn,gradIn)
    
    %% Bloch sim
    mag = BlochSimArb(rfIn,gradIn,maskedMaps);
    bloch_img=MakeFullImg(maskedMaps.mask_one_slice,mag.mz);

    %% Bloch sim for diff voltages
    Vol_Array = 100:10:250;
    
    for iDx=1:numel(Vol_Array)
        rf_2_sim = Vol_Array(iDx)*rfIn/max(abs(rfIn(:)));
        mag_temp = BlochSimArb(rf_2_sim,gradIn,maskedMaps);
        img_temp=MakeFullImg(maskedMaps.mask_one_slice,mag_temp.mz);
        sim_img(:,:,:,iDx) = img_temp;
    end
    
    %% Bloch sim for diff time offset
    Offset_array = -10:2:10;
    % raster time =  10 us
    
    %% Writiing zeros in RF and Gradient arrays to verify the performance of spoiler
    writeIniFile_ty(rfIn,0*gradIn)
    %% Writing a gaussian pulse to verify the performance of spoiler
    Gaussian_struct = tyMakeGaussian(600,5);
    rfIn_Gaussian = zeros(size(Gaussian_struct.RFOn));
    rfIn_Gaussian(Gaussian_struct.RFOn) = 120*[0 Gaussian_struct.RF_pulse 0]';
    rfIn_Gaussian = repmat(rfIn_Gaussian,[1 8]);
    gradIn_Gaussian = zeros(numel(Gaussian_struct.RFOn),3);
    gradIn_Gaussian(:,3) = 42.57E6*1e-5*Gaussian_struct.GradShape;  %mT/m to Gauss/cm
    writeIniFile_ty(rfIn_Gaussian,gradIn_Gaussian);
    
    %%
        %% Writing a sinc pulse to verify the performance of spoiler
    sinc_struct = tyMakeSinc(600,4,30);
    rfIn_sinc = zeros(size(sinc_struct.RFOn));
    rfIn_sinc(sinc_struct.RFOn) = 120*[0; sinc_struct.rfNormalized; 0];
    rfIn_sinc = repmat(rfIn_sinc,[1 8]);
    gradIn_sinc = zeros(numel(sinc_struct.RFOn),3);
    gradIn_sinc(:,1) = 42.57E6*1e-5*sinc_struct.GradShape;  %mT/m to Gauss/cm
    plot(10*(1:size(gradIn_sinc,1)),rfIn_sinc,10*(1:size(gradIn_sinc,1)),gradIn_sinc(:,1));
    writeIniFile_ty(rfIn_sinc,gradIn_sinc);
    
    %%
    rfIn_2us=zeros(size(rfIn,1)*2,size(rfIn,2));
    for iDx = 1:size(rfIn,2)
        rfIn_2us(:,iDx) = interp(rfIn(:,iDx),2);
    end
    
    %[RFTest,GradTest] = AdjustDelay(rfIn_2us, gradIn, -10);
    writeIniFile_ty(rfIn_2us,gradIn);
    
    %%    %%
    rfIn_10us=zeros(size(rfIn,1)*10,size(rfIn,2));
    for iDx = 1:size(rfIn,2)
        rfIn_10us(:,iDx) = interp(rfIn(:,iDx),10);
    end
    
    %[RFTest,GradTest] = AdjustDelay(rfIn_2us, gradIn, -10);
    writeIniFile_ty(rfIn_10us,gradIn)
    %% Change grad rf offset
    [RFTest,GradTest] = AdjustDelay(rfIn, gradIn, 10);
    writeIniFile_ty(RFTest,GradTest)

    
    %%
    figure(224)
    RF_pulse = reshape(OutCell{1}.bOut,[numel(OutCell{2}.bOut)/8 8]);
    set(gcf,'units','centimeters','position',[4 4 13 60],'paperunits','centimeters','paperposition',[4 4 13 60])

    RF_2_plot = [0; abs(RF_pulse(123:137,1)); 0];
    plot(RF_2_plot,'k','LineWidth',10)
    set(gcf,'color','w','InvertHardcopy','off')
    
    axis off
    %%
    sim_DiffV = sim_diffV(rfIn,gradIn,maskedMaps);
    %%
    delay_img = sim_delay(rfIn,gradIn,maskedMaps);
    
    %%
    OutCell_CP = calc_CP(AFull,tikhonovArray,param,maskedMaps);
    %%
        run_comparison_plot(CircleMaskFiltered,OutCell_CP{1}.finalMag,maskedMaps);
    %% Calculating the pulse for CP mode
    
   function OutCell = calc_CP(AFull,tikhonovArray,param,maskedMaps)
        [d1, d2] = size(AFull);
        AFull_reshaped = reshape(AFull,d1,d2/8,8);
        % voxels/RF time point/channels
        AFull_CP = sum(AFull_reshaped,3);
        ALambda_cell = cell(size(tikhonovArray));
            for iDx = 1:numel(tikhonovArray)
                ALambda_cell{iDx} = ((AFull_CP'*AFull_CP + 8*tikhonovArray(iDx)*...
                    eye(size(AFull_CP,2)))\sparse(eye(d2/8)))*AFull_CP';
            end
        OutCell = cell(size(tikhonovArray));
        for iDx = 1:numel(tikhonovArray)
            param_CP = param;
            param_CP.CGtikhonov = 8*tikhonovArray(iDx);
            OutCell{iDx} = run_variable_exchange(AFull_CP,param_CP,maskedMaps,ALambda_cell{iDx});
            %OutCell{3} = run_variable_exchange(AFull,param,maskedMaps,ALambda_cell{3});
        end            
    end
    
        %% Bloch sim for diff voltages
  function  sim_img = sim_diffV(rfIn,gradIn,maskedMaps)
    Vol_Array = 100:10:250;
    
    for iDx=1:numel(Vol_Array)
        rf_2_sim = Vol_Array(iDx)*rfIn/max(abs(rfIn(:)));
        mag_temp = BlochSimArb(rf_2_sim,gradIn,maskedMaps);
        img_temp=MakeFullImg(maskedMaps.mask_one_slice,mag_temp.mz);
        for jDx = 1:size(img_temp,3)
            sim_img(:,:,jDx,iDx) = img_temp(:,:,jDx).';
        end
    end
  end
    function  sim_img = sim_delay(rfIn,gradIn,maskedMaps)
    delay_array = -10:2:10;
    raster = 2;
    rf = interp1(0:size(rfIn,1),[zeros(1,size(rfIn,2)); rfIn],(1:5*size(rfIn,1))/5);
    grad = interp1(0:size(gradIn,1),[zeros(1,size(gradIn,2)); gradIn],(1:5*size(gradIn,1))/5);
    sample_2us = size(rf,1);
    for iDx=1:numel(delay_array)
        num_shift = delay_array(iDx)/raster;
        rf_2_sim = complex(zeros(size(rf,1)+5,size(rf,2)));
        grad_2_sim = zeros(size(grad,1)+5,size(grad,2));
        if num_shift<=0
            grad_2_sim(1:sample_2us,:) = grad;
            rf_2_sim((abs(num_shift)+1):(abs(num_shift)+sample_2us),:) = rf;
        else
            rf_2_sim(1:sample_2us,:) = rf;
            grad_2_sim((abs(num_shift)+1):(abs(num_shift)+sample_2us),:) = grad;
        end
        mag_temp = BlochSimArb(rf_2_sim,grad_2_sim,maskedMaps,raster);
        img_temp=MakeFullImg(maskedMaps.mask_one_slice,mag_temp.mz);
        for jDx = 1:size(img_temp,3)
            sim_img(:,:,jDx,iDx) = img_temp(:,:,jDx).';
        end
    end
  end
 function OutStruct = BlochSimArb(RF,Gradient,maskedMaps,varargin)
    %		mode= Bitmask mode:
    %		Bit 0:  0-Simulate from start or M0, 1-Steady State
    %		Bit 1:  1-Record m at time points.  0-just end time.
    if numel(varargin)==1
        dt = varargin{1}*1e-6;                 %10us
    else
        dt = 10e-6;
    end
    mode = 0;   
    %dt = 10e-6;                 %10us
    g = Gradient;
    sens = maskedMaps.b1SensMaskedHz;
    df = maskedMaps.b0MapMasked;       
    dp = maskedMaps.posVox;         % Assume isocentre
    t1 = inf;
    %t2 = 150e-3;                   % this value used for thesis
    [mx,my,mz] =  blochSim_mex_relax_multithread(RF,sens,g,dt,t1,t2,df,dp,mode);
    OutStruct = struct('mx',mx,'my',my,'mz',mz);
end
    %% plot the RF
    function plotRF(bOut)
    RF_pulse = reshape(bOut,[numel(bOut)/8 8]);
    figure(223)
    for iDx = 1:8
        set(gcf,'color','w','InvertHardcopy','off')
        set(gcf,'units','centimeters','position',[4 4 60 30],'paperunits','centimeters','paperposition',[4 4 60 30])
        subplot(2,4,iDx)
        plot((1:numel(bOut)/8)/100,abs(RF_pulse(:,iDx)),'LineWidth',1.2);ylim([0 280])
        title(sprintf('Channel %d',iDx))
        xlabel('Time (ms)')
        ylabel('Voltage (V)')
        set(gca,'FontSize',16)
    end
    end
    %%
    function [bSmooth,bSmooth_by_chan] = fileter_RF(bOut)
        b_by_chan = reshape(bOut,[numel(bOut)/8,8]);
        bSmooth_by_chan = zeros(size(b_by_chan));
        for iDx = 1:size(b_by_chan,2)
            bSmooth_by_chan(:,iDx) = smooth(b_by_chan(:,iDx));
        end
        bSmooth = reshape(bSmooth_by_chan,size(bOut));
    end
    %%
    function finalMag=MakeFullImg(mask_1_slice,ImgVec)
        finalMag = zeros(numel(mask_1_slice),1);
        finalMag(logical(mask_1_slice(:))) = ImgVec;
        finalMag = reshape(finalMag,size(mask_1_slice));
    end
       %%
       function T_init = makeTraj()
   options.Gmax         = 22.0;     % [mT/m]
   options.Smax         = 200.0;   	% YT[T/m/s]
   options.Vmax         = 150.0;    % [V]
   options.phi_bounds   = [   10.0000   10.0000   10.0000    2.0000    1.0000   -8.0000    0.5000 ; ...
                                     100.0000  100.0000  100.0000    7.5000   10.0000   +8.0000    2.0000 ].';
   options.GBF          =       struct(  'OnlyContinuousRanges',0,...
                                                        'GradSmooth_On',1,...
                                                        'GradSmooth_Lambda',1e-3,...
                                                        'GradSmooth_MaxDistFac',2,...
                                                        'verbosity',0,...
                                                        'UseMex',0,...
                                                        'SolverID','MinConF_SPG' );
   
   options.phi_init = [40  40  20   7    6    1    1 ].';
   [ kx_vec, ky_vec, kz_vec, R_mat, NoRFIndices ] = TrajectoryControlPoints_Shells(options.phi_init);
   NV = sum( R_mat(:,end) );
   G0 = zeros( NV,1 );
   T_init = PerformOptimization3D( kx_vec, ky_vec, kz_vec, R_mat, options.Gmax, options.Smax, G0, options.GBF );
       end
       
       function T_init = makeTraj_shells(kxy_extent,kz_extent,No_Shell,No_Rev)
   options.Gmax         = 22.0;     % [mT/m]
   options.Smax         = 200.0;   	% YT[T/m/s]
   options.Vmax         = 150.0;    % [V]
   options.phi_bounds   = [   10.0000   10.0000   10.0000    2.0000    1.0000   -8.0000    0.5000 ; ...
                                     100.0000  100.0000  100.0000    7.5000   10.0000   +8.0000    2.0000 ].';
   options.GBF          =       struct(  'OnlyContinuousRanges',0,...
                                                        'GradSmooth_On',1,...
                                                        'GradSmooth_Lambda',1e-3,...
                                                        'GradSmooth_MaxDistFac',2,...
                                                        'verbosity',0,...
                                                        'UseMex',0,...
                                                        'SolverID','MinConF_SPG' );
   
   options.phi_init = [kxy_extent  kxy_extent  kz_extent   No_Shell    No_Rev    1    1 ].';
   [ kx_vec, ky_vec, kz_vec, R_mat, NoRFIndices ] = TrajectoryControlPoints_Shells(options.phi_init);
   NV = sum( R_mat(:,end) );
   G0 = zeros( NV,1 );
   T_init = PerformOptimization3D( kx_vec, ky_vec, kz_vec, R_mat, options.Gmax, options.Smax, G0, options.GBF );
       end
    %%
    function run_comparison_plot(circle,result,maskedMaps)
        figure(222)
        set(gcf,'color','w','InvertHardcopy','off')
        set(gcf,'units','centimeters','position',[4 4 35 20],'paperunits','centimeters','paperposition',[4 4 35 20])
        %slice_array = 35:2:44;
        
        result_to_plot_array = 1:2:10;
        for iDx = 1:5
        ax1 = subplot(2,5,iDx);
        FAmap = rad2deg(abs(result(:,:,result_to_plot_array(iDx))));
        %FAmap([1:11 end-5:end],:) = [];
        imagesc(imrotate(flip(FAmap),-90),[0 100]);
        %imagesc(FAmap,[0 100]);
        %TT = title(sprintf('Slice %d',result_to_plot_array(iDx)));
        %TT.FontSize = 16;
        axis image; axis off;colormap(ax1, 'hot');
        text(2,4,sprintf('Slice %d FA',result_to_plot_array(iDx)),'Color','w','FontSize',16);
        if iDx == 5
            CLB1 = colorbar('FontSize',18);nudge(CLB1,[0.05 -0.048 0.01 0.095]);
            set(get(CLB1,'Title'),'String','deg')
            %ax1_pos = get(ax1,'position');
            %CLB1.Position(4) = ax1_pos(4);
            %CLB1.Position(2) = ax1_pos(2);
        end
        ax2 = subplot(2,5,iDx+5);
        target = squeeze(maskedMaps.mask_one_slice(:,:,result_to_plot_array(iDx))) - circle;
        diff = abs(rad2deg(result(:,:,result_to_plot_array(iDx))))-90*target;
        diff(diff == 0) = -inf;
        %diff([1:11 end-5:end],:) = [];
        imagesc(imrotate(flip(diff),-90),[-10 10]);
        %imagesc(diff,[-10 10]);
        %TT = title(sprintf('Slice %d diff',result_to_plot_array(iDx)));
        %TT.FontSize = 16;
        axis image; axis off; 
        nudge(ax2,[0 0.1 0 0]);
        text(2,4,sprintf('Slice %d Diff',result_to_plot_array(iDx)),'Color','w','FontSize',16);
        colormap(ax2, 'hot');
        if iDx == 5
            CLB2 = colorbar('FontSize',18);nudge(CLB2,[0.05 -0.048 0.01 0.095]);
            set(get(CLB2,'Title'),'String','deg')
        end

        end
        
        
    end
    %%
function out = run_variable_exchange(AFullFunc,param,maskedMaps,ALambda)
    MaxVEIter = 500;
    cost = inf;
    xCurr = zeros(size(AFullFunc,2),1);
    Tol = 0.001;    
    %targetFAInRad = maskedMaps.TargetMasked*deg2rad(param.targetFlipAngle);
    targetFAInRad = maskedMaps.target_masked*deg2rad(param.targetFlipAngle);
    
    phiTarget = angle(sum(maskedMaps.b1SensMaskedHz,2));
    z = exp(1i*phiTarget);
    mask_1_slice = maskedMaps.mask_one_slice;
    fullImage =  zeros(numel(mask_1_slice),1);
    CGtikhonov = param.CGtikhonov;
    
for iDx = 1:MaxVEIter
    targetFA = z .* targetFAInRad;

    xCurr = ALambda*targetFA;
    currFA = AFullFunc*xCurr;
    costNew = norm(currFA - targetFA) + CGtikhonov*(xCurr'*xCurr);
    
    costDiff = abs(cost-costNew)/cost;
    if costDiff < Tol
        break
    end
    cost = costNew;
    
    % Smooth phase
    fullImage(logical(mask_1_slice(:))) = angle(currFA);
    fullImage = reshape(fullImage,size(mask_1_slice));
    fullImage = imgaussfilt(fullImage,2,'FilterSize',5);
        
    phiTarget = fullImage(mask_1_slice(:));
    z = exp(1i*phiTarget);    
end
    bOut = xCurr;
    finalPwr = norm(xCurr);
    %finalCost = costNew;
    %finalCost = cost;
        finalMag = zeros(numel(mask_1_slice),1);
        finalMag(logical(mask_1_slice(:))) = AFullFunc*xCurr;
        finalMag = reshape(finalMag,size(mask_1_slice));

    %finalRMSE = norm(currFA - targetFA);
    finalRMSE = norm(abs(currFA) - abs(targetFA))/norm(abs(targetFA));
    out.bOut = bOut;
    out.finalPwr = finalPwr;
    out.finalMag = finalMag;
    out.finalRMSE = finalRMSE;
end
function out = run_SingVal(bVE,maskedMaps,param,AFull)
    bMatrix = reshape(bVE,[68 88]);
    [U,S,V] = svd(bMatrix,'econ');
    %U1 = U(:,1:30); S1 = S(1:30,1:30); V1 = V(:,1:30);
    %diagS1 = diag(S1);
    %xInitial = rand(30,1)*2*sum(diagS1)/30;
    xInitial = diag(S);
    xInitial = xInitial + 5*(rand(size(xInitial))-0.5);
    targetFA = maskedMaps.target_masked*deg2rad(param.targetFlipAngle);
    FunHandle = @(x) norm(abs(AFull*reshape(U*diag(x)*V',[5984 1]))-abs(targetFA))...
        /norm(abs(targetFA));
    A = ones(1,numel(xInitial));
    b = sum(xInitial);
    lb = zeros(1,numel(xInitial));
    optionsFMin = optimoptions(@fmincon);
    optionsFMin.Algorithm = 'sqp';
    optionsFMin.Display = 'iter';
    optionsFMin.MaxFunctionEvaluations = 1e5;
    optionsFMin.SpecifyConstraintGradient = false;
    optionsFMin.OptimalityTolerance = 1e-10;
    optionsFMin.StepTolerance = 1e-10;
    optionsFMin.DiffMinChange = 5;
    optionsFMin.FunctionTolerance = 1e-10;
    optionsFMin.FiniteDifferenceType = 'central';
    optionsFMin.MaxIterations = 500;
    optionsFMin.SubproblemAlgorithm = 'cg';
    [SingVal,out.fval,out.exitflag,out.output] = fmincon(FunHandle,xInitial,A,b,[],[],...
        lb.',[],[],optionsFMin);
    out.SingVal = diag(SingVal);
end
%%
function out = run_active_set(bVE,maskedMaps,param,AFull)
    %Setting upper and lower bounds
    maxV = 239; % Check where I got this.
    clear ub lb
    ub(1:numel(bVE)) = maxV;
    ub((1+numel(bVE)):(numel(bVE)*2)) = inf;

    lb(1:numel(bVE)) = 0;
    lb((1+numel(bVE)):(numel(bVE)*2)) = -inf;

    %Setting optimization options

    optionsFMin = optimoptions(@fmincon);
    optionsFMin.Algorithm = 'interior-point';
    optionsFMin.Display = 'iter';
    optionsFMin.MaxFunctionEvaluations = 1e5;
    optionsFMin.SpecifyConstraintGradient = false;
    optionsFMin.OptimalityTolerance = 1e-6;
    optionsFMin.StepTolerance = 1e-6;
    optionsFMin.FiniteDifferenceType = 'central';
    optionsFMin.MaxIterations = 20;
    optionsFMin.SubproblemAlgorithm = 'cg';
%     fmincon_options = optimoptions(	'fmincon', ...
%                         'MaxFunEvals',10^5,...
%                         'TolFun',1e-6,...
%                         'TolCon',1e-6,...
%                         'TolX',1e-6,...
%                         'MaxIter',20,...
%                         'Algorithm','interior-point',...
%                         'GradObj','on',...
%                         'GradConstr','on',...
%                         'SubproblemAlgorithm','cg',...
%                         'Display','iter'...
%                         );
    nonlincon = [];
    xInitial = [abs(bVE);angle(bVE)];
    targetFA = maskedMaps.target_masked*deg2rad(param.targetFlipAngle);
    FunHandle = @(x) norm(abs(AFull...
        *(x(1:numel(bVE)).*...
        exp(1i*x((numel(bVE)+1):numel(bVE)*2))))...
        - abs(targetFA))/norm(abs(targetFA));
    [bOutTemp,out.fval,out.exitflag,out.output] = fmincon(FunHandle,xInitial,[],[],[],[],...
        lb.',ub.',nonlincon,optionsFMin);
    
    out.bOut = bOutTemp(1:numel(bVE)).*exp(1i*bOutTemp((numel(bVE)+1):numel(bVE)*2));
    mask_1_slice = maskedMaps.mask_one_slice;
    finalMag = zeros(numel(mask_1_slice),1);
    finalMag(logical(mask_1_slice(:))) = AFull*out.bOut;
    finalMag = reshape(finalMag,size(mask_1_slice));
    out.finalMag = finalMag;
end
function [RFOut,GradOut] = AdjustDelay(RFIn, GradIn, Diff)
DiffIn10us = Diff/10;
% Positive Diff = RF starts earlier;
if DiffIn10us>0
    RFOut = [RFIn;zeros(DiffIn10us,8)];
    GradOut = [zeros(DiffIn10us,3);GradIn];
    
elseif DiffIn10us<0
    RFOut = [zeros(-DiffIn10us,8);RFIn];
    GradOut = [GradIn;zeros(-DiffIn10us,3)];
end
end
   %%
   function ImgNew = UnderSampleImage(Img,UndesampleFactor)
        % check image dimension and undersamplefactor size
        if ndims(Img) ~= numel(UndesampleFactor)
            disp('Image dimension and undersampling factor size do not match...')
            return
        end
        if ndims(Img) == 3
            ImgNew = Img(1:UndesampleFactor(1):end,1:UndesampleFactor(2):end,1:UndesampleFactor(3):end);
        elseif ndims(Img) == 4
            ImgNew = Img(1:UndesampleFactor(1):end,1:UndesampleFactor(2):end,...
                1:UndesampleFactor(3):end,1:UndesampleFactor(4):end);
        end
   end
    function PlotSingleSlice(ImgMtx,VerRes,HoriRes,Angle,bool_flip)
        if ~(Angle == 0 || Angle == 90 || Angle == 180)
            disp('Angle not recognised...')
            return
        end
        figure
        set(gcf,'color','w','InvertHardcopy','off')
        if bool_flip
            imagesc(imrotate(flip(ImgMtx,1),-Angle));
        else
            imagesc(imrotate(ImgMtx,-Angle));
        end
        if Angle == 0
            pbaspect([HoriRes*size(ImgMtx,2) VerRes*size(ImgMtx,1) 1]);
        elseif Angle == 90
            pbaspect([VerRes*size(ImgMtx,1) HoriRes*size(ImgMtx,2) 1]);
        elseif Angle == 180
            pbaspect([HoriRes*size(ImgMtx,2) VerRes*size(ImgMtx,1) 1]);
        end
        axis off; colormap gray;
    end
    %%    
    function h = tyPlotTrajectory( T, Gradient, OpenFigure )
    if nargin < 3
        OpenFigure = 1;
    end
    
    if OpenFigure
       figure
       hold on
       box on
       grid on 
    end

    if nargin < 2
       Gradient = []; 
    end

    %This method plots the trajectory, colorcoded by the velocity given by 
    % v = sqrt( Gx^2 + Gy^2 )
 
    %if vectors are very long, make them sparser
    AimLength   = 10000;
    Step        = numel(T.Sam.kx) / AimLength;
    kx_plot     = T.Sam.kx(1 : ceil(Step) : end);
    ky_plot     = T.Sam.ky(1 : ceil(Step) : end);
    kz_plot     = T.Sam.kz(1 : ceil(Step) : end);
    Gx_plot     = T.Sam.Gx(1 : ceil(Step) : end);
    Gy_plot     = T.Sam.Gy(1 : ceil(Step) : end);
    Gz_plot     = T.Sam.Gz(1 : ceil(Step) : end);
    if Step ~= round(Step)
        kx_plot = [ kx_plot, T.Sam.kx(end) ];
        ky_plot = [ ky_plot, T.Sam.ky(end) ];
        kz_plot = [ kz_plot, T.Sam.kz(end) ];
        Gx_plot = [ Gx_plot, T.Sam.Gx(end) ];
        Gy_plot = [ Gy_plot, T.Sam.Gy(end) ];
        Gz_plot = [ Gz_plot, T.Sam.Gz(end) ];
    end
    
    %create vector of velocitys
    V   = sqrt( (Gx_plot(:)).^2 + (Gy_plot(:)).^2 + (Gz_plot(:)).^2 );
    
    %plot colorcoded trajectory
    h = mesh([kx_plot(:), kx_plot(:)], [ky_plot(:), ky_plot(:)], [kz_plot(:), kz_plot(:)], [V V], ...
        'EdgeColor', 'interp', 'FaceColor', 'none', 'Linewidth', 3);
%     h = mesh([kx_plot(end-1000:end), kx_plot(end-1000:end)], [ky_plot(end-1000:end), ky_plot(end-1000:end)],...
%         [kz_plot(end-1000:end), kz_plot(end-1000:end)], [V(end-1000:end) V(end-1000:end)], ...
%         'EdgeColor', 'interp', 'FaceColor', 'none', 'Linewidth', 3);
    %set(gca,'Clim',[0, max(V(:))])
    set(gca,'Clim',[0, 20])
    
    kmax        = max(abs([kx_plot(:); ky_plot(:); kz_plot(:)]));
    RelLength   = 0.2;
    %plot gradient derivative
    if ~isempty( Gradient )
        if size( Gradient, 2 ) == 2
            Gradient = [Gradient, zeros( size(Gradient,1), 1 )];
        end
        %normalize
        Gradient = Gradient ./ max(abs(Gradient(:)));
        for k = 1 : size(Gradient,1)
           line([T.kx(k),T.kx(k) + kmax*RelLength*Gradient(k,1)], ...
                [T.ky(k),T.ky(k) + kmax*RelLength*Gradient(k,2)], ...
                [T.kz(k),T.kz(k) + kmax*RelLength*Gradient(k,3)], ...
               'Color','blue','Linewidth',1.5) 
        end
        
    end
    
    %plot control points
    hold on
    %hcp = PlotControlPoints( T );
    
    if OpenFigure
        %LegendStrings = {'k-space Trajectory','Unconstrained','Vector Constraint','Fixed Gradients'};
        %warning( 'off', 'MATLAB:legend:UnsupportedEdgeColor' )
        %legend([h,hcp(~isnan(hcp))],LegendStrings( find([1,~isnan(hcp)]) ))
        legend('k-space Trajectory')
        %warning( 'on', 'MATLAB:legend:UnsupportedEdgeColor')
    end
    
    %if z only contains zeros, assume that trajectory is only 2D, set view
    %accordingly
    
  	colormap(hot)
    t = colorbar;
    nudge(t,[0.04 0 0 0]);
    set(get(t,'ylabel'),'String', 'k-space Velocity $\|\textbf{G}\|_{2}$',...
        'Interpreter','latex','Fontsize',18,'Fontweight','bold');
    %set(gca,'Color',[0.8 0.8 0.8])
    xlabel('k_x (m^{-1})','Position',[-10 -50 -25]);ylabel('k_y (m^{-1})','Position',[-50 0 -25]);
    zlabel('k_z (m^{-1})','Rotation',0,'Position',[-45 40 28]);
    if all( T.Sam.kz == 0 )
        set(gca,'View',[0 90])
	else
        set(gca,'View',[-46 22])
    end
    zticks([-20 0 20]);
    %zlim([-30 30])
    %set figure size
    if OpenFigure
        set(gcf,'OuterPosition',[235   222   763   742])
        set(gcf,'Position',     [243   230   747   650])
        
        set(gca,'OuterPosition',[-0.1240   -0.0836    1.2193    1.1373])
        set(gca,'Position',     [ 0.0345+0.05    0.0415+0.015    0.8370-0.05    0.9269])
        
        
        
        set(gcf,'color','w','InvertHardcopy','off')
        set(gca,'FontSize',18)
        %movegui(gcf,'center')

        %axis equal 
    end
    end
    %%
 function writeIniFile_ty(rfIn,gradIn)
% From WTC writeIniFile.m
% Write an arbitrary pulse to the ini format that is expected by the
% SBBExcitationPTx class
% The pulse will be written to the local folder X
% If the scanner drives are mounted the pulse will be written to them.
% If the VM RFPulses shared folder is mounted it will also be written
% there.
% A matlab .mat file of the final pulse will also be written to Y

% Input:- rfIn: NxChannels array of complex transmit voltages
%       - gradIn: Nx3 array of gradient amplitudes. In Hz/cm
%       - sampleTime: Duration of each X step. In us.
% Comment
% Output:


%This code will apply the shift of gradient dimensions, make sure that this
%is not done before this function!
%% Checks
if size(gradIn,2) ~=3
    error('gradIn should be Nx3 array.')
end


if size(gradIn,1)==size(rfIn,1)
   %samples = size(gradIn,1);
   oversampling = 1;
elseif mod(size(rfIn,1),size(gradIn,1)) == 0
    % Integer scaling of RF samples
    %samples = size(gradIn,1);
    oversampling = size(rfIn,1)/size(gradIn,1);
    
    if oversampling>10
        error('Permitted range for oversampling of RF is 1-10.');
    end
else 
    error('The number of gradient steps must be an integer multiple of the rf samples.') % From the docs: Oversampling (optional): The RF oversamping factor. In other words, the number of RF samples per gradient sample. This parameter is optional and read as one if not present. Its valid range is from 1 to 10 (integer).
end
%% Do the gradient switch
% This is from the dicom coordinate system to an identity rotation matrix
% (HFS)
gradSwitched = gradIn(:,[2 1 3]);
gradSwitched(:,2) = -1*gradSwitched(:,2);
%% Units
GAMMA = 42.57E6; % Hz/T 

% Gradients from Hz/cm to mT/m
gmag_mTm = 1000*(gradSwitched*100)/GAMMA;
writeGrad_mTm = gmag_mTm;

% Normalise the RF
peakV = max(abs(rfIn(:)));
writeRfAmp = abs(rfIn)/peakV; %Scale the amplitude from 0 to 1;
writeRfPhase = angle(rfIn)+pi; % Scale phase from 0 to 2pi

numsamples = size(gmag_mTm,1);  %YT change
GradDelay = 0;                  %YT change
%% Write the file

fID = fopen('/Users/ytong/Documents/XP/pTxPulses/pTXArbitrary.ini','w');
fprintf(fID,'# Created in run_STA_KTPdesign.m \n');
fprintf(fID,'[pTXPulse]\n');
fprintf(fID,'\n');
fprintf(fID,'NUsedChannels        = %i\n',8);		
fprintf(fID,'MaxAbsRF             = %0.3f	# peak RF voltage [V]\n',peakV);
fprintf(fID,'Samples              = %i		# number of samples\n',numsamples);
fprintf(fID,'PulseName            = pTxArbitrary\n');
fprintf(fID,'Comment              = pTx pulse\n');
fprintf(fID,'Oversampling         = %i		# no oversampling\n',oversampling);
fprintf(fID,'NominalFlipAngle     = %i		# flip-angle of this pulse if played out with 502.385V peak voltage\n',10);
fprintf(fID,'SampleTime			  = 10\n');
fprintf(fID,'GradDelay			  = %i\n',GradDelay);
% Gradient
fprintf(fID,'\n');
fprintf(fID,'[Gradient]\n');
fprintf(fID,'\n');
for iDx = 1:size(writeGrad_mTm,1)
    fprintf(fID,'G[%i]= %f %f %f\n',iDx-1,writeGrad_mTm(iDx,1),writeGrad_mTm(iDx,2),writeGrad_mTm(iDx,3));
end
fprintf(fID,'\n');

% RF
for jDx = 1:size(writeRfAmp,2)
    fprintf(fID,'[pTXPulse_ch%i]\n',jDx-1);
    fprintf(fID,'\n');
    for iDx = 1:size(writeRfAmp,1)
        fprintf(fID,'RF[%i]=  %f\t%f\n',iDx-1,writeRfAmp(iDx,jDx),writeRfPhase(iDx,jDx));
    end
    fprintf(fID,'\n');
end

fclose(fID);

%% Write out matlab copy of parameters as well
save('/Users/ytong/Documents/XP/pTxPulses/pTXArbitrary.mat','peakV','writeRfAmp','writeRfPhase','writeGrad_mTm')

%% Copy to scanner
if isfolder('/Volumes/Disk_C')
    copyfile('/Users/ytong/Documents/XP/pTxPulses/pTXArbitrary.ini','/Volumes/Disk_C/MedCom/MriCustomer/seq/RFPulses/pTXArbitrary.ini')
end
if isfolder('/Volumes/Disk_C-1')
    copyfile('/Users/ytong/Documents/XP/pTxPulses/pTXArbitrary.ini','/Volumes/Disk_C-1/MedCom/MriCustomer/seq/RFPulses/pTXArbitrary.ini')
end
% %% Copy to VM
% YT not entirely sure what this is
% if isfolder('/Volumes/RFPulses')
%    copyfile('/Users/wclarke/Documents/VMShared/pTxArbitrary.ini','/Volumes/RFPulses/pTxArbitrary.ini')
% end 
end
%%
function RFShimWriteMLEV(RFToWrite)
%Write the file

    fID = fopen('/Users/ytong/Documents/MATLAB/tong-acptx/tySpiral/tyTRUST/MLEVShim.ini','w');
    fprintf(fID,'# Created in .m \n');
    fprintf(fID,'[pTXPulse]\n');
    fprintf(fID,'\n');
    fprintf(fID,'NUsedChannels        = %i\n',8);		
    fprintf(fID,'NShims               = %i\n',1);
    fprintf(fID,'PulseName            = MLEVShim\n');
    fprintf(fID,'Comment              = Shim for MLEV pulse\n');


    % RF
    for jDx = 1:8
        fprintf(fID,'[pTXPulse_ch%i]\n',jDx-1);
        fprintf(fID,'\n');
        fprintf(fID,'RF[%i]=  %f\t%f\n',0,RFToWrite(jDx,1),RFToWrite(jDx,2));
        fprintf(fID,'\n');
    end

    fclose(fID);
    if isfolder('/Volumes/Disk_C')
        copyfile('/Users/ytong/Documents/MATLAB/tong-acptx/tySpiral/tyTRUST/MLEVShim.ini',...
            '/Volumes/Disk_C/MedCom/MriCustomer/seq/RFPulses/pTXArbitrary.ini')
    end
    if isfolder('/Volumes/Disk_C-1')
        copyfile('/Users/ytong/Documents/MATLAB/tong-acptx/tySpiral/tyTRUST/MLEVShim.ini',...
            '/Volumes/Disk_C-1/MedCom/MriCustomer/seq/RFPulses/pTXArbitrary.ini')
    end
end
