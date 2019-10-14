    %% Load in the data and interporlate to B0
    pTxPath = '/Users/ytong/Documents/Data/7TDicom/20181105_F7T_2013_50_102';
    dt = Spectro.dicomTree('dir',pTxPath,'recursive',false);
    ptxFMObj = DicomFM.WTCpTxFieldmaps(dt,'B1String','dt_dream_wIce_60deg_90VRef__B1',...
        'B0String','fieldmap_ptx7t_iso4mm_trans_RL','LocaliserString','localizer_3D','InterpTarget',...
        'B0');
    %% Creating a mask using betted nii file
    mask_nii = niftiread('B0_mask.nii.gz');
    mask_nii = flip(mask_nii,2);
    
    %% Specify the slice number and make a mask
    Slice_Array = 35:44;
    ptxFMObj.createMask(@(x) DicomFM.maskFunctions.TRUST_mask(x,Slice_Array,3,mask_nii),true);
   
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
        TargetTemp(:,:,Slice_Array(iDx)) = squeeze(maskedMaps.mask(:,:,Slice_Array(iDx))) + CircleMaskFiltered;
    end
    %%
    %PlotSingleSlice(squeeze(TargetTemp(:,:,Slice_Array)),xz_res,y_res,90,true);
    maskedMaps.target_masked = TargetTemp(maskedMaps.mask)-1;


   %%

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
   %%
   options.phi_init = [40  40  40   7    6    1    1 ].';
   [ kx_vec, ky_vec, kz_vec, R_mat, NoRFIndices ] = TrajectoryControlPoints_Shells(options.phi_init);
   NV = sum( R_mat(:,end) );
   G0 = zeros( NV,1 );
   T_init = PerformOptimization3D( kx_vec, ky_vec, kz_vec, R_mat, options.Gmax, options.Smax, G0, options.GBF );
   %%
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
   %% Extending the RF pulse while keeping the number of points the same.
   Time_Vec_extended = Time_Vec*2;
   AFull_Extended = genAMatFull(1E-5*ones(numel(Time_Vec_extended),1),ones(numel(Time_Vec_extended),1),Gradient,maskedMaps.b1SensMaskedHz,...
            maskedMaps.b0MapMasked,maskedMaps.posVox);
   %% Building ALambda matrix. This is time-consuming. This is taken outside of the VE function
    param.targetFlipAngle = 180;
    param.numCh = 8;
    param.tol = 1e-5;
    %param.CGtikhonov = 1e-6;
    tikhonovArray = power(10,-6:1:-5);
    ALambda_cell = cell(size(tikhonovArray));
    %% 
    % run time: "\" < lsqminnorm < pinv
    for iDx = 1:numel(tikhonovArray)
        ALambda_cell{iDx} = ((AFull'*AFull + tikhonovArray(iDx)*eye(size(AFull,2)))\sparse(eye(5984)))...
            *AFull';
    end
    OutCell = cell(size(tikhonovArray));
    for iDx = 1:numel(tikhonovArray)
        param.CGtikhonov = tikhonovArray(iDx);
        OutCell{iDx} = run_variable_exchange(AFull,param,maskedMaps,ALambda_cell{iDx});
    end
    %% Plotting the target and final mag
    close all
    PlotSingleSlice(TargetTemp(:,:,Slice_Array),xz_res,xz_res,90,true);
    PlotSingleSlice(abs(finalMag),xz_res,xz_res,90,true);
    PlotSingleSlice(TargetTemp(:,:,Slice_Array)-1-abs(finalMag),xz_res,xz_res,90,true);
    %%
    [U1,S1,V1] = svd(AFull);
    %%
    SysMat_SVD.U = U1(:,1:600); SysMat_SVD.S = sparse(S1(1:600,1:600)); SysMat_SVD.V = V1(:,1:600);
    %%
    out = run_active_set(OutCell{3}.bOut,maskedMaps,param,SysMat_SVD);
    
    %%  
    out_singval = run_SingVal(OutCell{3}.bOut,maskedMaps,param,AFull);
    
    %%
    run_comparison_plot(CircleMaskFiltered,OutCell{2}.finalMag)
    %%
    RF_pulse = reshape(OutCell{2}.bOut,[748 8]);
    figure(223)
    for iDx = 1:8
        set(gcf,'color','w','InvertHardcopy','off')
        set(gcf,'units','centimeters','position',[4 4 60 30],'paperunits','centimeters','paperposition',[4 4 60 30])
        subplot(2,4,iDx)
        plot((1:748)/100,abs(RF_pulse(:,iDx)));ylim([0 300])
        title(sprintf('Channel %d',iDx))
        xlabel('Time (ms)')
        ylabel('Voltage (V)')
        set(gca,'FontSize',16)
    end
    %%
    function run_comparison_plot(target,result)
        figure(222)
        set(gcf,'color','w','InvertHardcopy','off')
        set(gcf,'units','centimeters','position',[4 4 50 20],'paperunits','centimeters','paperposition',[4 4 50 20])
        %slice_array = 35:2:44;
        result_to_plot_array = 1:2:10;
        for iDx = 1:5
        ax1 = subplot(2,5,iDx);
        imagesc(imrotate(rad2deg(abs(result(:,:,result_to_plot_array(iDx)))),-90),[0 180]);
        TT = title(sprintf('Slice %d',result_to_plot_array(iDx)));
        TT.FontSize = 16;
        axis off;colormap(ax1, 'hot');
        if iDx == 5
            CLB1 = colorbar('FontSize',18);nudge(CLB1,[0.04 0 0.01 0]);
        end
        ax2 = subplot(2,5,iDx+5);
        imagesc(imrotate(abs(rad2deg(result(:,:,result_to_plot_array(iDx))))-180*target,-90),[-20 20]);
        TT = title(sprintf('Slice %d diff',result_to_plot_array(iDx)));
        TT.FontSize = 16;
        axis off; colormap(ax2, 'parula');
        if iDx == 5
            CLB1 = colorbar('FontSize',18);nudge(CLB1,[0.04 0 0.01 0]);
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
function out = run_active_set(bVE,maskedMaps,param,SysMat_SVD)
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
    fmincon_options = optimoptions(	'fmincon', ...
                        'MaxFunEvals',10^5,...
                        'TolFun',1e-6,...
                        'TolCon',1e-6,...
                        'TolX',1e-6,...
                        'MaxIter',20,...
                        'Algorithm','interior-point',...
                        'GradObj','on',...
                        'GradConstr','on',...
                        'SubproblemAlgorithm','cg',...
                        'Display','iter'...
                        );
    nonlincon = [];
    xInitial = [abs(bVE);angle(bVE)];
    targetFA = maskedMaps.target_masked*deg2rad(param.targetFlipAngle);
    FunHandle = @(x) norm(abs(SysMat_SVD.U*(SysMat_SVD.S*(SysMat_SVD.V'...
        *(x( 1:numel(bVE)).*...
        exp(1i*x((numel(bVE)+1):numel(bVE)*2))))))...
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
    set(gca,'Clim',[0, max(V(:))])
    
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
    
  	colormap(parula)
    t = colorbar;
    
    set(get(t,'ylabel'),'String', 'k-space Velocity $\|\textbf{G}\|_{2}$','Interpreter','latex','Fontsize',14,'Fontweight','bold');
    %set(gca,'Color',[0.8 0.8 0.8])

    if all( T.Sam.kz == 0 )
        set(gca,'View',[0 90])
	else
        set(gca,'View',[-46 22])
    end
    
    %set figure size
    if OpenFigure
        set(gcf,'OuterPosition',[235   222   763   742])
        set(gcf,'Position',     [243   230   747   650])
        
        set(gca,'OuterPosition',[-0.1240   -0.0836    1.2193    1.1373])
        set(gca,'Position',     [ 0.0345    0.0415    0.8370    0.9269])
        
        set(gcf,'color','w','InvertHardcopy','off')
        set(gca,'FontSize',14)
        movegui(gcf,'center')

        axis equal 
    end
    end
    
