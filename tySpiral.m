%addpath(genpath('/Users/ytong/Documents/MATLAB/WTCpTx'));

%singleTxPath = '/Users/ytong/Documents/MATLAB/20170516_F7T_2017_PH_054';
    singleTxPath = '/Volumes/Data/DICOM/2017-07/20170704_F7T_2017_PH_069';
    singleTxObj = DicomFM.WTCSingleTxFieldmaps(singleTxPath);
    singleTxObj.interpolateTo('B1');
    slice = round(size(singleTxObj.getMask,3)/2);
    singleTxObj.createMask(@DicomFM.maskFunctions.multithreshMask,true);
    singleTxObj.setSlice(slice);
    singleTxObj.plotSlice(slice,3);
    mask = singleTxObj.getMask;
    b1 = singleTxObj.getB1PerV('Hz','delete');
    b0 = singleTxObj.getB0('Hz','delete');
    Positions = singleTxObj.getPositions('cm','delete'); %3 by Ns matrix
    Positions = Positions';%SysMtxSingle takes Ns by 3 matrix.
    Ns = size(b1,1);
%% Create the excitation pattern


  dim = size(mask);
  d = zeros(dim);
  d(round(0.35*dim(1)):round(0.5*dim(1)),round(0.25*dim(2)):round(0.65*dim(2))) = 1;
%     [xi,yi] = ndgrid(-dim(1)/2:dim(1)/2-1,-dim(2)/2:dim(2)/2-1);
%     xmi = round(mean(xi(logical(mask(:)))));
%     ymi = round(mean(yi(logical(mask(:)))));
%     d = circshift(d,[xmi ymi]); %Center the excitation pattern
  d = ift2((hamming(dim(1))*hamming(dim(2))').^2.*ft2(d));%Smooth it a bit.
  d = sind(90)*d;
  dFreqDomain = fftshift(ifft2(d));
  
  ImageRes = 0.4; % cm
  FDPixSize = 1./(dim*ImageRes); %cm
  FDRadius = 2; %number of points. 
  figure(11); imagesc([-FDPixSize(1)*dim(1)/2 FDPixSize(1)*dim(1)/2],...
      [-FDPixSize(2)*dim(2)/2 FDPixSize(2)*dim(2)/2],abs(dFreqDomain)); 
  title('2D FFT of the desired pattern'); 
  xlabel('kx 1/cm'); ylabel('ky 1/cm');
  TotalPwer = sum(abs(dFreqDomain(:)));
  PowerPercentage = 0;
  PowerThreshold = 0.9;
  while PowerPercentage < PowerThreshold && FDRadius <= floor(min(dim)/2)
      FDRadius = FDRadius + 1;
      PartialKSpace = dFreqDomain((dim(1)/2-(FDRadius-1)):(dim(1)/2+FDRadius), ...
          (dim(2)/2-(FDRadius-1)):(dim(2)/2+FDRadius));
      PartialPower = sum(abs(PartialKSpace(:)));
      PowerPercentage = PartialPower/TotalPwer;
  end
  fprintf('The power percentage with a %.4f cm(-1) by %.4f cm(-1) k-space is %.1f%%.\n', ...
      FDRadius*(FDPixSize(1)), FDRadius*(FDPixSize(2)),100*PowerPercentage);
  desiredRadius = max(FDPixSize)*FDRadius;
%% Get k-space trajectory
trajDensity = 'variable';
switch trajDensity
    case 'uniform'
        deltax = 15/50; % spatial resolution of trajectory % Modified. Originally 15/50.
        FOVsubsamp = 30; % cm, xfov of spiral
        densamp = 75; % duration of full density sampling (# of samples)
        dentrans = 75; % duration of transition from higher to lower
        nl = 1; % degree of undersampling outer part (any real #)
        gmax = 0.0004;  %T/cm ( = 0.4 mT/cm = 40 mT/m) % gmax = 4;  %g/cm
        % dgdtmax = 18000;      %g/cm/s; max slew rate.
        dgdtmax = 1.8;      %T/cm/s ( = 180 T/m/s) max slew rate. 
        dt = 10e-6;     %sec; sampling period in pulse sequence - 4 us
        MaxRadius = 1.87;
        while MaxRadius > desiredRadius
            deltax = deltax + 1/50;
            kdimxy = floor(FOVsubsamp/deltax); %#k-space lines
            [g,k,t,s,ddens,NN] = spiralgradlx6(FOVsubsamp,kdimxy,dt,dgdtmax/100,gmax,nl,densamp,dentrans);
            g = flip(g,2);
            g = g';
            gHzpercm = 42.577E6*[real(g(:)) imag(g(:))];% Hz/T
            [kTraj,MaxRadius] = tyTraj(gHzpercm,dt);
        end
        grad = gHzpercm(NN(2)+1:end,:);  
        gr = [grad zeros(size(grad,1),1)];%3 gradients instead of 2 for Bloch simulation
        gOut = gHzpercm/425.77; %from Hz/cm to mT/m
    case 'variable'
        vdMtxSize = 20;
        vdFOV = 32; 
        alpha = 3; %oversampling coefficient
        stVds = genVDSpirals(vdMtxSize, alpha, vdFOV, 1,'flagGrad',1);
        gHzpercm = 42.577E6*(1E-4)*[real(stVds.Grad') imag(stVds.Grad')]; %Hz/cm
        vdRadius = hypot(sum(gHzpercm(:,1)),sum(gHzpercm(:,2)))*dt;
        while vdRadius < desiredRadius
            vdFOV = vdFOV - 1;
            stVds = genVDSpirals(vdMtxSize, alpha, vdFOV, 1,'flagGrad',1);
            gHzpercm = 42.577E6*(1E-4)*[real(stVds.Grad') imag(stVds.Grad')];
            vdRadius = hypot(sum(gHzpercm(:,1)),sum(gHzpercm(:,2)))*dt;
        end
        %So for the gradient is spiral-out. Adding 2 parts for ramp-down: deceleration & returning to (0,0)
        decelPoints = 20;
        decel = [linspace(gHzpercm(end,1),0,decelPoints); linspace(gHzpercm(end,2),0,decelPoints)]';
        kFinal = (sum(gHzpercm) + sum(decel))*dt; %1/cm
        rampPoints = 75;
        % kFinal = 0.5*GrampMaxAmp*kFinalrampPoints*dt;
        GrampMaxAmp = (-kFinal)./(0.5*rampPoints*dt);
        rampDown = [linspace(0,GrampMaxAmp(1),(rampPoints+1)/2); linspace(0,GrampMaxAmp(2),(rampPoints+1)/2)]';
        rampDown = [rampDown; flip(rampDown(1:end-1,:),1)];
        grad = [gHzpercm; decel; rampDown]; %/Hz/cm
        grad = flip(grad,1);
        [kTraj,~] = tyTraj(grad,dt);
        gr = [grad zeros(size(grad,1),1)];%3 gradients instead of 2 for Bloch simulation
        gOut = grad/425.77; %from Hz/cm to mT/m
end
%%
figure(12)
plot(kTraj(:,1),kTraj(:,2))
xlabel('kx 1/cm')
ylabel('ky 1/cm')
title('k-space trajectory')

%% Build system matrix
rfOn = true(size(grad,1),1);
tp = dt * ones(size(grad,1),1);
sens = singleTxObj.getB1PerV('Hz','delete');
df = singleTxObj.getB0('Hz','delete');
dp = singleTxObj.getPositions('cm','delete').';
Nt = size(gr,1);
Nc = 1;
DA = genAMatFull(tp,rfOn,gr,sens,df,dp);
 
% Build RF waveform roughness penalty matrix 
    roughbeta = 10^2.25;
    R = sqrt(roughbeta)*spdiags([-ones(Nt,1) ones(Nt,1)],[0 1],Nt,Nt) + ...
        sqrt(roughbeta)*speye(Nt);
    Rfull = kron(speye(Nc),R); 
  % Create pulse
tikhonov = 1e-7;
[RFOut,finalCost,finalPwr,finalMag] = CG_SysMat(DA,d(mask(:)),tikhonov,mask);
fprintf('Max RF = %f\n',max(abs(RFOut)))
  % Plotting
%   RFOut = [0 ;RFOut];
%   RFOut(end) = [];
  magnetization = make_blochSim(RFOut,b1,b0,gOut*425.8,dt,Positions,mask);
  m = zeros(dim);m(mask) = DA*RFOut(:);
  err = norm(m-d)/norm(d);
  MagErr = norm((abs(m)-abs(d)))/norm(abs(d));
  roughness = norm(Rfull*RFOut(:));
  pow = sqrt(mean(abs(RFOut(:)).^2));
  fprintf('All-channels NRMSE: %0.2f%%. All-channels roughness: %f.\n',MagErr*100,roughness);
  figure(13)
  plot(0.01:0.01:numel(RFOut)*0.01,abs(RFOut))
  title('RF pulse shape')
  xlabel('Time (ms)')
  ylabel('RF (volts)')
  %% plotting

  maxamp = max(abs([d(:);m(:)]));
  figure(14)
  subplot(221)
  imagesc(imrotate(flip(abs(d),1),-90),[0 maxamp]);axis image;colorbar
  title 'Desired pattern'
  subplot(222)
  imagesc(imrotate(flip(abs(m),1),-90),[0 maxamp]);axis image;colorbar
  title(sprintf('Final pattern'));
  subplot(223)
  imagesc(imrotate(flip(abs(magnetization.mxy),1),-90),[0 maxamp]);axis image;colorbar
  title(sprintf('Bloch Simulation'));
  subplot(224)
  imagesc(imrotate(flip(abs(m)-abs(d),1),-90),[0 maxamp]);axis image;colorbar
  title(sprintf('Magnitude Error\nNRMSE = %0.2f%%',MagErr*100));
  fprintf('The slice position is %.4f mm.\n',10*Positions(1,3));%Positions in cm.

%% Writing into text file
  MakeDotH(RFOut,gOut,dt);