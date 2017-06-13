%addpath(genpath('/Users/ytong/Documents/MATLAB/WTCpTx'));

%singleTxPath = '/Users/ytong/Documents/MATLAB/20170516_F7T_2017_PH_054';
singleTxPath = '/Volumes/Data/DICOM/2017-06/20170613_F7T_2017_PH_054';
    singleTxObj = DicomFM.WTCSingleTxFieldmaps(singleTxPath);
    singleTxObj.interpolateTo('B1');
    singleTxObj.createMask(@DicomFM.maskFunctions.multithreshMask,false);
    slice = round(size(singleTxObj.getMask,3))/2;
    singleTxObj.setSlice(slice);
    mask = singleTxObj.getMask;
    b1 = singleTxObj.getB1PerV('Hz','delete');
    b0 = singleTxObj.getB0('Hz','delete');
    Positions = singleTxObj.getPositions('cm','delete'); %3 by Ns matrix
    Positions = Positions';%SysMtxSingle takes Ns by 3 matrix.
    Ns = size(b1,1);
   
%%
%Get k-space trajectory
deltax = 15/50; % spatial resolution of trajectory
FOVsubsamp = 10; % cm, xfov of spiral
densamp = 75; % duration of full density sampling (# of samples)
dentrans = 75; % duration of transition from higher to lower
nl = 1; % degree of undersampling outer part (any real #)

% gmax = 4;  %g/cm
gmax = 0.0004;  %T/cm ( = 0.4 mT/cm = 40 mT/m)
% dgdtmax = 18000;      %g/cm/s; max slew rate.
dgdtmax = 1.8;      %T/cm/s ( = 180 T/m/s) max slew rate. 

dt = 10e-6;     %sec; sampling period in pulse sequence - 4 us
kdimxy = floor(FOVsubsamp/deltax); %#k-space lines

[g,k,t,s,ddens,NN] = spiralgradlx6(FOVsubsamp,kdimxy,dt,dgdtmax/100,gmax,nl,densamp,dentrans);
g = flip(g,2);
g = g';
grad = [real(g(:)) imag(g(:))];
grad = grad(NN(2)+1:end,:);
grad = grad * 42.57E6; % Hz/T
gr = [grad zeros(size(grad,1),1)];

rfOn = true(size(grad,1),1);
tp = dt * ones(size(grad,1),1);
sens = singleTxObj.getB1PerV('Hz','delete');
df = singleTxObj.getB0('Hz','delete');
dp = singleTxObj.getPositions('cm','delete').';
Nt = size(gr,1);
Nc = 1;
roughbeta = 10^2.25;
DA = genAMatFull(tp,rfOn,gr,sens,df,dp);

lenghtRampUp = numel(g) - size(gr,1);
RampUp = [real(g(1:lenghtRampUp)) imag(g(1:lenghtRampUp))];
gOut = [RampUp;grad];%Adding gradient ramp.
gOut = gOut/425.8; %from Hz/cm to mT/m
%% Load excitation pattern and center it around the mask
  load /Users/ytong/Documents/MATLAB/Utility/Pattern/OffRect.mat
  dim = size(d);
  
  [xi,yi] = ndgrid(-dim(1)/2:dim(1)/2-1,-dim(2)/2:dim(2)/2-1);
  xmi = round(mean(xi(logical(mask(:)))));
  ymi = round(mean(yi(logical(mask(:)))));
  %d = circshift(d,[xmi ymi]); %Center the excitation pattern
  d = ift2((hamming(dim(1))*hamming(dim(2))').^2.*ft2(d));%Smooth it a bit.
  
% %% Build RF waveform roughness penalty matrix 
    R = sqrt(roughbeta)*spdiags([-ones(Nt,1) ones(Nt,1)],[0 1],Nt,Nt) + ...
        sqrt(roughbeta)*speye(Nt);
    Rfull = kron(speye(Nc),R);
% 
  %% Create pulse
tikhonov = 2e-5;
[RFOut,finalCost,finalPwr,finalMag] = CG_SysMat(DA,d(mask(:)),tikhonov,mask);
fprintf('Max RF = %f\n',max(abs(RFOut)))
  % Plotting
  magnetization = make_blochSim(RFOut,b1,b0,grad,dt,Positions,mask);
  m = zeros(dim);m(mask) = DA*RFOut(:);
  err = norm(m-d)/norm(d);
  roughness = norm(Rfull*RFOut(:));
  pow = sqrt(mean(abs(RFOut(:)).^2));
  fprintf('All-channels NRMSE: %0.2f%%. All-channels roughness: %f.\n',err*100,roughness);
  %%plotting
  maxamp = max(abs([d(:);m(:)]));
  figure(2)
  subplot(221)
  imagesc(abs(d),[0 maxamp]);axis image;colorbar
  title 'Desired pattern'
  subplot(222)
  imagesc(abs(m),[0 maxamp]);axis image;colorbar
  title(sprintf('Final pattern',Nc));
  subplot(223)
  imagesc(abs(magnetization.mxy),[0 maxamp]);axis image;colorbar
  title(sprintf('Bloch Simulation'));
  subplot(224)
  imagesc(abs(m)-abs(d),[0 maxamp]);axis image;colorbar
  title(sprintf('Magnitude Error\nComplex NRMSE = %0.2f%%',err*100));
  
%% Writing into text file
  MakeDotH(RFOut,gOut,dt);