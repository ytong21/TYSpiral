%addpath(genpath('/Users/ytong/Documents/MATLAB/WTCpTx'));

%singleTxPath = '/Users/ytong/Documents/MATLAB/20170516_F7T_2017_PH_054';
    singleTxPath = '/Volumes/Data/DICOM/2017-07/20170703_F7T_2017_PH_065';
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
    
    
   
%%
%Get k-space trajectory
deltax = 45/50; % spatial resolution of trajectory
% Modified. Originally 15/50.
FOVsubsamp = 45; % cm, xfov of spiral
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
gr = [grad zeros(size(grad,1),1)];%3 gradients instead of 2

rfOn = true(size(grad,1),1);
tp = dt * ones(size(grad,1),1);
sens = singleTxObj.getB1PerV('Hz','delete');
df = singleTxObj.getB0('Hz','delete');
dp = singleTxObj.getPositions('cm','delete').';
Nt = size(gr,1);
Nc = 1;
DA = genAMatFull(tp,rfOn,gr,sens,df,dp);

lenghtRampUp = numel(g) - size(gr,1);
RampUp = 42.57E6*[real(g(1:lenghtRampUp)) imag(g(1:lenghtRampUp))];% Hz/T
gOut = [RampUp;grad];%Adding gradient ramp.
gOut = gOut/425.8; %from Hz/cm to mT/m
gHzpercm = gOut*425.8;

kTraj = zeros(size(gHzpercm));
kTraj(1,:) = gHzpercm(1,:);
for ii = 2:size(gHzpercm,1)
    kTraj(ii,1) = kTraj(ii-1,1) + gHzpercm(ii,1);
    kTraj(ii,2) = kTraj(ii-1,2) + gHzpercm(ii,2);
end
kTraj = kTraj*dt; % to convert to 1/cm;  
figure(11)
plot(kTraj(:,1),kTraj(:,2))
xlabel('kx 1/cm')
ylabel('ky 1/cm')
title('k-space trajectory')
%% Load excitation pattern and center it around the mask
%   load OffRect.mat
%   dim = size(mask);
%   
%   [xi,yi] = ndgrid(-dim(1)/2:dim(1)/2-1,-dim(2)/2:dim(2)/2-1);
%   xmi = round(mean(xi(logical(mask(:)))));
%   ymi = round(mean(yi(logical(mask(:)))));
%   %d = circshift(d,[xmi ymi]); %Center the excitation pattern
  dim = size(mask);
  d = zeros(dim);
  d(round(0.35*dim(1)):round(0.5*dim(1)),round(0.25*dim(2)):round(0.65*dim(2))) = 1;
  d = ift2((hamming(dim(1))*hamming(dim(2))').^2.*ft2(d));%Smooth it a bit.
  dFreqDomain = fftshift(ifft2(d));figure(12);imagesc(abs(dFreqDomain));
  ImageRes = 0.4; % cm
  FDPixSize = 1./(dim*ImageRes); %cm
  FDRadius = 10; %number of points. 
  TotalEnergy = sum(abs(dFreqDomain(:)));
  PartialKSpace = dFreqDomain((dim(1)/2-(FDRadius-1)):(dim(1)/2+FDRadius), ...
      (dim(2)/2-(FDRadius-1)):(dim(2)/2+FDRadius));
  PartialEnergy = sum(abs(PartialKSpace(:)));
  PowerPercentage = PartialEnergy/TotalEnergy;
  fprintf('The power percentage with a k-space with size of %.4f cm(-1) by %.4f cm(-1) is %.1f%%.\n', ...
      FDRadius*(FDPixSize(1)), FDRadius*(FDPixSize(2)),100*PowerPercentage);
 %% Build RF waveform roughness penalty matrix 
    roughbeta = 10^2.25;
    R = sqrt(roughbeta)*spdiags([-ones(Nt,1) ones(Nt,1)],[0 1],Nt,Nt) + ...
        sqrt(roughbeta)*speye(Nt);
    Rfull = kron(speye(Nc),R); 
  % Create pulse
tikhonov = 1e-5;
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
  imagesc(imrotate(abs(d),90),[0 maxamp]);axis image;colorbar
  title 'Desired pattern'
  subplot(222)
  imagesc(imrotate(abs(m),90),[0 maxamp]);axis image;colorbar
  title(sprintf('Final pattern'));
  subplot(223)
  imagesc(imrotate(abs(magnetization.mxy),90),[0 maxamp]);axis image;colorbar
  title(sprintf('Bloch Simulation'));
  subplot(224)
  imagesc(imrotate(abs(m)-abs(d),90),[0 maxamp]);axis image;colorbar
  title(sprintf('Magnitude Error\nComplex NRMSE = %0.2f%%',MagErr*100));
  fprintf('The slice position is %.4f mm.\n',10*Positions(1,3));%Positions in cm.

%% Writing into text file
  MakeDotH(RFOut,gOut,dt);