%addpath(genpath('/Users/ytong/Documents/MATLAB/WTCpTx'));
singleTxPath = '/Users/ytong/Documents/MATLAB/20170516_F7T_2017_PH_038';
addpath ../spiral/util
%singleTxPath = '/Volumes/Data/DICOM/2017-05/20170516_F7T_2017_PH_038';
    singleTxObj = DicomFM.WTCSingleTxFieldmaps(singleTxPath);
    singleTxObj.interpolateTo('B1');
    singleTxObj.createMask(30,false);
    slice = round(size(singleTxObj.getMask([]),3))/2;
    mask = singleTxObj.getMask(slice);
    b1 = singleTxObj.getB1('Gauss','delete',slice);
    b0 = singleTxObj.getB1('Hz','delete',slice);
    Positions = singleTxObj.getPositions('cm','delete',slice);
    Ns = size(b1,1);
    fov = 25.6; %fov is hard coded at this momemnt. In cm.
%%
%Get k-space trajectory
  Nc = 1;
  roughbeta = 10^2.25; % nyu, xfov = 8
  traj = 'spiral';
  deltax = fov/50; % spatial resolution of trajectory
  FOVsubsamp = 7; % cm, xfov of spiral
  forwardspiral = 0;
  dorevspiralramp = 0;
  gmax = 4; % g/cm,
  dgdtmax = 8000; % 
  get_traj;
  k = k(NN(2)+1:end,:);
  Nt = size(k,1);
  dt = 10e-6;
  tb0 = (0:Nt-1)*dt;
%Load excitation pattern and center it around the mask
  load pattern_rect.mat
  dim = size(d);
  
  [xi,yi] = ndgrid(-dim(1)/2:dim(1)/2-1,-dim(2)/2:dim(2)/2-1);
  xmi = round(mean(xi(logical(mask(:)))));
  ymi = round(mean(yi(logical(mask(:)))));
  d = circshift(d,[xmi ymi]); %Center the excitation pattern
  d = ift2((hamming(dim(1))*hamming(dim(2))').^2.*ft2(d));%Smooth it a bit.
  
%% Build RF waveform roughness penalty matrix 
    R = sqrt(roughbeta)*spdiags([-ones(Nt,1) ones(Nt,1)],[0 1],Nt,Nt) + ...
        sqrt(roughbeta)*speye(Nt);
    Rfull = kron(speye(Nc),R);

    % Solve for rf with system matrix
    DA = SysMtxSingle(k,b0,b1,dt,tb0,Positions);
    ncgiters = floor(size(k,1));       % # CG Iterations per RF update
    tic
    rf = qpwls_pcg(zeros(size(k,1)*Nc,1),DA,1,d(mask),0,Rfull,1,ncgiters,mask); % CG
    toc
    rf = reshape(rf(:,end),[length(k) Nc]); % qpwls_pcg returns all iterates
    plot(abs(rf));
    
  %%
  m = zeros(dim);m(mask) = DA*rf(:);
  err = norm(m-d)/norm(d);
  roughness = norm(Rfull*rf(:));
  pow = sqrt(mean(abs(rf(:)).^2));
  fprintf('All-channels NRMSE: %0.2f%%. All-channels roughness: %f.\n',err*100,roughness);
  %%plotting
  maxamp = max(abs([d(:);m(:)]));
  figure
  subplot(221)
  imagesc(abs(d),[0 maxamp]);axis image;colorbar
  title 'Desired pattern'
  subplot(222)
  imagesc(abs(m),[0 maxamp]);axis image;colorbar
  title(sprintf('Final pattern',Nc));
  subplot(224)
  imagesc(abs(m-d));axis image;colorbar
  title(sprintf('Error\nNRMSE = %0.2f%%',err*100));
  