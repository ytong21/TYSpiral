%addpath(genpath('/Users/ytong/Documents/MATLAB/WTCpTx'));
singleTxPath = '/Users/ytong/Documents/MATLAB/20170516_F7T_2017_PH_038';
addpath ../spiral/util
%singleTxPath = '/Volumes/Data/DICOM/2017-05/20170516_F7T_2017_PH_038';
    singleTxObj = DicomFM.WTCSingleTxFieldmaps(singleTxPath);
    singleTxObj.interpolateTo('B1');
    singleTxObj.createMask(30,false);
    slice = round(size(singleTxObj.getMask([]),3))/2;
    mask = singleTxObj.getMask(slice);
    B1 = singleTxObj.getB1('Gauss','delete',slice);
    B0 = singleTxObj.getB1('Hz','delete',slice);
    Positions = singleTxObj.getPositions('cm','delete',slice);
    Ns = size(B1,1);
    fov = 25.6; %fov is hard coded at this momemnt. In cm.
%%
%Get k-space trajectory
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

%% Solve for rf with system matrix
DA = SysMtx(k,mask,Nc,b1,dt,tb0,b0,Positions);
  