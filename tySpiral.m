%addpath(genpath('/Users/ytong/Documents/MATLAB/WTCpTx'));
singleTxPath = '/Users/ytong/Documents/MATLAB/20170516_F7T_2017_PH_038';
%singleTxPath = '/Volumes/Data/DICOM/2017-05/20170516_F7T_2017_PH_038';
singleTxObj = DicomFM.WTCSingleTxFieldmaps(singleTxPath);
singleTxObj.interpolateTo('B1');
singleTxObj.createMask(30,false);
B1 = singleTxObj.getB1('Gauss','delete',[]);
B0 = singleTxObj.getB0('Hz','delete');
Loc = singleTxObj.getLoc('delete');
%fov = [25.6 25.6];
[dimb1(1),dimb1(2), ~] = size(WTCmask);
[DimX, DimY, NumSlices] = size(B0);clc
positions = singleTxObj.getPositions('cm','delete');
WTCmask = singleTxObj.getMask;
Ns = size(B1,1);
fov = 25.6; %fov is hard coded at this momemnt. In cm.
%%
%Get k-space trajectory
  traj = 'spiral';
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
  d = circshift(d,[xmi ymi]);
  d = ift2((hamming(dim(1))*hamming(dim(2))').^2.*ft2(d));
  