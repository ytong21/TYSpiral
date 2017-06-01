%addpath(genpath('/Users/ytong/Documents/MATLAB/WTCpTx'));
singleTxPath = '/Users/ytong/Documents/MATLAB/20170516_F7T_2017_PH_038';
%singleTxPath = '/Volumes/Data/DICOM/2017-05/20170516_F7T_2017_PH_038';
singleTxObj = DicomFM.WTCSingleTxFieldmaps(singleTxPath);
singleTxObj.interpolateTo('B1');
singleTxObj.createMask(30,false);
%%
B1 = singleTxObj.getB1('Gauss','delete',[]);
B0 = singleTxObj.getB0('Hz','delete');
Loc = singleTxObj.getLoc('delete');
fov = [25.6 25.6];
[DimX, DimY, NumSlices] = size(B0);clc
positions = singleTxObj.getPositions('cm','delete');
WTCmask = singleTxObj.getMask;

%%