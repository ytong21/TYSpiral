function [Cart,OutputStruct] = cartesian_run(MtxSize,param,maskedMaps,SINC)
%   Optimization of the relative k-space location and complex Sinc pusle weights 
%   of a two-spoke RF pulse based on initialization of N by N locations on a Cartesian grid.
%   The first step is via a variable exchange method with the CP mode 
%   B1 map phase as the target phase. The second step is via an Active-Set
%   algorithm where the k-space location and RF pulse weights are jointly
%   optimized.
%   Note that kIn/kOut are in units of (2*pi/FOX)
    FOX = param.FOX;
    kXAdj = linspace(-2.8,2.8,MtxSize);
    kYAdj = linspace(-2.8,2.8,MtxSize);
    kX = kXAdj*(2*pi)/FOX;
    kY = kYAdj*(2*pi)/FOX;
    KVec = zeros(2,numel(kX)*numel(kY));
    bOut = zeros(numel(kX)*numel(kY),2*param.numCh);    
    kOut = zeros(numel(kX)*numel(kY),2); 
    NumIt = 0;
for xDx = 1:numel(kX)
    for yDx = 1:numel(kY)
        NumIt = NumIt+1;
        KVec(:,NumIt) = [kX(xDx);kY(yDx)];
    end
end    
%%
NRMSE = zeros(size(KVec,2),1);
OptimType = 'Kb';
ExitFlag = zeros(size(KVec,2),1);
OutputStruct = cell(size(ExitFlag));
parfor iDx = 1:size(KVec,2)
    [bOut(iDx,:),NRMSE(iDx),kOut(iDx,:),ExitFlag(iDx),OutputStruct{iDx}] = VE_AS(OptimType,KVec(:,iDx),SINC,maskedMaps,param);
end
kOutAdj = kOut*FOX/(2*pi);
kOutAdj = kOutAdj';

Cart = struct('bOut',bOut,'kIn',KVec*FOX/(2*pi),'kOut',kOutAdj,'ExitFlag',ExitFlag,'NRMSE',NRMSE);
end
