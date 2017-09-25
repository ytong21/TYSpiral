function rand_run(NN,param,maskedMaps,SINC)
%   Optimization of the relative k-space location and complex Sinc pusle weights 
%   of a two-spoke RF pulse based on initialization of NN random locations.
%   The first step is via a variable exchange method with the CP mode 
%   B1 map phase as the target phase. The second step is via an Active-Set
%   algorithm, where the k-space location and RF pulse weights are jointly
%   optimized. 

FOX = param.FOX;                %Field of excitation in cm.
kVecAdj = rand(2,NN)*6-3;      %randomly distributed in [-3,3].
kVec = kVecAdj * (2*pi)/FOX;    %kVec in 1/cm.
OptimType = 'Kb';
% Runing the optimization
parfor iDx = 1:size(kVec,2)
    [bOut(iDx,:),kOut(iDx,:)] = VE_AS(OptimType,kVec(:,iDx),SINC,maskedMaps,param);
end