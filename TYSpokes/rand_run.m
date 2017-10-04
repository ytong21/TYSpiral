function NRMSE_Vec = rand_run(NumIt,Numk,param,maskedMaps,SINC)
%   Optimization of the relative k-space location and complex Sinc pusle weights 
%   of a two-spoke RF pulse based on initialization of NN random locations.
%   The first step is via a variable exchange method with the CP mode 
%   B1 map phase as the target phase. The second step is via an Active-Set
%   algorithm, where the k-space location and RF pulse weights are jointly
%   optimized. 
    function minNRMSE = rand_run_internal(Numk,param,maskedMaps,SINC)
        FOX = param.FOX;                %Field of excitation in cm.
        kVecAdj = rand(2,Numk)*6-3;      %randomly distributed in [-3,3].
        kVec = kVecAdj * (2*pi)/FOX;    %kVec in 1/cm.
        OptimType = 'Kb';
        % Runing the optimization
        %kOut = zeros(size(kVecAdj));
        %bOut = zeros(size(kVecAdj));
        NRMSE = zeros(size(kVecAdj,2),1);
        for iDx = 1:size(kVec,2)
            [~,NRMSE(iDx,:),~] = VE_AS(OptimType,kVec(:,iDx),SINC,maskedMaps,param);
        end
        minNRMSE = min(NRMSE(:));
    end
FuncHandle = @rand_run_internal;
NRMSE_Vec = zeros(NumIt,1);
parfor kk = 1:NumIt
    NRMSE_Vec(kk) = feval(FuncHandle,Numk,param,maskedMaps,SINC);
end
    
end

