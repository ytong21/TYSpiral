function Landscape = landscape_run(StepSize,param,maskedMaps,SINC)
%   Optimization of complex Sinc pusle weights  of a two-spoke RF pulse based on
%   initialization of N by N locations on a Cartesian grid.
%   The first step is via a variable exchange method with the CP mode 
%   B1 map phase as the target phase. The second step is via an Active-Set
%   algorithm. 

    deltaXAdj = -3:StepSize:3;
    deltaYAdj = -3:StepSize:3;
    FOX = 25; 
    deltaX = deltaXAdj*(2*pi)/FOX;
    deltaY = deltaYAdj*(2*pi)/FOX;
    

     NumIt = 0;
     deltaKVec = zeros(2,numel(deltaX)*numel(deltaY));
     NRMSE = zeros(numel(deltaX)*numel(deltaY),1);
     bOut = zeros(numel(deltaX)*numel(deltaY),2*param.numCh);
for xDx = 1:numel(deltaX)
    for yDx = 1:numel(deltaY)
        NumIt = NumIt+1;
        deltaKVec(:,NumIt) = [deltaX(xDx);deltaY(yDx)];
    end
end
OptimType = 'b';
tic
parfor iDx = 1:size(deltaKVec,2)
    [bOut(iDx,:),NRMSE(iDx,:),~] = VE_AS(OptimType,deltaKVec(:,iDx),SINC,maskedMaps,param);
end
toc
    Landscape = struct('kVec',deltaKVec);
    Landscape.ComputationTime = toc;
    bOut = reshape(bOut,[numel(deltaX),numel(deltaY),2*param.numCh]);
    NRMSE = reshape(NRMSE,[numel(deltaX),numel(deltaY)]);
    Landscape.bOut = bOut;   Landscape.NRMSE = NRMSE;
    