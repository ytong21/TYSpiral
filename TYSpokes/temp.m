

    kXAdj = linspace(-3,3,11);
    kYAdj = linspace(-3,3,11);
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
parfor iDx = 1:size(KVec,2)
    [bOut(iDx,:),kOut(iDx,:)] = VE_AS(OptimType,KVec(:,iDx),SINC,maskedMaps,param);
end
kOut = kOut*FOX/(2*pi);
kOut = kOut';
%%        
figure(102)
clf
imagesc(deltaXAdj,deltaYAdj,StructToSave.NRMSE)
colorbar

hold on
plot(KVec(1,:)*FOX/(2*pi),KVec(2,:)*FOX/(2*pi),'.k', 'MarkerSize', 18);
hold on
plot(kOut(1,:),kOut(2,:),'+r', 'MarkerSize', 10);

