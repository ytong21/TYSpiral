function [bOut,finalRMSE,finalPwr,finalMag] =runVE(AFullFunc,param,maskedMaps)

    MaxVEIter = 50;
    cost = inf;
    xCurr = zeros(param.numCh,1);
    Tol = 0.001;    
    %targetFAInRad = maskedMaps.TargetMasked*deg2rad(param.targetFlipAngle);
    targetFAInRad = ones(size(maskedMaps.b0MapMasked))*deg2rad(param.targetFlipAngle);
    
    phiTarget = angle(sum(maskedMaps.b1SensMaskedHz,2));
    z = exp(1i*phiTarget);
    fullImage =  zeros(numel(maskedMaps.mask),1);
    CGtikhonov = param.CGtikhonov;
    
ALambda = pinv(AFullFunc'*AFullFunc + CGtikhonov*eye(numel(xCurr)))*AFullFunc';

for iDx = 1:MaxVEIter
    
    targetFA = z .* targetFAInRad;
    
%     [xs,~ ] = qpwls_pcg(xCurr, AFullFunc,1, targetM, 0, CGtikhonov, 1, cgIterations);
%     xCurr = xs(:,end);
    xCurr = ALambda*targetFA;

    
    currFA = AFullFunc*xCurr;
    costNew = norm(currFA - targetFA) + CGtikhonov*(xCurr'*xCurr);
    
    costDiff = abs(cost-costNew)/cost;
    if costDiff < Tol
        break
    end
    cost = costNew;
    
    % Smooth phase
    fullImage(logical(maskedMaps.mask(:))) = angle(currFA);
    fullImage = reshape(fullImage,size(maskedMaps.mask));
    fullImage = imgaussfilt3(fullImage,2,'FilterSize',5);
        
    phiTarget = fullImage(maskedMaps.mask(:));
    z = exp(1i*phiTarget);    
end

bOut = xCurr;
finalPwr = norm(xCurr);
%finalCost = costNew;
%finalCost = cost;
    finalMag = zeros(numel(maskedMaps.mask),1);
    finalMag(logical(maskedMaps.mask(:))) = AFullFunc*xCurr;
    finalMag = reshape(finalMag,size(maskedMaps.mask));

%finalRMSE = norm(currFA - targetFA);
finalRMSE = norm(abs(currFA) - abs(targetFA))/norm(abs(targetFA));

% fprintf('\n Tikhonov = %f \n',CGtikhonov)
% fprintf('Power = %f \n',finalPwr)
% fprintf('Cost = %f \n',finalCost))