function [] = STASpokes(AFullSpokes,param,maskedMaps)
%%
% Define some limits
% RF limits / parameters
maxV = 239; % Max voltage permitted by scanner on each channel.

%  The following limits are taken from the coil file of the Nova Medical
%  8Tx/32Rx head coil.
% Component protection limits
powerLimit.maxPerChannel = 1250; % Instantaneous maximum power of any one channel.

% Protected mode constraints
powerLimit.sum6Min = 8; % W. Limit on the average power of the magnitude sum of all channels over 6 minutes.
powerLimit.sum10Sec = 8*2; % W. Limit over 10 seconds

powerLimit.perChannel6Min = 1.5; % W. Limit on the average per-channel power over 6 minutes. 
powerLimit.perChannel10Sec = 1.5*2; % W. Limit over 10 seconds.

% Gradient limits
gradLimit.maxAmp = 4; %mT/m
gradLimit.maxSlew = 100; % T/m/s

%%
    tikhonovVec = [ 0 10.^(-7:0.5:1)];
    for tDx = 1:numel(tikhonovVec)
        param.CGtikhonov = tikhonovVec(tDx);
        [bOut(:,tDx),costOut(tDx),powerOut(tDx),magOut(:,tDx)] = designSpokes(AFullSpokes,param,maskedMaps);
    end
 %% Check we are within power constraints
% Calculate the power for each
addpath('/Users/ytong/Documents/MATLAB/WTCpTx/pulseSAR/ProtectedMode');
impedance = 50; %ohm
power = @(x) abs(x).^2/impedance; % Watts

for iDx = 1:size(bOut,2)
    pulseStruct.points = reshape(bOut(:,iDx),[designParams.numCh designParams.points]);
    pulseStruct.timestep = designParams.subPulseDuration;
    [perChannel10s(:,iDx), perChannel6min(:,iDx), sum10s(iDx), sum6min(iDx)] = calculateAverageSAR(pulseStruct,repetitionTime);
    
    % Max instantaneous power
    maxInstPower(:,iDx) = max(power(pulseStruct.points),[],2);
    
end

powerLimtCheck.maxPerChannel = maxInstPower>powerLimit.maxPerChannel;
powerLimtCheck.perChannel10Sec = perChannel10s>powerLimit.perChannel10Sec;
powerLimtCheck.perChannel6Min = perChannel6min>powerLimit.perChannel6Min;
powerLimtCheck.sum10Sec = sum10s>powerLimit.sum10Sec;
powerLimtCheck.sum6Min = sum6min>powerLimit.sum6Min;


powerLimtCheck.all = any([powerLimtCheck.maxPerChannel;powerLimtCheck.perChannel10Sec;...
    powerLimtCheck.perChannel6Min;powerLimtCheck.sum10Sec;powerLimtCheck.sum6Min],1);

%% plot l-curve
figure(150)
clf
plot(powerOut,costOut,'x--')
hold on
for iDx = 1:numel(powerLimtCheck.all)
    if powerLimtCheck.all(iDx)==1
        plot(powerOut(iDx),costOut(iDx),'rx')
    else
        plot(powerOut(iDx),costOut(iDx),'gx')
    end
end
% Select the point to use from the l-curve analysis and the restrictions on
% power
pointToUse = find(~powerLimtCheck.all,1,'first');
bFinal = bOut(:,pointToUse);