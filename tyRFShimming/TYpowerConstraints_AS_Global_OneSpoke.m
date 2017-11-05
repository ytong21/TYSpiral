function [c,ceq] = TYpowerConstraints_AS_Global_OneSpoke(x,currentConstraints,pulseTR,RF_pulse)
ceq = [];
% c<=0
% ceq=0
% See https://uk.mathworks.com/help/optim/ug/nonlinear-constraints.html for
% the form of g and geq

% x is the vector of optimisation parameters (typically RF mag, RF phase, gradient??)
% currentConstraints are the constraints defined in the coilfile, modified suitabily for any other pulses using some of the limit
% pulseTR is how often the pulse gets repeated. All calculations currently assume that the sequence will run for a full 10s/6minutes uninterupted.
% conversionFunction    is a function handle which converts from the
%                       optimisation parameters to a full description of the RF (and gradients in
%                       the future), the second output is an magnitude integral
%                       corresponding to each optimisation point.

% equality constraints are unused but must be returned.


% Initial calculations

%First Nc*2 elements =  abs parts of complex per channel weights for 2 spokes
%Second Nc*2 elements = angle parts of complex per channel weights for 2 spokes
%For case 'Kb', last two elements = [deltaKx;deltaKy];
    % Define
    impedance = 50; % ohms matched
    WTCpower =@(x) x.^2./impedance;
    %voltage = @(p) sqrt(p*impedance);
    if numel(x) == 16   %   Full, amp and phase of 8 channels
        RFAmp = x(1:8);
    elseif numel(x) == 1    %   CP mode. 1 amp value only.
        RFAmp = ones(8,1)*x;
    elseif numel(x) == 8    %   Phase-only. 1 amp and 8 phases. 
        RFAmp = ones(8,1)*x(1);
    else
        disp('Please check the size of input vector')
    end
    cPeakPwrPerChan = WTCpower(RFAmp) - currentConstraints.flPTxCoilPerChPeakPwrLimit;
    
    RFAmpFull = zeros(8,numel(RF_pulse));
    for iDx = 1:8
        RFAmpFull(iDx,:) = RFAmp(iDx)*RF_pulse;
    end
    PowerIntPerSpoke = sum(WTCpower(RFAmpFull)*(10E-6),2);
    AvgPowerPerTR = PowerIntPerSpoke/pulseTR;
    %VoltsPerSpoke = sum(abs(RF_pulse))*(10e-6);  % s   tyMakeSinc has a 10us raster
    %PowerPerTR = (WTCpower(RFAmp)).*VoltsPerSpoke;    %W per TR
    RelaxFactor = 3;
    cAvgPwrPerChan10s = AvgPowerPerTR/RelaxFactor - currentConstraints.flPTxCoilPerChAvPwrLimit_10s;
    cAvgPwrPerChan6min = AvgPowerPerTR/RelaxFactor - currentConstraints.flPTxCoilPerChAvPwrLimit_LT;
    SumPwrPerTR = sum(AvgPowerPerTR,1);    %W per TR
    cSumPwrPerTR10s = SumPwrPerTR/RelaxFactor - currentConstraints.flPTxCoilSumAvPwrLimit_10s; 
    cSumPwrPerTR6min = SumPwrPerTR/RelaxFactor - currentConstraints.flPTxCoilSumAvPwrLimit_LT; 
 
 c = [cPeakPwrPerChan; cAvgPwrPerChan10s; cAvgPwrPerChan6min; cSumPwrPerTR10s; cSumPwrPerTR6min];


