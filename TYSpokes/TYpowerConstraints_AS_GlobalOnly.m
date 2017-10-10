function c = TYpowerConstraints_AS_GlobalOnly(x,currentConstraints,pulseTR,sinc_pulse)

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
power =@(x) x.^2./impedance;
    %voltage = @(p) sqrt(p*impedance);
    RFAmp = x(1:16);
    cPeakPwrPerChan = power(RFAmp) - currentConstraints.flPTxCoilPerChPeakPwrLimit;
    VoltsPerSpoke = sum(abs(sinc_pulse))*(10e-6);  % s   tyMakeSinc has a 10us raster
    PowerPerTR = (power(RFAmp(1:8)) + power(RFAmp(9:16))).*VoltsPerSpoke;    %W per TR
    cAvgPwrPerChan10s = PowerPerTR*(10/pulseTR)/10 - currentConstraints.flPTxCoilPerChAvPwrLimit_10s;
    cAvgPwrPerChan6min = PowerPerTR*(6*60/pulseTR)/(6*60) - currentConstraints.flPTxCoilPerChAvPwrLimit_LT;
    SumPwrPerTR = sum(PowerPerTR,1);    %W per TR
    cSumPwrPerTR10s = SumPwrPerTR.*VoltsPerSpoke*(10/pulseTR)/10 - currentConstraints.flPTxCoilSumAvPwrLimit_10s; 
    cSumPwrPerTR6min = SumPwrPerTR.*VoltsPerSpoke*(6*60/pulseTR)/(6*60) - currentConstraints.flPTxCoilSumAvPwrLimit_LT; 
 
 c = [cPeakPwrPerChan;cAvgPwrPerChan10s;cAvgPwrPerChan6min;cSumPwrPerTR10s;cSumPwrPerTR6min];


