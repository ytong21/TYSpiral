function [Efficiency,FinalMag] = LabelEff(FAinDeg,RFStruct, RF_sep_sec)
Efficiency = zeros(size(FAinDeg));
FinalMag = zeros(size(FAinDeg));
    FAinRad = deg2rad(FAinDeg);
    RFAmpInuT = FAinRad/RFStruct.FArBloch;
    RFAmpInmT = RFAmpInuT/1000;
    
    for iDx = 1:numel(FAinRad)
        [M, ~, ~, ~, ~] = test_regPCASL_seq_meanGz(0.3, 0.1, 0, 0.8, 6, 'gaussian_dsv', ...
            [], RFAmpInmT(iDx), RFStruct.DurationInSec, RF_sep_sec, [0 0], [1 1], 1, [0 0], 20e-6, inf, 0.3, 0);
        FinalMag(iDx) = M(3,end);
        Efficiency(iDx) = (1-FinalMag(iDx))/2;
    end