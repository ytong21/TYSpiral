function bPhase = runPhaseOnly(AFull,param,RFStruct)
    maxV = 239; % Check where I got this. (WTC)
    ub = [maxV; ones(8,1)*inf];
    lb = [0; -ones(8,1)*inf];
    
    
    optionsFMin = optimoptions(@fmincon);
    optionsFMin.Algorithm = 'active-set';
    optionsFMin.Display = 'none';
    optionsFMin.MaxFunctionEvaluations = 20000;
    optionsFMin.SpecifyConstraintGradient = false;
    optionsFMin.OptimalityTolerance = param.tol;
    optionsFMin.FiniteDifferenceType = 'central';
    
    protectedModeConstraints = CoilConstraints.novaCoil( true );
    nonlincon = @(x) TYpowerConstraints_AS_Global_OneSpoke(x,protectedModeConstraints,...
        param.TR,RFStruct.RF_pulse);

    TargetFA = ones(size(AFull,1),1)*deg2rad(param.targetFlipAngle);
    FunHandle = @(x) -goodnessOfFit(abs(AFull*(x(1)*exp(1i*x(2:9)))),TargetFA,'NRMSE');
    %VoltageArray = 10:5:80;
    %bSingle = zeros(size(VoltageArray));
    for iDx = 1:numel(VoltageArray)
        [bPhase,~,~,~] = fmincon(FunHandle,xInitial,[],[],[],[],lb.',ub.',nonlincon,optionsFMin);
    end
end