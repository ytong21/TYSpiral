function bCP = runCP(AFull,param,RFStruct)
    maxV = 239; % Check where I got this. (WTC)
    ub = maxV;
    lb = 0;
    
    AFullCP = AFull*ones(8,1);
    
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

    TargetFA = ones(size(AFullCP,1),1)*deg2rad(param.targetFlipAngle);
    xInitial = abs((TargetFA')/(AFullCP'));
    FunHandle = @(x) -goodnessOfFit(abs(AFullCP*x),TargetFA,'NRMSE');
    %VoltageArray = 10:5:80;
    %bSingle = zeros(size(VoltageArray));
    %for iDx = 1:numel(VoltageArray)
        [bSingle,~,~,~] = fmincon(FunHandle,xInitial,[],[],[],[],lb.',ub.',nonlincon,optionsFMin);
    %end
    bCP = ones(8,1)*bSingle;
end