function [bOut, ErrorOut] = runPhaseOnly(AFull,param,RFStruct)
    maxV = 239; % Check where I got this. (WTC)
    ub = [maxV; ones(8,1)*inf];
    lb = [0; -ones(8,1)*inf];
    
    optionsFMin = optimoptions(@fmincon);
    optionsFMin.Algorithm = 'active-set';
    optionsFMin.Display = 'none';
    optionsFMin.MaxFunctionEvaluations = 5000;
    optionsFMin.SpecifyConstraintGradient = false;
    optionsFMin.OptimalityTolerance = param.tol;
    optionsFMin.FiniteDifferenceType = 'central';
    
    protectedModeConstraints = CoilConstraints.novaCoil( true );
    nonlincon = @(x) TYpowerConstraints_AS_Global_OneSpoke(x,protectedModeConstraints,...
        param.TR,RFStruct.RF_pulse);

    TargetFA = ones(size(AFull,1),1)*deg2rad(param.targetFlipAngle);
    AFullCP = AFull*ones(8,1);
    VoltCP = abs((TargetFA')/(AFullCP'));
    FunHandle = @(x) -goodnessOfFit(abs(AFull*(x(1)*exp(1i*x(2:9)))),TargetFA,'NRMSE');
    Iterations = 20;
    ErrorOut = zeros(1,Iterations);
    bOut = zeros(9,Iterations);
    for iDx = 1:Iterations
        xInitial = [VoltCP; (rand(8,1)*(2*pi)-pi)];
        [bOut(:,iDx),ErrorOut(iDx),~,~] = fmincon(FunHandle,xInitial,[],[],[],[],lb.',ub.',nonlincon,optionsFMin);
    end
    %bPhase = bOut(1)*exp(1i*bOut(2:9));
end