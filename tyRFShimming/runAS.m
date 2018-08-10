function [bOut,fval,exitflag,output] = runAS(bVE,RFStruct,maskedMaps,param,AFull)


% sens = maskedMaps.b1SensMasked;
% dp = maskedMaps.posVox;
    %df = maskedMaps.b0MapMasked;


    %Setting upper and lower bounds
    maxV = 239; % Check where I got this.
    clear ub lb
    ub(1:param.numCh) = maxV;
    ub((1+param.numCh):(param.numCh*2)) = inf;

    lb(1:param.numCh) = 0;
    lb((1+param.numCh):(param.numCh*2)) = -inf;

    %Setting optimization options

    optionsFMin = optimoptions(@fmincon);
    optionsFMin.Algorithm = 'active-set';
    optionsFMin.Display = 'none';
    optionsFMin.MaxFunctionEvaluations = param.MaxEvaluation;
    optionsFMin.SpecifyConstraintGradient = false;
    optionsFMin.OptimalityTolerance = param.tol;
    optionsFMin.FiniteDifferenceType = 'central';
    optionsFMin.MaxIterations = 40000;
    protectedModeConstraints = CoilConstraints.novaCoil( true );
    nonlincon = @(x) TYpowerConstraints_AS_Global_OneSpoke(x,protectedModeConstraints,...
        param.RFSep,RFStruct.RF_pulse);
    
    % Creating a function handle based on getAMatSimp
    xInitial = [abs(bVE);angle(bVE)];
    %TargetFA = maskedMaps.TargetMasked*deg2rad(param.targetFlipAngle);
    TargetFA = ones(size(maskedMaps.b0MapMasked))*deg2rad(param.targetFlipAngle);
    FunHandle = @(x) (norm(abs(AFull*(x(1:8).*exp(1i*x(9:16)))) - TargetFA)/norm(TargetFA));
    %FunHandle = @(x) sum(abs(abs(AFull*(x(1:8).*exp(1i*x(9:16))))-TargetFA))/sum(TargetFA);
    %FunHandle = @(x) -goodnessOfFit(abs(AFull*(x(1:8).*exp(1i*x(9:16)))),TargetFA,'NRMSE');
    %options = optimoptions('fmincon','MaxIterations',4000);
    [bOut,fval,exitflag,output] = fmincon(FunHandle,xInitial,[],[],[],[],lb.',ub.',nonlincon,optionsFMin);
    
    
end