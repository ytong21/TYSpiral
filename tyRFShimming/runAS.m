function [bOut,fval,exitflag,output] = runAS(bVE,RFStruct,maskedMaps,param,AFull)


% sens = maskedMaps.b1SensMasked;
% dp = maskedMaps.posVox;
df = maskedMaps.b0MapMasked;


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
    
    protectedModeConstraints = CoilConstraints.novaCoil( true );
    nonlincon = @(x) TYpowerConstraints_AS_Global_OneSpoke(x,protectedModeConstraints,...
        param.TR,RFStruct.RF_pulse);
    
    % Creating a function handle based on getAMatSimp
    xInitial = [abs(bVE);angle(bVE)];
    TargetFA = ones(size(df))*deg2rad(param.targetFlipAngle);
    %FunHandle = @(x) norm(abs(AFull*(x(1:8).*exp(1i*x(9:16)))) - TargetFA)/norm(TargetFA);
    FunHandle = @(x) -goodnessOfFit(abs(AFull*(x(1:8).*exp(1i*x(9:16)))),TargetFA,'NRMSE');
    
    [bOut,fval,exitflag,output] = fmincon(FunHandle,xInitial,[],[],[],[],lb.',ub.',nonlincon,optionsFMin);
    
    
end