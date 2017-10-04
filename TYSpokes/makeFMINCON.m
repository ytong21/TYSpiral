function [bOut,fval,exitflag,output] = makeFMINCON(OptimType,deltaK,bVE,SINC,maskedMaps,param)


VecInInitial = bVE;

%Setting upper and lower bounds
maxV = 239; % Check where I got this.
clear ub lb
ub(1:param.numCh*2) = maxV;
ub((1+param.numCh*2):(param.numCh*4)) = inf;

lb(1:param.numCh*2) = 0;
lb((1+param.numCh*2):(param.numCh*4)) = -inf;

%Setting optimization options

optionsFMin = optimoptions(@fmincon);
optionsFMin.Algorithm = 'active-set';
optionsFMin.Display = 'none';
optionsFMin.MaxFunctionEvaluations = param.MaxEvaluation;
optionsFMin.SpecifyConstraintGradient = true;
optionsFMin.OptimalityTolerance = param.tol;

switch OptimType 
    case 'Kb'
        FOX = 25; %25cm 
        ub(param.numCh*4+1:param.numCh*4+2) = [3.5*(2*pi)/FOX 3.5*(2*pi)/FOX];   
        lb(param.numCh*4+1:param.numCh*4+2) = [-3.5*(2*pi)/FOX -3.5*(2*pi)/FOX];
        VecInInitial(param.numCh*4+1:param.numCh*4+2,1) = deltaK;
        FunFinal = getAMatSpokesKb(SINC,maskedMaps,param);
    case 'b'
        FunFinal = getAMatSpokesb(SINC,maskedMaps,deltaK,param);
end

[bOut,fval,exitflag,output] = fmincon(FunFinal,VecInInitial,[],[],[],[],lb.',ub.',[],optionsFMin);