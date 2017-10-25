function [bOut,NRMSE,kOut,exitflag,output] = VE_AS(OptimType,deltaK,SINC,maskedMaps,param)       

        AFullSpokes = genAMatTwoSpokes(SINC,maskedMaps.b1SensMasked,maskedMaps.b0MapMasked,...
            maskedMaps.posVox,deltaK); %Nv-(2*Nc) sys mtx, rad/V
        [bVE,~,~,~] = designSpokes(AFullSpokes,param,maskedMaps);
switch OptimType
    case 'b'
        [OutAll,NRMSE,exitflag,output] = makeFMINCON(OptimType,deltaK,[abs(bVE); angle(bVE)],SINC,maskedMaps,param);
        %bOut = complex(bOutReal(1:2*param.numCh),bOutReal(1+2*param.numCh:4*param.numCh));
        bOut = OutAll(1:2*param.numCh).*exp(OutAll(1+2*param.numCh:4*param.numCh));
        kOut = [];
        
    case 'Kb';
        [OutAll,NRMSE,exitflag,output] = makeFMINCON(OptimType,deltaK,[abs(bVE); angle(bVE);deltaK],SINC,maskedMaps,param);
        %bOut = complex(OutReal(1:2*param.numCh),OutReal(1+2*param.numCh:4*param.numCh));
        %bOut = OutAll(1:4*param.numCh);
        bOut = OutAll(1:2*param.numCh).*exp(OutAll(1+2*param.numCh:4*param.numCh));
        kOut = OutAll(4*param.numCh+1:end);
end
end