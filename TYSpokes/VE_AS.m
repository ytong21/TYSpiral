function [bOut,NRMSE,kOut] = VE_AS(OptimType,deltaK,SINC,maskedMaps,param)       

        AFullSpokes = genAMatTwoSpokes(SINC,maskedMaps.b1SensMasked,maskedMaps.b0MapMasked,...
            maskedMaps.posVox,deltaK); %Nv-(2*Nc) sys mtx, rad/V
        [bVE,~,~,~] = designSpokes(AFullSpokes,param,maskedMaps);
switch OptimType
    case 'b'
        [bOutReal,NRMSE] = makeFMINCON(OptimType,deltaK,[abs(bVE); angle(bVE)],SINC,maskedMaps,param);
        bOut = complex(bOutReal(1:2*param.numCh),bOutReal(1+2*param.numCh:4*param.numCh));
        kOut = [];
        
    case 'Kb';
        [OutReal,NRMSE] = makeFMINCON(OptimType,deltaK,[abs(bVE); angle(bVE);deltaK],SINC,maskedMaps,param);
        bOut = complex(OutReal(1:2*param.numCh),OutReal(1+2*param.numCh:4*param.numCh));
        kOut = OutReal(4*param.numCh+1:end);
end