function NRMSE = getAMatSpokesb(SINC,maskedMaps,deltaK,param)
%(Nc*2)*2+2 by 1 Vec
%First Nc*2 elements =  abs of complex per channel weights for 2 spokes
%Second Nc*2 elements = angle of complex per channel weights for 2 spokes
%Duration in seconds
%sense in uT/Volt

sens = maskedMaps.b1SensMasked;
dp = maskedMaps.posVox;
df = maskedMaps.b0MapMasked;
Nv = size(df,1);

if size(sens,1) ~= size(dp,1)
    error('Sensitivity first dimension must be the same size as the position first dimension.');
end

    NCha = size(sens,2);
    TWOPI = 2*pi;
    totalDuration = SINC.DurationInSec;


%deltaK = @(VecIn) VecIn(end-1:end,:);
bInAbs = @(VecIn) VecIn(1:NCha*2,:);
bInAngle = @(VecIn) VecIn(NCha*2+1:NCha*4,:);
bInComplex = @(VecIn) bInAbs(VecIn).*(cos(bInAngle(VecIn))+1i*sin(bInAngle(VecIn)));

    b0TermSpoke1 = exp(-1i*TWOPI*df*(totalDuration*1.5));
    ASpoke1 = b0TermSpoke1;
    b0TermSpoke2 = exp(-1i*TWOPI*df*(totalDuration*0.5));
    PhaseBlip = exp(1i*(deltaK(1)*dp(:,1)+deltaK(2)*dp(:,2)));
    ASpoke2 = bsxfun(@times,b0TermSpoke2,PhaseBlip);
    AFullSpoke1 = zeros(size(ASpoke1,1),NCha);
    AFullSpoke2 = zeros(size(ASpoke2,1),NCha);
     
for cDx = 1:NCha
    AFullSpoke1(:,cDx) = bsxfun(@times,ASpoke1,sens(:,cDx));
    AFullSpoke2(:,cDx) = bsxfun(@times,ASpoke2,sens(:,cDx));
end

AFullSpokes = 1i*SINC.FAr*[AFullSpoke1 AFullSpoke2];

FinalFA = @(VecIn) AFullSpokes*bInComplex(VecIn);
TargetFA = ones(Nv,1)*deg2rad(param.targetFlipAngle);
NRMSE = @(VecIn) norm(abs(FinalFA(VecIn)) - TargetFA)/norm(TargetFA);
    