function NRMSE = getAMatSpokesKb(SINC,maskedMaps,param)
%(Nc*2)*2+2 by 1 Vec
%First Nc*2 elements =  real parts of complex per channel weights for 2 spokes
%Second Nc*2 elements = real parts of complex per channel weights for 2 spokes
%Last two elements = [deltaKx;deltaKy];
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


deltaK = @(VecIn) VecIn(end-1:end,:);
bInAbs = @(VecIn) VecIn(1:NCha*2,:);
bInAngle = @(VecIn) VecIn(NCha*2+1:NCha*4,:);
bInComplex = @(VecIn) bInAbs(VecIn).*(cos(bInAngle(VecIn))+1i*sin(bInAngle(VecIn)));

%Calculate A for spoke 1
    b0TermSpoke1 = exp(-1i*TWOPI*df*(totalDuration*1.5));
    ASpoke1 = b0TermSpoke1;
    AFullSpoke1 = zeros(size(ASpoke1,1),NCha);
    for cDx = 1:NCha
        AFullSpoke1(:,cDx) = bsxfun(@times,ASpoke1,sens(:,cDx));
    end

%Calculate A for spoke 2 with anonymous function
    b0TermSpoke2 = exp(-1i*TWOPI*df*(totalDuration*0.5));

PhaseBlip = @(VecIn) exp(1i*(dp(:,1:2)*deltaK(VecIn)));
ASpoke2 = @(VecIn) PhaseBlip(VecIn).*b0TermSpoke2;
AFullSpoke2 = @(VecIn) repmat(ASpoke2(VecIn),1,NCha).*sens;
AFullSpokes = @(VecIn) 1i*SINC.FAr*[AFullSpoke1 AFullSpoke2(VecIn)];

FinalFA = @(VecIn) AFullSpokes(VecIn)*bInComplex(VecIn);
TargetFA = ones(Nv,1)*deg2rad(param.targetFlipAngle);
NRMSE = @(VecIn) norm(abs(FinalFA(VecIn)) - TargetFA)/norm(TargetFA);
    