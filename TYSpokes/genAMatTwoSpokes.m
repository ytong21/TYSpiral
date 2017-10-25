function AFullSpokes = genAMatTwoSpokes(SINC,sens,df,dp,deltaK)

 
if size(sens,1) ~= size(dp,1)
    error('Sensitivity first dimension must be the same size as the position first dimension.');
end

%Duration in seconds
%sense in uT/Volt

TWOPI = 2*pi;
totalDuration = SINC.DurationInSec;
%PhaseBlip = zeros(size(df));
    b0TermSpoke1 = exp(-1i*TWOPI*df*(totalDuration*1.5));
    ASpoke1 = b0TermSpoke1;
    b0TermSpoke2 = exp(-1i*TWOPI*df*(totalDuration*0.5));
    PhaseBlip = exp(1i*(deltaK(1)*dp(:,1)+deltaK(2)*dp(:,2)));
    ASpoke2 = bsxfun(@times,b0TermSpoke2,PhaseBlip);
%     AFullSpoke1 = zeros(size(ASpoke1,1),NCha);
%     AFullSpoke2 = zeros(size(ASpoke2,1),NCha);
     
% for cDx = 1:NCha
    AFullSpoke1= bsxfun(@times,ASpoke1,sens);
    AFullSpoke2 = bsxfun(@times,ASpoke2,sens);
% end

AFullSpokes = 1i*SINC.FArBloch*[AFullSpoke1 AFullSpoke2];
