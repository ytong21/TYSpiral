function AFullSpokes = getAMatSimp(RFStruct,sens,df,dp)

%RF_pulse = RF_pulse_mag.*exp(1i*RF_pulse_phase);
 
if size(sens,1) ~= size(dp,1)
    error('Sensitivity first dimension must be the same size as the position first dimension.');
end

%Duration in seconds
%sense in uT/Volt
TWOPI = 2*pi;
totalDuration = RFStruct.DurationInSec;
%PhaseBlip = zeros(size(df));
    b0TermSpoke1 = exp(-1i*TWOPI*df*(totalDuration*1.5));
    ASpoke1 = b0TermSpoke1;
    AFullSpoke1= bsxfun(@times,ASpoke1,sens);

AFullSpokes = 1i*RFStruct.FArBloch*AFullSpoke1;
