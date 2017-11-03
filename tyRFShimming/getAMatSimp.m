function AFullSpokes = getAMatSimp(RFStruct,sens,df)

%RF_pulse = RF_pulse_mag.*exp(1i*RF_pulse_phase);
 


%Duration in seconds
%sense in uT/Volt
TWOPI = 2*pi;
totalDuration = RFStruct.DurationInSec;
%PhaseBlip = zeros(size(df));
    b0TermSpoke1 = exp(-1i*TWOPI*df*(totalDuration*1.5));
    ASpoke1 = b0TermSpoke1;
    AFullSpoke1= bsxfun(@times,ASpoke1,sens);

AFullSpokes = 1i*RFStruct.FArBloch*AFullSpoke1;
