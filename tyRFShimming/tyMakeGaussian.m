function RFStruct = tyMakeGaussian(Duration,SliceThickness)

%Duration "us"
%SliceThickness "mm"
    DurationInSec = Duration*(1e-6);
    SliceThicknessInMeter = SliceThickness*(1e-3);
    %GyroRatio = 42.57e6; %Hz/T
    T1 = 1500e-3;
    T2 = 70e-3;
    
dt = 10; % in us.
RF_pulse = gaussian_RF_dsv(dt, Duration, 1);
RF_pulse = RF_pulse(1,:);
% Note that gaussian_RF_dsv, the resulting RF pulse shape always has an odd
% number of points. 
%Grissom seems to prefer g/cm and cm.
% Building the gradient shape, ramp up and ramp down of 1000us included 
    G_shape = [linspace(0,1,25),ones(1,numel(RF_pulse)),linspace(1,0,25)]';
    RFOn = false(size(G_shape));
    RFOn(G_shape == 1) = true;

    OffResonance = linspace(-5000,5000,2000);
    % In the dupas paper, FAr is the FA produced by a pulse whose peak
    % amplutude is 1uT, which is 42.57Hz.
    RF_pulse_bloch = 42.57*RF_pulse;
    [mxss1,myss1_mt,~] =  bloch_CTR_Hz...
            (complex(real(RF_pulse_bloch),imag(RF_pulse_bloch)),zeros(numel(RF_pulse_bloch),3),...
            10E-6,T1,T2,OffResonance,[0,0,0],0);
    mxy = abs(mxss1+1i*myss1_mt);
    MagOnRes = max(mxy(:));
    FArBloch = asin(MagOnRes);
    
    %Calculating required gradient amplitutde

DeltaF = fwhm(OffResonance,mxy);
GyroRatioPerMT = 42.57e3; %Hz/T
G_amp = (2*pi*DeltaF)/(2*pi*GyroRatioPerMT*SliceThicknessInMeter); %mT/m



RFStruct = struct('RF_pulse',RF_pulse,'GradShape',G_shape*G_amp,...
    'RFOn',RFOn,'DurationInSec',DurationInSec,'DeltaF',DeltaF,'mxy',mxy,...
    'FArBloch',FArBloch,'G_amp',G_amp);
RFStruct.pulseInt = trapz(RF_pulse)*dt*1e-6;
RFStruct.unitEnergy = trapz(RF_pulse.^2/50)*dt*1e-6;