function RFStruct = tyMakeHanning(Duration,SliceThickness)
%Duration "us"
%NumZeros unitless
%SliceThickness "mm"
    NumZeros = 2;
    T1 = 1500e-3;
    T2 = 70e-3;
%FAr in rad/uT
    DurationInSec = Duration*(1e-6);
    SliceThicknessInMeter = SliceThickness*(1e-3);
%Downsample by 10. Assuming peak amp is 1uT.
    DownSamp = 10;
    NumPoints = Duration/DownSamp;%10 us raster
    
RF_pulse = Hanning_window(NumPoints);
GyroRatio = 42.57e6; %Hz/T

%Calculating required gradient amplitutde
t0 = (Duration*(1e-6))/NumZeros;
DeltaF = 1/t0;
GyroRatioPerMT = 42.57e3; %Hz/T
G_amp = (2*pi*DeltaF)/(2*pi*GyroRatioPerMT*SliceThicknessInMeter); %mT/m

% Building k-space trajectory
%g = ones();

k = zeros(NumPoints,2);
k(:,2) = linspace(0,GyroRatio*G_amp*(Duration*(1e-6)),NumPoints); %1/m
k = k/100; %1/cm. For simulation use. 
%Grissom seems to prefer g/cm and cm.
% Building the gradient shape, ramp up and ramp down of 1000us included 
    G_shape = [linspace(0,1,25),ones(1,NumPoints),linspace(1,0,25)]';
    RFOn = false(size(G_shape));
    RFOn(G_shape == 1) = true;

    OffResonance = linspace(-5000,5000,2000);
    RF_pulse_bloch = 42.57*RF_pulse;
    [mxss1,myss1_mt,~] =  bloch_CTR_Hz...
            (complex(real(RF_pulse_bloch),imag(RF_pulse_bloch)),zeros(numel(RF_pulse_bloch),3),...
            10E-6,T1,T2,OffResonance,[0,0,0],0);
    mxy = abs(mxss1+1i*myss1_mt);
    MagOnRes = max(mxy(:));
    FArBloch = asin(MagOnRes);
%plotting
% figure(21)
% plot(linspace(0,Duration,NumPoints),sinc_pulse)
% xlim([0,Duration])

RFStruct = struct('kTraj',k,'RF_pulse',RF_pulse,'GradShape',G_shape*G_amp,...
    'RFOn',RFOn,'DurationInSec',DurationInSec,'DeltaF',DeltaF,'mxy',mxy,...
    'FArBloch',FArBloch);

