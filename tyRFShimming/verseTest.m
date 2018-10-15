RFStruct_Example = tyMakeGaussian(600,5);
b = [0 0 RFStruct_Example.RF_pulse 0 0];  % arbitrary unit, peak=1
g = [0 0 6*ones(size(RFStruct_Example.RF_pulse)) 0 0];            % in mT/m
dt = 10e-3;                     % in ms
bmax = 0.5;
gmax = 22;                      % in mT/m
smax = 200;                     % in mT/(m*ms)
dtout = 10e-3;                  % in ms
emax = sum(b.^2)*dt*0.5;        % in units of b*b*dt


[bv,gv] = mintverse(b,g,dt,bmax,gmax,smax,dtout,emax);
bv_fine = interp1(1:numel(bv),bv,1:0.2:numel(bv));
gv_fine = interp1(1:numel(gv),gv,1:0.2:numel(gv));

% write into a .h file
phase_vector_1mm = -0.001*cumsum(gv_fine*42577)*2E-6*2*pi;
MakeDotH_VERSE(bv_fine,gv,dt_fine,phase_vector_1mm);

%% Temporary values. Need a bit more work.
numCh = 8;
peakV = 60;
numsamples = numel(bv);
targetFlipAngle = 20;


%%
fID = fopen('/Users/ytong/Documents/XP/GaussianVERSE.ini','w');
writeGrad_mTm = [zeros(numel(gv),2) gv];
writeRfAmp = abs(bv);
writeRfPhase = angle(bv);
fprintf(fID,'# Created in verseTest.m \n');
fprintf(fID,'[pTXPulse]\n');
fprintf(fID,'\n');
fprintf(fID,'NUsedChannels        = %i\n',numCh);		
fprintf(fID,'MaxAbsRF             = %0.3f	# peak RF voltage [V]\n',peakV);
fprintf(fID,'Samples              = %i		# number of samples\n',numsamples);
fprintf(fID,'PulseName            = GaussianVERSE\n');
fprintf(fID,'Comment              = pTx pulse\n');
fprintf(fID,'Oversampling         = %i		# no oversampling\n',1);
fprintf(fID,'NominalFlipAngle     = %i		# flip-angle of this pulse if played out with ??V peak voltage\n',targetFlipAngle);
fprintf(fID,'SampleTime			 = 10\n');

% Gradient
fprintf(fID,'\n');
fprintf(fID,'[Gradient]\n');
fprintf(fID,'\n');
for iDx = 1:size(writeGrad_mTm,1)
    fprintf(fID,'G[%i]= %f %f %f\n',iDx-1,writeGrad_mTm(iDx,1),writeGrad_mTm(iDx,2),writeGrad_mTm(iDx,3));
end
fprintf(fID,'\n');

% RF
for jDx = 1:size(writeRfAmp,2)
    fprintf(fID,'[pTXPulse_ch%i]\n',jDx-1);
    fprintf(fID,'\n');
    for iDx = 1:size(writeRfAmp,1)
        fprintf(fID,'RF[%i]=  %f\t%f\n',iDx-1,writeRfAmp(iDx,jDx),writeRfPhase(iDx,jDx));
    end
    fprintf(fID,'\n');
end

fclose(fID);
%%
figure(100)
clf
FontSize = 16;
lgdFontSize = 14;
subplot(1,3,1)
plot(10*(1:numel(b)),b,'b')
hold on
plot(10*(1:numel(bv)),bv,'r--')
title('RF waveform','FontSize',FontSize);xlabel('Time (us)','FontSize',FontSize,'FontWeight','bold');
ylabel('Arbitrary unit','FontSize',FontSize,'FontWeight','bold');
lgd = legend('Gaussian','VERSE');lgd.FontSize = lgdFontSize;

subplot(1,3,2)
plot(10*(1:numel(g)),g,'b')
hold on
plot(10*(1:numel(gv)),gv,'r--')
title('Gradient waveform','FontSize',FontSize);xlabel('Time (us)','FontSize',FontSize,'FontWeight','bold');
ylabel('mT/m','FontSize',FontSize,'FontWeight','bold');
lgd = legend('Gaussian','VERSE');lgd.FontSize = lgdFontSize;
dp = -0.5:0.005:0.5;
[mx,my,mz] = bloch_CTR_Hz(200*bv,((circshift(gv,0))/10)*42.57E3,10e-6,inf,inf,0,dp,0);

[mx1,my1,mz1] = bloch_CTR_Hz(200*b,g*42.57E3/10,10e-6,inf,inf,0,dp,0);
subplot(1,3,3)
plot(dp,mz1,'b')
hold on
plot(dp,mz,'r--')
title('Slice profile','FontSize',FontSize)
xlabel('Position (cm)','FontSize',FontSize,'FontWeight','bold')
ylabel('Magnetisation','FontSize',FontSize,'FontWeight','bold')
lgd=legend('Gaussian','VERSE');lgd.FontSize = lgdFontSize;
%%
dp = -5:0.01:5;
TimeArray = 0:1e-6:((numel(b)-1)*10e-6);
%FreqShift = (2*6*42.57E3/10);
FreqShift = (2*6*42.57E3/10);
PhaseArray = TimeArray*FreqShift;
New_b = interp1(0:10e-6:((numel(b)-1)*10e-6),b,0:1e-6:((numel(b)-1)*10e-6));
New_g = interp1(0:10e-6:((numel(b)-1)*10e-6),g,0:1e-6:((numel(b)-1)*10e-6));
[mx2,my2,mz2] = bloch_CTR_Hz(200*New_b.*exp(-1i*PhaseArray*pi*2),New_g*42.57E3/10,1e-6,3,0.1,0,dp,0);
figure(600)
clf
plot(dp,mz2)

%%
dp = -5:0.01:5;
TimeArray = 0:1e-6:((numel(bv)-1)*10e-6);
%FreqShift = (2*6*42.57E3/10);
FreqShift = (0.2*6*42.57E3/10);
PhaseArray = TimeArray*FreqShift;
New_bv = interp1(0:10e-6:((numel(bv)-1)*10e-6),bv,0:1e-6:((numel(bv)-1)*10e-6));
New_gv = interp1(0:10e-6:((numel(bv)-1)*10e-6),gv,0:1e-6:((numel(bv)-1)*10e-6));
[mx3,my3,mz3] = bloch_CTR_Hz(200*New_bv.*exp(-1i*PhaseArray*pi*2),New_gv*42.57E3/10,1e-6,3,0.1,0,dp,0);
figure(601)
clf
plot(dp,mz3)

%%
dp = -5:0.01:5;
gv_Tm = gv*1E-3;
gv_Hzperm = gv_Tm*42.57E6;
gvCUMSUM = cumsum(gv_Hzperm)*10E-6;
distance_m = 0.02;
phase_vector = -gvCUMSUM*distance_m*pi*2;


[mx4,my4,mz4] = bloch_CTR_Hz(200*bv.*exp(-1i*phase_vector),gv_Hzperm/100,10e-6,inf,inf,0,dp,0);
figure(60)
clf
subplot(1,3,1)
yyaxis left
plot(10*(1:numel(bv)),200*bv); title('RF waveform');ylabel('Amplitude (Hz)');xlabel('time (us)')
yyaxis right
plot(10*(1:numel(bv)),phase_vector);ylabel('Phase (radians)');xlabel('Time (us)')

subplot(1,3,2)
plot(10*(1:numel(bv)),gv,'m');title('Gradient waveform');ylabel('mT/m');xlabel('Time (us)')
subplot(1,3,3)

plot(dp,mz4,'b')
title('Slice profile')
xlabel('Position (cm)')
ylabel('Magnetisation')




%% Making an rf array with finer sampling
% tyMakeGaussian temporarily changed
% 2us raster time. 300 points/600 us for rf/gradient flat top
% 15 points/30 us for ramp up/down
RFStruct_Fine = tyMakeGaussian(600,5);
b_fine = RFStruct_Fine.RF_pulse(2:end);
g_fine = ones(size(b_fine))*6;
StepSize = 6/15;
g_fine = [0:StepSize:StepSize*14 g_fine StepSize*14:-StepSize:0];
b_fine = [zeros(1,15) b_fine zeros(1,15)];
%b_fine = circshift(b_fine,5);
phase_vector_1mm_gaussian_fine = -0.001*cumsum(g_fine*42577)*2E-6*2*pi;
g_coarse = [2 4 6 6*ones(1,60) 4 2 0];
dt_fine.rf = 2E-6; dt_fine.grad = 10E-6;
MakeDotH_VERSE(b_fine,g_coarse,dt_fine,phase_vector_1mm_gaussian_fine);
%% Making a gaussian into a .h file
b2write = [0 0 RFStruct_Example.RF_pulse 0 0 0];
g2write = [0 2 4 6*ones(1,numel(RFStruct_Example.RF_pulse)-1) 4 2 0];
phase_vector_1mm_gaussian = -0.001*cumsum(g2write*42577)*10E-6*2*pi;
MakeDotH_VERSE(b2write,g2write,dt*1E-3,phase_vector_1mm_gaussian);