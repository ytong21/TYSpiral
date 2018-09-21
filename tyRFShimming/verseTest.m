RFStruct_Example = tyMakeGaussian(600,5);
b = [0 0 RFStruct_Example.RF_pulse 0 0];  % arbitrary unit, peak=1
g = [0 0 6*ones(size(RFStruct_Example.RF_pulse)) 0 0];            % in mT/m
dt = 10e-3;                     % in ms
bmax = 0.5;
gmax = 22;                      % in mT/m
smax = 200;                     % in mT(m*ms)
dtout = 10e-3;                  % in ms
emax = sum(b.^2)*dt*0.5;        % in units of b*b*dt


[bv,gv] = mintverse(b,g,dt,bmax,gmax,smax,dtout,emax);




%% Temporary values. Need a bit more work.
numCh = 8;
peakV = 60;
numsamples = numel(bv);
targetFlipAngle = 20;

%% write into a .h file

MakeDotH_VERSE(bv,gv,dt*1E-3);
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
figure(60)
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
figure(60)
clf
plot(dp,mz3)

%%
dp = -5:0.01:5;
gv_Tm = gv*1E-3;
gv_Hzperm = gv_Tm*42.57E6;
gvCUMSUM = cumsum(gv_Hzperm)*10E-6;
distance_m = -0.02;
phase_vector = gvCUMSUM*distance_m;

[mx4,my4,mz4] = bloch_CTR_Hz(200*bv.*exp(-1i*phase_vector*pi*2),gv_Hzperm/100,10e-6,inf,inf,0,dp,0);
figure(60)
clf
plot(dp,mz4)