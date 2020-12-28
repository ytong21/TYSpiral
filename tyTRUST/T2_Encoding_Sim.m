% a script to simulate different T2 encoding schemes


% Scheme 1: simple hard pulses 
% (90x)-(180y)-(180y)-(180-y)-(180-y)-(90-x)
% time resolution 100us
CO_RF_duration = 600e-6;     %us
time_step = 10e-6;    %us

CO_RF_amp = (pi/2)/(2*pi*CO_RF_duration);

CO_unit_hard_pulse = create_ones(CO_RF_duration,time_step);
CO_gap_90_180 = create_zeros((2500-600)*1e-6,time_step);
CO_gap_180_180 = create_zeros((5e3-600)*1e-6,time_step);

complex_plus_x = complex(1,0);
complex_minus_x = complex(-1,0);
complex_plus_y = complex(0,1);
complex_minus_y = complex(0,-1);

%CO_MLEV_single = [CO_gap_90_180;2*CO_unit_hard_pulse*complex_plus_y;...     % 180y
CO_MLEV_single = [2*CO_unit_hard_pulse*complex_plus_y;...     % 180y
    CO_gap_180_180;2*CO_unit_hard_pulse*complex_plus_y;...    % 180y
    CO_gap_180_180;2*CO_unit_hard_pulse*complex_minus_y;...   % 180-y
    CO_gap_180_180;2*CO_unit_hard_pulse*complex_minus_y;...  % 180-y
    CO_gap_180_180];
    eTE = 40;
MLEV_train = make_MLEV_train(eTE,CO_MLEV_single,CO_unit_hard_pulse,CO_gap_90_180,CO_gap_180_180);
% CO_RF_array = [CO_unit_hard_pulse*complex_plus_x;...       % 90x
%     CO_gap_90_180;2*CO_unit_hard_pulse*complex_plus_y;...     % 180y
%     CO_gap_180_180;2*CO_unit_hard_pulse*complex_plus_y;...    % 180y
%     CO_gap_180_180;2*CO_unit_hard_pulse*complex_minus_y;...   % 180-y
%     CO_gap_180_180;2*CO_unit_hard_pulse*complex_minus_y;...   % 180-y
%     CO_gap_90_180;CO_unit_hard_pulse*complex_minus_x];        % 90-x
%%
make_seq_diag(CO_RF_array,21)
%%
make_seq_diag(MLEV_train,42)
%% Bloch Sim
B1_variation = 0.2:0.05:1.8;
%%
%CO_mag = BlochSimT2Prep(MLEV_train,B1_variation*CO_RF_amp);
%%
% Scheme 2: composite hard pulses 
% (90x)-(90x-180y-90x)-(90x-180y-90x)
% (90-x-180-y-90-x)-(90-x-180-y-90-x)-(90-x)
% time resolution 100us
Lu_RF_duration_90 = 200e-6;
Lu_RF_duration_180 = 400e-6;
Lu_unit_hard_pulse_90 = create_ones(Lu_RF_duration_90,time_step);
Lu_unit_hard_pulse_180 = create_ones(Lu_RF_duration_180,time_step);

Lu_composite_180_plus = ...
    [Lu_unit_hard_pulse_90*complex_plus_x;...    % 90x
    Lu_unit_hard_pulse_180*complex_plus_y;...    % 180y
    Lu_unit_hard_pulse_90*complex_plus_x];       % 90x

Lu_composite_180_minus = ...
    [Lu_unit_hard_pulse_90*complex_minus_x;...    % 90-x
    Lu_unit_hard_pulse_180*complex_minus_y;...    % 180-y
    Lu_unit_hard_pulse_90*complex_minus_x];       % 90-x
Lu_gap_90_180 = create_zeros((2500-100-400)*1e-6,time_step);
Lu_gap_180_180 = create_zeros((5000-800)*1e-6,time_step);

Lu_RF_array = [0;Lu_unit_hard_pulse_90*complex_plus_x;...   %90x
    Lu_gap_90_180;Lu_composite_180_plus;...                 %90x-180y-90x
    Lu_gap_180_180;Lu_composite_180_plus;...                %90x-180y-90x
    Lu_gap_180_180;Lu_composite_180_minus;...               %90-x-180-y-90-x
    Lu_gap_180_180;Lu_composite_180_minus;...               %90-x-180-y-90-x
    Lu_gap_90_180;Lu_unit_hard_pulse_90*complex_minus_x;0]; %90-x

Lu_RF_amp = (pi/2)/(2*pi*Lu_RF_duration_90);

%%
make_seq_diag(Lu_RF_array,22)
%%
%Lu_mag = BlochSimT2Prep(Lu_RF_array,B1_variation*Lu_RF_amp);
OffRess = -300:25:300;
%Lu_mag = BlochSimT2Prep_OffRes(Lu_RF_array,B1_variation*Lu_RF_amp,OffRess);
CO_mag = BlochSimT2Prep_OffRes(CO_RF_array,B1_variation*Lu_RF_amp,OffRess);
%%
Comp180_mag = BlochSimT2Prep(Lu_composite_180_minus,B1_variation*Lu_RF_amp);
%%
figure(44)
    set(gcf,'color','w','InvertHardcopy','off')
    set(gcf,'units','centimeters','position',[4 4 40 20],'paperunits','centimeters','paperposition',[4 4 40 20])
plot(B1_variation,CO_mag.mz,'LineWidth',1.5)
hold on
plot(B1_variation,Lu_mag.mz,'LineWidth',1.5)
legend('O Brien','Lu')
ylabel('Mz')
xlabel('B1 inhomogeneity')
set(gca,'FontSize',14)
%%
Mag_Scatter(CO_mag)
%%
Mag_Scatter(Lu_mag)
%%
PlotMzSingle(CO_mag,B1_variation,OffRess)
%%
PlotMz(CO_mag,Lu_mag,B1_variation,OffRess,b0,b1CP)
 %%   
function vec = create_ones(duration, time_step)
    vec = ones(round(duration/time_step),1);
end 
function vec = create_zeros(duration, time_step)
    vec = zeros(round(duration/time_step),1);
end 
%%
function make_seq_diag(CO_RF_array,time_limit)
    figure
    %figure(555)
    %CO_FA_array = 90*CO_RF_array;
    set(gcf,'color','w','InvertHardcopy','off')
    set(gcf,'units','centimeters','position',[4 4 40 20],'paperunits','centimeters','paperposition',[4 4 40 20])
    hold on
    CO_time_vec = (1:numel(CO_RF_array))/100;
%     index_plus_x = find((real(CO_FA_array)>0));
%     index_minus_x = find((real(CO_FA_array)<0));
%     index_plus_y = find(real(CO_FA_array)>0);
%     index_minus_y = find(real(CO_FA_array)<0);
    
%     full_plus_x = abs(CO_FA_array(index_plus_x)) + 
    %plot(CO_time_vec((real(CO_FA_array)>0)),abs(real(CO_FA_array((real(CO_FA_array)>0)))));
%     plot(CO_time_vec(index_plus_x),abs(real(CO_FA_array(index_plus_x))),...
%         CO_time_vec(index_minus_x),abs(real(CO_FA_array(index_minus_x))),...
%         CO_time_vec(index_plus_y),abs(imag(CO_FA_array(index_plus_y))),...
%         CO_time_vec(index_minus_y),abs(imag(CO_FA_array(index_minus_y))),...
%         'LineWidth',2);
%     plot(CO_time_vec,abs(real(CO_FA_array)).*(real(CO_FA_array)>0),...
%         CO_time_vec,abs(real(CO_FA_array)).*(real(CO_FA_array)<0),...
%         CO_time_vec,abs(imag(CO_FA_array)).*(imag(CO_FA_array)>0),...
%         CO_time_vec,abs(imag(CO_FA_array)).*(imag(CO_FA_array)<0),...
%         'LineWidth',2);
    
    area(CO_time_vec,abs(real(CO_RF_array)).*(real(CO_RF_array)>0));
    area(CO_time_vec,abs(real(CO_RF_array)).*(real(CO_RF_array)<0));
    area(CO_time_vec,abs(imag(CO_RF_array)).*(imag(CO_RF_array)>0));
    area(CO_time_vec,abs(imag(CO_RF_array)).*(imag(CO_RF_array)<0));
%         CO_time_vec,abs(real(CO_FA_array)).*(real(CO_FA_array)<0),...
%         CO_time_vec,abs(imag(CO_FA_array)).*(imag(CO_FA_array)>0),...
%         CO_time_vec,abs(imag(CO_FA_array)).*(imag(CO_FA_array)<0));    
    
     %plot(1:time_limit,zeros(time_limit,1),'k','LineWidth',1.5);
    %plot((1:numel(CO_RF_array))/100,zeros(size(CO_RF_array)),'k','LineWidth',2);
    %plot((1:numel(CO_RF_array))/100,abs(CO_RF_array));
    xlim([0 time_limit])
    %ylim([0 2.2])
    %yticks(0:45:90);
    xlabel('Time (ms)');
    ylabel('RF Amplitude (a.u.)')
    legend('+x','-x','+y','-y');
    %ldg.Title.FontSize = 20;
    set(gca,'FontSize',20,'linewidth',1.5,'ytick',[]);
end
%%
function OutStruct = BlochSimT2Prep(RF_array,RF_amp_array_Hz)
mode = 2;
dt = 10e-6;
g = zeros(numel(RF_array),3);
sens = complex(RF_amp_array_Hz);
df = zeros(size(sens));         % Assume there is no off-resonannce at the moment
dp = zeros(numel(sens),3);         % Assume isocentre
t1 = inf;
t2 = 150e-3;
[mx,my,mz] =  blochSim_mex_relax_multithread(RF_array,sens,g,dt,t1,t2,df,dp,mode);
OutStruct = struct('mx',mx,'my',my,'mz',mz);
end
%%
function MLEV_train = make_MLEV_train(eTE,MLEV_single,hard_pulse,gap1,gap2)
remainder = mod(eTE,20);
if eTE < 20 || remainder~=0
    error('eTE has to be a multiple of 20 ms.')
end
number_of_MLEVs = eTE/20;
all_180s = repmat(MLEV_single,number_of_MLEVs,1);
RF_90x_plus = hard_pulse*complex(1,0);
RF_90x_minus = hard_pulse*complex(-1,0);
% gap2 is the gap between 2 180 pulses.
gap_180s_length = numel(gap2);
all_180s(end+1-gap_180s_length:end) =[];
all_180s = [gap1;all_180s;gap1];
MLEV_train = [RF_90x_plus;all_180s;RF_90x_minus];
end
%%
function OutStruct = BlochSimT2Prep_OffRes(RF_array,RF_amp_array_Hz,OffResArray)
mode = 0;
dt = 10e-6;
g = zeros(numel(RF_array),3);
sens = complex(RF_amp_array_Hz);
%df = zeros(size(sens));         % Assume there is no off-resonannce at the moment
dp = zeros(numel(sens),3);         % Assume isocentre
t1 = 2587e-3;
t2 = 150e-3;
mx = zeros(numel(OffResArray),numel(sens));
my = mx;    mz = mx;
for iDx = 1:numel(OffResArray)
    df = ones(size(sens))*OffResArray(iDx);
    [mx(iDx,:),my(iDx,:),mz(iDx,:)] =  blochSim_mex_relax_multithread(RF_array,sens,g,dt,t1,t2,df,dp,mode);
end
    OutStruct = struct('mx',mx,'my',my,'mz',mz);
end

function OutStruct = Bloch_CTR_Hz_T2Prep(RF_array,RF_amp_array_Hz)
mode = 2;
dt = 10e-6;
g = zeros(numel(RF_array),3);
sens = complex(RF_amp_array_Hz);
df = zeros(size(sens));         % Assume there is no off-resonannce at the moment
dp = zeros(numel(sens),3);         % Assume isocentre
t1 = inf;
t2 = 150e-3;
[mx,my,mz] =  bloch_CTR_Hz(RF_array,sens,g,dt,t1,t2,df,dp,mode);
OutStruct = struct('mx',mx,'my',my,'mz',mz);
end
%%
function PlotMzSingle(CO_mag,B1_variation,OffRess)
figure(557)
    %CO_FA_array = 90*CO_RF_array;
    set(gcf,'color','w','InvertHardcopy','off')
    set(gcf,'units','centimeters','position',[4 4 30 20],'paperunits','centimeters','paperposition',[4 4 30 20])
    level = 0.875173319042947;      %ideal Mz after T2 prep
    %plot_Range = [0.4 max([CO_mag.mz(:);Lu_mag.mz(:)])];
    plot_Range = [0 1];
    %OffRess_flipped = flip(OffRess);
    CO = (CO_mag.mz);
    %Lu = (Lu_mag.mz);
    x = [B1_variation(1)  B1_variation(end)];
    y = [OffRess(1) OffRess(end)];
    %subplot(1,2,1)
    imagesc(x,y,CO,plot_Range)
    c = colorbar;
    c.Label.String = 'Mz';
%     line('parent',c,'xdata',[-5 5],...
%      'ydata',[level level],'color','k','LineWidth',3)
    %rectangle('Position',[0.8 -100 0.4 200],'LineWidth',3,'EdgeColor',[0.3010 0.7450 0.9330]);
%     rectangle('Position',[b1_pct(1) b0_pct(1) b1_pct(2)-b1_pct(1) b0_pct(2)-b0_pct(1)],...
%         'LineWidth',3,'EdgeColor','b');
    set(gca,'FontSize',16);
    colormap hot
    xlabel('B_{1} inhomogeneity')
    ylabel('Off-resonance (Hz)')
%     title('Simple hard pulse')
%     subplot(1,2,2)
%     imagesc(x,y,Lu,plot_Range)
%     c = colorbar;
%     c.Label.String = 'Mz';
% 
%     line('parent',hCbar,'xdata',[-5 5],...
%      'ydata',[level level],'color','k','LineWidth',3)
%     colormap hot
%     set(gca,'FontSize',16);
%     xlabel('B1 inhomogeneity')
%     ylabel('Off-resonance (Hz)')
%     title('Composite pulse')
%     rectangle('Position',[b1_pct(1) b0_pct(1) b1_pct(2)-b1_pct(1) b0_pct(2)-b0_pct(1)],...
%         'LineWidth',3,'EdgeColor','b');
end
%%
function PlotMz(CO_mag,Lu_mag,B1_variation,OffRess,b0,b1CP)
    figure(556)
    %CO_FA_array = 90*CO_RF_array;
    set(gcf,'color','w','InvertHardcopy','off')
    set(gcf,'units','centimeters','position',[4 4 45 18],'paperunits','centimeters','paperposition',[4 4 45 18])
    level = 0.875173319042947;      %ideal Mz after T2 prep
    b0_pct = zeros(2,1);
    b0_pct(1) = prctile(b0,5);
    b0_pct(2) = prctile(b0,95);
    b1_pct = zeros(2,1);
    b1_pct(1) = prctile(b1CP,5)/mean(b1CP);
    b1_pct(2) = prctile(b1CP,95)/mean(b1CP); 
    %plot_Range = [0.4 max([CO_mag.mz(:);Lu_mag.mz(:)])];
    plot_Range = [0 1];
    %OffRess_flipped = flip(OffRess);
    CO = (CO_mag.mz);
    Lu = (Lu_mag.mz);
    x = [B1_variation(1)  B1_variation(end)];
    y = [OffRess(1) OffRess(end)];
    subplot(1,2,1)
    imagesc(x,y,CO,plot_Range)
    c = colorbar;
    c.Label.String = 'Mz';
    line('parent',c,'xdata',[-5 5],...
     'ydata',[level level],'color','k','LineWidth',3)
    %rectangle('Position',[0.8 -100 0.4 200],'LineWidth',3,'EdgeColor',[0.3010 0.7450 0.9330]);
    rectangle('Position',[b1_pct(1) b0_pct(1) b1_pct(2)-b1_pct(1) b0_pct(2)-b0_pct(1)],...
        'LineWidth',3,'EdgeColor','b');
    set(gca,'FontSize',16);
    colormap hot
    xlabel('B_{1} inhomogeneity')
    ylabel('Off-resonance (Hz)')
    title('Simple hard pulse')
    subplot(1,2,2)
    imagesc(x,y,Lu,plot_Range)
    c = colorbar;
    c.Label.String = 'Mz';

    line('parent',hCbar,'xdata',[-5 5],...
     'ydata',[level level],'color','k','LineWidth',3)
    colormap hot
    set(gca,'FontSize',16);
    xlabel('B1 inhomogeneity')
    ylabel('Off-resonance (Hz)')
    title('Composite pulse')
    rectangle('Position',[b1_pct(1) b0_pct(1) b1_pct(2)-b1_pct(1) b0_pct(2)-b0_pct(1)],...
        'LineWidth',3,'EdgeColor','b');
end
function Mag_Scatter(mag)
figure

[x,y,z] = sphere ;
surf(x,y,z,'FaceColor','none');
num_pts = size(mag.my,1);
C = parula(num_pts);
hold on
mx = mag.mx(:,17);my = mag.my(:,17);mz = mag.mz(:,17);
xlabel('Mx');ylabel('My');zlabel('Mz');
scatter3(mx,my,mz,36,C,'filled');
colormap( jet(numel(mag.mx)) );
set(gcf,'color','w','InvertHardcopy','off');
set(gca,'FontSize',14);
% for iDx = 1:num_pts
%     scatter3(mx(iDx),my(iDx),mz(iDx),36,C(iDx),'filled')
%     pause(0.02)
% end
end