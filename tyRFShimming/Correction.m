filed_strength_comp
%%
plot_signal_T1
%%
function filed_strength_comp()
% setting parameters
TE = 13;            %ms
PLD = 1800;         %ms recommended value from ASL white paper
%Alsop et al. 2015 (10.1002/mrm.25197)
T2_star = [55, 32]; %ms 3T/7T
%van der Zwaag et al. 2009 (10.1016/j.neuroimage.2009.05.015)
T1 = [1649, 2087];  %ms 3T/7T
%Zhang et al. 2013 (10.1002/mrm.24550)

% relaxation
T2s_decay = @(t,T2) exp(-t/T2);
T1_decay = @(t,T1) 1-2*exp(-t/T1);vg
T2s_values = [T2s_decay(TE,T2_star(1)) T2s_decay(TE,T2_star(2))];
T2s_ratio = T2s_values(2)/T2s_values(1);
T1_values = [T1_decay(PLD,T1(1)) T2s_decay(PLD,T1(2))];
T1_ratio = T1_values(2)/T1_values(1);

% B1
Imaging_B1_ratio = 1-0.11;          % in the grey matter
% see Tong et al. (10.1002/mrm.28173)

% labelling efficiency
FA_CP = [7.679 8.671 7.183 7.566 6.159];
FA_verse = [14.567 15.477 15.146 12.426 13.323];             
% mean FA in 5 subj.
load('/Users/ytong/Documents/MATLAB/tong-acptx/tySpiral/tyRFShimming/GaussEfficiency.mat')
% "Efficiency" contains simulated efficieny for different FAs
% See
% /Users/ytong/Documents/MATLAB/tong-acptx/tySpiral/tyRFShimming/LabelEff
% for details
FAToSim = 0:0.05:25;
CalcEff = @(x) mean(interp1(FAToSim,Efficiency,x,'spline'));
Eff_20deg = Efficiency(find(0:0.05:25==20));
Eff_ratio_CP = CalcEff(FA_CP)/Eff_20deg;
Eff_ratio_verse = CalcEff(FA_verse)/Eff_20deg;

% instrinsic SNR assuming SNR proportional to B0
SNR_ratio = 7/3;

ratio_array = [Eff_ratio_CP,T2s_ratio,Imaging_B1_ratio,T1_ratio,SNR_ratio];
X = categorical({'Labelling B_1^+','T_2^*','Imaging B_1^+','T_1','Intrinsic SNR'});
X = reordercats(X,{'Labelling B_1^+','T_2^*','Imaging B_1^+','T_1','Intrinsic SNR'});
bar(X,[Eff_ratio_verse 0 0 0 0],0.5,...
    'FaceColor',[0.8500 0.3250 0.0980])
hold on
b = bar(X,ratio_array,0.5,'FaceColor',[0.9290 0.6940 0.1250]);
b.FaceColor = 'flat';
b.CData(4,:) = [0.4660 0.6740 0.1880];
b.CData(5,:) = [0.4660 0.6740 0.1880];
yline(1,'--k','LineWidth',3);
ylabel('SNR ratio')
title('7 T vs. 3 T')

gain = [prod(ratio_array), prod([Eff_ratio_verse ratio_array(2:end)])];
fprintf('CP gain: %0.2f\nVERSE gain: %0.2f\n',gain(1),gain(2));
box off
set(gca,'FontSize',20)
    plot_dim = [22 20];
    set(gcf,'color','w','InvertHardcopy','off')
    set(gcf,'units','centimeters','position',[4 4 plot_dim(1) plot_dim(2)],...
        'paperunits','centimeters','paperposition',[4 4 plot_dim(1) plot_dim(2)])
end

function plot_signal_T1
% a function that plots the signal for tag/control/diff images in ASL
% T1 values form Zhang et al. 2013 (10.1002/mrm.24550)

    T1 = [1649, 2087];
    T1_decay = @(t,T1,Mz) 1-(1-Mz)*exp(-t/T1);
    time_array = 0:10:2000;
    tag_curve = [T1_decay(time_array,T1(1),-1);T1_decay(time_array,T1(2),-1)];
    control_curve = ones(size(tag_curve));
    diff_curve = control_curve-tag_curve;
    color_cell = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980]};
    for iDx = 1:2
        subplot(2,1,iDx)
        if iDx ==1
            % 3T tag
            plot(time_array,tag_curve(1,:),'-','MarkerSize',15,'LineWidth',3,'Color',color_cell{1})
            hold on
            % 7T tag
            plot(time_array,tag_curve(2,:),'-','MarkerSize',15,'LineWidth',3,'Color',color_cell{2})
            % 3T control
            plot(time_array,control_curve(1,:),'--','MarkerSize',15,'LineWidth',3,'Color',color_cell{1})
            % 7T control
            plot(time_array,control_curve(2,:),'--','MarkerSize',15,'LineWidth',3,'Color',color_cell{2})
            legend('3 T label','7 T label','3 T control','7 T control','Location','southeast')
            title('Label & control signal')
        elseif iDx ==2
            plot(time_array,diff_curve(1,:),'-','MarkerSize',15,'LineWidth',3,'Color',color_cell{1})
            hold on
            plot(time_array,diff_curve(2,:),'-','MarkerSize',15,'LineWidth',3,'Color',color_cell{2})
            legend('3 T','7 T','Location','northeast')
            title('Difference signal')            
        end
        set(gca,'FontSize',20)
        xlabel('Time (ms)'); ylabel('Signal (a.u.)')
    end
    
    plot_dim = [25 22];
    set(gcf,'color','w','InvertHardcopy','off')
    set(gcf,'units','centimeters','position',[4 4 plot_dim(1) plot_dim(2)],...
        'paperunits','centimeters','paperposition',[4 4 plot_dim(1) plot_dim(2)])
end