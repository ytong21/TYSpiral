%%
AA = make_mosaic(s012,[7,8]);
%%
PathName = '/Users/ytong/Documents/Data/TRUST/20200723_Phantom_23_07';
PlotRange = [3 1200];
%%
[~,AllSlices180V] = load_dcm(PathName,15:24);
AllSlices180V_chop = trim3D_vert(AllSlices180V,6);
AllSlices180V_mosaic = make_mosaic(AllSlices180V_chop,[2,5]);
plot_mosaic(AllSlices180V_mosaic,[35 20])

%%

[~,AllSlices223V] = load_dcm(PathName,29:38);
AllSlices223V_chop = trim3D_vert(AllSlices223V,6);
AllSlices223V_mosaic = make_mosaic(AllSlices223V_chop,[2,5]);
plot_mosaic(AllSlices223V_mosaic,[35 20],PlotRange)


%%
[~,Top_DiffV] = load_dcm(PathName,39:54);
[~,Mid_DiffV] = load_dcm(PathName,55:70);
[~,Low_DiffV] = load_dcm(PathName,72:87);
[~,Mid_Delay_223] = load_dcm(PathName,[89:98 100]);
[~,Mid_Delay_180] = load_dcm(PathName,101:111);
[~,Spoils_180] = load_dcm(PathName,113:121);
[~,Spoils_223] = load_dcm(PathName,122:130);
[~,FullSlab223V] = load_dcm(PathName,29:38);
[~,FullSlab180V] = load_dcm(PathName,15:24);
%%
plot_diffV(Top_DiffV,PlotRange)
%%
plot_delay(Mid_Delay_223,PlotRange)
%%
plot_spoil(Spoils_223,PlotRange)
%%
plot_full_slab(FullSlab180V,PlotRange)
%%
fraction = find_fraction(AA, seg_new);
%%
chart_spoil(Spoils_180, Spoils_223,seg_new)
%%
chart_diffV({Low_DiffV,Mid_DiffV,Top_DiffV},seg)
%%
BB = seg;
BB(DD>0)=0;
CC = BB.*AA;
figure
imagesc(CC)

%%
DD = zeros(size(seg));
DD(seg==3)=1;
DD = imgaussfilt(DD,0.5);
imagesc(DD)

%%
seg_new = seg;
seg_new(DD>0)=3;
imagesc(seg_new)
%%
out = find_mean_back(FullSlab180V,FullSlab223V,seg_new);
%%
fraction = zeros(1,size(Mid_DiffV,3));
for iDx = 1:size(Mid_DiffV,3)
    fraction(iDx) = find_fraction(Mid_DiffV(:,:,iDx), seg_new);
end
%%
sim_conta_volume(FullSlab180V,seg_new)
%%
    run_comparison_plot(CircleMaskFiltered,OutCell{1}.finalMag,maskedMaps);

    
    %%
    plot_diffV(abs(permute(squeeze(sim_DiffV(:,:,5,:)),[2 1 3])),[0 1])
    
    %%
    plot_diffV(abs(sim_DiffV(:,:,5,:)),[0 1])
    %%
    plot_delay(abs(delay_img(:,:,5,:)),[0 1])
    
    %%
    plot_sat_profile(Spoils_180)
    %%
    [Sz_0, s_perc] = calc_Sz_t0(100.1,1229,80,1370)
%%
plot_T1_field
    %%
    function plot_T1_field
    T1_3T = 1649;   %ms
    T1_7T = 2087;   %ms
    % this is for Hct level of 42%
    Hct = 0.42;  Y = 0.6;
    T2_3T = calc_T2_from_Hct_3T(Hct,Y);                     %From Lu 2012 10.1002/mrm.22970
    T2_7T = 1/(dot([14.6, -31.2, 223.5],[1,1-Y,(1-Y)^2]));  %From Krishnamurthy 2014 10.1002/mrm.24868
    T2_decay = @(t,T2) exp(-t/T2);
    eTE20_Mz = [T2_decay(20e-3,T2_3T), T2_decay(20e-3,T2_7T)];
    time_array = 0:10:1500;
    T1_decay = @(t,T1,Mz) 1-(1-Mz)*exp(-t/T1);
    tag_curve = [T1_decay(time_array,T1_3T,-1);T1_decay(time_array,T1_7T,-1);
        T1_decay(time_array,T1_3T,-eTE20_Mz(1));T1_decay(time_array,T1_7T,-eTE20_Mz(2))];
    control_curve = [ones(2,numel(time_array));...
        T1_decay(time_array,T1_3T,eTE20_Mz(1));T1_decay(time_array,T1_7T,eTE20_Mz(2))];
    color_cell = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980]};
    title_cell = {'eTE = 0','eTE = 20 ms'};
    for iDx = 1:2
        subplot(2,1,iDx)
            for jDx = 1:2
            plot(time_array,tag_curve((iDx-1)*2+jDx,:),'-d','MarkerSize',15,'LineWidth',4,...
                'MarkerIndices',121,'Color',color_cell{jDx})
            hold on
            plot(time_array,control_curve((iDx-1)*2+jDx,:),'--d','MarkerSize',15,'LineWidth',4,...
                'MarkerIndices',121,'Color',color_cell{jDx})    
            end
            % 121 is for a TI = 1200ms
            xlabel('Time (ms)');ylabel('M_z');
            ylim([-1 1]);title(title_cell{iDx})
            set(gca,'FontSize',20)
    end
    plot_dim = [30 30];
    set(gcf,'color','w','InvertHardcopy','off')
    set(gcf,'units','centimeters','position',[4 4 plot_dim(1) plot_dim(2)],...
        'paperunits','centimeters','paperposition',[4 4 plot_dim(1) plot_dim(2)])
    legend('3 T label','3 T control','7 T label','7 T control','Location','southeast')
    
    end
    %%
    function T2 = calc_T2_from_Hct_3T(Hct,Y)
    % calibration data from 10.1002/mrm.22970
    % This is for tau_cmpg = 5ms
    a_array = [-4.4 39.1 -33.5];
    b_array = [1.5 4.7];
    c_array = 167.8;
    A = dot(a_array,[1 Hct Hct^2]);
    B = dot(b_array,[Hct Hct^2]);
    C = c_array*Hct*(1-Hct);
    R2 = dot([A B C],[1 1-Y (1-Y)^2]);
    T2 = 1/R2;
    end
function plot_sat_profile(Img3D)
Img = Img3D(:,:,4:6);   %middle slice
% Y = 21;
%line_index = 15:37;
line_profile = squeeze(Img(:,25,:));
specification = {'-','--','-.'};
figure
hold on
for iDx = 1:size(line_profile,2)
    plot(line_profile(:,iDx),specification{iDx},'LineWidth',4);
    width(iDx) = fwhm(1:size(line_profile,1), line_profile(:,iDx));
end
xlim([1 size(line_profile,1)])
xtick off;
ylabel('Signal intensity (a.u.)')
%ytick off
plot_dim = [30 20];
set(gcf,'color','w','InvertHardcopy','off')
set(gcf,'units','centimeters','position',[4 4 plot_dim(1) plot_dim(2)],...
    'paperunits','centimeters','paperposition',[4 4 plot_dim(1) plot_dim(2)])
set(gca,'FontSize',30)
legend('1 sat','2 sat','3 sat')
disp(width)
end
function [Sz_0,percentage] = calc_Sz_t0(Sz_t,S0,t,T1)
    Sz_0 = (Sz_t-S0*(1-exp(-t/T1)))/exp(-t/T1);
    percentage = Sz_0/S0;
end
    function run_comparison_plot(circle,result,maskedMaps)
    % taken from main for phantom
        %figure(222)
        set(gcf,'color','w','InvertHardcopy','off')
        set(gcf,'units','centimeters','position',[4 4 35 20],'paperunits','centimeters','paperposition',[4 4 35 20])
        %slice_array = 35:2:44;        
        result_to_plot_array = 1:2:10;  toChop = 6;
        FAmap = rad2deg(abs(result(:,:,result_to_plot_array)));
        FAmap = trim3D_vert(FAmap,toChop);
        % calc difference
        target = squeeze(maskedMaps.mask_one_slice(:,:,result_to_plot_array)) - circle;
        diff = abs(rad2deg(result(:,:,result_to_plot_array)))-90*target;
        diff(diff == 0) = -inf;
        diff = trim3D_vert(diff,toChop);
        for iDx = 1:size(FAmap,3)
            FA_plot(:,:,iDx) = FAmap(:,:,iDx)';
            diff_plot(:,:,iDx) = diff(:,:,iDx)';
        end
        FA_mosaic = make_mosaic(FAmap,[1 5]);
        diff_mosaic = make_mosaic(diff,[1 5]);
        figure
        FigDim = [42 18];
        set(gcf,'color','w','InvertHardcopy','off')
        set(gcf,'units','centimeters','position',[4 4 FigDim(1) FigDim(2)],...
            'paperunits','centimeters','paperposition',[4 4 FigDim(1) FigDim(2)])
        subplot(2,1,1)
        plot_mosaic(FA_mosaic,[40 18],[0 100],false,true,'deg');
        write_text(1:5,1:5,result_to_plot_array,FAmap);
        title('Achieved FA','FontSize',20);
        subplot(2,1,2)
        plot_mosaic(diff_mosaic,[40 18],[-10 10],false,true,'deg')
        write_text(1:5,1:5,result_to_plot_array,FAmap);
        title('Difference','FontSize',20)
        %set(gca,'FontSize',24)
    end
    function write_text(index_array,index_matrix,text_array,Img)
    [d1,d2,~] = size(Img);
        for iDx = 1:numel(index_array)
            [row,col] = find(index_matrix == index_array(iDx));
            location_row = (row-1)*d1+3;
            location_col = (col-1)*d2+1;
            string = sprintf('Slice %d',text_array(iDx));
            text(location_col,location_row,string,'Color','w',...
                'FontSize',18)
        end
    end
function sim_conta_volume(Img,Seg)
seg_3D = repmat(Seg,[1 1 size(Img,3)]);
ROI = Img(seg_3D==3);
Back = Img(seg_3D==2);
mean_back = calc_Sz_t0(mean(Back),1229,80,1370);
mean_ROI = mean(ROI);
ratio_array = 0:40;
contamination = (mean_back*(1-ratio_array/100))./...
    (mean_back*(1-ratio_array/100)+mean_ROI*(ratio_array/100));
figure(131)
plot(ratio_array,contamination*100,'LineWidth',5)
xlabel('ROI volume (%)');ylabel('Contamination (%)')
            set(gcf,'color','w','InvertHardcopy','off')
            set(gcf,'units','centimeters','position',[4 4 30 20],...
                'paperunits','centimeters','paperposition',[4 4 30 20])
            set(gca,'FontSize',24)
end
function Yv_sim
%Y= 40%, 65%
% See Lu et al. 2008 MRM equations 4&5.
T1b = 2100; %ms
T2b = [calc_T2_fromY(0.4,2),calc_T2_fromY(0.65,2)];  %ms
C = 1/(T1b*1e-3)-1./(T2b*1e-3);
eTE = 0:160;    decay_curve = zeros(2,numel(eTE));

for iDx = 1:numel(C)
    decay_curve(iDx,:) = exp(eTE/C(iDx));
end
bool_plot = [true,true,false, false];
if bool_plot(1) == true
        figure(121)
        plot(eTE,decay_curve','LineWidth',5);
        xlabel('eTE (ms)');ylabel('MR signal (a.u.)')
            set(gcf,'color','w','InvertHardcopy','off')
            set(gcf,'units','centimeters','position',[4 4 35 25],...
                'paperunits','centimeters','paperposition',[4 4 35 25])
            set(gca,'FontSize',24)
        legend('Y = 40%','Y = 65%')    
end
        c_contami = 0:0.05:0.4;
        decay_bi_exp = zeros(numel(c_contami),numel(eTE));
        for iDx = 1:numel(c_contami)
            decay_bi_exp(iDx,:) = c_contami(iDx)*decay_curve(2,:)+(1-c_contami(iDx))*decay_curve(1,:);
            legend_cell{iDx} = strcat('Contamination=',num2str(c_contami(iDx)*100),'%');
        end
        
if bool_plot(2) == true
        figure(122)
        plot(eTE,decay_bi_exp','LineWidth',3);
        xlabel('eTE (ms)');ylabel('MR signal (a.u.)')
            set(gcf,'color','w','InvertHardcopy','off')
            set(gcf,'units','centimeters','position',[4 4 35 25],...
                'paperunits','centimeters','paperposition',[4 4 35 25])
            set(gca,'FontSize',24)
            legend(legend_cell)
end
% Expotential fit, perfect
fitted = cell(size(decay_bi_exp,1),1);
T2_fitted = zeros(size(decay_bi_exp,1),1);
Y_fitted = zeros(size(decay_bi_exp,1),1);
% Fit with 0, 20, 40 ms for eTE with no noise
eTE_new = eTE(1:20:41); decay_bi_exp_new = decay_bi_exp(:,1:20:41);
fitted_wo_noise = cell(size(decay_bi_exp,1),1);
T2_fitted_wo_noise = zeros(size(decay_bi_exp,1),1);
% Fit with 0, 20, 40 ms for eTE with noise
fitted_w_noise = cell(size(decay_bi_exp,1),1);
T2_fitted_w_noise = zeros(size(decay_bi_exp,1),1);
syms x
for iDx = 1:size(decay_bi_exp,1)
    fitted{iDx} = fit(eTE',decay_bi_exp(iDx,:)','exp1');
    T2_fitted(iDx) = 1000*1/(1/(T1b*1e-3)-1/fitted{iDx}.b);
    % Find corresponding Y
%     eqn = 264*x^2 - 17.6*x + 14.6 - 1/(1e-3*T2_fitted(iDx)) == 0;
%     solx = solve(eqn, x);
    Y_fitted(iDx) = calc_Y_fromT2(T2_fitted(iDx),2);
    % Fit with 0, 20, 40 ms for eTE with noise
    fitted_wo_noise{iDx} = fit(eTE_new',decay_bi_exp_new(iDx,:)','exp1');
    T2_fitted_wo_noise(iDx) = 1000*1/(1/(T1b*1e-3)-1/fitted_wo_noise{iDx}.b);
    % add some noise
    exp_noise = decay_bi_exp_new(iDx,:)'+0.02*randn(size(decay_bi_exp_new(iDx,:)'));
    fitted_w_noise{iDx} = fit(eTE_new',exp_noise,'exp1');
    T2_fitted_w_noise(iDx) = 1000*1/(1/(T1b*1e-3)-1/fitted_w_noise{iDx}.b);
end
if bool_plot(3) == true
        figure(123)
        specification = {'-','--','-.';'-o','--d','-.^'};
        fit_2_plot = [T2_fitted,T2_fitted_wo_noise,T2_fitted_w_noise];
        hold on
        for iDx = 1:size(fit_2_plot,2)
            plot(c_contami*100,fit_2_plot(:,iDx)',specification{2,iDx},...
                'LineWidth',5,'MarkerSize',20);
        end
        xlabel('Contamination (%)');ylabel('T_2 (ms)')
            set(gcf,'color','w','InvertHardcopy','off')
            set(gcf,'units','centimeters','position',[4 4 35 25],...
                'paperunits','centimeters','paperposition',[4 4 35 25])
            set(gca,'FontSize',24)
            legend('Perfect fit','3 eTEs wo/ noise','3 eTEs w/ noise',...
                'Location','NorthWest')
end
if bool_plot(4) == true
        figure(124)
        hold on
        yyaxis left
        plot(c_contami*100,T2_fitted','-o',...
                'LineWidth',5,'MarkerSize',20);
        xlabel('Contamination (%)');ylabel('T_2 (ms)')
        
        yyaxis right
        plot(c_contami*100,100*Y_fitted','--d',...
                'LineWidth',5,'MarkerSize',20);
        xlabel('Contamination (%)');ylabel('Y (%)')        
            set(gcf,'color','w','InvertHardcopy','off')
            set(gcf,'units','centimeters','position',[4 4 35 20],...
                'paperunits','centimeters','paperposition',[4 4 35 20])
            set(gca,'FontSize',24)
             legend('T_2 vs. contamination','Y vs. contamination',...
                 'Location','NorthWest')
end
end
function Y = calc_Y_fromT2(T2,HctMode)
if sum(HctMode==[1,2,3])==0
    error('HctMode must be 1, 2, or 3.')
end
% 1 2 3 corresponds to 34%, 42%, 54%
Hct_ABC = [14.6,    -31.2,  223.5;
           14.9,    -17.6,  264.0;
           16.7,     3.7,   240.9];
syms x
eqn = Hct_ABC(HctMode,3)*x^2 + Hct_ABC(HctMode,2)*x + Hct_ABC(HctMode,1)...
        - 1/(1e-3*T2) == 0;
    solx = solve(eqn, x);
    Y = 1-double(solx(2));
end
function T2 = calc_T2_fromY(Y,HctMode)
if sum(HctMode==[1,2,3])==0
    error('HctMode must be 1, 2, or 3.')
end
% 1 2 3 corresponds to 34%, 42%, 54%
Hct_ABC = [14.6,    -31.2,  223.5;
           14.9,    -17.6,  264.0;
           16.7,     3.7,   240.9];
T2 =  1000/(Hct_ABC(HctMode,:)*[1; 1-Y; (1-Y)^2]);
end
function sim_T2_Y
%Refer to Krishnamurthy et al. MRM 2014
%10.1002/mrm.24868 table 1
Hct_ABC = [14.6,    -31.2,  223.5;
           14.9,    -17.6,  264.0;
           16.7,     3.7,   240.9];
c1 = 0:0.05:0.4;
c2 = 1-c1;
T2_concentration = zeros(3,numel(c1));
Y1 = 0.65;   Y2 = 0.4;  Y = 0.3:0.01:1;
T2_Y = zeros(3,numel(Y));
specification = {'-','--','-.';'-o','--d','-.x'};
lgd = {'Hct=34%','Hct=42%','Hct=54%'};
for jDx = 1:3
    for iDx = 1:numel(c1)
        Y_temp = c1(iDx)*Y1+c2(iDx)*Y2;
        T2_concentration(jDx,iDx) = 1000/(Hct_ABC(jDx,:)*[1; 1-Y_temp; (1-Y_temp)^2]);
        % convert to ms.
    end
    for iDx = 1:numel(Y)
        T2_Y(jDx,iDx) = 1000/(Hct_ABC(jDx,:)*[1; 1-Y(iDx); (1-Y(iDx))^2]);
        % convert to ms.
    end
end
figure(111)
hold on
for jDx = 1:3
    plot(100*Y,T2_Y(jDx,:),specification{1,jDx},'LineWidth',5);
end
    set(gcf,'color','w','InvertHardcopy','off')
    set(gcf,'units','centimeters','position',[4 4 35 25],...
        'paperunits','centimeters','paperposition',[4 4 35 25])
    set(gca,'FontSize',24)
    xlabel('Oxygenation Y(%)');ylabel('T_2 (ms)');
    legend(lgd,'Location','NorthWest');
figure(112)
hold on
for jDx = 1:3
    plot(100*c1,T2_concentration(jDx,:),specification{2,jDx},...
        'LineWidth',5,'MarkerSize',15);
end
    set(gcf,'color','w','InvertHardcopy','off')
    set(gcf,'units','centimeters','position',[4 4 35 25],...
        'paperunits','centimeters','paperposition',[4 4 35 25])
    set(gca,'FontSize',24)
    xlabel('Singal contamination (%)');ylabel('T_2 (ms)');
    legend(lgd,'Location','NorthWest');
end
function out = find_mean_back(FullSlab180V,FullSlab223V,seg_new)
%     voxels_180 = [];
%     voxels_223 = [];
    mask = repmat(seg_new,[1,1,10]);
    back_180 = FullSlab180V(mask==2);
    back_223 = FullSlab223V(mask==2);
    ROI_180 = FullSlab180V(mask==3);
    ROI_223 = FullSlab223V(mask==3); 
    [Sz_0_180, ~] = calc_Sz_t0(mean(back_180),1229,80,1370);
    [Sz_0_223, ~] = calc_Sz_t0(mean(back_223),1229,80,1370);
    contam(1) = Sz_0_180*numel(back_180)/(Sz_0_180*numel(back_180)+...
        sum(ROI_180));
    contam(2) = Sz_0_223*numel(back_223)/(Sz_0_223*numel(back_223)+...
        sum(ROI_223));
    out = {[back_180,back_223],[ROI_180,ROI_223],contam};
end
function chart_diffV(Img_Cell,seg)

    diffV_metric = zeros(numel(Img_Cell),size(Img_Cell{1},3));
    % dim1: (low, mid, high) dim2: voltage dim3: mean/std
    for iDx = 1:numel(Img_Cell)
        for jDx = 1:size(Img_Cell{1},3)
            [~,mean_std_temp,~]=find_fraction(Img_Cell{iDx}(:,:,jDx),seg);
            diffV_metric(iDx,jDx,1) = mean_std_temp(1);
            diffV_metric(iDx,jDx,2) = mean_std_temp(2);
        end
    end
    voltage = 100:10:250;
    color_cell = {[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4660 0.6740 0.1880]};
    marker_cell = {'-o','-x','-s'};
    hold on
    for iDx = 1:numel(Img_Cell)
        errorbar(voltage,diffV_metric(iDx,:,1),diffV_metric(iDx,:,2),...
            marker_cell{iDx},'LineWidth',3,'MarkerSize',17,'Color',color_cell{iDx},'CapSize',15)
    end
    legend('Low','Mid','High')
    xlabel('Pulse voltage');ylabel('Backgroud pixel intensity (a.u.)')
    set(gcf,'color','w','InvertHardcopy','off')
    set(gcf,'units','centimeters','position',[4 4 35 25],...
        'paperunits','centimeters','paperposition',[4 4 35 25])
    set(gca,'FontSize',22)
end
function chart_spoil(Spoils_180, Spoils_223,seg)
    plot_mode = 'ROI';
    mean_back_180 = zeros(1,size(Spoils_180,3));
    mean_back_223 = zeros(1,size(Spoils_223,3));
    mean_ROI_180 = zeros(1,size(Spoils_180,3));
    mean_ROI_223 = zeros(1,size(Spoils_223,3));
    for iDx = 1:numel(mean_back_180)
        [~,back_temp,ROI_temp]=find_fraction(Spoils_180(:,:,iDx),seg);
        mean_back_180(iDx) = back_temp(1);
        mean_ROI_180(iDx) = ROI_temp(1);
        [~,back_temp,ROI_temp]=find_fraction(Spoils_223(:,:,iDx),seg);
        mean_back_223(iDx) = back_temp(1);
        mean_ROI_223(iDx) = ROI_temp(1);
    end
    if     strcmp(plot_mode, 'ROI')
        m_180 = reshape(mean_ROI_180,[3,3]);
        m_223 = reshape(mean_ROI_223,[3,3]);
    elseif strcmp(plot_mode, 'backgroud')
        m_180 = reshape(mean_back_180,[3,3]);
        m_223 = reshape(mean_back_223,[3,3]);
    end
    x_data = 1:3;       % 1,2,3 spoils
    color_cell = {[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4660 0.6740 0.1880]};
    hold on
    for iDx = 1:3
        plot(x_data,m_180(:,iDx),'-o','LineWidth',3,'MarkerSize',17,'Color',color_cell{iDx})
    end
    for iDx = 1:3
        plot(x_data,m_223(:,iDx),'--d','LineWidth',3,'MarkerSize',17,'Color',color_cell{iDx})
    end
    xlim([0.8 3.2]);xticks(1:3);%yticks(100:100:500)
    xlabel('No. of sat pulses');
    ylabel(sprintf('Mean %s pixel intensity (a.u.)',plot_mode))
    legend('180V low','180V mid','180V high','223V low','223V mid','223V high')
    set(gcf,'color','w','InvertHardcopy','off')
    set(gcf,'units','centimeters','position',[4 4 30 30],...
        'paperunits','centimeters','paperposition',[4 4 30 30])
    set(gca,'FontSize',22)
end
function plot_diffV(Img,PlotRange)
mode = 'sim';
if strcmp(mode,'experiment')
    Img_trim = trim3D_vert(Img,5);
elseif strcmp(mode,'sim')
    Img_trim = trim3D_hori(Img,5);
end
[d1,d2,~] = size(Img_trim);
Img_mosaic = make_mosaic(Img_trim,[4,4]);
plot_mosaic(Img_mosaic,[30 35],PlotRange,true,true,'')
% Trying to add some text to the mosaic images
% First try to find the correct locations to add texts
voltage_array = 100:10:250;
index_array = 1:numel(voltage_array);
index_matrix = reshape(index_array,[4,4])';
for iDx = 1:numel(voltage_array)
    [row,col] = find(index_matrix == index_array(iDx));
    location_row = (row-1)*d1+4;
    location_col = (col-1)*d2+1;
    string = sprintf('%d V',voltage_array(iDx));
    text(location_col,location_row,string,'Color','w',...
        'FontSize',16)
end
end
function plot_nii
    nii_slice_array = 32:2:40;
    for iDxx = 1:5
        FileName = strcat('/Users/ytong/Documents/Data/TRUST/2013_50_123_nii/',num2str(nii_slice_array(iDxx)),'_masked.nii.gz');
        nii_tmp = niftiread(FileName);
        nii_tmp = imrotate(nii_tmp,90);
        TRUST_nii(:,:,iDxx) = nii_tmp;
    end
    TRUST_nii(:,:,iDxx+1) = zeros(size(TRUST_nii,1),size(TRUST_nii,2));
    Img_trim_1 = trim3D_vert(TRUST_nii,8);
    Img_trim = trim3D_hori(Img_trim_1,8);
    Img_mosaic = make_mosaic(Img_trim,[2,3]);
    [d1,d2,~] = size(Img_trim);
    index_array = 1:6;
    index_matrix = reshape(index_array,[3,2])';
    plot_mosaic(Img_mosaic,[30 35])
    for iDx = 1:numel(index_array)-1
        [row,col] = find(index_matrix == index_array(iDx));
        location_row = (row-1)*d1+4;
        location_col = (col-1)*d2+1;
        string = sprintf('Slice %d',(iDx-1)*2+1);
        text(location_col,location_row,string,'Color','w',...
            'FontSize',16)
    end
end
function plot_delay(Img,PlotRange)
%Img_trim = trim3D_vert(Img,5);
mode = 'sim';
if strcmp(mode,'experiment')
    Img_trim = trim3D_vert(Img,5);
elseif strcmp(mode,'sim')
    Img_trim = trim3D_hori(Img,5);
end
[d1,d2,~] = size(Img_trim);
% add an extra empty image
Img_trim_new = cat(3,Img_trim,zeros(d1,d2));
Img_mosaic = make_mosaic(Img_trim_new,[3,4]);
plot_mosaic(Img_mosaic,[28 38],PlotRange,true,true,'')
% Trying to add some text to the mosaic images
% First try to find the correct locations to add texts
delay_array = [-10:2:10 nan];
index_array = 1:numel(delay_array);
index_matrix = reshape(index_array,[4,3])';
for iDx = 1:numel(delay_array)-1
    [row,col] = find(index_matrix == index_array(iDx));
    location_row = (row-1)*d1+4;
    location_col = (col-1)*d2+1;
    string = strcat(num2str(delay_array(iDx)),'\mus');
    text(location_col,location_row,string,'Color','w',...
        'FontSize',16)
end
end
function plot_full_slab(Img,PlotRange)
Img_trim = trim3D_vert(Img,5);
[d1,d2,~] = size(Img_trim);
% add an extra empty image
Img_trim_new = cat(3,Img_trim,zeros(d1,d2),zeros(d1,d2));
Img_mosaic = make_mosaic(Img_trim_new,[3,4]);
plot_mosaic(Img_mosaic,[28 38],PlotRange)
% Trying to add some text to the mosaic images
% First try to find the correct locations to add texts
index_array = 1:size(Img_trim_new,3);
index_matrix = reshape(index_array,[4,3])';
for iDx = 1:numel(index_array)-2
    [row,col] = find(index_matrix == iDx);
    location_row = (row-1)*d1+4;
    location_col = (col-1)*d2+1;
    string = sprintf('Slice %d',iDx);
    text(location_col,location_row,string,'Color','w',...
        'FontSize',16)
end
end
function plot_spoil(Img,PlotRange)
Img_trim = trim3D_vert(Img,5);
[d1,d2,d3] = size(Img_trim);
Img_trim_new = Img_trim;
Img_trim_new(:,:,1:3) = Img_trim(:,:,7:9);
Img_trim_new(:,:,7:9) = Img_trim(:,:,1:3);
Img_mosaic = make_mosaic(Img_trim_new,[3,3]);
plot_mosaic(Img_mosaic,[30 35],PlotRange)
% Trying to add some text to the mosaic images
% First try to find the correct locations to add texts
pos_cell = {'high','mid','low'};
spoil_cell = {'1 spoil','2 spoils','3 spoils'};
index_array = 1:d3;
index_matrix = reshape(index_array,[3,3])';
for iDx = 1:d3
    [row,col] = find(index_matrix == iDx);
    location_row = (row-1)*d1+4;
    location_col = (col-1)*d2+1;
    string = strcat(spoil_cell{col},{' '},pos_cell{row});
    text(location_col,location_row,string,'Color','w',...
        'FontSize',16)
end
end

function [fraction,out_back,out_ROI] = find_fraction(Img, Segmentation)
% The value in segmentation can be 1,2,3
Background = Img(Segmentation==2);
ROI = Img(Segmentation==3);
fraction = sum(ROI)/(sum(ROI)+sum(Background));
mean_background = mean(Background);
mean_ROI = mean(ROI);
std_background = std(Background);
std_ROI = std(ROI);
out_back = [mean_background,std_background];
out_ROI = [mean_ROI,std_ROI];
end
function new_img = make_mosaic(input_img,mosaic_mtx)
    if ndims(input_img) ~= 3
        error('The input image has to be 3D.')
    end
    [d1,d2,d3] = size(input_img);
    new_img = zeros(d1*mosaic_mtx(1),d2*mosaic_mtx(2));
    if d3 ~= mosaic_mtx(1)*mosaic_mtx(2)
        error('Cannot make mosaic images of given size.')
    end
    for iDx = 1:d3
        rem = mod(iDx,mosaic_mtx(2));
        if rem == 0
            loc_2 = mosaic_mtx(2);
        else
            loc_2 = rem;
        end
        
        loc_1 = (iDx-loc_2)/mosaic_mtx(2)+1;
        d1_array = d1*(loc_1-1)+1:d1*loc_1;
        d2_array = d2*(loc_2-1)+1:d2*loc_2;
        new_img(d1_array,d2_array) = input_img(:,:,iDx);
    end

end
function [DicomImgCell,DicomImgMtx] = load_dcm(PathName,SeriesArray)
% A script to load dicom files
    dt = Spectro.dicomTree('dir',PathName,'recursive',false);
    DicomImgCell = cell(numel(SeriesArray),1);
    for iDx = 1:numel(SeriesArray)
        matched = dt.search('target','instance','query',...
            @(inst,ser,stu) ser.SeriesNumber==SeriesArray(iDx));       
        for jDx = 1:numel(matched)
            DicomImg_temp =  Spectro.dicomImage({matched(jDx).Filename});
            ImgRaw(:,:,jDx) = DicomImg_temp.image';
            %ImgVol4D(:,:,:,iDx) = Chop(ImgRaw(:,:,iDx));
        end
        DicomImgCell{iDx} = ImgRaw;
    end
    [img_dim_1,img_dim_2,~] = size(ImgRaw);
    DicomImgMtx = zeros(img_dim_1,img_dim_2,numel(SeriesArray));
    for iDx = 1:numel(DicomImgCell)
        DicomImgMtx(:,:,iDx) = DicomImgCell{iDx}(:,:,1);
    end
end
function OutputImg = trim3D_vert(InputImg,voxel2chop)
[dim1,dim2,dim3] = size(InputImg);
OutputImg = zeros(dim1,dim2-voxel2chop*2,dim3);
    for iDx = 1:dim3
        Temp = InputImg(:,:,iDx);
        Temp(:,1:voxel2chop) = [];
        Temp(:,(end-voxel2chop+1):end) = [];
        OutputImg(:,:,iDx) = Temp;
    end
end
function OutputImg = trim3D_hori(InputImg,voxel2chop)
[dim1,dim2,dim3] = size(InputImg);
OutputImg = zeros(dim1-voxel2chop*2,dim2,dim3);
    for iDx = 1:dim3
        Temp = InputImg(:,:,iDx);
        Temp(1:voxel2chop,:) = [];
        Temp((end-voxel2chop+1):end,:) = [];
        OutputImg(:,:,iDx) = Temp;
    end
end
%a function to make image collages
function plot_mosaic(Img,FigDim,varargin)
    %figure()
    if numel(varargin) == 0
        imagesc(Img)
    elseif numel(varargin) == 1
        PlotRange = varargin{1};
        imagesc(Img,PlotRange)
    elseif numel(varargin) == 2
        PlotRange = varargin{1};
        imagesc(Img,PlotRange)
        bool_newfig = varargin{2};
    elseif numel(varargin) == 3
        PlotRange = varargin{1};
        imagesc(Img,PlotRange)
        bool_newfig = varargin{2};
        bool_colorbar = varargin{3};
    elseif numel(varargin) == 4
        PlotRange = varargin{1};
        imagesc(Img,PlotRange)
        bool_newfig = varargin{2};
        bool_colorbar = varargin{3};
        str_colorbar = varargin{4};
    else
        error('Only 2, 3, 4, 5 or 6 inputs are allowed.');
    end
    if (numel(varargin)>=3) && bool_newfig        
        set(gcf,'color','w','InvertHardcopy','off')
        set(gcf,'units','centimeters','position',[4 4 FigDim(1) FigDim(2)],...
            'paperunits','centimeters','paperposition',[4 4 FigDim(1) FigDim(2)])
    end
    axis image        
    axis off
    colormap hot
    if (numel(varargin)>=3) && bool_colorbar
        colorbar
        CLB2 = colorbar('FontSize',20);%nudge(CLB2,[0.05 -0.048 0.01 0.095]);
        set(get(CLB2,'Title'),'String',str_colorbar)
        %nudge(CLB2,[0 0 0.1 0]);
    end
end