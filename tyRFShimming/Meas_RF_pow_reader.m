%filetext = fileread('/Users/ytong/Documents/Data/RFSWDHistoryTXArray2018-08-13-10.07.03.log');
filetext = fileread('/Users/ytong/Documents/Data/RFSWDHistoryTXArray2018-08-14-14.47.54.log');
expr = struct;
matches = struct;
extract = struct;
expr.max_pow = '[^\n]*max peak[^\n]*';
expr.max_vol = '[^\n]*sTXSPEC\.aRFPULSE\[0\]\.flAmplitude[^\n]*';
expr.ref_vol = '[^\n]*sTXSPEC\.asNucleusInfo\[0\]\.flReferenceAmplitude[^\n]*';
expr.six_min_forw = '[^\n]*max 6min average forw[^\n]*';
expr.ten_sec_forw = '[^\n]*max 10s  average forw[^\n]*';
expr.seq_name = '[^\n]*tSequenceFileName[^\n]*';
% Converting cell arrays to strings
    FieldNameStr = fieldnames(expr);
for iDx = 1:numel(fieldnames(expr))
    Temp = FieldNameStr{iDx};
    matches.(Temp) = regexp(filetext,expr.(Temp),'match');
    if strcmp(Temp,'seq_name') 
        extract.(Temp) = cell(1,numel(matches.(Temp)));
    else
        extract.(Temp) = zeros(1,numel(matches.(Temp)));
    end
    
    for iDy = 1:numel(matches.(Temp))
        if strcmp(Temp,'six_min_forw')||strcmp(Temp,'ten_sec_forw')||strcmp(Temp,'max_pow')
            TempStr = regexp(matches.(Temp){iDy}, '\S\w\S[:]\s[\d\S]+', 'match');
            TempStr = TempStr{1};
            TempStr(1:5) = [];
            extract.(Temp)(1,iDy) = str2double(TempStr);
        elseif ~strcmp(Temp,'seq_name')
            TempStr = regexp(matches.(Temp){iDy}, '=\s\d+', 'match');
            TempStr = TempStr{1};
            TempStr(1:2) = [];
            extract.(Temp)(1,iDy) = str2double(TempStr);
        else
            TempStr = regexp(matches.(Temp){iDy}, '=\s\S+', 'match');
            TempStr = TempStr{1};
            TempStr(1:4) = []; TempStr(end) = [];
            extract.(Temp){1,iDy} = TempStr;
        end
    end
end

%%  Looking for a new way to read in the data. Separate the string based on a delimiter
filetext = fileread('/Users/ytong/Documents/Data/RFSWDHistoryTXArray2018-08-14-14.47.54.log');
string_separated = strsplit(filetext,'-----------------------------------');
% expr.max_vol = '[^\n]*sTXSPEC\.aRFPULSE\[0\]\.flAmplitude[^\n]*';
% expr.ref_vol = '[^\n]*sTXSPEC\.asNucleusInfo\[0\]\.flReferenceAmplitude[^\n]*';
% expr.seq_name = '[^\n]*tSequenceFileName[^\n]*';
% expr.prot_name = '[^\n]*tProtocolName[^\n]*';

expression.power.max_pow = '[^\n]*max peak[^\n]*';
expression.power.six_min_forw = '[^\n]*max 6min average forw[^\n]*';
expression.power.ten_sec_forw = '[^\n]*max 10s  average forw[^\n]*';

expression.RF.max_vol = '[^\n]*sTXSPEC\.aRFPULSE\[0\]\.flAmplitude[^\n]*';
expression.RF.ref_vol = '[^\n]*sTXSPEC\.asNucleusInfo\[0\]\.flReferenceAmplitude[^\n]*';
expression.seq.seq_name = '[^\n]*tSequenceFileName[^\n]*';
expression.seq.prot_name = '[^\n]*tProtocolName[^\n]*';
info_extracted = struct;
%%
IdxCount = 0;
for iDx = 1:numel(string_separated)-1
    %   First find out the max peak power and rthe 6min and 10s avg power
    if strfind(string_separated{iDx+1},'PALI MESSAGE: power limit for sum channel exceeded:')
        FieldNameStr = fieldnames(expression);
        for iDx2 = 1:numel(FieldNameStr)
            TempField = FieldNameStr{iDx2};
            info_extracted(iDx).(TempField) = [];
        end
    else
        FieldNameStr = fieldnames(expression.power);
        for iDx2 = 1:numel(FieldNameStr)
            TempField = FieldNameStr{iDx2};
            The_Line = regexp(string_separated{iDx+1},expression.power.(TempField),'match');
            TempStr = regexp(The_Line{1}, '\S\w\S[:]\s[\d\S]+', 'match');
            TempStr = TempStr{1};
            TempStr(1:5) = [];
            info_extracted(iDx).power.(TempField) = str2double(TempStr);
        end
        %   Find out if the string contains "### ASCCONV BEGIN ###"
        if strfind(string_separated{iDx+1},'### ASCCONV BEGIN ###')
        %   Find sequence name and protocol name
            FieldNameStr = fieldnames(expression.seq);
            for iDx2 = 1:numel(FieldNameStr)
                TempField = FieldNameStr{iDx2};
                The_Line = regexp(string_separated{iDx+1},expression.seq.(TempField),'match');
                TempStr = regexp(The_Line{1}, '=\s\S+', 'match');
                TempStr = TempStr{1};
                TempStr(1:4) = []; TempStr(end) = [];
                info_extracted(iDx).seq.(TempField) = TempStr;
            end
        %else
        end
    end
end

%% figure(100)
% clf
% plot(8*AA.^2/50,[2103 3723 5895],'x')
% hold on
% plot([0 6000],1.64723661485319*[0 6000]-46.3678756476679,'k:')
% 
% %%
%  Expected_Pow = 8*extract.max_vol(12:24).^2/50;
%  Meas_Pow = extract.max_pow(13:25);
% figure(128)
% clf
% set(gcf,'color','w','InvertHardcopy','off')
% set(gcf,'units','centimeters','position',[4 4 25 20],'paperunits','centimeters','paperposition',[0 0 25 20])
% plot(Expected_Pow,Meas_Pow,'or','MarkerSize',10)
% pbaspect([1 1 1])
% hold on
% p_CP = polyfit(Expected_Pow,Meas_Pow,1);
% plot([0 3000],p_CP(1)*[0 3000]+p_CP(2),'k:')
% pbaspect([1 1 1])
% title('PALI MC vs expected coil plug power (CP mode)','FontSize',14)
% xlabel('Expected coil pug power (W)','FontSize',14)
% ylabel('Measured PALI MC power (W)','FontSize',14)
% if sign(p_CP(2)) == 1
%     signstr = '+';
% elseif sign(p_CP(2)) == -1
%     signstr = '-';
% end
% fit_CP = sprintf('Linear fit:  y = %.1fx %s %.1f',p_CP(1),signstr, abs(p_CP(2)));
% lgd = legend('Measurements',fit_CP);
% lgd.FontSize = 12;
% 
% %%
% 
% Expected_Pow_onechan = extract.max_vol(25:end).^2/50;
% Meas_Pow_onechan = extract.max_pow(27:end);
% p_shim = polyfit(Expected_Pow_onechan,Meas_Pow_onechan,1);
% figure(129)
% clf
% set(gcf,'color','w','InvertHardcopy','off')
% set(gcf,'units','centimeters','position',[4 4 25 20],'paperunits','centimeters','paperposition',[0 0 25 20])
% plot(Expected_Pow_onechan,Meas_Pow_onechan,'or','MarkerSize',10)
% pbaspect([1 1 1])
% hold on
% plot([0 350],p_shim(1)*[0 350]+p_shim(2),'k:')
% pbaspect([1 1 1])
% title('PALI MC vs expected coil plug power (one chan on)','FontSize',14)
% xlabel('Expected coil pug power (W)','FontSize',14)
% ylabel('Measured PALI MC power (W)','FontSize',14)
% if sign(p_shim(2)) == 1
%     signstr = '+';
% elseif sign(p_shim(2)) == -1
%     signstr = '-';
% end
% fit_shim = sprintf('Linear fit:  y = %.1fx %s %.1f',p_shim(1),signstr, abs(p_shim(2)));
% lgd = legend('Measurements',fit_shim);
% lgd.FontSize = 12;

%% Comparison between CP avg and one master and one slave channel on
Expected_Pow_onechan = extract.max_vol(25:end).^2/50;
Meas_Pow_onechan = extract.max_pow(27:end);

p_130 = zeros(3,2);
p_130(1,:) = polyfit(Expected_Pow,Meas_Pow,1);
p_130(2,:) = polyfit(Expected_Pow_onechan,Meas_Pow_onechan,1);
p_130(3,:) = polyfit(slave_ramp(1,:),slave_ramp(2,:),1);
%%
figure(130)
clf
set(gcf,'color','w','InvertHardcopy','off')
set(gcf,'units','centimeters','position',[4 4 25 20],'paperunits','centimeters','paperposition',[0 0 25 20])
plot(Expected_Pow/8,Meas_Pow/8,'*b','MarkerSize',9)             % CP
pbaspect([1 1 1])
hold on
plot(Expected_Pow_onechan,Meas_Pow_onechan,'dm','MarkerSize',9) % master
plot(slave_ramp(1,:),slave_ramp(2,:),'or','MarkerSize',9)       % slave 1
title('PALI MC vs coil plug power','FontSize',14)
xlabel('Expected coil plug power (W)','FontSize',14)
ylabel('Measured PALI MC power (W)','FontSize',14)
lgd = legend('CP mode avg','Master chan 0','Slave chan 1');

plot([0 350],p_130(1,1)*[0 350]+p_130(1,2)/8,'b:')
plot([0 350],p_130(2,1)*[0 350]+p_130(2,2),'m:')
plot([0 350],p_130(3,1)*[0 350]+p_130(3,2),'r:')
%%
Multi_chan_plug = ones(1,9)*120^2/50;
Multi_chan_PALI = extract.max_pow(22:29);
Multi_chan_PALI = [extract.max_pow(12)/8 Multi_chan_PALI];

%%
figure(131)
clf
set(gcf,'color','w','InvertHardcopy','off')
set(gcf,'units','centimeters','position',[4 4 18 20],'paperunits','centimeters','paperposition',[0 0 18 20])
plot(Multi_chan_plug(1),Multi_chan_PALI(1),'*b','MarkerSize',10)
%pbaspect([1 1 1])
xlim([278 298])
hold on
%plot(Expected_Pow_onechan,Meas_Pow_onechan,'ob','MarkerSize',10)

title('PALI MC vs coil plug power','FontSize',14)
xlabel('Expected coil plug power (W)','FontSize',14)
ylabel('Measured PALI MC power (W)','FontSize',14)
plot(Multi_chan_plug(9),Multi_chan_PALI(9),'dm','MarkerSize',10)
lgdstr = cell(1,9);
lgdstr{1} = 'CP mode avg';
lgdstr{2} = 'Master 0';

for iDx = 2:8
    if iDx == 2
        plot(Multi_chan_plug(iDx),Multi_chan_PALI(iDx),'or','MarkerSize',10)
    else
        plot(Multi_chan_plug(iDx),Multi_chan_PALI(iDx),'o','MarkerSize',10)
    end
    lgdstr{iDx+1} = ['Slave ',num2str(iDx-1)];
end
set(gca,'XTick',Multi_chan_plug(1))
lgdd = legend(lgdstr);
lgdd.FontSize = 12;
 %%
% figure(132)
% clf
% 
% set(gcf,'color','w','InvertHardcopy','off')
% set(gcf,'units','centimeters','position',[4 4 25 20],'paperunits','centimeters','paperposition',[0 0 25 20])
% 
%%
CP_ramp = zeros(2,6);
CP_ramp(1,:) = extract.max_vol(7:12).^2/50;
CP_ramp(2,:) = extract.max_pow(7:12)/8;

slave_ramp = zeros(2,6);
slave_ramp(1,:) = CP_ramp(1,:);
slave_ramp(2,:) = extract.max_pow(14:19);
slave_ramp(2,2) = extract.max_pow(end);

% 
% plot(CP_ramp(1,:),CP_ramp(2,:),'o','MarkerSize',10)
% ylim([0 600])
% hold on
% plot(slave_ramp(1,:),slave_ramp(2,:),'d','MarkerSize',10)
% title('PALI MC vs expected coil plug power','FontSize',14)
% xlabel('Expected coil pug power (W)','FontSize',14)
% ylabel('Measured PALI MC power (W)','FontSize',14)
% lgdd = legend('CP mode avg','Slave 1');
% lgdd.FontSize = 12;
% 
