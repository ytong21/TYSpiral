function info_extracted = read_RFSWD_pTx(filename)
%  Looking for a new way to read in the data. Separate the string based on a delimiter
filetext = fileread(filename);
string_separated = strsplit(filetext,'-----------------------------------');
% expr.max_vol = '[^\n]*sTXSPEC\.aRFPULSE\[0\]\.flAmplitude[^\n]*';
% expr.ref_vol = '[^\n]*sTXSPEC\.asNucleusInfo\[0\]\.flReferenceAmplitude[^\n]*';
% expr.seq_name = '[^\n]*tSequenceFileName[^\n]*';
% expr.prot_name = '[^\n]*tProtocolName[^\n]*';

expression.power.max_pow = '[^\n]*max peak[^\n]*';
expression.power.six_min_forw = '[^\n]*max 6min average forw[^\n]*';
expression.power.ten_sec_forw = '[^\n]*max 10s  average forw[^\n]*';

%expression.RF.max_vol = '[^\n]*sTXSPEC\.aRFPULSE\[0\]\.flAmplitude[^\n]*';
%expression.RF.ref_vol = '[^\n]*sTXSPEC\.asNucleusInfo\[0\]\.flReferenceAmplitude[^\n]*';
expression.RF.gen = '[^\n]*sTXSPEC\.aRFPULSE\[\d\]\.[^\n]*';
%expression.RF.tName = '[^\n]*sTXSPEC\.aRFPULSE\[\d\]\.[^\n]*';
expression.seq.seq_name = '[^\n]*tSequenceFileName[^\n]*';
expression.seq.prot_name = '[^\n]*tProtocolName[^\n]*';
%expression.RF.parts.tName = '[^\n]*sTXSPEC\.aRFPULSE\[\d\]\.tName[^\n]*';
%expression.RF.parts.tName = '[^\n]*sTXSPEC\.aRFPULSE\[\d\]\.tName[^\n]*';
expression.seq.TR = '[^\n]*alTR\[0\][^\n]*';
expression.seq.tag_dur = '[^\n]*sWiPMemBlock\.adFree\[11\][^\n]*';
expression.seq.gaussian_dur = '[^\n]*sWiPMemBlock\.alFree\[7\][^\n]*';
expression.seq.rf_sep = '[^\n]*sWiPMemBlock\.alFree\[8\][^\n]*';
info_extracted = struct;

IdxCount = 0;
ExceededIdx = [];
for iDx = 1:numel(string_separated)-1
    if ~isempty(strfind(string_separated{iDx+1},'PALI MESSAGE: power limit for sum channel exceeded:'))
        ExceededIdx = [ExceededIdx IdxCount+1];
    elseif isempty(strfind(string_separated{iDx+1},'PALI MESSAGE: power limit for sum channel exceeded:')) ...
            && ~isempty(strfind(string_separated{iDx+1},'### ASCCONV BEGIN ###')) ...
            && isempty(strfind(string_separated{iDx+1},'%SERV_SEQ%\txa_adjust')) ...
            && isempty(strfind(string_separated{iDx+1},'localizer'))...
            && isempty(strfind(string_separated{iDx+1},'Coil ID                               : Body'))...
            && isempty(strfind(string_separated{iDx+1},'sj_gre_sample'))
        IdxCount = IdxCount+1;
    %   First find out the max peak power and rthe 6min and 10s avg power
        FieldNameStr = fieldnames(expression.power);
        for iDx2 = 1:numel(FieldNameStr)
            TempField = FieldNameStr{iDx2};
            The_Line = regexp(string_separated{iDx+1},expression.power.(TempField),'match');
            TempStr = regexp(The_Line{1}, '\S\w\S[:]\s[\d\S]+', 'match');
            TempStr = TempStr{1};
            TempStr(1:5) = [];
            info_extracted(IdxCount).power.(TempField) = str2double(TempStr);
        end

        %   Find sequence name and protocol name, tag_dur, gauss_dur and rf
        %   sep. (if it is VEPCASL)
        if isempty(strfind(string_separated{iDx+1},'to_VEPCASL_perf_ep2d_pasl_BGS_Arb_mPLD'))
            FieldNameStr = {'seq_name','prot_name','TR'};
        else
            FieldNameStr = fieldnames(expression.seq);
        end
            for iDx2 = 1:numel(FieldNameStr)
                TempField = FieldNameStr{iDx2};
                The_Line = regexp(string_separated{iDx+1},expression.seq.(TempField),'match');
                TempStr = regexp(The_Line{1}, '=\s\S+', 'match');
                TempStr = TempStr{1};
                if strcmp(TempField,'seq_name')||strcmp(TempField,'prot_name')
                    TempStr(1:3) = []; TempStr(end) = [];
                else
                    TempStr(1:2) = []; TempStr = str2double(TempStr);
                end
                info_extracted(IdxCount).seq.(TempField) = TempStr;
            end

        
        %   Find the info about RF pulses. 
        RFTempAll = regexp(string_separated{iDx+1},expression.RF.gen,'match');
        %   There should be 3 matches for each rf pulse: tName
        %   bAmplitudeValid and flAmplitude/bAmplitudeValid
        num_of_pulses = regexp(RFTempAll{end},'\[\d*\]','match'); 
        num_of_pulses = num_of_pulses{1};
        num_of_pulses([1 end]) = [];
        num_of_pulses = str2double(num_of_pulses)+1;
        info_extracted(IdxCount).RF.tName = cell(1,num_of_pulses);
        info_extracted(IdxCount).RF.pulse_amp = zeros(1,num_of_pulses);
        
        for iDx2 = 1:num_of_pulses
            % Constructing the expressions
            temp_expression.tName = ['[^\n]*sTXSPEC\.aRFPULSE\[' num2str(iDx2-1) '\]\.tName[^\n]*'];
            temp_expression.flAmplitude = ['[^\n]*sTXSPEC\.aRFPULSE\[' num2str(iDx2-1) '\]\.flAmplitude[^\n]*'];
            temp_expression.bManualAmplitude = ['[^\n]*sTXSPEC\.aRFPULSE\[' num2str(iDx2-1) '\]\.bManualAmplitude[^\n]*'];
            
            % Finding tName
            The_Line = regexp(string_separated{iDx+1}, temp_expression.tName , 'match');
            TempStr = regexp(The_Line{1}, '=\s\S+', 'match');   TempStr = TempStr{1};
            TempStr([1:3 end]) = [];
            info_extracted(IdxCount).RF.tName{iDx2} = TempStr;
            
            % Finding pulse amp
            % If the amp is changed mannually, and flAmplitude does not
            % exisit, amp = 0;
            The_Line = regexp(string_separated{iDx+1},temp_expression.flAmplitude, 'match');
            if ~isempty(The_Line)
                TempStr = regexp(The_Line, '=\s\S+', 'match');
                TempStr = TempStr{1}{1};
                TempStr(1:2) = [];
                info_extracted(IdxCount).RF.pulse_amp(iDx2) = str2double(TempStr);
            end
        end
        
%         for iDx2 = 1:numel(RFTempAll)
%             remainder = mod(iDx2,3);
%             RF_index = (iDx2-remainder)/3+1;
%             if remainder == 1
%             	TempStr = regexp(RFTempAll{iDx2}, '=\s\S+', 'match'); 
%                 TempStr = TempStr{1};
%                 TempStr(1:3) = []; TempStr(end) = [];
%                 info_extracted(IdxCount).RF.tName{RF_index} = TempStr;
%             elseif remainder == 0
%             	TempStr = regexp(RFTempAll{iDx2}, '=\s\d+', 'match');
%                 TempStr = TempStr{1};
%                 TempStr(1:2) = [];
%                 info_extracted(IdxCount).RF.pulse_amp(RF_index) = str2double(TempStr);
%             end
% 
%         end
    end
end