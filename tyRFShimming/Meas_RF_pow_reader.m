filetext = fileread('/Users/ytong/Documents/Data/RFSWDHistoryTXArray2018-07-11-14.57.33.log');
expr_max_pow = '[^\n]*max peak[^\n]*';
matches_max_pow = regexp(filetext,expr_max_pow,'match');
%%
expr_max_vol = '[^\n]*sTXSPEC\.aRFPULSE\[0\]\.flAmplitude[^\n]*';
matches_max_vol = regexp(filetext,expr_max_vol,'match');

%%
meas_pow = zeros(1,numel(matches_max_pow));
for iDx = 1:numel(matches_max_pow)
     temp = regexp(matches_max_pow{1,iDx},'\d*[.]\d*','match');
     meas_pow(1,iDx) = str2double(temp{1});
end
%%
expected_vol = zeros(1,numel(matches_max_vol));
for iDx = 1:numel(matches_max_vol)
    temp = regexp(matches_max_vol{1,iDx},'=\s\d+','match');
    temp = temp{1};
    temp(1:2) = [];
    expected_vol(1,iDx) = str2double(temp);
end