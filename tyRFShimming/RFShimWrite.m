function RFShimWrite(RFToWrite)
%Write the file

fID = fopen('/Users/ytong/Documents/MATLAB/tong-acptx/tySpiral/tyRFShimming/pTxVEPCASLShim.ini','w');
fprintf(fID,'# Created in run_STA_KTPdesign.m \n');
fprintf(fID,'[pTXPulse]\n');
fprintf(fID,'\n');
fprintf(fID,'NUsedChannels        = %i\n',8);		
fprintf(fID,'NShims               = %i\n',1);
fprintf(fID,'PulseName            = TagShim\n');
fprintf(fID,'Comment              = Shim for VEPCASL tag pulse\n');


% RF
for jDx = 1:8
    fprintf(fID,'[pTXPulse_ch%i]\n',jDx-1);
    fprintf(fID,'\n');
    fprintf(fID,'RF[%i]=  %f\t%f\n',0,RFToWrite(jDx,1),RFToWrite(jDx,2));
    fprintf(fID,'\n');
end

fclose(fID);