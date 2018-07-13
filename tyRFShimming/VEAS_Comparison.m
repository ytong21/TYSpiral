function VEAS_Comparison(AFull,bVE,bAS,minIndex)

if size(bVE,2) ~= size(bAS,2)
    disp('variable-exchange and active-set solution size mismatch');
else
    %FA_VE = zeros(size(AFull,1),1);
    %FA_AS = zeros(size(AFull,1),1);
    %for iDx = 1:size(bVE,2)
   
        FA_VE = AFull*(bVE(1:8,minIndex));
        FA_AS = AFull*(bAS(1:8,minIndex).*exp(1i*bAS(9:16,minIndex)));
    %end
    fprintf('FA_VE is %d of FA_AS.\n',mean(FA_VE)/mean(FA_AS));
end