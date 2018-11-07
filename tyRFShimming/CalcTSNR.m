function CalcTSNR(dt,SeriesArray)   
DicomImg = cell(size(SeriesArray));
ImgRaw = cell(size(SeriesArray));
ImgVol4D = cell(size(SeriesArray));
DiffImgTimeSeries = cell(size(SeriesArray));
TagArray = 1:2:31;
ControlArray = 2:2:32;
for iDx = 1:numel(SeriesArray)
    [DicomImg{iDx},ImgRaw{iDx},ImgVol4D{iDx}] = LoadNSort(dt,SeriesArray(iDx));
    DiffImgTimeSeries{iDx} = ImgVol4D{iDx}(:,:,:,ControlArray) - ImgVol4D{iDx}(:,:,:,TagArray) ;
end
DiffImgSum = squeeze(sum(DiffImgTimeSeries{1}+DiffImgTimeSeries{2}+DiffImgTimeSeries{3},4));
SliceIdx = 10;

%
figure(197);
imagesc(DiffImgSum(:,:,SliceIdx),[0 400])
        hTmp = imellipse;
        wait(hTmp);
level = prctile(reshape(DiffImgSum(:,:,SliceIdx),[numel(DiffImgSum(:,:,SliceIdx)),1]),70);
mask = createMask(hTmp).*(DiffImgSum(:,:,SliceIdx)>level);
close(197)
figure(198)
imagesc(mask.*DiffImgSum(:,:,SliceIdx))
pause(1);
%ImgToCalcTSNR = mask.*DiffImgSum(:,:,SliceIdx);
TimeSeriesArray = zeros(nnz(mask),size(DiffImgTimeSeries{1},4),3);
count = 0;
for iDx = 1:size(DiffImgTimeSeries{1},1)
    for jDx = 1:size(DiffImgTimeSeries{1},2)
        if mask(iDx,jDx)
            count = count+1;
            TimeSeriesArray(count,:,1) = squeeze(DiffImgTimeSeries{1}(iDx,jDx,SliceIdx,:));
            TimeSeriesArray(count,:,2) = squeeze(DiffImgTimeSeries{2}(iDx,jDx,SliceIdx,:));
            TimeSeriesArray(count,:,3) = squeeze(DiffImgTimeSeries{3}(iDx,jDx,SliceIdx,:));
        end
    end
end
AA = DiffImgTimeSeries;

end
function [DicomImg,ImgRaw,ImgVol4D] = LoadNSort(dt,SeriesNum)
    matched = dt.search('target','instance','query',...
        @(inst,ser,stu) ser.SeriesNumber==SeriesNum);
    DicomImg = cell(numel(matched),1);
    ImgRaw = zeros(550,550,numel(matched));
    for iDx = 1:numel(matched)
        DicomImg{iDx} =  Spectro.dicomImage({matched(iDx).Filename});
        ImgRaw(:,:,iDx) = DicomImg{iDx}.image';
        ImgVol4D(:,:,:,iDx) = Chop(ImgRaw(:,:,iDx));
    end
end
function ImgTmp = Chop(ImgRaw)
        ImgTmp = zeros(110,110,25);
        for ii = 1:size(ImgRaw,1)
            for jj = 1:size(ImgRaw,2)
                vert_img_num = floor((ii-0.5)/110)+1;
                hori_img_num = floor((jj-0.5)/110)+1;
                slice_num = (vert_img_num-1)*5+hori_img_num;
                vert_idx = rem(ii,110);
                hori_idx = rem(jj,110);
                if vert_idx == 0
                    vert_idx = 110;
                end
                if hori_idx == 0;
                    hori_idx = 110;
                end
                ImgTmp(vert_idx,hori_idx,slice_num) = ImgRaw(ii,jj);
            end
        end
end