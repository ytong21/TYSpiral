%% Testing on Will's data
Path_3T = '/Users/ytong/Documents/Data/3T/F3T_2013_50_107';
SeriesArray_3T = [8,11,14,17,20];
% element no 1: Matching sequence with 20 deg FA
% element no 2: Matching sequence with matched FA
% element no 3: Matching TR with optimised labelling (no double inversion)
% element no 4: Optmised 3T (white paper protocol see Alsop et al. 2015)
% element no 5: M0
TimeSeires_3T = cell(numel(SeriesArray_3T),1);
for iDx = 1:numel(SeriesArray_3T)
    TimeSeires_3T{iDx} = get_time_series(Path_3T,SeriesArray_3T(iDx));
end
%%
Path_7T = '/Users/ytong/Documents/Data/ForISMRMabstract/20181105_F7T_2013_50_102';
SeriesArray_7T = [18,22,28,32,26];
% element no 1: VERSE shimmed
% element no 2: Gaussian shimmed
% element no 3: VERSE CP
% element no 4: Gaussian CP
% element no 5: M0
TimeSeires_7T = cell(numel(SeriesArray_7T),1);
for iDx = 1:numel(SeriesArray_7T)
    TimeSeires_7T{iDx} = get_time_series(Path_7T,SeriesArray_7T(iDx));
end

%%
SliceIndex = 10;
mask_3T = Make_GM_Mask_One_Slice(TimeSeires_3T,SliceIndex);
mask_7T = Make_GM_Mask_One_Slice(TimeSeires_7T,SliceIndex);

%%
tSNR_3T = cell(numel(TimeSeires_3T)-1,1);
for iDx = 1:numel(tSNR_3T)
    tSNR_3T{iDx} = calc_tSNR(squeeze(TimeSeires_3T{iDx}(:,:,SliceIndex,:)),mask_3T);
end
tSNR_7T = cell(numel(TimeSeires_7T)-1,1);
for iDx = 1:numel(tSNR_7T)
    tSNR_7T{iDx} = calc_tSNR(squeeze(TimeSeires_7T{iDx}(:,:,SliceIndex,:)),mask_7T);
end

%%
figure(140)
UpperLim = 4;

for iDx = 1:4
    subplot(2,4,iDx)
    imagesc(tSNR_7T{iDx},[0 UpperLim]);
    axis equal;axis off;colormap('hot')
    subplot(2,4,iDx+4)
    imagesc(tSNR_3T{iDx},[0 UpperLim]);
    axis equal;axis off;colormap('hot')
end

%%
function time_series = get_time_series(path,series_num)
    ObjTemp = ASL_analyse(path,series_num);
    time_series = ObjTemp.Get_DiffTimeSeries();
end

function mask = Make_GM_Mask_One_Slice(TimeSeires,SliceIndex)
    SumImg = zeros(110,110);
    % Might need to change this
    for iDx = 1:numel(TimeSeires)
        SumImg = SumImg + squeeze(sum(TimeSeires{iDx}(:,:,SliceIndex,:),4));
    end
    %
    figure(197);
    imagesc(SumImg,[0 400])
            hTmp = imellipse;
            wait(hTmp);
    ImgArray = reshape(SumImg,[numel(SumImg) 1]);
    level = prctile(ImgArray,70);
    mask = createMask(hTmp).*(SumImg>level);
    close(197)
end

function tSNR_map = calc_tSNR(Img2D_times_series,mask)
    count = 0;
    TimeSeriesArray = zeros(nnz(mask),size(Img2D_times_series,3));

for iDx = 1:size(Img2D_times_series,1)
    for jDx = 1:size(Img2D_times_series,2)
        if mask(iDx,jDx)
            count = count+1;
            TimeSeriesArray(count,:) = squeeze(Img2D_times_series(iDx,jDx,:));
        end
    end
end
    MeanVec = mean(TimeSeriesArray,2);
    StdVec = std(TimeSeriesArray,0,2);
    tSNR = MeanVec./StdVec;
    tSNR_map = zeros(size(mask));
    tSNR_map(logical(mask)) = tSNR;
end