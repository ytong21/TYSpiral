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
%%
[~,Spoils_180] = load_dcm(PathName,113:121);
[~,Spoils_223] = load_dcm(PathName,122:130);
%%
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
fraction = zeros(1,size(Mid_DiffV,3));
for iDx = 1:size(Mid_DiffV,3)
    fraction(iDx) = find_fraction(Mid_DiffV(:,:,iDx), seg_new);
end
%% 
function chart_spoil(Spoils_180, Spoils_223,seg)
    mean_back_180 = zeros(1,size(Spoils_180,3));
    mean_back_223 = zeros(1,size(Spoils_223,3));
    for iDx = 1:numel(mean_back_180)
        [~,mean_std_temp,~]=find_fraction(Spoils_180(:,:,iDx),seg);
        mean_back_180(iDx) = mean_std_temp(1);
        [~,mean_std_temp,~]=find_fraction(Spoils_223(:,:,iDx),seg);
        mean_back_223(iDx) = mean_std_temp(1);
    end
    m_180 = reshape(mean_back_180,[3,3]);
    m_223 = reshape(mean_back_223,[3,3]);
    x_data = 1:3;       % 1,2,3 spoils
    color_cell = {[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4660 0.6740 0.1880]};
    hold on
    for iDx = 1:3
        plot(x_data,m_180(:,iDx),'-o','LineWidth',3,'MarkerSize',17,'Color',color_cell{iDx})
    end
    for iDx = 1:3
            plot(x_data,m_223(:,iDx),'--d','LineWidth',3,'MarkerSize',17,'Color',color_cell{iDx})
    end
    xlim([0.8 3.2]);xticks(1:3);yticks(100:100:500)
    xlabel('No. of sat pulses');ylabel('Mean backgroud pixel intensity (a.u.)')
    legend('180V low','180V mid','180V high','223V low','223V mid','223V high')
    set(gcf,'color','w','InvertHardcopy','off')
    set(gcf,'units','centimeters','position',[4 4 30 30],...
        'paperunits','centimeters','paperposition',[4 4 30 30])
    set(gca,'FontSize',22)
end
function plot_diffV(Img,PlotRange)
Img_trim = trim3D_vert(Img,5);
[d1,d2,~] = size(Img_trim);
Img_mosaic = make_mosaic(Img_trim,[4,4]);
plot_mosaic(Img_mosaic,[30 35],PlotRange)
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
Img_trim = trim3D_vert(Img,5);
[d1,d2,~] = size(Img_trim);
% add an extra empty image
Img_trim_new = cat(3,Img_trim,zeros(d1,d2));
Img_mosaic = make_mosaic(Img_trim_new,[3,4]);
plot_mosaic(Img_mosaic,[28 38],PlotRange)
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
    figure()
    if numel(varargin) == 0
        imagesc(Img)
    elseif numel(varargin) == 1
        PlotRange = varargin{1};
        imagesc(Img,PlotRange)
    else
        error('Only 2 or 3 inputs are allowed.');
    end
    set(gcf,'color','w','InvertHardcopy','off')
    set(gcf,'units','centimeters','position',[4 4 FigDim(1) FigDim(2)],...
        'paperunits','centimeters','paperposition',[4 4 FigDim(1) FigDim(2)])
    axis image        
    axis off
    colormap hot
end