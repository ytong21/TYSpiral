function PlotPerfusionNew(dt,SerieNumberArray,LabEff,EPI_series)

%   SerieNumberArray: 4 element array 
%   1. VShimed 2.G_Shimmed 3.G_CP 4.M0
matched = cell(size(SerieNumberArray));
DicomImg = cell(size(SerieNumberArray));
ImgRaw = cell(size(SerieNumberArray));
ImgVol = cell(size(SerieNumberArray));
Img2Plot = cell(size(SerieNumberArray));
ImgRawSub = cell(size(SerieNumberArray)-1);
plot_array = 6:17;
CBF = cell(1,3);
for iDx = 1:numel(SerieNumberArray)
    matched{iDx} = dt.search('target','instance','query',...
        @(inst,ser,stu) ser.SeriesNumber==SerieNumberArray(iDx));
    if isempty(matched{iDx})
        disp('Please check series number. Not images found');
        return
    end
    
    DicomImg{iDx} = Spectro.dicomImage({matched{iDx}(1).Filename});
    ImgRaw{iDx} = DicomImg{iDx}.image';
    
    %Slope = DicomImg{iDx}.info{1}.RescaleSlope;
    %Intercept = DicomImg{iDx}.info{1}.RescaleIntercept;
    ImgTmp = Chop(ImgRaw{iDx});
    
    

%     if iDx < 4
%          AA = -ImgTmp*Slope+Intercept;
%          BB = AA(:,:,plot_array);
%          CC = reshape(BB,[110,1320]);
%          ImgRawSub{iDx} = [CC(:,1:440);CC(:,441:880);CC(:,881:1320)];
%     end
    if iDx ~= 4
        ImgVol{iDx} = 4095-ImgTmp;
    else
        ImgVol{iDx} = ImgTmp;
    end
    ToPlotTemp = ImgVol{iDx}(:,:,plot_array);
    ToPlotTemp = reshape(ToPlotTemp,[110,1320]);
    Img2Plot{iDx} = [ToPlotTemp(:,1:440);ToPlotTemp(:,441:880);ToPlotTemp(:,881:1320)];
end


    matched{5} = dt.search('target','instance','query',...
        @(inst,ser,stu) ser.SeriesNumber==EPI_series);
    if isempty(matched{5})
        disp('Please check series number. Not images found');
        return
    end
    DicomImg{5} = cell(numel(matched{5}),1);
    for iDx = 1:numel(matched{5})
            DicomImg{5}{iDx} = Spectro.dicomImage({matched{5}(iDx).Filename});
            ImgRaw{5}(:,:,iDx) = DicomImg{5}{iDx}.image';
            ImgTmp = Chop(ImgRaw{5}(:,:,iDx));
            ImgVol{5}(:,:,:,iDx) = ImgTmp;
            ToPlotTemp = ImgVol{5}(:,:,plot_array,iDx);
            ToPlotTemp = reshape(ToPlotTemp,[110,1320]);
            Img2Plot{5}(:,:,iDx) = [ToPlotTemp(:,1:440);ToPlotTemp(:,441:880);ToPlotTemp(:,881:1320)];
    end
    TagArray = 1:2:31;
    ControlArray = 2:2:32;
    AvgDiff = mean(Img2Plot{5}(:,:,ControlArray) - Img2Plot{5}(:,:,TagArray),3);
    AvgDiff(AvgDiff<=0) = 0;
    %TotDiff = sum(Img2Plot{5}(:,:,ControlArray) - Img2Plot{5}(:,:,TagArray),3);
    Thresh = 200;
    M0 = Img2Plot{4};
    M0(Img2Plot{4}<Thresh) = inf;
    T1b = 2587e-3;
    up = 6000*0.9*AvgDiff*exp(1/T1b);
    down = 2*LabEff*T1b*M0*(1-exp(-1/T1b));
    CBF{1} = up./down;
    CBF{1}(CBF{1}>500) = 0;
    %%
figure(199)
    set(gcf,'color','w','InvertHardcopy','off')
    Width = 35; Length = 30;LowerBound = 1800;FontSizeTmp = 13;
    set(gcf,'units','centimeters','position',[4 4 Width Length],'paperunits',...
        'centimeters','paperposition',[0 0 Width Length]);
    clf
    All_Points = [Img2Plot{1}(:); Img2Plot{2}(:); Img2Plot{3}(:)];
    UpperBound = prctile(All_Points,99);
    Title = {'VERSE Shimmed','Gaussian Shimmed','Gaussian CP mode','CBF Map from VERSE'};
    SubPlotsCell = cell(4,1);
for iDx = 1:3
        SubPlotsCell{iDx} = subplot(2,2,iDx);
        imagesc(Img2Plot{iDx},[LowerBound UpperBound]);colormap('gray');axis equal;axis off
        title(Title{iDx},'FontSize',FontSizeTmp)
end
SubPlotsCell{4} = subplot(2,2,4);

imagesc(CBF{1},[0 130]);colormap(gca,'jet');axis equal;axis off
        title(Title{4},'FontSize',FontSizeTmp)
clb_obj = colorbar('FontSize',13);
clb_obj.XLabel.String = 'ml/100g/min';
nudge(clb_obj,[0.05 0 0 0]);
%Amplify = 1.15;
% for iDx = 1:3
%         SubPlotsCell{iDx}.Title.Units = 'normalized';
%         SubPlotsCell{iDx}.Position = [SubPlotsCell{iDx}.Position(1)-0.05 SubPlotsCell{iDx}.Position(2)-0.05...
%             SubPlotsCell{iDx}.Position(3)*Amplify SubPlotsCell{iDx}.Position(4)*Amplify];
% 
%         SubPlotsCell{iDx}.Title.String = Title{iDx};
%         SubPlotsCell{iDx}.Title.FontWeight = 'Bold';
%         SubPlotsCell{iDx}.Title.FontSize = FontSizeTmp;
%         
%         SubPlotsCell{iDx}.Title.Position(2) = SubPlotsCell{iDx}.Title.Position(2)-0.1;
%         SubPlotsCell{iDx}.Title.Position(1) = SubPlotsCell{iDx}.Position(1)...
%             +SubPlotsCell{iDx}.Position(3)/2;%-SubPlotsCell{iDx}.Title.Position(2)/2;

% end
%%
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
end