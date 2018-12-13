classdef ASL_analyse < handle
  % pTx_design:  a pTx pulse design class
   properties (Access = public)
 
   end
   properties (Access = protected)
        DataPathName;
        DicomTree;
        FieldStrength;
        SeriesFolderName3T;
        SeriesNum;
        DicomImg;
        ImgRaw;
        Img3D_TimeSeries;
        DiffImgTimeSeries;
        DiffImgTimeSeriesTotal;
   end
   methods
       function obj = ASL_analyse(DataPathName,SeriesNum)
           % Constructor
           obj.DataPathName = DataPathName;
           obj.DetermineFieldStrength();
           obj.SeriesNum = SeriesNum;
           if obj.FieldStrength == 3
               obj.Find3tSeriesFolders();
           end
           obj.SetDicomTree();
           obj.load();
           obj.Calc_Difference_Img();
       end
       
       function obj = DetermineFieldStrength(obj)
           [~,FileName,~] = fileparts(obj.DataPathName);
           PathNameSplitted = strsplit(FileName,'_');
           if (numel(PathNameSplitted) == 5) && (strcmp(PathNameSplitted{2},'F7T'))
               obj.FieldStrength = 7;
           elseif (numel(PathNameSplitted) == 4) && (strcmp(PathNameSplitted{1},'F3T'))
               obj.FieldStrength = 3;
           else
               disp('Path format unregocnised.');
               obj.FieldStrength = NaN;
           end
       end
       
       function obj = Find3tSeriesFolders(obj)
           if obj.FieldStrength == 3
               DirStruct = dir(obj.DataPathName);
               SeriesFolderName = cell(numel(DirStruct)-2,1);
               % Excluding '.' and '..'
               for iDx = 1:numel(DirStruct)
                   if (~strcmp(DirStruct(iDx).name,'.')) && (~strcmp(DirStruct(iDx).name,'..'))...
                           && (~strcmp(DirStruct(iDx).name,'.DS_Store'))
                       TempCell = strsplit(DirStruct(iDx).name,'_');
                       SeriesIdxTemp = int8(str2double(TempCell{end-1}));
                       SeriesFolderName{SeriesIdxTemp} = DirStruct(iDx).name;
                   else
                       disp('Path format unregocnised.');
                   end
               end
               obj.SeriesFolderName3T = SeriesFolderName;
           else
               disp('This function intended for 3T data only.')
           end
       end
       
       function obj = SetDicomTree(obj)
           if ~isprop(obj,'FieldStrength')
               disp('FieldStrength has to be defined first.')
               return
           else
               switch obj.FieldStrength
                   case 7
                       obj.DicomTree = Spectro.dicomTree('dir',obj.DataPathName,'recursive',false);
                   case 3
                       PathTemp = obj.SeriesFolderName3T{obj.SeriesNum};
                       PathFull = fullfile(obj.DataPathName,PathTemp);
                       obj.DicomTree = Spectro.dicomTree('dir',PathFull,'recursive',false);
               end
           end
       end
       function ImgRawSorted = SortGluedImage(~,ImgRaw2D)
           SortedIndicies = [1:110;111:220;221:330;331:440;441:550];
           % 5 by 110 matrix
           Rows = cell(5,1);
           ImgCell = cell(5,5);
           for iDx = 1:5
               Rows{iDx} = ImgRaw2D(SortedIndicies(iDx,:),:);
               TempRow = Rows{iDx};
               for jDx = 1:5
                   ImgCell{iDx,jDx} = TempRow(:,SortedIndicies(iDx,:));
               end
           end
           % One step at a time
           ImgRawSorted = zeros(110,110,25);
           for iDx = 1:5
               for jDx = 1:5
                   ImgRawSorted(:,:,(iDx-1)*5+jDx) = ImgCell{iDx,jDx};
               end
           end
       end
       function obj = load(obj)
           if ~isprop(obj,'FieldStrength')
               disp('FieldStrength has to be defined first.')
               return
           else
               switch obj.FieldStrength
                   case 7
                        matched = obj.DicomTree.search('target','instance','query',...
                            @(inst,ser,stu) ser.SeriesNumber==obj.SeriesNum);
                        %   Load dicom images based on series number
                        obj.DicomImg = cell(numel(matched),1);
                        %   Load dicom images into cells
                        obj.ImgRaw = zeros(550,550,numel(matched));
                        %   These are collated images
                        obj.Img3D_TimeSeries = zeros(110,110,25,numel(matched));
                        % Might need to check if the size is indeed 550 by 550
                        for iDx = 1:numel(matched)
                            obj.DicomImg{iDx} =  Spectro.dicomImage({matched(iDx).Filename});
                            obj.ImgRaw(:,:,iDx) = obj.DicomImg{iDx}.image';
                            obj.Img3D_TimeSeries(:,:,:,iDx) = obj.SortGluedImage(obj.DicomImg{iDx}.image');
                        end
                        obj.Img3D_TimeSeries = squeeze(obj.Img3D_TimeSeries);
                        % In case of loading a M0 image. 
                   case 3
                       matched = obj.DicomTree.search('target','series','query',@(varargin) true);
                       % In the 3T case, the DicomTree includes only one
                       % folder.
                       obj.DicomImg = cell(numel(matched.instance),1);
                       %   Load dicom images into cells
                       obj.ImgRaw = zeros(550,550,numel(matched.instance));
                       %   These are collated images
                       obj.Img3D_TimeSeries = zeros(110,110,25,numel(matched.instance));
                       for iDx = 1:numel(matched.instance)
                            obj.DicomImg{iDx} =  Spectro.dicomImage({matched.instance(iDx).Filename});
                            obj.ImgRaw(:,:,iDx) = obj.DicomImg{iDx}.image';
                            obj.Img3D_TimeSeries(:,:,:,iDx) = obj.SortGluedImage(obj.DicomImg{iDx}.image');
                       end
                       obj.Img3D_TimeSeries = squeeze(obj.Img3D_TimeSeries);
                       % In case of loading a M0 image. 
               end
           end
       end
       function obj = Calc_Difference_Img(obj)
           Num_of_Meas = size(obj.Img3D_TimeSeries,4);
           Tag_Indices = 1:2:Num_of_Meas;
           Control_Indices = 2:2:Num_of_Meas;
           obj.DiffImgTimeSeries = obj.Img3D_TimeSeries(:,:,:,Control_Indices) -...
               obj.Img3D_TimeSeries(:,:,:,Tag_Indices);
       end      
       
       function obj = Calc_Difference_Img_Total(obj)
            obj.DiffImgTimeSeriesTotal = squeeze(sum(obj.DiffImgTimeSeries,4));
       end
       
       function DiffImgTimeSeries = Get_DiffTimeSeries(obj)
           DiffImgTimeSeries = obj.DiffImgTimeSeries;
       end
       
       function Get_Brain_Mask(obj)
           Num_of_Meas = size(obj.Img3D_TimeSeries,4);
           Control_Indices = 2:2:Num_of_Meas;
           RawImg_Sum = squeeze(sum(obj.Img3D_TimeSeries(:,:,Control_Indices),4));
           % sum of all the control images
           input = '/Users/ytong/Documents/Data/MaskTemp/EPI.nii';
           niftiwrite(RawImg_Sum,input);
           output = '/Users/ytong/Documents/Data/MaskTemp/EPImasked.nii';
           command_str = sprintf('bet %s %s -m',input,output);
           unix(command_str);
       end

   end
end