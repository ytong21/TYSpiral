function run_fsl_asl(pathname_in)
PATH = getenv('PATH');
setenv('PATH', [PATH ':/usr/local/bin:/usr/local/fsl/bin']);
setenv('FSLOUTPUTTYPE','NIFTI_GZ')
setenv('FSLMULTIFILEQUIT','TRUE')
setenv('FSLTCLSH','/usr/local/fsl/bin/fsltclsh')
setenv('FSLWISH','/usr/local/fsl/bin/fslwish')
setenv('FSLDIR','/usr/local/fsl')
% Adding path for fsl

    [~,FileName,~] = fileparts(pathname_in);
    PathNameSplitted = strsplit(FileName,'_');
    if (numel(PathNameSplitted) == 5) && (strcmp(PathNameSplitted{2},'F7T'))
        FieldStrength = 7;
    elseif (numel(PathNameSplitted) == 4) && (strcmp(PathNameSplitted{1},'F3T'))
    	FieldStrength = 3;
    else
    	disp('Path format unregocnised.');
        return
    end
%% cd to the folder containing the data.
%system(sprintf('cd %s',pathname_in));
cd(pathname_in);

%% Get subject info
SubjectInfoSingle = FindSubject(FileName, FieldStrength);
%% Running fslreorient2std on all .nii.gz images
% nii_all = dir(strcat(pathname_in,'/*nii.gz'));
% disp('Running fslreorient2std on all .nii.gz images...');
% for iDx = 1:numel(nii_all)
%     reorient_cmd = sprintf('fslreorient2std %s %s',nii_all(iDx).name,nii_all(iDx).name);
%     RunCommand(reorient_cmd);
% end
M0struct = dir(strcat(pathname_in,'/*M0.nii.gz'));
if FieldStrength == 3
    split_cmd = sprintf('fslsplit %s %s',M0struct(1).name,'M0_');
    RunCommand(split_cmd)
    M0_filename = 'M0_0000.nii.gz';
elseif FieldStrength == 7
    M0_filename = M0struct(1).name;
end
%% Running bet and fast for T1 images (3T only)
% if FieldStrength == 3
%     disp('Generating GM mask from the T1 image...');
%     bet_cmd = 'bet images_003_t1mprax1mmiso32chv21001 t1_brain -m -o -R -f 0.4';
%     RunCommand(bet_cmd);
%     fast_cmd = 'fast -B -t 1 -o t1_brain t1_brain';
%     RunCommand(fast_cmd);
% end
% %% Prepare field maps (7T only)
% % While the dicom images from the field mapping sequence has 5 series 
% % (2 amp, 2 phase, and 1 phase diff)
% % Nifti images have only 3 series (1 amp with 2 vols, 1 phase with 2 vols, and 1 phase diff)
% % First use fslsplit to separate the 2 amp images.
% 
% if FieldStrength == 7
%     struct_field_map_amp_phase = dir(strcat(pathname_in,'/*dtfieldmapptx7tiso4mmtransR1001*nii.gz'));  
%     struct_field_map_phasediff = dir(strcat(pathname_in,'/*dtfieldmapptx7tiso4mmtransR1001*nii.gz'));
%     name_field_map_amp_2vols = struct_field_map_amp_phase(1).name;
%     split_cmd = strcat('fslsplit ',name_field_map_amp_2vols, ' fieldmap_');
%     % according to fsl website, "pick one with better contrast".
%     % The second volume (vol0001) has better contrast
% end
%% Creat GM Mask
 PCASL = dir(strcat(pathname_in,'/*PCASL*nii.gz'));
% Register the M0 to the t1 structural
if FieldStrength == 7
    t1_file = sprintf('/Users/ytong/Documents/Data/3T/F3T_2013_50_%s/%s',...
        num2str(SubjectInfoSingle.num_3T),'images_003_t1mprax1mmiso32chv21001.nii.gz');
elseif FieldStrength == 3
    t1_file = 'images_003_t1mprax1mmiso32chv21001.nii.gz';
end

% epi_reg_cmd = sprintf(strcat('epi_reg --epi=%s --t1=',...
%     '%s --t1brain=t1_brain.nii.gz --wmseg=t1_brain_pve_2.nii.gz --out=EPI_T1 -v'),...
%     M0_filename,t1_file);
% RunCommand(epi_reg_cmd);
% convert_cmd = 'convert_xfm -omat T1_to_M0.mat -inverse EPI_T1.mat';
% RunCommand(convert_cmd);
% flirt_cmd = sprintf('flirt -in t1_brain_pve_1 -ref %s -applyxfm -init T1_to_M0.mat -out GM_pve_ASL_space',...
%      M0_filename);
% RunCommand(flirt_cmd);

%% Generating perfusion-weighted images...
mcflirt_cmd = cell(numel(PCASL),1);
substraction_cmd = cell(numel(PCASL),1);
copygeom1 = cell(numel(PCASL),1);   %copygeom2 = cell(numel(PCASL),1);
std_cmd = cell(numel(PCASL),1);     quant_cmd = cell(numel(PCASL),1);
tSNR_cmd = cell(numel(PCASL),1);    mean_cmd = cell(numel(PCASL),1);

% Creating GM mask in ASL space
% disp('Creating GM mask in ASL space ...')
%     threshold_cmd = strcat('fslmaths GM_pve_ASL_space -thr',...
%         ' 0.4 GM_pve_ASL_space_thresh');
%     RunCommand(threshold_cmd);
%     createmask_cmd = strcat('fslmaths GM_pve_ASL_space_thresh',...
%         ' -bin GM_pve_ASL_space_thresh_mask');
%     RunCommand(createmask_cmd);

disp('Generating perfusion-weighted images...');
for iDx = 1:numel(PCASL)
    cell_temp = strsplit(PCASL(iDx).name,'.');
    filename_wo_ext = cell_temp{1};
    % Make a subdirectory for each ASL series
    %mkdir(filename_wo_ext);
    %copyfile(PCASL(iDx).name,filename_wo_ext);
    cd(filename_wo_ext);
%     mcflirt_cmd{iDx} = sprintf('mcflirt -in %s -refvol %s',...
%         PCASL(iDx).name,strcat('../',M0_filename));
%     RunCommand(mcflirt_cmd{iDx});
%     substraction_cmd{iDx} = sprintf(['asl_file --data=%s --ntis=1 --ibf=tis --iaf=tc '...
%         '--diff --out=%s'],strcat(filename_wo_ext,'_mcf'),strcat(filename_wo_ext,'_diff'));
%     RunCommand(substraction_cmd{iDx})
%     copygeom1{iDx} = sprintf('fslcpgeom %s %s',PCASL(iDx).name,strcat(filename_wo_ext,'_diff'));
    %copygeom2{iDx} = sprintf('fslcpgeom %s %s',PCASL(iDx).name,strcat(filename_wo_ext,'_diff_mean'));
%     RunCommand(copygeom1{iDx});%RunCommand(copygeom2{iDx});
%     mean_cmd{iDx} = sprintf('fslmaths %s -Tmean %s',strcat(filename_wo_ext,'_diff'),...
%         strcat(filename_wo_ext,'_mean'));
%     RunCommand(mean_cmd{iDx});
%     std_cmd{iDx} = sprintf('fslmaths %s -Tstd %s',strcat(filename_wo_ext,'_diff'),...
%         strcat(filename_wo_ext,'_std'));
%     RunCommand(std_cmd{iDx});
%     tSNR_cmd{iDx} = sprintf('fslmaths %s -div %s %s',strcat(filename_wo_ext,'_mean'),...
%         strcat(filename_wo_ext,'_std'), strcat(filename_wo_ext,'_tSNR'));
%     RunCommand(tSNR_cmd{iDx});
    % Assign different analysis parameters based on filename
    [label_dur, TI, LabEff] = DetermineParams(filename_wo_ext,FieldStrength,SubjectInfoSingle);
    % Using oxford_asl for quantification. M0 image used as reference
    quant_cmd{iDx} = sprintf(['oxford_asl -i %s '...
        '--tis=%s --bolus=%s --casl --slicedt=0.0425 --cmethod=voxel '...
        '-c %s --alpha %s --tr 3.2'], ...
        strcat(filename_wo_ext,'_diff'),num2str(TI),num2str(label_dur),strcat('../',M0_filename),...
        num2str(LabEff));
    RunCommand(quant_cmd{iDx});
%         applymask_cmd = sprintf('fslmaths %s -div %s %s',strcat(filename_wo_ext,'_tSNR'),...
%         strcat('../','GM_pve_ASL_space_thresh_mask'), strcat(filename_wo_ext,'_tSNR_masked'));
%     RunCommand(applymask_cmd);
    cd ..
end



end
%%
function RunCommand(command)
    [status,cmdout] = system(command,'-echo');
    if status~=0
       error('Status = %i.\nCmdout=\n%s',status,cmdout)
    end
end
function SubjectInfoSingle = FindSubject(FileName, FieldStrength)
    load SubjectInfo.mat SubjectInfo
    PathNameSplitted = strsplit(FileName,'_');
    SubjectNum = str2double(PathNameSplitted{end});

    for iDx = 1:numel(SubjectInfo)
        if FieldStrength == 3
            NumFromStorage = SubjectInfo(iDx).num_3T;
        elseif FieldStrength == 7
            NumFromStorage = SubjectInfo(iDx).num_7T;
        end
        if SubjectNum == NumFromStorage
            SubjectInfoSingle = SubjectInfo(iDx);
            return
        end
    end
end
function [label_dur, TI, LabEff] = DetermineParams(filename_wo_ext,FieldStrength,SubjectInfoSingle)
    cell_temp2 = strsplit(filename_wo_ext,'_');
    filename_wo_prefix = cell_temp2{end};
    PLD = 1;
    if FieldStrength == 3
        switch filename_wo_prefix
            case 'toep2dPCASLmatch7T'
                label_dur = 1;      LabEff = 0.8662;
            case 'toep2dPCASLmatch7TMatchFlip'
                label_dur = 1;      LabEff = SubjectInfoSingle.Prisma_matched_7TFA_eff;
            case {'toep2dPCASL3Ttagging','toep2dPCASL3Toptimized'}
                label_dur = 1.8;    LabEff = 0.9370;
        end
    elseif FieldStrength == 7
        if  contains(filename_wo_prefix,'43')
                label_dur = 1;      LabEff = SubjectInfoSingle.Magnetom_eff(1);
        elseif (contains(filename_wo_prefix,'53') ||    contains(filename_wo_prefix,'52'))
                label_dur = 1;      LabEff = SubjectInfoSingle.Magnetom_eff(2);
        elseif contains(filename_wo_prefix,'23')   
                label_dur = 1;      LabEff = SubjectInfoSingle.Magnetom_eff(3);
        elseif contains(filename_wo_prefix,'28')   
                label_dur = 1;      LabEff = SubjectInfoSingle.Magnetom_eff(4); 
        end
    else
        disp('Field strength not recognised.')
        return
    end
    TI = label_dur + PLD;
end