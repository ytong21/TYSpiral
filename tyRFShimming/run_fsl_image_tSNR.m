function run_fsl_image_tSNR(pathname_in)
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
    if (numel(PathNameSplitted) == 4) && (strcmp(PathNameSplitted{1},'F7T'))
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
%SubjectInfoSingle = FindSubject(FileName, FieldStrength);
M0 = 'M0.nii.gz';

%% Build the command for fslmaths
VERSE_dir = dir(strcat(pathname_in,'/*VERSE28'));
cd(VERSE_dir.name);
control_series = dir('*odd_masked.nii.gz');
std_command = sprintf('fslmaths %s -Tstd %s',control_series.name,...
    'VERSE_control_std');
[~,~] = RunCommand(std_command);
copy_geom_command = sprintf('fslcpgeom %s %s -d',strcat('..',filesep,'EPI_unwarp',filesep,M0),...
    'VERSE_control_std');
[~,~] = RunCommand(copy_geom_command);
% threshold_cmd = strcat('fslmaths GM_pve_ASL_space -thr',...
%         ' 0.4 GM_pve_ASL_space_thresh');
% RunCommand(threshold_cmd);
% createmask_cmd = strcat('fslmaths GM_pve_ASL_space_thresh',...
%         ' -bin GM_pve_ASL_space_thresh_mask');
%     RunCommand(createmask_cmd);
div_command = sprintf('fslmaths %s -div %s %s',strcat('..',filesep,'EPI_unwarp',filesep,M0),...
    'VERSE_control_std','VERSE_img_tSNR');
[~,~] = RunCommand(div_command);
% mean_command = sprintf('fslstats %s -M','VERSE_img_tSNR');
% [~,tSNR_mean] = RunCommand(mean_command);
% std_command = sprintf('fslstats %s -S','VERSE_img_tSNR');
% [~,tSNR_std] = RunCommand(std_command);
% voxel_command = sprintf('fslstats %s -M','VERSE_img_tSNR');
% [~,voxels] = RunCommand(Voxel_command);
%%
function [status,cmdout] = RunCommand(command)
    [status,cmdout] = system(command,'-echo');
    if status~=0
       error('Status = %i.\nCmdout=\n%s',status,cmdout)
    end
end
end