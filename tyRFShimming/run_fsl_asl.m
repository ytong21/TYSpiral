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
system(sprintf('cd %s',pathname_in));

%% Running fslreorient2std on all .nii.gz images
nii_all = dir(strcat(pathname_in,'/*nii.gz'));
disp('Running fslreorient2std on all .nii.gz images...');
for iDx = 1:numel(nii_all)
    reorient_cmd = sprintf('fslreorient2std %s %s',nii_all(iDx).name,nii_all(iDx).name);
    RunCommand(reorient_cmd);
end
%% Running bet and fast for T1 images (3T only)
if FieldStrength == 3
    disp('Generating GM mask from the T1 image...');
    bet_cmd = 'bet images_003_t1mprax1mmiso32chv21001 t1_brain -m -o -R -f 0.4';
    RunCommand(bet_cmd);
    fast_cmd = 'fast -B -t 1 -o t1_brain t1_brain';
    RunCommand(fast_cmd);
    threshold_cmd = strcat('fslmaths t1_brain_pve_1 -thr',...
        ' 0.5 t1_brain_GM_thresh');
    RunCommand(threshold_cmd);
    createmask_cmd = strcat('fslmaths t1_brain_GM_thresh ',...
        '-bin t1_brain_GM_thresh_mask');
    RunCommand(createmask_cmd);
end
%% Prepare field maps (7T only)
% While the dicom images from the field mapping sequence has 5 series 
% (2 amp, 2 phase, and 1 phase diff)
% Nifti images have only 3 series (1 amp with 2 vols, 1 phase with 2 vols, and 1 phase diff)
% First use fslsplit to separate the 2 amp images.

if FieldStrength == 7
    struct_field_map_amp_phase = dir(strcat(pathname_in,'/*dtfieldmapptx7tiso4mmtransR1001*nii.gz'));  
    struct_field_map_phasediff = dir(strcat(pathname_in,'/*dtfieldmapptx7tiso4mmtransR1001*nii.gz'));
    name_field_map_amp_2vols = struct_field_map_amp_phase(1).name;
    split_cmd = strcat('fslsplit ',name_field_map_amp_2vols, ' fieldmap_');
    % according to fsl website, "pick one with better contrast".
    % The second volume (vol0001) has better contrast
end
%% Generating perfusion-weighted images...
PCASL = dir(strcat(pathname_in,'/*PCASL*nii.gz'));
mcflirt_cmd = cell(numel(PCASL),1);
substraction_cmd = cell(numel(PCASL),1);
copygeom1 = cell(numel(PCASL),1);   copygeom2 = cell(numel(PCASL),1);
std_cmd = cell(numel(PCASL),1);     quant_cmd = cell(numel(PCASL),1);
disp('Generating perfusion-weighted images...');
M0struct = dir(strcat(pathname_in,'/*M0.nii.gz'));
for iDx = 1:numel(PCASL)
    cell_temp = strsplit(PCASL(iDx).name,'.');
    filename_wo_ext = cell_temp{1};
    mcflirt_cmd{iDx} = sprintf('mcflirt -in %s -refvol %s',PCASL(iDx).name,M0struct(1).name);
    RunCommand(mcflirt_cmd{iDx});
    substraction_cmd{iDx} = sprintf(['asl_file --data=%s --ntis=1 --ibf=tis --iaf=tc '...
        '--diff --out=%s mean=%s'],strcat(filename_wo_ext,'_mcf'),strcat(filename_wo_ext,'_diff'),...
        strcat(filename_wo_ext,'_diff_mean'));
    RunCommand(substraction_cmd{iDx})
    copygeom1{iDx} = sprintf('fslcpgeom %s %s',PCASL(iDx).name,strcat(filename_wo_ext,'_diff'));
    copygeom2{iDx} = sprintf('fslcpgeom %s %s',PCASL(iDx).name,strcat(filename_wo_ext,'_diff_mean'));
    RunCommand(copygeom1{iDx});RunCommand(copygeom2{iDx});
    std_cmd{iDx} = sprintf('fslmaths %s -Tstd %s',strcat(filename_wo_ext,'_diff'),...
        strcat(filename_wo_ext,'_std'));
    RunCommand(std_cmd{iDx});
    
    % Assign different analysis parameters based on filename
    cell_temp2 = strsplit(filename_wo_ext,'_');
    filename_wo_prefix = cell_temp2{end};
    switch filename_wo_prefix
        case {'toep2dPCASLmatch7T','toep2dPCASLmatch7TMatchFlip'}
            label_dur = 1;
        case {'toep2dPCASL3Ttagging','toep2dPCASL3Toptimized'}
            label_dur = 1.8;
    end
    PLD = 1;
    TI = label_dur + PLD;        
    
    % Using oxford_asl for quantification. M0 image used as reference
    quant_cmd{iDx} = sprintf(['oxford_asl -i %s --ntis=1 --ibf=tis --iaf=tc '...
        '--tis=%s --bolus=%s --casl --slicedt=0.0425 --cmethod=voxel '...
        '-c %s'], ...
        strcat(filename_wo_ext,'_mcf'),num2str(TI),num2str(label_dur),M0struct(1).name);
end


end
function RunCommand(command)
    [status,cmdout] = system(command,'-echo');
    if status~=0
       error('Status = %i.\nCmdout=\n%s',status,cmdout)
    end
end