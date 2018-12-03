cd('/Users/44737483/Dropbox/MEG_data/1412'); load('trans_matrix.mat');
load('sens.mat'); load('mri_realigned.mat'); load('headmodel_singleshell.mat');

fif_file = '/Users/44737483/Documents/scripts_mcq/MEMES/elekta/rs_asd_1413_aliens_quat_tsss.fif';

mesh = ft_read_headshape({[path_to_MRI_library '917255/MEG/anatomy/917255.L.midthickness.4k_fs_LR.surf.gii'],[path_to_MRI_library '917255/MEG/anatomy/917255.R.midthickness.4k_fs_LR.surf.gii']});
mesh = ft_convert_units(mesh,'cm');

% Get original MRI file
mri_file = [path_to_MRI_library '917255' '/MEG/anatomy/T1w_acpc_dc_restore.nii.gz'];
mri_orig = ft_read_mri(mri_file); mri_orig = ft_convert_units(mri_orig,'cm');

% Check original MRI and mesh match
ft_determine_coordsys(mri_orig,'interactive','no'); hold on;
ft_plot_mesh(mesh,'facealpha',0.8); camlight; hold on;
%ft_plot_mesh(mesh,'facealpha',0.5); hold on;  

%% Transform 1 (MESH --> coreg via manual marking of fiducial points)

mesh.pos = ft_warp_apply(inv(mri_orig.transform),mesh.pos);
mesh.pos = ft_warp_apply(mri_realigned.transform,mesh.pos);


ft_determine_coordsys(mri_realigned,'interactive','no'); hold on;
ft_plot_mesh(mesh,'facealpha',0.8); camlight; hold on;

%% Transform 2 (MESH --> coreg via ICP adjustment)

%mesh.pos = ft_warp_apply(trans_matrix,mesh.pos);

%mri_realigned2 = ft_transform_geometry(trans_matrix,mri_realigned);

%ft_determine_coordsys(mri_realigned2,'interactive','no'); hold on;
%ft_plot_mesh(mesh,'facealpha',0.8); camlight; hold on;

%% Final figure

headshape = ft_read_headshape(fif_file); headshape = ft_convert_units(headshape,'cm');

%ft_determine_coordsys(mri_realigned,'interactive','no'); hold on;
figure;
ft_plot_sens(sens); 
ft_plot_headshape(headshape);
%ft_plot_vol(headmodel_singleshell,'facealpha',0.8); camlight; hold on;

ft_plot_mesh(mesh,'facealpha',0.8); camlight; hold on;


cd('/Users/44737483/Dropbox/MEG_data/1412'); save mesh mesh



