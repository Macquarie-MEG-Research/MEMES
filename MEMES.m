%% Create a library of MRIs

subject = {'100307';'102816';'104012';'105923';'106521';'108323';'109123';...
    '111514';'112920';'113922';'116524';'116726';'125525';'133019';'140117';...
    '146129';'149741';'151526';'153732';'154532'};

mri_realigned_library = [];
scalp_library = [];
mesh_library = [];

for i = 1:length(subject)

mri_file = ['/Users/44737483/Documents/scripts_mcq/mri_less/' subject{i} '/MEG/anatomy/T1w_acpc_dc_restore.nii.gz'];

% Load in MRI
mri_orig                    = ft_read_mri(mri_file); % in mm, read in mri from DICOM
mri_orig = ft_convert_units(mri_orig,'cm'); mri_orig.coordsys = 'neuromag';

% Give rough estimate of fiducial points
cfg                         = [];
cfg.method                  = 'interactive';
cfg.viewmode                = 'ortho';
cfg.coordsys                = 'neuromag';
[mri_realigned]             = ft_volumerealign(cfg, mri_orig);

mri_realigned_library{i} = mri_realigned;

% check that the MRI is consistent after realignment
ft_determine_coordsys(mri_realigned, 'interactive', 'no');
hold on; % add the subsequent objects to the figure
drawnow; % workaround to prevent some MATLAB versions (2012b and 2014b) from crashing
ft_plot_headshape(headshape_downsampled);
drawnow;



cfg = [];
cfg.output    = 'scalp';
cfg.scalpsmooth = 5;
cfg.scalpthreshold = 0.08;
scalp  = ft_volumesegment(cfg, mri_realigned);

scalp_library{i} = scalp; 

%% Create mesh out of scalp surface
cfg = [];
cfg.method = 'projectmesh';
cfg.numvertices = 50000;
mesh = ft_prepare_mesh(cfg,scalp);
mesh = ft_convert_units(mesh,'cm');

figure;ft_plot_mesh(mesh,'facecolor',[238,206,179]./255,'EdgeColor','none','facealpha',0.8); hold on;
camlight; hold on; view([-180,-10]); drawnow;

mesh_library{i} = mesh;

clear mesh scalp mri_realigned

end

cd('/Users/44737483/Documents/scripts_mcq/mri_less');

save mri_realigned_library scalp_library mesh_library








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% MRI Estimation for MEG Sourcespace (MEMES)

% Add functions from coreg (include as subfunctions later)
addpath('/Users/44737483/Documents/scripts_mcq/alien');

elpfile = '/Users/44737483/Documents/mcq_data/2704/meg/2704_AT_ME160_2017_09_27.elp';
hspfile = '/Users/44737483/Documents/mcq_data/2704/meg/2704_AT_ME160_2017_09_27.hsp';
confile = '/Users/44737483/Documents/mcq_data/2704/meg/run-alien/2704_AT_ME160_2017_09_27_aliens.con';
mrkfile = '/Users/44737483/Documents/mcq_data/2704/meg/run-rs/2704_AT_ME160_2017_09_27_rs_PRE.mrk';

% Get Polhemus Points
[shape] = parsePolhemus(elpfile,hspfile);

% Read the grads from the con file
grad_con                    = ft_read_sens(confile); %in cm, load grads

% Read the mrk file
mrk                         = ft_read_headshape(mrkfile,'format','yokogawa_mrk');
markers                     = mrk.fid.pos([2 3 1 4 5],:);%reorder mrk to match order in shape
[R,T,Yf,Err]                = rot3dfit(markers,shape.fid.pnt(4:end,:));%calc rotation transform
meg2head_transm             = [[R;T]'; 0 0 0 1];%reorganise and make 4*4 transformation matrix

% Transform sensors based on the MRKfile
grad_trans      = ft_transform_geometry_PFS_hacked(meg2head_transm,grad_con); %Use my hacked version of the ft function - accuracy checking removed not sure if this is good or not
grad_trans.fid  = shape; %add in the head information
% Rotate about z-axis
rot180mat       = rotate_about_z(180);
grad_trans      = ft_transform_geometry(rot180mat,grad_trans);

% Get headshape downsampled to 100 points with facial info removed
headshape_downsampled = downsample_headshape_noface(hspfile,200,grad_trans);
% Rotate about z-axis
headshape_downsampled = ft_transform_geometry(rot180mat,headshape_downsampled);

%% Perform ICP

% Error term variable
error_term = zeros(1,length(scalp_library)); 
% Variable to hold the transformation matrices
trans_matrix_library = [];

for m = 1:length(scalp_library)

numiter = 50; disp(m);

% Perform ICP
[R, t, err, dummy, info] = icp(mesh_library{m}.pos', headshape_downsampled.pos', numiter, 'Minimize', 'plane', 'Extrapolation', true,'WorstRejection', 0.05);

% Add error to error_term
error_term(m) = err(end);

% Add transformation matrix to trans_matrix_library
trans_matrix_library{m} = inv([real(R) real(t);0 0 0 1]);

end

% Create figure to summarise the fits
figure; 

for i = 1:length(scalp_library)

mesh_spare = mesh_library{i};
mesh_spare.pos = ft_warp_apply(trans_matrix_library{i}, mesh_spare.pos);

subplot(4,5,i)
ft_plot_mesh(mesh_spare,'facecolor',[238,206,179]./255,'EdgeColor','none','facealpha',0.8); hold on;
camlight; hold on; view([-180,-10]);
title(error_term(i));
ft_plot_headshape(headshape_downsampled);
end

print('2704_with_hcp','-dpdf','-r200');

%% Use the best for to create a source model for MEG source analysis

winner = find(error_term == min(error_term));

% Create figure to show ICP fit
mesh_spare = mesh_library{winner};
mesh_spare.pos = ft_warp_apply(trans_matrix_library{winner}, mesh_spare.pos);

figure;ft_plot_mesh(mesh_spare,'facecolor',[238,206,179]./255,'EdgeColor','none','facealpha',0.8); hold on;
camlight; hold on; view([-180,-10]);
title(error_term(winner));
ft_plot_headshape(headshape_downsampled);

print('winning_sourcemodel','-dpdf','-r200');

% % Make fancy video
c = datestr(clock); %time and date

figure;
ft_plot_mesh(mesh_spare,'facecolor',[238,206,179]./255,'EdgeColor','none','facealpha',0.8); hold on;
camlight; hold on;
ft_plot_headshape(headshape_downsampled); title(sprintf('%s.   Error of ICP fit = %d' , c, error_term(winner)));
OptionZ.FrameRate=15;OptionZ.Duration=5.5;OptionZ.Periodic=true;
CaptureFigVid([0,0; 360,0], 'ICP_quality',OptionZ)

% Segment
cfg           = [];
cfg.output    = 'brain';
mri_segmented  = ft_volumesegment(cfg, mri_realigned_library{winner});

% Create singleshell headmodel
cfg = [];
cfg.method='singleshell';

headmodel_singleshell = ft_prepare_headmodel(cfg, mri_segmented); % in cm, create headmodel

% Apply transformation matrix
headmodel_singleshell.bnd.pos = ft_warp_apply(trans_matrix_library{winner},headmodel_singleshell.bnd.pos);

figure;ft_plot_headshape(headshape_downsampled) %plot headshape
ft_plot_sens(grad_trans, 'style', 'k*')
ft_plot_vol(headmodel_singleshell,  'facecolor', 'cortex', 'edgecolor', 'none');alpha 1; camlight
view([90,0]); title('After Coreg');
%print('headmodel_quality','-dpdf');


