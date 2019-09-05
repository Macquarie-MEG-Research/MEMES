function child_MEMES_manual(dir_name,grad_trans,headshape_downsampled,...
    path_to_MRI_library, age_list)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MRI Estimation for MEG Sourcespace (MEMES) cutomised for child MEG.
%
% Written by Robert Seymour (Macquarie Univ Dept of Cognitive Science, July
% 2018). Some sub-functions written by Associate Professor Paul Sowman.
%
%%%%%%%%%%%
% Inputs:
%%%%%%%%%%
%
% - dir_name            = directory for saving
% - path_to_MRI_library = path to HCP MRI library

%%%%%%%%%%%%%%%%%%
% Variable Inputs:
%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%
% Outputs:
%%%%%%%%%%%
%
% - trans_matrix            = transformation matrix applied to headmodel
%                           and sourcemodel
% - sourcemodel3d           = 8mm sourcemodel warped to MNI space
% - headmodel               = singleshell headmodel (10000 vertices)

%%%%%%%%%%%%%%%%%%%%%
% Other Information:
%%%%%%%%%%%%%%%%%%%%%

% Example function call:
% child_MEMES(dir_name,elpfile,hspfile,confile,mrkfile,...
% path_to_MRI_library,'no','rot3dfit')

% The template
% MRIs (all ages 2.5 - 10.5) can be downloaded from:
% http://jerlab.psych.sc.edu/neurodevelopmentalmridatabase/
%
% John Richards (USC, USA) retains all copyrights to the templates.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(['\nThis is MEMES for child MEG data v0.1\n\nMake sure you have ',...
    'asked Robert for the mesh, headmodel and sourcemodel library\n\n']);
% Check inputs
disp('Performing input check');
assert(path_to_MRI_library(end) == '/','path_to_MRI_library needs to end with /\n');

% CD to right place
cd(dir_name); fprintf('\nCDd to the right place\n');

% Load the mesh
load([path_to_MRI_library age_list '/mesh.mat'],'mesh');

% Load fiducial information

load([path_to_MRI_library age_list '/fiducials.mat']);

[R,T,Yf,Err]                = rot3dfit(headshape_downsampled.fid.pos...
    ,fiducials);%calc rotation transform
meg2head_transm             = [[R;T]'; 0 0 0 1];%reorganise and make 4*4 transformation matrix


% [TR, TT, ER, t, info] = icp(fiducials',headshape_downsampled.fid.pos');
% head2head = [[TR TT]; 0 0 0 1];
% 
% % Realign mesh based on fids
% mesh2 = mesh;
% mesh2.pos = ft_warp_apply(meg2head_transm,mesh.pos);

cfg = [];
cfg.method = 'fiducial';
cfg.coordsys = 'ctf';
%cfg.fiducial.nas    = fiducials(:,1); 
%cfg.fiducial.lpa    = fiducials(:,2);
%cfg.fiducial.rpa    = fiducials(:,3);
mesh2 = ft_meshrealign(cfg, mesh);

mesh2 = mesh;
mesh2.pos = ft_warp_apply(inv(meg2head_transm),mesh.pos)

figure; 
ft_plot_mesh(mesh2,'facecolor','r','facealpha',0.3);
ft_plot_mesh(ft_warp_apply(inv(meg2head_transm),fiducials),'vertexcolor','y','vertexsize',20);
ft_plot_headshape(headshape_downsampled);





%% Now refine coreg with ICP
% Create scalp mesh
% Segment the mri preserving scalp info
cfg             = [];
cfg.output      = 'scalp';
cfg.scalpsmooth = 5;
cfg.scalpthreshold = 0.08;
scalp           = ft_volumesegment(cfg, mri_realigned);

% Create mesh out of scalp surface
cfg = [];
cfg.method = 'isosurface';
cfg.numvertices = 10000;
mesh = ft_prepare_mesh(cfg,scalp);
mesh = ft_convert_units(mesh,'mm');

% Number of iterations for the ICP algorithm
numiter = 60;

% Perform ICP
[R, t, err, dummy, ~] = icp(mesh.pos', headshape_downsampled.pos', ...
    numiter, 'Minimize', 'plane', 'Extrapolation', true,'WorstRejection', 0.1);

% Get the transformation matrix
trans_matrix = inv([real(R) real(t);0 0 0 1]);

%% Load 

%% Load the segmented brain
fprintf('Loading the segmented MRI\n');

mri_file_without_head = ['/Users/44737483/Documents/scripts_mcq/'...
    'MRIDataBase_JohnRichards_USC/Children/Head/ANTS' age_list...
    'Years_head_brain.nii.gz'];

mri_orig2 = ft_read_mri(mri_file_without_head);
mri_orig2 = ft_convert_units(mri_orig2,'mm'); % in mm
mri_orig2.coordsys = 'ctf';

% Bit of jiggery-pokery to trick FT into thinking you've segmented the MRI
% using SPM...
mri_segmented = mri_orig2;
mri_segmented.brain = mri_segmented.anatomy;
mri_segmented.brain(find(mri_segmented.brain>0)) = 1;
mri_segmented.brain = logical(mri_segmented.brain);

%% Create singleshell headmodel
fprintf('Creating singleshell headmodel\n');
cfg = [];
cfg.tissue = 'brain';
cfg.method='singleshell';
headmodel = ft_prepare_headmodel(cfg, mri_segmented); % in mm, create headmodel







%%
% Create figure to check headodel and sourcemodel match
figure;
ft_plot_vol(headmodel,  'facecolor', 'cortex', 'edgecolor', 'none');
alpha 0.4; camlight;
ft_plot_mesh(sourcemodel3d.pos(sourcemodel3d.inside,:),'vertexsize',5);
view([0 0]);

view_angle = [0 90 180 270];

% Create figure to show final coregiration
figure; hold on;
for rep = 1:4
    subplot(2,2,rep);
    ft_plot_vol(headmodel,  'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.6; camlight;
    ft_plot_mesh(sourcemodel3d.pos(sourcemodel3d.inside,:),'vertexsize',3);
    ft_plot_sens(grad_trans, 'style', 'r*')
    ft_plot_headshape(headshape_downsampled) %plot headshape
    view([view_angle(rep),0]);
    ft_plot_mesh(mesh_spare,'facecolor',[238,206,179]./255,'EdgeColor','none','facealpha',0.5);
    camlight; lighting phong; material dull;
end

print('coregistration_volumetric_quality_check','-dpng','-r100');

end

