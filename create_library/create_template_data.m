%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create meshes, headmodels and sourcemodels for MEG sourcespace analysis
% from template MRIs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set paths to MRI library and path to new libary
path_to_MRI_library = '/Users/44737483/Documents/scripts_mcq/MRI_database/';
path_to_new_library = '/Users/44737483/Documents/scripts_mcq/new_HCP_library_for_MEMES/';

% Get subject names from MRI database
try
    cd(path_to_MRI_library);
    % Get a list of all files and folders in this folder.
    files = dir(path_to_MRI_library);
    files(1:2) = [];
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    
    % Now these names to a variable called subject
    subject = [];
    
    for sub = 1 : length(subFolders)
        subject{sub} = subFolders(sub).name;
    end
    
    fprintf('%d subjects found in the MRI library: from %s to %s\n',...
        length(subject),subject{1}, subject{end});
    
catch
    warning('Something is wrong with your MRI library... Check the path!\n');
end


%%

% For every subject
for sub = 1:length(subject)
    
    %% Set up directory to hold info for each sub
    dir_for_loop = [path_to_new_library subject{sub}];
    
    mkdir(dir_for_loop);
    cd(dir_for_loop);
    
    %% Load MRI (THIS WILL BE MRI LIBRARY SPECIFIC)
    mri_to_use = [path_to_MRI_library subject{sub}...
        '/MEG/anatomy/T1w_acpc_dc_restore.nii.gz'];
    
    mri_orig = ft_read_mri(mri_to_use);
    mri_orig            = ft_convert_units(mri_orig,'mm'); % in mm
    mri_orig.coordsys = 'ras';
    
    %% Transform MRI based on fiducials
    % Manually mark fiducials
    cfg                         = [];
    cfg.headshape               = headshape_downsampled;
    cfg.viewmode                = 'ortho';
    cfg.method                  = 'interactive';
    cfg.coordsys                = 'bti';
    [mri_realigned]             = ft_volumerealign(cfg, mri_orig);
    
    % Also realign using ICP if necessary
    cfg                         = [];
    cfg.headshape = headshape_downsampled;
    cfg.method                  = 'headshape';
    cfg.interactive = 'no';
    cfg.coordsys                = 'bti';
    [mri_realigned]             = ft_volumerealign(cfg, mri_realigned);
    
    % Save realigned MRI to file
    save mri_realigned mri_realigned
    
    % Make a figure to check you've marked LPA and RPA the right way round(!)
    ft_determine_coordsys(mri_realigned, 'interactive', 'no');
    hold on; % add the subsequent objects to the figure
    drawnow; % workaround to prevent some MATLAB versions (2012b and 2014b) from crashing
    ft_plot_headshape(headshape_downsampled); view([180 0]);
    print('qc_mri_realigned','-dpng','-r100');
    
%     % Get marked fiducials from mri structure
%     vox_Nas = mri_realigned.cfg.fiducial.nas;
%     vox_Lpa = mri_realigned.cfg.fiducial.lpa;
%     vox_Rpa = mri_realigned.cfg.fiducial.rpa;
%     
%     % transformation matrix of individual MRI
%     vox2head = (mri_orig.transform);
%     
%     % transform voxel indices to MRI head coordinates
%     head_Nas          = ft_warp_apply(vox2head, vox_Nas, 'homogenous'); % nasion
%     head_Lpa          = ft_warp_apply(vox2head, vox_Lpa, 'homogenous'); % Left preauricular
%     head_Rpa          = ft_warp_apply(vox2head, vox_Rpa, 'homogenous'); % Right preauricular
%     
%     % Save marked fiducials for later
%     fids = [head_Nas;head_Lpa;head_Rpa];
%     save('fiducials.txt', 'fids','-ascii', '-double', '-tabs')
    
    %% Create scalp mesh
    fprintf('Creating facial mesh');
    % Segment the mri preserving scalp info
    cfg = [];
    cfg.output    = 'scalp';
    cfg.scalpsmooth = 5;
    cfg.scalpthreshold = 0.08;
    scalp  = ft_volumesegment(cfg, mri_realigned);
    
    % Create mesh out of scalp surface
    cfg = [];
    cfg.method = 'isosurface';
    cfg.numvertices = 10000;
    mesh = ft_prepare_mesh(cfg,scalp);
    mesh = ft_convert_units(mesh,'mm');
    
    % Plot figure of mesh
    figure;ft_plot_mesh(mesh); alpha 0.3;
    ft_plot_headshape(headshape_downsampled); view([0,0]);
    print('qc_mesh','-dpng','-r100');
    
    % Save to mesh library variable
    fprintf('Saving to the mesh\n')
    save mesh mesh
    
    % House-keeping
    close all force
    clear mesh mri_orig mri_realigned scalp
    
end






%% Load the segmented brain
fprintf('Loading the segmented MRI\n');

mri_file_without_head = ['/Users/44737483/Documents/scripts_mcq/'...
    'MRIDataBase_JohnRichards_USC/Children/Head/ANTS' subject{sub}...
    'Years_head_brain.nii.gz'];

%% Segment
cfg           = [];
cfg.output    = 'brain';
mri_segmented  = ft_volumesegment(cfg, mri_realigned);

%% Create singleshell headmodel
fprintf('Creating singleshell headmodel\n');
cfg = [];
cfg.tissue = 'brain';
cfg.method='singleshell';
headmodel = ft_prepare_headmodel(cfg, mri_segmented); % in mm, create headmodel

% Create Figure
figure;ft_plot_vol(headmodel);
ft_plot_mesh(mesh); alpha 0.3; view([0,0]);
print('qc_headmodel','-dpng','-r100');

% Save headmodel
save headmodel headmodel;

%% Create sourcemodel

sourcemodel_mm = [10 8 5];

for size = 1:length(sourcemodel_mm);
    
    fprintf('Creating %dmm sourcemodel\n');
    
    % create the subject specific grid, using the template grid that has just been created
    cfg                = [];
    cfg.grid.warpmni   = 'yes';
    cfg.grid.resolution = sourcemodel_mm(size);
    cfg.grid.nonlinear = 'yes'; % use non-linear normalization
    cfg.mri            = mri_realigned;
    cfg.grid.unit      ='mm';
    cfg.inwardshift = '-1.5';
    grid               = ft_prepare_sourcemodel(cfg);
    
    % Transform based on fiducials
    %         grid.pos = ft_warp_apply(inv(mri_orig.transform),grid.pos);
    %         grid.pos = ft_warp_apply(mri_realigned.transform,grid.pos);
    
    % Save
    sourcemodel3d = grid;
    save(sprintf('sourcemodel3d_%dmm',sourcemodel_mm(size)),'sourcemodel3d');
    
    % Create figure and save
    figure;ft_plot_mesh(grid.pos(grid.inside,:)); view([0,0]);
    print(sprintf('qc_sourcemodel3d_shape_%dmm',sourcemodel_mm(size)),'-dpng','-r100');
    
    figure;
    figure; ft_plot_vol(headmodel_singleshell); alpha 0.3;
    ft_plot_mesh(mesh); alpha 0.3;
    ft_plot_mesh(grid.pos(grid.inside,:)); view([0,0]);
    
    print(sprintf('qc_sourcemodel3d_%dmm',sourcemodel_mm(size)),'-dpng','-r100');
    
    % Clear for next loop
    clear grid sourcemodel3d
end

% Clear for next loop
clear mesh headmodel mri_orig mri_realigned headmodel_singleshell scalp

close all

fprintf('Finished... CHECK for quality control\n');

end


%% Example call to child_MEMES

dir_name    = '/Users/44737483/Documents/scripts_mcq/child_test/2913/resting_state'
elpfile     = '/Users/44737483/Documents/scripts_mcq/child_test/2913/2913_ES_ME125_2018_02_24.elp';
hspfile     = '/Users/44737483/Documents/scripts_mcq/child_test/2913/2913_ES_ME125_2018_02_24.hsp';
confile     = '/Users/44737483/Documents/scripts_mcq/child_test/2913/resting_state/2913_ES_ME125_2018_02_24_B2.con';
mrkfile     = '/Users/44737483/Documents/scripts_mcq/child_test/2913/2913_ES_ME125_2018_02_24_INI.mrk';

path_to_MRI_library = '/Users/44737483/Documents/scripts_mcq/MRIDataBase_JohnRichards_USC/database_for_MEMES/';

child_MEMES(dir_name,elpfile,hspfile,confile,mrkfile,path_to_MRI_library,'')

%%

