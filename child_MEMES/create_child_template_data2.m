%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create meshes, headmodels and sourcemodels for MEG sourcespace analysis 
% from child templates (obtained from John Richards at USC).
%
% John Richards (USC, USA) retains all copyrights to the templates.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

list_of_ages = {'2-5','3-0','3-5','4-0','4-5','5-0'}

mesh_library = {};

for age = 1:length(list_of_ages)
    
    %% Set up directory to hold info for each age
    dir_for_loop = ['/Users/44737483/Documents/scripts_mcq/'...
        'MRIDataBase_JohnRichards_USC/database_for_MEMES_child/' list_of_ages{age}]
    
    mkdir(dir_for_loop);
    cd(dir_for_loop);
    
    %% Load MRI with face
    mri_file_with_head = ['/Users/44737483/Documents/scripts_mcq/'...
        'MRIDataBase_JohnRichards_USC/Children/Head/ANTS' list_of_ages{age}...
        'Years_head.nii'];
    
    mri_orig            = ft_read_mri(mri_file_with_head);
    mri_orig            = ft_convert_units(mri_orig,'mm'); % in mm
    
%     cfg         = [];
%     mri_orig    = ft_volumereslice(cfg, mri_orig);
    
    cfg         = [];
    cfg.method  = 'interactive';
    cfg.coordsys = 'ctf';
    mri_orig    = ft_volumerealign(cfg, mri_orig);
    
    % Fiducials on MRI
    fid_pnt = mri_orig.cfg.fiducial
    fid_vox = [fid_pnt.nas; fid_pnt.lpa; fid_pnt.rpa];
    
    % Represent fiducial points in the head coordinate system for MRI
    fiducials = ft_warp_apply(mri_orig.transform, fid_vox, 'homogeneous');
    
    % Make a figure to check you've marked LPA and RPA the right way round(!)
    ft_determine_coordsys(mri_orig, 'interactive', 'no');
    hold on; % add the subsequent objects to the figure
    drawnow; % workaround to prevent some MATLAB versions (2012b and 2014b) from crashing
    ft_plot_headshape(headshape_downsampled); view([0,0]);
    ft_plot_mesh(fiducials,'vertexcolor','y','vertexsize',20);
    
    % Save data for later
    disp('Saving mri_realigned and fiducial information');
    mri_realigned = mri_orig;
    save mri_realigned mri_realigned
    
    % Save the fiducials for later
    save fiducials fiducials
    
    %% Create scalp mesh
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
    fprintf('Saving to the mesh libary\n')
    save mesh mesh
    
    %% Load the segmented brain
    fprintf('Loading the segmented MRI\n');
    
    mri_file_without_head = ['/Users/44737483/Documents/scripts_mcq/'...
        'MRIDataBase_JohnRichards_USC/Children/Head/ANTS' list_of_ages{age}...
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
    
    % Transform based on the fiducials you marked
    headmodel = ft_transform_vol(inv(mri_realigned.transformorig),headmodel);
    headmodel = ft_transform_vol(mri_realigned.transform,headmodel);
    
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
        cfg                 = [];
        cfg.grid.warpmni    = 'yes';
        cfg.grid.resolution = sourcemodel_mm(size);
        cfg.grid.nonlinear  = 'yes'; % use non-linear normalization
        cfg.mri             = mri_realigned;
        cfg.grid.unit       ='mm';
        cfg.inwardshift     = '-1.5';
        grid                = ft_prepare_sourcemodel(cfg);
        
        % Save
        sourcemodel3d = grid;
        save(sprintf('sourcemodel3d_%dmm',sourcemodel_mm(size)),'sourcemodel3d');
        
        % Create figure and save
        figure;ft_plot_mesh(grid.pos(grid.inside,:)); view([0,0]);
        print(sprintf('qc_sourcemodel3d_shape_%dmm',sourcemodel_mm(size)),'-dpng','-r100');
        
        figure;
        figure; ft_plot_vol(headmodel); alpha 0.3;
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

