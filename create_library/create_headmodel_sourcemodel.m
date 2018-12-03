
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create meshes, headmodels and sourcemodels for MEG sourcespace analysis 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path_to_MRI_library = '/Users/44737483/Documents/scripts_mcq/new_HCP_library_for_MEMES/';

% Get a list of all files and folders in this folder.
files = dir(path_to_MRI_library)
files(1:2) = [];
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);

% Now these names to a variable called subject
list_of_subs = [];

for k = 1 : length(subFolders)
    list_of_subs{k} = subFolders(k).name;
end

%list_of_subs = listFolders('D:\Judy\PACE\SLIM-database-marked');

for sub = 1:length(list_of_subs)
    try
        %% Go to subject directory
        dir_for_loop = [path_to_MRI_library...
            list_of_subs{sub}];
        
        cd(dir_for_loop)
        
        %% Load the brain
        load('mri_realigned.mat');
        load('mesh.mat');
        
        %% Segment
        cfg           = [];
        cfg.output    = 'brain';
        mri_segmented  = ft_volumesegment(cfg, mri_realigned);
        
        %% Create singleshell headmodel
        cfg = [];
        cfg.tissue = 'brain';
        cfg.method='singleshell';
        headmodel_singleshell = ft_prepare_headmodel(cfg, mri_segmented);
        
        % Create Figure
        figure;ft_plot_vol(headmodel_singleshell);
        ft_plot_mesh(mesh); alpha 0.3; view([0,0]);
        print('qc_headmodel','-dpng','-r100');
        
        % Save headmodel
        headmodel = headmodel_singleshell;
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
            
            % Save
            sourcemodel3d = grid;
            save(sprintf('sourcemodel3d_%dmm',sourcemodel_mm(size)),'sourcemodel3d');
            
            % Create figure and save
            figure;ft_plot_mesh(grid.pos(grid.inside,:)); view([0,0]);
            print(sprintf('qc_sourcemodel3d_shape_%dmm',sourcemodel_mm(size)),'-dpng','-r100');
            
            figure; ft_plot_vol(headmodel_singleshell); alpha 0.3;
            ft_plot_mesh(mesh); alpha 0.3;
            ft_plot_mesh(grid.pos(grid.inside,:)); view([0,0]);
            
            print(sprintf('qc_sourcemodel3d_%dmm',sourcemodel_mm(size)),'-dpng','-r100');
            
            % Clear for next loop
            clear grid sourcemodel3d
        end
        
        % Clear for next loop
        clear mesh headmodel mri_realigned headmodel_singleshell mri_segmented
        close all
        
        fprintf('Finished... CHECK for quality control\n');
        
    catch
        cd(path_to_MRI_library)
        txt_to_save = [list_of_subs{sub} ' could not be saved'];
        save(sprintf('%s',list_of_subs{sub}),'txt_to_save');
        clear txt_to_save
    end
end


