function MEMES3(dir_name,grad_trans,headshape_downsampled,...
    path_to_MRI_library,method,scaling,varargin)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MRI Estimation for MEG Sourcespace (MEMES)
%
%%%%%%%%%%%
% Inputs:
%%%%%%%%%%%
%
% - dir_name              = directory for saving
% - grad_trans            = MEG sensors realigned based on elp/mrk files
% - headshape_downsampled = headshape downsampled to 100-200 scalp points
% - path_to_MRI_library   = path to HCP MRI library
% - method                = method for creating pseudo head- and
%                         source-model: 'best' or 'average'
% - scaling               = scaling factor range applied to MRIs
%
%%%%%%%%%%%%%%%%%%
% Variable Inputs:
%%%%%%%%%%%%%%%%%%
%
% sourcemodel_size      = size of sourcemodel grid (5,8 or 10mm)
%
%%%%%%%%%%%
% Outputs:
%%%%%%%%%%%
%
% - grad_trans              = sensors transformed to correct
% - shape                   = headshape and fiducial information
% - headshape_downsampled   = headshape downsampled to 100 points
% - trans_matrix            = transformation matrix applied to headmodel
%                           and sourcemodel
% - sourcemodel3d           = sourcemodel warped to MNI space
% - headmodel               = singleshell headmodel (10000 vertices)
%
%%%%%%%%%%%%%%%%%%%%%
% Other Information:
%%%%%%%%%%%%%%%%%%%%%
%
% Example function call:
% MEMES3(dir_name,grad_trans,headshape_downsampled,...
%    path_to_MRI_library,method,[0.98:1.02],8)

% I have introduced a variable scaling parameter for the MRIs to
% help with coregistration. For example to apply -2% to +2% scaling to
% every MRI specify: scaling = [0.98:0.01:1.2].
%
% However NOTE: the more scaling factors you apply the longer it will take
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(['\nThis is MEMES v0.3\n\nMake sure you have asked Robert'...
    'for an MRI library\n\n']);
warning('on')

%% Check inputs
disp('Performing input check');
% If Path to MRI library doesn't end with / or \ throw up and error
if ismember(path_to_MRI_library(end),['/','\']) == 0
    error('!!! Path to MRI library must end with / or \ !!!');
end

%assert(method == 'average','method = average is not yet supported. Use best\n');

if length(scaling) == 1
    scaling = 1;
end

% If variable inputs are empty use defaults
if isempty(varargin)
    sourcemodel_size    = 8;
else
    sourcemodel_size    = varargin{1};
end

%% Extract subject names from your MRI library

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

% Now try to load relevent information from the first subject
fprintf('Now checking the MRI library is organised correctly...\n');
try
    load([path_to_MRI_library subject{1} '/mesh.mat']);
    load([path_to_MRI_library subject{1} '/headmodel.mat']);
    load([path_to_MRI_library subject{1} '/mri_realigned.mat']);
    load([path_to_MRI_library subject{1} '/sourcemodel3d_8mm.mat']);
    clear mesh headmodel mri_realigned sourcemodel3d
    fprintf('...Subject %s is organised correctly!\n',subject{1});
    
catch
    warning('Your MRI library is not organised correctly');
    disp('Each folder should contain: mesh.mat, headmodel.mat, mri_realigned.mat, sourcemodel3d_8mm.mat');
end

%% CD to the right place

% CD to right place
cd(dir_name); fprintf('\nCDd to the right place\n');

%% Perform ICP

% Error term variable - MEMES will crash here if your MRI library path is
% wrong..
error_term = zeros(1,length(subject));
% Variable to hold the transformation matrices
trans_matrix_library = [];
scaling_factor_all = zeros(1,length(subject));
count = 1;

% For each subject...
for m = 1:length(subject)
    
    % Load the mesh
    load([path_to_MRI_library subject{m} '/mesh.mat'])
    
    numiter = 30; count2 = 1;
    
    trans_matrix_temp = []; error_2 = [];
    
    % Perform ICP fit with different scaling factors
    for scale = scaling
        fprintf('Completed iteration %d of %d ; %d of %d MRIs\n',count2,length(scaling),m,length(subject));
        mesh_coord_scaled = ft_warp_apply([scale 0 0 0;0 scale 0 0; 0 0 scale 0; 0 0 0 1],mesh.pos);
        % Perform ICP
        [R, t, err, dummy, info] = icp(mesh_coord_scaled', ...
            headshape_downsampled.pos', numiter, 'Minimize', 'plane',...
            'Extrapolation', true,'WorstRejection', 0.05);
        error_2(count2) = err(end);
        trans_matrix_temp{count2} = inv([real(R) real(t);0 0 0 1]);
        count2 = count2+1;
    end
    
    % Find scaling factor with smallest error
    min_error = min(error_2);
    % Add error to error_term
    error_term(m) = min_error;
    
    % Add transformation matrix to trans_matrix_library
    trans_matrix_library{m} = trans_matrix_temp{find(error_2==min_error)};
    % Add scaling factor
    scaling_factor_all(m) = scaling(find(error_2==min_error));
    
    fprintf('Best scaling factor is %.2f\n',scaling(find(error_2==min_error)));
    
    % Clear mesh for next loop
    clear mesh
end

fprintf(' Finished the iterations\n');

%% Make pretty figure
fprintf('\n Finding good, OK and bad examples\n');

error_term_sorted = sort(error_term, 'ascend');
middle_num = length(error_term_sorted)/2;
winners = find(ismember(error_term,error_term_sorted(1:3)));
middles = find(ismember(error_term,error_term_sorted(middle_num-1:middle_num+1)));
losers = find(ismember(error_term,error_term_sorted(end-2:end)));

concat = [winners middles losers];

% Create figure to summarise the losers,middles and winners
figure;
for i = 1:9
    load([path_to_MRI_library subject{(concat(i))} '/mesh.mat'])
    mesh_spare = mesh;
    mesh_spare.pos = ft_warp_apply([scaling_factor_all(concat(i)) 0 0 0;...
        0 scaling_factor_all(concat(i)) 0 0; ...
        0 0 scaling_factor_all(concat(i)) 0; 0 0 0 1],mesh_spare.pos);
    mesh_spare.pos = ft_warp_apply(trans_matrix_library{(concat(i))}, mesh_spare.pos);
    
    subplot(3,3,i)
    ft_plot_mesh(mesh_spare,'facecolor',[238,206,179]./255,'EdgeColor','none','facealpha',0.8); hold on;
    camlight; hold on; view([-270,-10]);
    if ismember(i,1:3)
        title(sprintf('BEST: %d', error_term((concat(i)))));
    elseif ismember(i,4:6)
        title(sprintf('MIDDLE: %d', error_term((concat(i)))));
    elseif ismember(i,7:9)
        title(sprintf('WORST: %d', error_term((concat(i)))));
    end
    
    ft_plot_headshape(headshape_downsampled);
    
    clear mesh mesh_spare
    
    if i == 9
        print('best_middle_worst_examples','-dpng','-r100');
    end
end

%% Create figure to show different scaling factors

if length(scaling) > 1
    try
        figure;hist(scaling_factor_all,length(scaling));
        ylabel('Count');
        xlabel('Scaling Parameter');
        
        % Get information about the same
        % histogram by returning arguments
        [n,x] = hist(scaling_factor_all,5);
        % Create strings for each bar count
        barstrings = num2str(n');
        
        barstrings2 = num2str(scaling');
        
        % Create text objects at each location
        ylim([0 max(n)+5]);
        text(x,n,barstrings,'horizontalalignment','center','verticalalignment','bottom');
        
        xticks(scaling);
        xTick = get(gca,'xtick');
        
        h = findobj(gca,'Type','patch');
        h.FaceColor = [0 0.5 0.5];
        h.EdgeColor = 'w';
        set(gca,'FontSize',15);
        print('scaling_factor_distribution','-dpng','-r100');
    catch
        disp('Cannot Display scaling factors (?)');
    end
end


%% Use the best for to create a source model for MEG source analysis

% winner = find(error_term == min(min(error_term)));
% fprintf('\nThe winning MRI is number %d of %d\n',winner,length(mesh_library));
% trans_matrix = trans_matrix_library{winner};
%
% % Create figure to show ICP fit
% mesh_spare = mesh_library{winner};
% mesh_spare.pos = ft_warp_apply(trans_matrix, mesh_spare.pos);
%
% figure;ft_plot_mesh(mesh_spare,'facecolor',[238,206,179]./255,'EdgeColor','none','facealpha',0.8); hold on;
% camlight; hold on; view([-180,-10]);
% title(error_term(winner));
% ft_plot_headshape(headshape_downsampled);
%
% % print('winning_sourcemodel','-dpng','-r100');
%
% try
%     % % Make fancy video
%     c = datestr(clock); %time and date
%
%     figure;
%     ft_plot_mesh(mesh_spare,'facecolor',[238,206,179]./255,'EdgeColor','none','facealpha',0.8); hold on;
%     camlight; hold on;
%     ft_plot_headshape(headshape_downsampled); title(sprintf('%s.   Error of ICP fit = %d' , c, error_term(winner)));
%     OptionZ.FrameRate=15;OptionZ.Duration=5.5;OptionZ.Periodic=true;
%     CaptureFigVid([0,0; 360,0], 'ICP_quality',OptionZ)
%
% catch
%     fprintf('You need CaptureFigVid in your path for fancy videos\n');
% end

fprintf('\n Constructing the headmodel and sourcemodel \n');

switch method
    case 'average'
        fprintf('USE WITH CAUTION - Still testing \n');
        
        % Average over how many? N=20 the best?
        average_over_n = 20;
        % Variable to hold average sourcemodel .pos
        average_sourcemodel_all = [];
        % Variable to hold average sourcemodel .pos
        average_headmodel_all = [];
        
        for rep = 1:average_over_n
            % Find the number of the nth MRI
            winner_rep = find(ismember(error_term,error_term_sorted(rep)));
            
            % Update the user
            fprintf('Loaded MRI %d of %d : %s ... Scaling factor: %.2f\n',...
                rep,average_over_n,subject{winner_rep},...
                scaling_factor_all(winner_rep));
            
            % Get the transformation matrix of the winner
            trans_matrix = trans_matrix_library{winner_rep};
            
            %% Get mesh
            % Get facial mesh of 1st winner
            if rep == 1
                load([path_to_MRI_library subject{winner_rep} '/mesh.mat'])
                mesh.pos = ft_warp_apply([scaling_factor_all(winner_rep) 0 0 0;0 ...
                    scaling_factor_all(winner_rep) 0 0; 0 0 scaling_factor_all(winner_rep) 0;...
                    0 0 0 1],mesh.pos);
                mesh.pos = ft_warp_apply(trans_matrix, mesh.pos);
                mesh_spare = mesh;
            end
            
            clear mesh
            
            %% Create Headmodel (in mm)
            load([path_to_MRI_library subject{winner_rep} '/headmodel.mat']);
            
            % Scale
            headmodel.bnd.pos = ft_warp_apply([scaling_factor_all(winner_rep) 0 0 0;0 ...
                scaling_factor_all(winner_rep) 0 0; 0 0 scaling_factor_all(winner_rep) 0; 0 0 0 1],...
                headmodel.bnd.pos);
            
            % Transform (MESH --> coreg via ICP adjustment)
            headmodel.bnd.pos = ft_warp_apply(trans_matrix,headmodel.bnd.pos);
            
            % Add the pos field to the array outside the loop
            average_headmodel_all(rep,:,:) = headmodel.bnd.pos(:,:);
            
            % Reserve the first headmodel for later
            if rep == 1
                headmodel_for_outside_loop = headmodel;
            end
            
            clear headmodel
            
            %% Create Sourcemodel (in mm)
            
            % Load specified sized sourcemodel
            load([path_to_MRI_library ...
                subject{winner_rep} '/sourcemodel3d_' num2str(sourcemodel_size) 'mm.mat']);
            
            % Scale
            sourcemodel3d.pos = ft_warp_apply([scaling_factor_all(winner_rep)...
                0 0 0;0 scaling_factor_all(winner_rep) 0 0; 0 0 ...
                scaling_factor_all(winner_rep) 0; 0 0 0 1],sourcemodel3d.pos);
            
            % Transform (MESH --> coreg via ICP adjustment)
            sourcemodel3d.pos = ft_warp_apply(trans_matrix,sourcemodel3d.pos);
            
            average_sourcemodel_all(rep,:,:) = sourcemodel3d.pos;
            
            % Reserve the first headmodel for later
            if rep == 1
                sourcemodel_for_outside_loop = sourcemodel3d;
            end
            
            clear trans_matrix sourcemodel3d winner_rep
            
        end
        
        % Average Headmodel
        fprintf('Averaging Headmodel\n');
        headmodel = headmodel_for_outside_loop;
        headmodel.bnd.pos = squeeze(mean(average_headmodel_all,1));
        
        % Average Sourcemodel
        fprintf('Averaging Sourcemodel\n');
        sourcemodel3d = sourcemodel_for_outside_loop;
        sourcemodel3d.pos = squeeze(mean(average_sourcemodel_all,1));
        
        % Create figure to check headodel and sourcemodel match
        figure;
        ft_plot_vol(headmodel,  'facecolor', 'cortex', 'edgecolor', 'none');
        alpha 0.4; camlight;
        ft_plot_mesh(sourcemodel3d.pos(sourcemodel3d.inside,:),'vertexsize',5);
        view([0 0]);
        
        view_angle = [0 90 180 270];
        
        % Create figure to show final coregiration (with mesh of 1st place
        % MRI)
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
        
        %% SAVE
        fprintf('\nSaving the necessary data\n');
        
        save headmodel headmodel
        %save trans_matrix trans_matrix
        save grad_trans grad_trans
        save sourcemodel3d sourcemodel3d
        %save mri_realigned_MEMES mri_realigned_MEMES
        
        fprintf('\nCOMPLETED - check the output for quality control\n');
        
        
    case 'best'
        
        % Find the MRI with the lowest ICP error between Polhemus points
        % and 3D scalp mesh
        winner = find(error_term == min(min(error_term)));
        fprintf('\nThe winning MRI is number %d of %d : %s\n',winner,length(subject),subject{winner});
        
        % Get the transformation matrix of the winner
        trans_matrix = trans_matrix_library{winner};
        
        % Get facial mesh of winner
        load([path_to_MRI_library subject{winner} '/mesh.mat'])
        mesh_spare = mesh;
        mesh_spare.pos = ft_warp_apply([scaling_factor_all(winner) 0 0 0;0 ...
            scaling_factor_all(winner) 0 0; 0 0 scaling_factor_all(winner) 0;...
            0 0 0 1],mesh_spare.pos);
        mesh_spare.pos = ft_warp_apply(trans_matrix, mesh_spare.pos);
        
        % Get MRI of winning subject
        fprintf('Transforming the MRI\n');
        load([path_to_MRI_library subject{winner} '/mri_realigned.mat'],'mri_realigned');
        disp('done loading');
        mri_realigned_MEMES = ft_transform_geometry(trans_matrix,...
            mri_realigned);
        
        %% Create Headmodel (in mm)
        fprintf(' Creating Headmodel in mm\n');
        
        load([path_to_MRI_library subject{winner} '/headmodel.mat']);
        
        % Scale
        headmodel.bnd.pos = ft_warp_apply([scaling_factor_all(winner) 0 0 0;0 ...
            scaling_factor_all(winner) 0 0; 0 0 scaling_factor_all(winner) 0; 0 0 0 1],...
            headmodel.bnd.pos);
        
        % Transform (MESH --> coreg via ICP adjustment)
        headmodel.bnd.pos = ft_warp_apply(trans_matrix,headmodel.bnd.pos);
        
        figure;
        ft_plot_vol(headmodel);
        ft_plot_headshape(headshape_downsampled);
        
        %% Create Sourcemodel (in mm)
        fprintf('Creating an %dmm Sourcemodel in mm\n',sourcemodel_size);
        
        % Load specified sized sourcemodel
        load([path_to_MRI_library ...
            subject{winner} '/sourcemodel3d_' num2str(sourcemodel_size) 'mm.mat']);
        
        % Scale
        sourcemodel3d.pos = ft_warp_apply([scaling_factor_all(winner) 0 0 0;0 scaling_factor_all(winner) 0 0; 0 0 scaling_factor_all(winner) 0; 0 0 0 1],sourcemodel3d.pos);
        
        % Transform (MESH --> coreg via ICP adjustment)
        sourcemodel3d.pos = ft_warp_apply(trans_matrix,sourcemodel3d.pos);
        
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
        
        
        %         %% Create coregistered 3D cortical mesh
        %         mesh = ft_read_headshape({[path_to_MRI_library ...
        %             subject{winner} '/MEG/anatomy/' subject{winner} '.L.midthickness.4k_fs_LR.surf.gii'],...
        %             [path_to_MRI_library subject{winner} '/MEG/anatomy/' subject{winner} ...
        %             '.R.midthickness.4k_fs_LR.surf.gii']});
        %
        %         mesh = ft_convert_units(mesh,'mm');
        %
        %         % Transform 1 (MESH --> coreg via manual marking of fiducial points)
        %         mesh.pos = ft_warp_apply(rmatx,mesh.pos);
        %         mesh.pos = ft_warp_apply(initial_mri_realign{winner},mesh.pos);
        %
        %         % Scale
        %         mesh.pos = ft_warp_apply([scaling_factor_all(winner) 0 0 0;0 scaling_factor_all(winner) 0 0; 0 0 scaling_factor_all(winner) 0; 0 0 0 1],mesh.pos);
        %
        %         %transform 2 (MESH --> coreg via ICP adjustment)
        %         mesh.pos = ft_warp_apply(trans_matrix,mesh.pos);
        %
        %         %ft_determine_coordsys(mri_realigned2,'interactive','no'); hold on;
        %         figure;ft_plot_sens(grad_trans);
        %         ft_plot_headshape(headshape_downsampled) %plot headshape
        %         ft_plot_mesh(mesh,'facealpha',0.8); camlight; hold on; view([100 4]);
        %         print('headmodel_3D_cortical_mesh_quality','-dpng');
        
        %% SAVE
        fprintf('\nSaving the necessary data\n');
        
        save headmodel headmodel
        save trans_matrix trans_matrix
        save grad_trans grad_trans
        save sourcemodel3d sourcemodel3d
        save mri_realigned_MEMES mri_realigned_MEMES
        %        save mesh mesh
        
        
        fprintf('\nCOMPLETED - check the output for quality control\n');
        
    otherwise
        fprintf('Something went wrong - did you specify *average* or *best*')
end
end