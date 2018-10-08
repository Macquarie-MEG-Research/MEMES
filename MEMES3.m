
function MEMES3(dir_name,elpfile,hspfile,confile,mrkfile,path_to_MRI_library,bad_coil,method,scaling,varargin)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MRI Estimation for MEG Sourcespace (MEMES)
%
%%%%%%%%%%%
% Inputs:
%%%%%%%%%%%
%
% - dir_name            = directory for saving
% - elpfile             = path to elp file
% - hspfile             = path to hsp file
% - confile             = path to con file
% - mrkfile             = path to mrk file
% - path_to_MRI_library = path to HCP MRI library
% - mesh_library        = mesh library (in mm) created from HCP MRIs
% - initial_mri_realign = transform for initial realigning estimate
% - bad_coil            = list of bad coils (up to length of 2). Enter as:
%                         {LPAred','RPAyel','PFblue','LPFwh','RPFblack'}
% - method              = method for creating pseudo head- and
%                       source-model: 'best' or 'average'
% - scaling             = scaling factor range applied to MRIs
%
%%%%%%%%%%%%%%%%%%
% Variable Inputs:
%%%%%%%%%%%%%%%%%%
%
% sourcemodel_size      = size of sourcemodel grid (5,8 or 10mm)
% include_face          = inclue facial points acquired during head
%                       digitisation ('yes' = default)
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

%%%%%%%%%%%%%%%%%%%%%
% Other Information:
%%%%%%%%%%%%%%%%%%%%%

% Example function call:
% MEMES3(dir_name,elpfile,hspfile,confile,mrkfile,path_to_MRI_library,...
% '','best',[0.98:0.01:1.02],8,'no')

% I have introduced a variable scaling parameter for the MRIs to
% help with coregistration. For example to apply -2% to +2% scaling to
% every MRI specify: scaling = [0.98:0.01:1.2].
%
% However NOTE: the more scaling factors you apply the longer it will take
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nThis is MEMES v0.3\n\nMake sure you have asked Robert for an MRI library\n\n');
warning('on')

%% Check inputs
disp('Performing input check');
assert(length(bad_coil)<3,'You need at least 3 good coils for accurate alignment. Also make sure you enter bad_coil strings in curly brackets {}');
% If Path to MRI library doesn't end with / or \ throw up and error
if ismember(path_to_MRI_library(end),['/','\']) == 0
    error('!!! Path to MRI library must end with / or \ !!!');
end

% Check if bad_coils are entered correctly
if strcmp(bad_coil,'')
    disp('No bad coils marked');
else
    for check1 = 1:length(bad_coil)
        if ismember(bad_coil{check1},{'','LPAred','RPAyel','PFblue','LPFwh','RPFblack'}) == 0
            error('!!! Please enter bad_coils correctly in the form {LPAred,RPAyel,PFblue,LPFwh,RPFblack} !!!');
        end
    end
end


%assert(method == 'average','method = average is not yet supported. Use best\n');

if length(scaling) == 1
    scaling = 1;
end

% If variable inputs are empty use defaults
if isempty(varargin)
    sourcemodel_size    = 8;
    include_face        = 'yes';
else
    sourcemodel_size    = varargin{1};
    include_face        = varargin{2};
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

%% Read the relevent subject data from file

% CD to right place
cd(dir_name); fprintf('\nCDd to the right place\n');

% Get Polhemus Points
disp('Reading elp and hspfile');
[shape]  = parsePolhemus(elpfile,hspfile);
shape   = ft_convert_units(shape,'mm');

% Read the grads from the con file
disp('Reading Sensors');
grad_con = ft_read_sens(confile); %in cm, load grads
grad_con = ft_convert_units(grad_con,'mm'); %in mm

% Read mrk_file
disp('Reading the mrk file');
mrk      = ft_read_headshape(mrkfile,'format','yokogawa_mrk');
mrk      = ft_convert_units(mrk,'mm'); %in mm

%% Perform Realighment Using Paul's Fancy Functions
if strcmp(bad_coil,'')
    disp('NO BAD MARKERS');
    markers                     = mrk.fid.pos([2 3 1 4 5],:);%reorder mrk to match order in shape
    [R,T,Yf,Err]                = rot3dfit(markers,shape.fid.pnt(4:end,:));%calc rotation transform
    meg2head_transm             = [[R;T]'; 0 0 0 1];%reorganise and make 4*4 transformation matrix

    disp('Performing re-alignment');
    grad_trans                  = ft_transform_geometry_PFS_hacked(meg2head_transm,grad_con); %Use my hacked version of the ft function - accuracy checking removed not sure if this is good or not
    grad_trans.fid              = shape; %add in the head information
    save grad_trans grad_trans

    % Else if there is a bad marker
else
    fprintf(''); disp('TAKING OUT BAD MARKER(S)');

    badcoilpos = [];

    % Identify the bad coil
    for num_bad_coil = 1:length(bad_coil)
        pos_of_bad_coil = find(ismember(shape.fid.label,bad_coil{num_bad_coil}))-3;
        badcoilpos(num_bad_coil) = pos_of_bad_coil;
    end

    % Re-order mrk file to match elp file
    markers               = mrk.fid.pos([2 3 1 4 5],:);%reorder mrk to match order in shape
    % Now take out the bad marker(s) when you realign
    markers(badcoilpos,:) = [];

    % Get marker positions from elp file
    fids_2_use = shape.fid.pnt(4:end,:);
    % Now take out the bad marker(s) when you realign
    fids_2_use(badcoilpos,:) = [];

    % If there are two bad coils use the ICP method, if only one use
    % rot3dfit as usual
    disp('Performing re-alignment');

    if length(bad_coil) == 2
        [R, T, err, dummy, info]    = icp(fids_2_use', markers','Minimize', 'point');
        meg2head_transm             = [[R T]; 0 0 0 1];%reorganise and make 4*4 transformation matrix
        grad_trans                  = ft_transform_geometry_PFS_hacked(meg2head_transm,grad_con); %Use my hacked version of the ft function - accuracy checking removed not sure if this is good or not
        grad_trans.fid              = shape; %add in the head information
    else
        [R,T,Yf,Err]                = rot3dfit(markers,fids_2_use);%calc rotation transform
        meg2head_transm             = [[R;T]'; 0 0 0 1];%reorganise and make 4*4 transformation matrix
        grad_trans                  = ft_transform_geometry_PFS_hacked(meg2head_transm,grad_con); %Use my hacked version of the ft function - accuracy checking removed not sure if this is good or not
        grad_trans.fid              = shape; %add in the head information
    end
end

% Create figure to view relignment
hfig = figure;
subplot(2,2,1);ft_plot_headshape(shape);
hold on; ft_plot_sens(grad_trans); view([180, 0]);
subplot(2,2,2);ft_plot_headshape(shape);
hold on; ft_plot_sens(grad_trans); view([-90, 0]);
subplot(2,2,3);ft_plot_headshape(shape);
hold on; ft_plot_sens(grad_trans); view([0, 0]);
hax = subplot(2,2,4);ft_plot_headshape(shape);
hold on; ft_plot_sens(grad_trans); view([90, 0]);

% Downsample the headshape information with or without facial points
fprintf('Downsampling headshape information to %d points\n'...
    ,100);
headshape_downsampled = downsample_headshape(hspfile,100,include_face);
headshape_downsampled = ft_convert_units(headshape_downsampled,'mm'); %in mm
disp('Saving headshape downsampled');
save headshape_downsampled headshape_downsampled

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
        [R, t, err, dummy, info] = icp(mesh_coord_scaled', headshape_downsampled.pos', numiter, 'Minimize', 'plane', 'Extrapolation', true,'WorstRejection', 0.05);
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

    clear mesh

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
        fprintf('NOT SUPPORTED YET \n');
        %         %% Calculate headmodel and sourcemodel over first n MRIs
        %
        %         average_over_n = 10; average_sourcemodel_all = []; average_mri_all = [];
        %
        %         average_n = find(ismember(error_term,error_term_sorted(1:average_over_n)));
        %
        %
        %         for rep = 1:average_over_n
        %
        %             fprintf(' Getting sourcemodel %d of %d \n',rep,average_over_n);
        %
        %             % Get MRI of winning subject
        %             mri_file = [path_to_MRI_library subject{average_n(rep)} '/MEG/anatomy/T1w_acpc_dc_restore.nii.gz'];
        %
        %             mri_orig                    = ft_read_mri(mri_file); % in mm, read in mri from DICOM
        %             mri_orig = ft_convert_units(mri_orig,'cm'); mri_orig.coordsys = 'neuromag';
        %
        %             mri_orig.transform = initial_mri_realign{average_n(rep)};
        %             mri_realigned = mri_orig;
        %
        %             mri_realigned = ft_transform_geometry(trans_matrix_library{average_n(rep)},mri_realigned);
        %
        %             average_mri_all{rep} = mri_realigned;
        %
        %             %% Create Sourcemodel (in cm)
        %             fprintf(' Creating Sourcemodel in cm\n');
        %
        %             % Get transformation matrix to convert BTI to SPM
        %             % (N.B. subject 1 not working?)
        %             path_to_transform   = [path_to_MRI_library subject{average_n(rep)} '/MEG/anatomy/' subject{average_n(rep)} '_MEG_anatomy_transform.txt'];
        %
        %             txt = textscan(fopen(path_to_transform),'%s%s%s%s\n'); %open text file
        %             transform.bti2spm = zeros(4);
        %             fclose('all');
        %
        %             for reps = 1:4 %some jiggery-pokery to extract transform matrix from text file
        %                 transform.bti2spm(:,reps) = str2num(strjoin(txt{1,reps}(32:35)))';
        %             end
        %
        %             clear reps
        %
        %             path_to_sourcemodel = [path_to_MRI_library subject{average_n(rep)} '/MEG/anatomy/' subject{average_n(rep)} '_MEG_anatomy_sourcemodel_3d8mm.mat'];
        %
        %             % Load sourcemodel and convert to 'mm'
        %             load(path_to_sourcemodel); sourcemodel3d = ft_convert_units(sourcemodel3d,'mm');
        %
        %             % Transform BTI --> MRI-space
        %             sourcemodel3d.pos = ft_warp_apply(transform.bti2spm,sourcemodel3d.pos);
        %
        %             % Convert to cm
        %             sourcemodel3d = ft_convert_units(sourcemodel3d,'cm');
        %
        %             % Transform 1 (MESH --> coreg via manual marking of fiducial points)
        %             sourcemodel3d.pos = ft_warp_apply(inv(mri_orig.transform),sourcemodel3d.pos);
        %             sourcemodel3d.pos = ft_warp_apply(initial_mri_realign{average_n(rep)},sourcemodel3d.pos);
        %
        %             %transform 2 (MESH --> coreg via ICP adjustment)
        %             sourcemodel3d.pos = ft_warp_apply(trans_matrix_library{average_n(rep)},sourcemodel3d.pos);
        %
        %             average_sourcemodel_all{rep} = sourcemodel3d;
        %
        %         end
        %
        %         pos_headmodel = []; pos_sourcemodel = []; tri_headmodel = [];
        %
        %         for rep = 1:average_over_n
        %             pos_headmodel(rep,:,:,:) = average_headmodel_all{1,rep}.bnd.pnt;
        %             tri_headmodel(rep,:,:,:) = average_headmodel_all{1,rep}.bnd.tri;
        %             pos_sourcemodel(rep,:,:,:) = average_sourcemodel_all{1,rep}.pos;
        %         end
        %
        %         average_headmodel = average_headmodel_all{1};
        %         average_headmodel.bnd.pnt = squeeze(mean(pos_headmodel));
        %         average_headmodel.bnd.tri = squeeze(mean(tri_headmodel));
        %
        %         average_sourcemodel = average_sourcemodel_all{1};
        %         average_headmodel.pos = squeeze(mean(pos_sourcemodel));
        %
        %         figure;ft_plot_vol(average_headmodel,  ...
        %                 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.6; camlight;
        %             ft_plot_mesh(average_sourcemodel.pos(average_sourcemodel.inside,:),'vertexsize',2);
        %             %ft_plot_sens(grad_trans, 'style', 'r*')
        %             %ft_plot_headshape(headshape_downsampled) %plot headshape
        %
        %         volume(average_headmodel_all{1,1}.bnd.pnt)
        %
        %
        %         view_angle = [0 90 180 270];
        %         figure; hold on;
        %         for rep = 1:4
        %             subplot(2,2,rep);%ft_plot_vol(average_headmodel,  ...
        %                 %'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.6; camlight;
        %             ft_plot_mesh(average_sourcemodel.pos(average_sourcemodel.inside,:),'vertexsize',2);
        %             ft_plot_sens(grad_trans, 'style', 'r*')
        %             ft_plot_headshape(headshape_downsampled) %plot headshape
        %             view([view_angle(rep),0]);
        %         end
        %
        %         %% SAVE
        %         fprintf('\nSaving the necessary data\n');
        %
        %         save average_headmodel average_headmodel
        %         save grad_trans grad_trans
        %         save average_sourcemodel average_sourcemodel
        %
        %
        %         fprintf('\nCOMPLETED - check the output for quality control\n');


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
        mri_realigned = load([path_to_MRI_library...
            subject{winner} '/mri_realigned.mat']);

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
        %        save mesh mesh


        fprintf('\nCOMPLETED - check the output for quality control\n');

    otherwise
        fprintf('Something went wrong - did you specify *average* or *best*')
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Subfunctions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [shape] = parsePolhemus(elpfile,hspfile)

        fid1 = fopen(elpfile);
        C = fscanf(fid1,'%c');
        fclose(fid1);

        E = regexprep(C,'\r','xx');
        E = regexprep(E,'\t','yy');

        returnsi = strfind(E,'xx');
        tabsi = strfind(E,'yy');
        sensornamesi = strfind(E,'%N');
        fiducialsstarti = strfind(E,'%F');
        lastfidendi = strfind(E(fiducialsstarti(3):fiducialsstarti(length(fiducialsstarti))+100),'xx');
        fiducialsendi = fiducialsstarti(1)+strfind(E(fiducialsstarti(1):fiducialsstarti(length(fiducialsstarti))+lastfidendi(1)),'xx');

        NASION = E(fiducialsstarti(1)+4:fiducialsendi(1)-2);
        NASION = regexprep(NASION,'yy','\t');
        NASION = str2num(NASION);

        LPA = E(fiducialsstarti(2)+4:fiducialsendi(2)-2);
        LPA = regexprep(LPA,'yy','\t');
        LPA = str2num(LPA);

        RPA = E(fiducialsstarti(3)+4:fiducialsendi(3)-2);
        RPA = regexprep(RPA,'yy','\t');
        RPA = str2num(RPA);

        LPAredstarti = strfind(E,'LPAred');
        LPAredendi = strfind(E(LPAredstarti(1):LPAredstarti(length(LPAredstarti))+45),'xx');
        LPAred = E(LPAredstarti(1)+11:LPAredstarti(1)+LPAredendi(2)-2);
        LPAred = regexprep(LPAred,'yy','\t');
        LPAred = str2num(LPAred);

        RPAyelstarti = strfind(E,'RPAyel');
        RPAyelendi = strfind(E(RPAyelstarti(1):RPAyelstarti(length(RPAyelstarti))+45),'xx');
        RPAyel = E(RPAyelstarti(1)+11:RPAyelstarti(1)+RPAyelendi(2)-2);
        RPAyel = regexprep(RPAyel,'yy','\t');
        RPAyel = str2num(RPAyel);

        PFbluestarti = strfind(E,'PFblue');
        PFblueendi = strfind(E(PFbluestarti(1):PFbluestarti(length(PFbluestarti))+45),'xx');
        PFblue = E(PFbluestarti(1)+11:PFbluestarti(1)+PFblueendi(2)-2);
        PFblue = regexprep(PFblue,'yy','\t');
        PFblue = str2num(PFblue);

        LPFwhstarti = strfind(E,'LPFwh');
        LPFwhendi = strfind(E(LPFwhstarti(1):LPFwhstarti(length(LPFwhstarti))+45),'xx');
        LPFwh = E(LPFwhstarti(1)+11:LPFwhstarti(1)+LPFwhendi(2)-2);
        LPFwh = regexprep(LPFwh,'yy','\t');
        LPFwh = str2num(LPFwh);

        RPFblackstarti = strfind(E,'RPFblack');
        RPFblackendi = strfind(E(RPFblackstarti(1):end),'xx');
        RPFblack = E(RPFblackstarti(1)+11:RPFblackstarti(1)+RPFblackendi(2)-2);
        RPFblack = regexprep(RPFblack,'yy','\t');
        RPFblack = str2num(RPFblack);

        allfids = [NASION;LPA;RPA;LPAred;RPAyel;PFblue;LPFwh;RPFblack];
        fidslabels = {'NASION';'LPA';'RPA';'LPAred';'RPAyel';'PFblue';'LPFwh';'RPFblack'};

        fid2 = fopen(hspfile);
        C = fscanf(fid2,'%c');
        fclose(fid2);
        E = regexprep(C,'\r','xx'); %replace returns with "xx"
        E = regexprep(E,'\t','yy'); %replace tabs with "yy"
        returnsi = strfind(E,'xx');
        tabsi = strfind(E,'yy');

        headshapestarti = strfind(E,'position of digitized points');
        headshapestartii = strfind(E(headshapestarti(1):end),'xx');
        headshape = E(headshapestarti(1)+headshapestartii(2)+2:end);
        headshape = regexprep(headshape,'yy','\t');
        headshape = regexprep(headshape,'xx','');
        headshape = str2num(headshape);

        shape.pnt = headshape;
        shape.fid.pnt = allfids;
        shape.fid.label = fidslabels;

        %convert to BESA style coordinates so can use the .pos file or sensor
        %config from .con
        %         shape.pnt = cat(2,fliplr(shape.pnt(:,1:2)),shape.pnt(:,3)).*1000;
        %         %shape.pnt = shape.pnt(1:length(shape.pnt)-15,:); % get rid of nose points may want to alter or comment this depending on your digitisation
        %         %shape.pnt = shape.pnt*1000;
        %         neg = shape.pnt(:,2)*-1;
        %         shape.pnt(:,2) = neg;
        %
        %         shape.fid.pnt = cat(2,fliplr(shape.fid.pnt(:,1:2)),shape.fid.pnt(:,3)).*1000;
        %         %shape.fid.pnt = shape.fid.pnt*1000;
        %         neg2 = shape.fid.pnt(:,2)*-1;
        %         shape.fid.pnt(:,2) = neg2;
        shape.unit='m';
        %        shape = ft_convert_units(shape,'cm');

        new_name2 = ['shape.mat'];
        save (new_name2,'shape');
    end

    function [R,T,Yf,Err] = rot3dfit(X,Y)
        %ROT3DFIT Determine least-square rigid rotation and translation.
        % [R,T,Yf] = ROT3DFIT(X,Y) permforms a least-square fit for the
        % linear form
        %
        % Y = X*R + T
        %
        % where R is a 3 x 3 orthogonal rotation matrix, T is a 1 x 3
        % translation vector, and X and Y are 3D points sets defined as
        % N x 3 matrices. Yf is the best-fit matrix.
        %
        % See also SVD, NORM.
        %
        % rot3dfit: Frank Evans, NHLBI/NIH, 30 November 2001
        %

        % ROT3DFIT uses the method described by K. S. Arun, T. S. Huang,and
        % S. D. Blostein, "Least-Squares Fitting of Two 3-D Point Sets",
        % IEEE Transactions on Pattern Analysis and Machine Intelligence,
        % PAMI-9(5): 698 - 700, 1987.
        %
        % A better theoretical development is found in B. K. P. Horn,
        % H. M. Hilden, and S. Negahdaripour, "Closed-form solution of
        % absolute orientation using orthonormal matrices", Journal of the
        % Optical Society of America A, 5(7): 1127 - 1135, 1988.
        %
        % Special cases, e.g. colinear and coplanar points, are not
        % implemented.

        %error(nargchk(2,2,nargin));
        narginchk(2,2); %PFS Change to update
        if size(X,2) ~= 3, error('X must be N x 3'); end;
        if size(Y,2) ~= 3, error('Y must be N x 3'); end;
        if size(X,1) ~= size(Y,1), error('X and Y must be the same size'); end;

        % mean correct

        Xm = mean(X,1); X1 = X - ones(size(X,1),1)*Xm;
        Ym = mean(Y,1); Y1 = Y - ones(size(Y,1),1)*Ym;

        % calculate best rotation using algorithm 12.4.1 from
        % G. H. Golub and C. F. van Loan, "Matrix Computations"
        % 2nd Edition, Baltimore: Johns Hopkins, 1989, p. 582.

        XtY = (X1')*Y1;
        [U,S,V] = svd(XtY);
        R = U*(V');

        % solve for the translation vector

        T = Ym - Xm*R;

        % calculate fit points

        Yf = X*R + ones(size(X,1),1)*T;

        % calculate the error

        dY = Y - Yf;
        Err = norm(dY,'fro'); % must use Frobenius norm
    end

    function [output] = ft_transform_geometry_PFS_hacked(transform, input)

        % FT_TRANSFORM_GEOMETRY applies a homogeneous coordinate transformation to
        % a structure with geometric information, for example a volume conduction model
        % for the head, gradiometer of electrode structure containing EEG or MEG
        % sensor positions and MEG coil orientations, a head shape or a source model.
        %
        % The units in which the transformation matrix is expressed are assumed to
        % be the same units as the units in which the geometric object is
        % expressed. Depending on the input object, the homogeneous transformation
        % matrix should be limited to a rigid-body translation plus rotation
        % (MEG-gradiometer array), or to a rigid-body translation plus rotation
        % plus a global rescaling (volume conductor geometry).
        %
        % Use as
        %   output = ft_transform_geometry(transform, input)
        %
        % See also FT_WARP_APPLY, FT_HEADCOORDINATES

        % Copyright (C) 2011, Jan-Mathijs Schoffelen
        %
        % This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
        % for the documentation and details.
        %
        %    FieldTrip is free software: you can redistribute it and/or modify
        %    it under the terms of the GNU General Public License as published by
        %    the Free Software Foundation, either version 3 of the License, or
        %    (at your option) any later version.
        %
        %    FieldTrip is distributed in the hope that it will be useful,
        %    but WITHOUT ANY WARRANTY; without even the implied warranty of
        %    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        %    GNU General Public License for more details.
        %
        %    You should have received a copy of the GNU General Public License
        %    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
        %
        % $Id: ft_transform_geometry.m$

        % flg rescaling check
        allowscaling = ~ft_senstype(input, 'meg');

        % determine the rotation matrix
        rotation = eye(4);
        rotation(1:3,1:3) = transform(1:3,1:3);

        if any(abs(transform(4,:)-[0 0 0 1])>100*eps)
            error('invalid transformation matrix');
        end

        %%### get rid of this accuracy checking below as some of the transformation
        %%matricies will be a bit hairy###
        if ~allowscaling
            % allow for some numerical imprecision
            %if abs(det(rotation)-1)>1e-6%100*eps
            %if abs(det(rotation)-1)>100*eps  % allow for some numerical imprecision
            %error('only a rigid body transformation without rescaling is allowed');
            %end
        end

        if allowscaling
            % FIXME build in a check for uniform rescaling probably do svd or so
            % FIXME insert check for nonuniform scaling, should give an error
        end

        tfields   = {'pos' 'pnt' 'o' 'coilpos' 'chanpos' 'chanposold' 'chanposorg' 'elecpos', 'nas', 'lpa', 'rpa', 'zpoint'}; % apply rotation plus translation
        rfields   = {'ori' 'nrm'     'coilori' 'chanori' 'chanoriold' 'chanoriorg'};                                          % only apply rotation
        mfields   = {'transform'};           % plain matrix multiplication
        recfields = {'fid' 'bnd' 'orig'};    % recurse into these fields
        % the field 'r' is not included here, because it applies to a volume
        % conductor model, and scaling is not allowed, so r will not change.

        fnames    = fieldnames(input);
        for k = 1:numel(fnames)
            if ~isempty(input.(fnames{k}))
                if any(strcmp(fnames{k}, tfields))
                    input.(fnames{k}) = apply(transform, input.(fnames{k}));
                elseif any(strcmp(fnames{k}, rfields))
                    input.(fnames{k}) = apply(rotation, input.(fnames{k}));
                elseif any(strcmp(fnames{k}, mfields))
                    input.(fnames{k}) = transform*input.(fnames{k});
                elseif any(strcmp(fnames{k}, recfields))
                    for j = 1:numel(input.(fnames{k}))
                        input.(fnames{k})(j) = ft_transform_geometry(transform, input.(fnames{k})(j));
                    end
                else
                    % do nothing
                end
            end
        end
        output = input;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that applies the homogeneous transformation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [new] = apply(transform, old)
        old(:,4) = 1;
        new = old * transform';
        new = new(:,1:3);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rotate_about_z - make a rotation matix for arbitrary rotation in degrees
% around z axis
%
% Written by Paul Sowman Oct 2017 (http://web.iitd.ac.in/~hegde/cad/lecture/L6_3dtrans.pdf - page 4)
%
% INPUTS:
% - deg        = degrees of rotation required
%
% OUTPUTS:
% - rmatx      = a 4*4 rotation matrix for deg degrees about z
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function rmatx=rotate_about_z(deg)

        deg   = deg2rad(deg);
        rmatx = [cos(deg) sin(deg) 0 0;-sin(deg) cos(deg) 0 0;0 0 1 0;0 0 0 1];
    end

%     function [headshape_downsampled] = downsample_headshape_noface(path_to_headshape,numvertices,sensors)
%         % Get headshape
%         headshape = ft_read_headshape(path_to_headshape);
%         % Convert to cm
%         headshape = ft_convert_units(headshape,'cm');
%         % Convert to BESA co-ordinates
%         headshape.pos = cat(2,fliplr(headshape.pos(:,1:2)),headshape.pos(:,3));
%         headshape.pos(:,2) = headshape.pos(:,2).*-1;
%
%         % Get indices of facial points (up to 4cm above nasion)
%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % Is 4cm the correct distance?
%         % Possibly different for child system?
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%         count_facialpoints = find(headshape.pos(:,3)<4);
%         if isempty(count_facialpoints)
%             disp('CANNOT FIND ANY FACIAL POINTS');
%         else
%             facialpoints = headshape.pos(count_facialpoints,:,:);
%             rrr = 1:4:length(facialpoints);
%             facialpoints = facialpoints(rrr,:); clear rrr;
%         end
%
%         % Remove facial points for now
%         headshape.pos(count_facialpoints,:) = [];
%
%         % Create mesh out of headshape downsampled to x points specified in the
%         % function call
%         cfg.numvertices = numvertices;
%         cfg.method = 'headshape';
%         cfg.headshape = headshape.pos;
%         mesh = ft_prepare_mesh(cfg, headshape);
%
%         % Replace the headshape info with the mesh points
%         headshape.pos = mesh.pos;
%
%         % Create figure for quality checking
%         figure; subplot(2,2,1);ft_plot_mesh(mesh); hold on;
%         title('Downsampled Mesh');
%         view(0,0);
%         subplot(2,2,2);ft_plot_mesh(headshape); hold on;
%         title('Downsampled Headshape View 1');
%         view(0,0);
%         subplot(2,2,3);ft_plot_mesh(headshape); hold on;
%         title('Downsampled Headshape View 2');
%         view(90,0);
%         subplot(2,2,4);ft_plot_mesh(headshape); hold on;
%         title('Downsampled Headshape View 3');
%         view(180,0);
%         print('headshape_quality','-dpdf');
%
%         % Add in names of the fiducials from the sensor
%         headshape.fid.label = {'NASION','LPA','RPA'};
%
%         % Convert fiducial points to BESA
%         headshape.fid.pos = cat(2,fliplr(headshape.fid.pos(:,1:2)),headshape.fid.pos(:,3));
%         headshape.fid.pos(:,2) = headshape.fid.pos(:,2).*-1;
%
%         % Plot for quality checking
%         figure;ft_plot_sens(sensors) %plot channel position : between the 1st and 2nd coils
%         ft_plot_headshape(headshape) %plot headshape
%         view(0,0);
%         print('headshape_quality2','-dpdf');
%
%         % Export filename
%         headshape_downsampled = headshape;
%
%     end

    function [headshape_downsampled] = downsample_headshape(path_to_headshape,...
            numvertices,varargin)

        % If not specified include the facial points
        if isempty(varargin)
            include_facial_points = 'yes';

        else
            include_facial_points = varargin{1};
        end


        % Get headshape
        headshape = ft_read_headshape(path_to_headshape);
        % Convert to cm
        headshape = ft_convert_units(headshape,'cm');
        % Convert to BESA co-ordinates
        %         headshape.pos = cat(2,fliplr(headshape.pos(:,1:2)),headshape.pos(:,3));
        %         headshape.pos(:,2) = headshape.pos(:,2).*-1;

        % Get indices of facial points (up to 3cm above nasion)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Is 3cm the correct distance?
        % Possibly different for child system?
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        count_facialpoints = find(headshape.pos(:,3)<3);
        if isempty(count_facialpoints)
            disp('CANNOT FIND ANY FACIAL POINTS - COREG BY ICP MAY BE INACCURATE');
        else
            facialpoints = headshape.pos(count_facialpoints,:,:);
            rrr = 1:4:length(facialpoints);
            facialpoints = facialpoints(rrr,:); clear rrr;
        end

        % Remove facial points for now
        headshape.pos(count_facialpoints,:) = [];

        % Create mesh out of headshape downsampled to x points specified in the
        % function call
        cfg.numvertices = numvertices;
        cfg.method = 'headshape';
        cfg.headshape = headshape.pos;
        mesh = ft_prepare_mesh(cfg, headshape);

        % Replace the headshape info with the mesh points
        headshape.pos = mesh.pos;

        % Create figure for quality checking
        figure; subplot(2,2,1);ft_plot_mesh(mesh); hold on;
        title('Downsampled Mesh');
        view(0,0);
        subplot(2,2,2);ft_plot_mesh(headshape); hold on;
        title('Downsampled Headshape View 1');
        view(0,0);
        subplot(2,2,3);ft_plot_mesh(headshape); hold on;
        title('Downsampled Headshape View 2');
        view(90,0);
        subplot(2,2,4);ft_plot_mesh(headshape); hold on;
        title('Downsampled Headshape View 3');
        view(180,0);
        print('headshape_quality','-dpdf');

        % Add the facial points back in (default) or leave out if user specified
        % 'no' in function call
        if strcmp(include_facial_points,'yes')
            try
                % Add the facial info back in
                headshape.pos = vertcat(headshape.pos,facialpoints);
            catch
                disp('Cannot add facial info back into headshape');
            end
        else
            headshape.pos = headshape.pos;
            disp('Not adding facial points back into headshape');
        end

        %Add in names of the fiducials from the sensor
        headshape.fid.label = {'NASION','LPA','RPA'};

        % Convert fiducial points to BESA
        %         headshape.fid.pos = cat(2,fliplr(headshape.fid.pos(:,1:2)),headshape.fid.pos(:,3));
        %         headshape.fid.pos(:,2) = headshape.fid.pos(:,2).*-1;

        % Plot for quality checking
        figure;%ft_plot_sens(sensors) %plot channel position : between the 1st and 2nd coils
        ft_plot_headshape(headshape) %plot headshape
        view(0,0);
        print('headshape_quality2','-dpdf');

        % Export filename
        headshape_downsampled = headshape;

    end

end
