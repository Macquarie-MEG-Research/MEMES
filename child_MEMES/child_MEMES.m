function child_MEMES(dir_name,grad_trans,headshape_downsampled,...
    path_to_MRI_library)
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
% - elpfile             = path to elp file
% - hspfile             = path to hsp file
% - confile             = path to con file
% - mrkfile             = path to mrk file
% - path_to_MRI_library = path to HCP MRI library
% - bad_coil            = list of bad coils (up to length of 2). Enter as:
%                         {'LPAred','RPAyel','PFblue','LPFwh','RPFblack'}
%
%%%%%%%%%%%%%%%%%%
% Variable Inputs:
%%%%%%%%%%%%%%%%%%
%
% - transform_sensors = 'yes' or 'no' (default = 'no')
% - include_face      = inclue facial points acquired during head
%                       digitisation ('yes' = default)
% - sens_coreg_method = method used to realign MEG sensors based on 5
%                       marker coils. Use 'rot3dfit' or 'icp'. For some
%                       reason the usual rot3dfit method seems to fail
%                       sometimes. Try using 'icp' in this case...
%
%%%%%%%%%%%
% Outputs:
%%%%%%%%%%%
%
% - grad_trans              = sensors transformed to correct
% - shape                   = headshape and fiducial information
% - headshape_downsampled   = headshape downsampled to 100 points with facial
%                           information preserved
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

% This script estimates headmodel and sourcemodel for MEG sourcespace
% analysis by matching polhemus points to a databsse of developmental
% template MRIs, via an iterative closest point alogorithm. The template
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


%% ICP loop

age_list = {'2-5','3-0','4-0','4-5','5-0','5-5','6-0','6-5','7-0',...
    '7-5','8-0','8-5','9-0','9-5','10-0','10-5'};

% Error term variable
error_term = zeros(1,length(age_list));
% Variable to hold the transformation matrices
trans_matrix_library = [];

% For every member of the age list..
for m = 1:length(age_list)
    
    % Load the mesh
    load([path_to_MRI_library age_list{m} '/mesh.mat'],'mesh');
    
    % Number of iterations for the ICP algorithm
    numiter = 30;
    
    % Perform ICP
    [R, t, err, dummy, ~] = icp(mesh.pos', headshape_downsampled.pos', ...
        numiter, 'Minimize', 'plane', 'Extrapolation', true,'WorstRejection', 0.1);
    
    % Add error to error_term list
    error_term(m) = err(end);
    
    % Get the transformation matrix
    trans_matrix = inv([real(R) real(t);0 0 0 1]);
    
    % Add transformation matrix to trans_matrix_library
    trans_matrix_library{m} = trans_matrix;
    
    % Clear mesh for next loop and display message
    clear mesh R t err trans_matrix
    fprintf('Completed age %s\n',age_list{m});
    
end

fprintf('\nFinished the iterations\n');

%% Create figure to show good, OK and bad examples of the iterations
fprintf('\nFinding good, OK and bad examples\n');

error_term_sorted = sort(error_term, 'ascend');
winners = find(ismember(error_term,error_term_sorted(1)));
middles = find(ismember(error_term,error_term_sorted(round(length(error_term)/2))));
losers = find(ismember(error_term,error_term_sorted(end)));

concat = [winners middles losers];

% Create figure to summarise the losers,middles and winners
figure;
for i = 1:3
    
    load([path_to_MRI_library age_list{concat(i)} '/mesh.mat'])
    
    mesh_spare = mesh;
    mesh_spare.pos = ft_warp_apply(trans_matrix_library{(concat(i))}, mesh_spare.pos);
    
    subplot(1,3,i)
    ft_plot_mesh(mesh_spare,'facecolor',[238,206,179]./255,'EdgeColor','none','facealpha',0.8); hold on;
    camlight; hold on; view([-270,-10]);
    if ismember(i,1)
        title(sprintf('BEST: %.2f', error_term((concat(i)))));
    elseif ismember(i,2)
        title(sprintf('MIDDLE: %.2f', error_term((concat(i)))));
    elseif ismember(i,3)
        title(sprintf('WORST: %.2f', error_term((concat(i)))));
    end
    
    ft_plot_headshape(headshape_downsampled);
end

%% Create graph of error over age
figure;
plot(1:length(age_list),error_term,'LineWidth',3);
xticks([1:1:16]);
xticklabels(age_list);
ylabel('Error','FontSize',15);
xlabel('Age Template','FontSize',15);
print('error_age','-dpng','-r100');

%% Determine the winning MRI and load the facial mesh
winner = find(error_term == min(min(error_term)));
fprintf('\nThe winning MRI is number %d of %d\n',winner,length(age_list));
trans_matrix = trans_matrix_library{winner};

% Save the trans matrix to disk
save trans_matrix trans_matrix

% Get facial mesh of winner
load([path_to_MRI_library age_list{winner} '/mesh.mat'])
mesh_spare = mesh;
mesh_spare.pos = ft_warp_apply(trans_matrix, mesh_spare.pos);

%% Create Headmodel (in mm)
fprintf('Creating Headmodel in mm\n');

path_to_headmodel = [path_to_MRI_library age_list{winner} '/headmodel.mat'];

% Load headmodel
load(path_to_headmodel);

% Transform (MESH --> coreg via ICP adjustment)
headmodel = ft_transform_vol(trans_matrix,headmodel);

% Save
save headmodel headmodel

figure;
ft_plot_vol(headmodel);
ft_plot_headshape(headshape_downsampled); view([0,0]);

%% Create Sourcemodel (in mm)
fprintf('Creating an 8mm Sourcemodel in mm\n');

% This loads the 8mm one, but you can change to 5mm
path_to_sourcemodel = [path_to_MRI_library age_list{winner} '/sourcemodel3d_8mm.mat'];

% Load sourcemodel and convert to 'mm'
load(path_to_sourcemodel);
sourcemodel3d = ft_convert_units(sourcemodel3d,'mm');

% Transform (MESH --> coreg via ICP adjustment)
sourcemodel3d.pos = ft_warp_apply(trans_matrix,sourcemodel3d.pos);

% Save
save sourcemodel3d sourcemodel3d

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

