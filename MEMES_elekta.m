%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MRI Estimation for MEG Sourcespace (MEMES) using Elekta rather than
% Yokogawa data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MEMES(dir_name,fif_file,path_to_MRI_library,mesh_library,initial_mri_realign)

% list of HCP subjects
subject = {'100307';'102816';'104012';'105923';'106521';'108323';...
    '109123';'111514';'112920';'113922';'116524';'116726';'125525';...
    '133019';'140117';'146129';'149741';'151526';'153732';'154532';...
    '156334';'158136';'162026';'162935';'164636';'166438';'169040';...
    '172029';'174841';'175237';'175540';'177746';'179245';'181232';...
    '182840';'185442';'187547';'189349';'191033';'191437';'191841';...
    '192641';'195041';'198653';'200109';'204521';'205119';'212318';...
    '212823';'214524';'221319';'223929';'233326';'248339';'250427';...
    '255639';'257845';'283543';'287248';'293748';'352132';'352738';...
    '353740';'358144';'406836';'433839';'500222';'512835';'555348';...
    '559053';'568963';'581450';'599671';'601127';'660951';'662551';...
    '665254';'667056';'679770';'680957';'706040';'707749';'715950';...
    '725751';'735148';'783462';'814649';'825048';'872764';'877168';...
    '891667';'898176';'912447';'917255';'990366'}

% CD to right place
cd(dir_name); fprintf('\n CDd to the right place\n');

%% Read headshape and get sensor information
% If you have too many headshape points you may need to downsample a little
% - guess it wouldn't be too difficult to write a script...
headshape  = ft_read_headshape(fif_file); figure;ft_plot_headshape(headshape);
grad_trans = ft_read_sens(fif_file);

%% Perform ICP

% Error term variable
error_term = zeros(1,length(mesh_library));
% Variable to hold the transformation matrices
trans_matrix_library = []; count = 1;

for m = 1:length(mesh_library)
    
    fprintf('Completed iteration %d of %d\n',m,length(mesh_library));
    
    numiter = 50;
    
    % Perform ICP
    [R, t, err, dummy, info] = icp(mesh_library{m}.pos', headshape.pos', numiter, 'Minimize', 'plane', 'Extrapolation', true,'WorstRejection', 0.05);
    
    % Add error to error_term
    error_term(m) = err(end);
    
    % Add transformation matrix to trans_matrix_library
    trans_matrix_library{m} = inv([real(R) real(t);0 0 0 1]);
end


%% Make pretty figure
error_term_sorted = sort(error_term, 'descend');
losers = find(ismember(error_term,error_term_sorted(1:3)));
middles = find(ismember(error_term,error_term_sorted(46:48)));
winners = find(ismember(error_term,error_term_sorted(end-2:end)));

concat = [winners middles losers];

% Create figure to summarise the losers,middles and winners
figure;
for i = 1:9
    
    mesh_spare = mesh_library{(concat(i))};
    mesh_spare.pos = ft_warp_apply(trans_matrix_library{(concat(i))}, mesh_spare.pos);
    
    subplot(3,3,i)
    ft_plot_mesh(mesh_spare,'facecolor',[238,206,179]./255,'EdgeColor','none','facealpha',0.8); hold on;
    camlight; hold on; view([-180,-10]);
    if ismember(i,1:3)
        title(sprintf('BEST: %d', error_term((concat(i)))));
    elseif ismember(i,4:6)
        title(sprintf('MIDDLE: %d', error_term((concat(i)))));
    elseif ismember(i,7:9)
        title(sprintf('WORST: %d', error_term((concat(i)))));
    end
    
    ft_plot_headshape(headshape);
    
    if i == 9
        print('best_middle_worst_examples','-dpdf','-r200');
    end
end

%% Use the best for to create a source model for MEG source analysis

winner = find(error_term == min(min(error_term)));
fprintf('\nThe winning MRI is number %d of %d\n',winner,length(mesh_library));
trans_matrix = trans_matrix_library{winner};

% Create figure to show ICP fit
mesh_spare = mesh_library{winner};
mesh_spare.pos = ft_warp_apply(trans_matrix, mesh_spare.pos);

figure;ft_plot_mesh(mesh_spare,'facecolor',[238,206,179]./255,'EdgeColor','none','facealpha',0.8); hold on;
camlight; hold on; view([-180,-10]);
title(error_term(winner));
ft_plot_headshape(headshape);

print('winning_sourcemodel','-dpdf','-r200');

try
    % % Make fancy video
    c = datestr(clock); %time and date
    
    figure;
    ft_plot_mesh(mesh_spare,'facecolor',[238,206,179]./255,'EdgeColor','none','facealpha',0.8); hold on;
    camlight; hold on;
    ft_plot_headshape(headshape_downsampled); title(sprintf('%s.   Error of ICP fit = %d' , c, error_term(winner)));
    OptionZ.FrameRate=15;OptionZ.Duration=5.5;OptionZ.Periodic=true;
    CaptureFigVid([0,0; 360,0], 'ICP_quality',OptionZ)
    
catch
    fprintf('You need CaptureFigVid in your path for fancy videos\n');
end

% Get MRI of winning subject
mri_file = [path_to_MRI_library subject{winner} '/MEG/anatomy/T1w_acpc_dc_restore.nii.gz'];

mri_orig                    = ft_read_mri(mri_file); % in mm, read in mri from DICOM
mri_orig = ft_convert_units(mri_orig,'cm'); mri_orig.coordsys = 'neuromag';

mri_orig.transform = initial_mri_realign{winner};
mri_realigned = mri_orig; clear mri_orig;

% Apply transformation matrix
mri_realigned = ft_transform_geometry(trans_matrix,mri_realigned);

% Create figure for quality checking
ft_determine_coordsys(mri_realigned,'interactive','no'); hold on;
ft_plot_headshape(headshape); view([102 5]);

% Segment
fprintf('\nSegmenting the MRI... This may take a while...\n');
cfg           = [];
cfg.output    = 'brain';
mri_segmented  = ft_volumesegment(cfg, mri_realigned);

% Create singleshell headmodel
cfg = [];
cfg.method='singleshell';

headmodel_singleshell = ft_prepare_headmodel(cfg, mri_segmented); % in cm, create headmodel

figure;ft_plot_headshape(headshape) %plot headshape
ft_plot_sens(grad_trans, 'style', 'k*')
ft_plot_vol(headmodel_singleshell,  'facecolor', 'cortex', 'edgecolor', 'none');alpha 1; camlight
view([90,0]); title('After Coreg');
print('headmodel_quality','-dpdf');

fprintf('\nSaving the necessary data\n');

save headmodel_singleshell headmodel_singleshell
save mri_realigned mri_realigned
save trans_matrix trans_matrix
save grad_trans grad_trans

fprintf('\nCOMPLETED - check the output for quality control\n');

end