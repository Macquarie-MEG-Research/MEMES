
function MEMES2(dir_name,elpfile,hspfile,confile,mrkfile,path_to_MRI_library,mesh_library,initial_mri_realign,bad_coil,method,scaling)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MRI Estimation for MEG Sourcespace (MEMES)
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
% - mesh_library        = mesh library (in mm) created from HCP MRIs
% - initial_mri_realign = transform for initial realigning estimate
% - bad_coil            = list of bad coils 
% - method              = method for creating pseudo head- and
%                       source-model: 'best' or 'average'
% - scaling             = scaling factor range applied to MRIs 
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
% - mesh                    = coregistered 3D cortical mesh (4000 vertices
%                           per hemisphere
% - sourcemodel3d           = 8mm sourcemodel warped to MNI space
% - headmodel               = singleshell headmodel (10000 vertices) 

%%%%%%%%%%%%%%%%%%%%%
% Other Information:
%%%%%%%%%%%%%%%%%%%%%

% Example function call:
% MEMES2(dir_name,elpfile,hspfile,confile,mrkfile,path_to_MRI_library,...
% mesh_library,initial_mri_realign,'','best',[0.98:0.01:1.02])

% Here I have introduced a variable scaling parameter for the MRIs to 
% help with coregistration. For example to apply -2% to +2% scaling to 
% every MRI specify: scaling = [0.98:0.01:1.2]. 

% However NOTE: the more scaling factors you apply the longer it will take

% This script is specific for adult MEG-HCP data (95 MRIs, headmodels and
% sourcemodels), but could easily be adapted for other datasets. 

% The transforms are a bit ad-hoc (BTI --> SPM ; rotate 90deg ; initial
% ICP realign; subject-specific scaling ; final subject-specific ICP).

% To adapt this be careful about the first three transforms or create
% example head and sourcemodels from initial realign rather than faff 
% about later on...
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nThis is MEMES v0.2\n\nMake sure you have asked Robert for the corrected HCP library\n\n');
% Check inputs
disp('Performing input check');
assert(length(bad_coil)<3,'You need at least 3 good coils for accurate alignment\n');
assert(path_to_MRI_library(end) == '/','path_to_MRI_library needs to end with /\n');
%assert(method == 'average','method = average is not yet supported. Use best\n');

if length(scaling) == 1
    scaling = 1;
end


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
    '891667';'898176';'912447';'917255';'990366'};

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
if isempty(bad_coil)
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
    
    % Identify the bad coil
    badcoilpos = find(ismember(shape.fid.label,bad_coil{1}));
    
    % Take away the bad marker
    marker_order = [2 3 1 4 5];
    markers                     = mrk.fid.pos(marker_order,:);%reorder mrk to match order in shape
    % Now take out the bad marker when you realign
    markers(badcoilpos-3,:) = [];
    
    fids_2_use = shape.fid.pnt(4:end,:); fids_2_use(badcoilpos-3,:) = [];
    
    % If there is a second bad coil remove this now
    if length(bad_coil) == 2
        badcoilpos2 = find(ismember(shape.fid.label,bad_coil{2}));
        markers(badcoilpos2-4,:) = [];
        fids_2_use(badcoilpos2-4,:) = [];
    end
    
    [R,T,Yf,Err]                = rot3dfit(markers,fids_2_use);%calc rotation transform
    meg2head_transm             = [[R;T]'; 0 0 0 1];%reorganise and make 4*4 transformation matrix
    
    disp('Performing re-alignment');
    grad_trans                  = ft_transform_geometry_PFS_hacked(meg2head_transm,grad_con); %Use my hacked version of the ft function - accuracy checking removed not sure if this is good or not
    grad_trans.fid              = shape; %add in the head information
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

fprintf('Downsampling headshape information to %d points whilst preserving facial information\n'...
    ,100); 
headshape_downsampled = downsample_headshape(hspfile,100);
headshape_downsampled = ft_convert_units(headshape_downsampled,'mm'); %in mm
disp('Saving headshape downsampled');
save headshape_downsampled headshape_downsampled

figure;ft_plot_mesh(mesh_library{1,1}); ft_plot_headshape(headshape_downsampled);

%% Perform ICP

% Error term variable
error_term = zeros(1,length(mesh_library));
% Variable to hold the transformation matrices
trans_matrix_library = []; 
scaling_factor_all = zeros(1,length(mesh_library));
count = 1;

for m = 1:length(mesh_library)
    
    numiter = 30; count2 = 1;
    
    trans_matrix_temp = []; error_2 = [];
    
    for scale = scaling
        fprintf('Completed iteration %d of %d ; %d of %d MRIs\n',count2,length(scaling),m,length(mesh_library));
        mesh_coord_scaled = ft_warp_apply([scale 0 0 0;0 scale 0 0; 0 0 scale 0; 0 0 0 1],mesh_library{m}.pos);
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

end

fprintf(' Finished the iterations\n');

%% Make pretty figure
fprintf('\n Finding good, OK and bad examples\n');

error_term_sorted = sort(error_term, 'ascend');
winners = find(ismember(error_term,error_term_sorted(1:3)));
middles = find(ismember(error_term,error_term_sorted(46:48)));
losers = find(ismember(error_term,error_term_sorted(end-2:end)));

concat = [winners middles losers];

% Create figure to summarise the losers,middles and winners
figure;
for i = 1:9
    
    mesh_spare = mesh_library{(concat(i))};
    mesh_spare.pos = ft_warp_apply([scaling_factor_all(concat(i)) 0 0 0;0 scaling_factor_all(concat(i)) 0 0; 0 0 scaling_factor_all(concat(i)) 0; 0 0 0 1],mesh_spare.pos);

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
    
    if i == 9
        print('best_middle_worst_examples','-dpng','-r100');
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
        
        %% Calculate headmodel and sourcemodel over first n MRIs
        
        average_over_n = 10; average_sourcemodel_all = []; average_mri_all = [];
        
        average_n = find(ismember(error_term,error_term_sorted(1:average_over_n)));
        
        
        for rep = 1:average_over_n
            
            fprintf(' Getting sourcemodel %d of %d \n',rep,average_over_n);
            
            % Get MRI of winning subject
            mri_file = [path_to_MRI_library subject{average_n(rep)} '/MEG/anatomy/T1w_acpc_dc_restore.nii.gz'];
            
            mri_orig                    = ft_read_mri(mri_file); % in mm, read in mri from DICOM
            mri_orig = ft_convert_units(mri_orig,'cm'); mri_orig.coordsys = 'neuromag';
            
            mri_orig.transform = initial_mri_realign{average_n(rep)};
            mri_realigned = mri_orig;
            
            mri_realigned = ft_transform_geometry(trans_matrix_library{average_n(rep)},mri_realigned);
            
            average_mri_all{rep} = mri_realigned;
            
            %% Create Sourcemodel (in cm)
            fprintf(' Creating Sourcemodel in cm\n');
            
            % Get transformation matrix to convert BTI to SPM
            % (N.B. subject 1 not working?)
            path_to_transform   = [path_to_MRI_library subject{average_n(rep)} '/MEG/anatomy/' subject{average_n(rep)} '_MEG_anatomy_transform.txt'];
            
            txt = textscan(fopen(path_to_transform),'%s%s%s%s\n'); %open text file
            transform.bti2spm = zeros(4);
            fclose('all');
            
            for reps = 1:4 %some jiggery-pokery to extract transform matrix from text file
                transform.bti2spm(:,reps) = str2num(strjoin(txt{1,reps}(32:35)))';
            end
            
            clear reps
            
            path_to_sourcemodel = [path_to_MRI_library subject{average_n(rep)} '/MEG/anatomy/' subject{average_n(rep)} '_MEG_anatomy_sourcemodel_3d8mm.mat'];
            
            % Load sourcemodel and convert to 'mm'
            load(path_to_sourcemodel); sourcemodel3d = ft_convert_units(sourcemodel3d,'mm');
            
            % Transform BTI --> MRI-space
            sourcemodel3d.pos = ft_warp_apply(transform.bti2spm,sourcemodel3d.pos);
            
            % Convert to cm
            sourcemodel3d = ft_convert_units(sourcemodel3d,'cm');
            
            % Transform 1 (MESH --> coreg via manual marking of fiducial points)
            sourcemodel3d.pos = ft_warp_apply(inv(mri_orig.transform),sourcemodel3d.pos);
            sourcemodel3d.pos = ft_warp_apply(initial_mri_realign{average_n(rep)},sourcemodel3d.pos);
            
            %transform 2 (MESH --> coreg via ICP adjustment)
            sourcemodel3d.pos = ft_warp_apply(trans_matrix_library{average_n(rep)},sourcemodel3d.pos);
            
            average_sourcemodel_all{rep} = sourcemodel3d;
            
        end
        
        pos_headmodel = []; pos_sourcemodel = []; tri_headmodel = [];
        
        for rep = 1:average_over_n
            pos_headmodel(rep,:,:,:) = average_headmodel_all{1,rep}.bnd.pnt;
            tri_headmodel(rep,:,:,:) = average_headmodel_all{1,rep}.bnd.tri;
            pos_sourcemodel(rep,:,:,:) = average_sourcemodel_all{1,rep}.pos;
        end
        
        average_headmodel = average_headmodel_all{1};
        average_headmodel.bnd.pnt = squeeze(mean(pos_headmodel));
        average_headmodel.bnd.tri = squeeze(mean(tri_headmodel));
        
        average_sourcemodel = average_sourcemodel_all{1};
        average_headmodel.pos = squeeze(mean(pos_sourcemodel));
        
        figure;ft_plot_vol(average_headmodel,  ...
                'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.6; camlight;
            ft_plot_mesh(average_sourcemodel.pos(average_sourcemodel.inside,:),'vertexsize',2);
            %ft_plot_sens(grad_trans, 'style', 'r*')
            %ft_plot_headshape(headshape_downsampled) %plot headshape

        volume(average_headmodel_all{1,1}.bnd.pnt)
        
        
        view_angle = [0 90 180 270];
        figure; hold on;
        for rep = 1:4
            subplot(2,2,rep);%ft_plot_vol(average_headmodel,  ...
                %'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.6; camlight;
            ft_plot_mesh(average_sourcemodel.pos(average_sourcemodel.inside,:),'vertexsize',2);
            ft_plot_sens(grad_trans, 'style', 'r*')
            ft_plot_headshape(headshape_downsampled) %plot headshape
            view([view_angle(rep),0]);
        end
        
        %% SAVE
        fprintf('\nSaving the necessary data\n');
        
        save average_headmodel average_headmodel
        save grad_trans grad_trans
        save average_sourcemodel average_sourcemodel
        
        
        fprintf('\nCOMPLETED - check the output for quality control\n');
        
        
        
    case 'best'
        
        winner = find(error_term == min(min(error_term)));
        fprintf('\nThe winning MRI is number %d of %d\n',winner,length(mesh_library));
        trans_matrix = trans_matrix_library{winner};
        
        % Get facial mesh of winner
        mesh_spare = mesh_library{1,winner};
        mesh_spare.pos = ft_warp_apply([scaling_factor_all(winner) 0 0 0;0 ...
            scaling_factor_all(winner) 0 0; 0 0 scaling_factor_all(winner) 0;...
            0 0 0 1],mesh_spare.pos);
        mesh_spare.pos = ft_warp_apply(trans_matrix, mesh_spare.pos);
        
        % Get MRI of winning subject
        mri_file = [path_to_MRI_library subject{winner} '/MEG/anatomy/T1w_acpc_dc_restore.nii.gz'];
        
        mri_orig                    = ft_read_mri(mri_file); % in mm, read in mri from DICOM
        mri_orig = ft_convert_units(mri_orig,'mm'); mri_orig.coordsys = 'neuromag';
        
        %% Create Headmodel (in cm)
        fprintf(' Creating Headmodel in mm\n');
        
        path_to_headmodel = [path_to_MRI_library subject{winner} '/MEG/anatomy/' subject{winner} '_MEG_anatomy_headmodel.mat'];
        
        % Load headmodel
        load(path_to_headmodel);
        
        % Get transformation matrix to convert BTI to SPM
        % This requires some re-ordering of the original HCP-MEG dataset
        % Essentially you just need to apply the transform.bti2spm to the
        % headmodel and sourcemodel (later)
        path_to_transform   = [path_to_MRI_library subject{winner} '/MEG/anatomy/' subject{winner} '_MEG_anatomy_transform.txt'];
        
        txt = textscan(fopen(path_to_transform),'%s%s%s%s\n'); %open text file
        indx_txt = find(contains(txt{1,1},'transform.bti2spm')); %find index of bti2spm
        transform.bti2spm = zeros(4);
        
        for rep = 1:4 %some jiggery-pokery to extract transform matrix from text file
            transform.bti2spm(:,rep) = str2num(strjoin(txt{1,rep}(indx_txt+1:indx_txt+4)))';
        end
        
        % Transform BTI --> MRI-space
        headmodel.bnd.pnt = ft_warp_apply(transform.bti2spm,headmodel.bnd.pnt);
        
        % Convert to mm
        headmodel = ft_convert_units(headmodel,'mm');
        
        % 90deg rotation matrix
        rmatx = [6.12323399573677e-17,1,0,0;-1,6.12323399573677e-17,...
            0,0;0,0,1,0;0,0,0,1];
        
        % Transform 1 (MESH --> coreg via manual marking of fiducial points)        
        headmodel.bnd.pnt = ft_warp_apply(rmatx,headmodel.bnd.pnt);
        headmodel.bnd.pnt = ft_warp_apply(initial_mri_realign{winner},headmodel.bnd.pnt);
        
        % Scale
        headmodel.bnd.pnt = ft_warp_apply([scaling_factor_all(winner) 0 0 0;0 ...
            scaling_factor_all(winner) 0 0; 0 0 scaling_factor_all(winner) 0; 0 0 0 1],...
            headmodel.bnd.pnt);

        %transform 2 (MESH --> coreg via ICP adjustment)
        headmodel.bnd.pnt = ft_warp_apply(trans_matrix,headmodel.bnd.pnt);
        
        figure;
        ft_plot_vol(headmodel);
        ft_plot_headshape(headshape_downsampled);
        
        %% Create Sourcemodel (in mm)
        fprintf(' Creating an 8mm Sourcemodel in mm\n');
        
        % This loads the 8mm one, but you can change to 5mm
        path_to_sourcemodel = [path_to_MRI_library subject{winner} '/MEG/anatomy/' subject{winner} '_MEG_anatomy_sourcemodel_3d8mm.mat'];
        
        % Load sourcemodel and convert to 'mm'
        load(path_to_sourcemodel); sourcemodel3d = ft_convert_units(sourcemodel3d,'mm');
        
        % Transform BTI --> MRI-space
        sourcemodel3d.pos = ft_warp_apply(transform.bti2spm,sourcemodel3d.pos);
        
        % Convert to mm
        sourcemodel3d = ft_convert_units(sourcemodel3d,'mm');
        
        % Transform 1 (MESH --> coreg via manual marking of fiducial points)
        sourcemodel3d.pos = ft_warp_apply(rmatx,sourcemodel3d.pos);
        sourcemodel3d.pos = ft_warp_apply(initial_mri_realign{winner},sourcemodel3d.pos);
        
        % Scale
        sourcemodel3d.pos = ft_warp_apply([scaling_factor_all(winner) 0 0 0;0 scaling_factor_all(winner) 0 0; 0 0 scaling_factor_all(winner) 0; 0 0 0 1],sourcemodel3d.pos);
        
        %transform 2 (MESH --> coreg via ICP adjustment)
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
        
        %% Create coregistered 3D cortical mesh
        mesh = ft_read_headshape({[path_to_MRI_library ...
            subject{winner} '/MEG/anatomy/' subject{winner} '.L.midthickness.4k_fs_LR.surf.gii'],...
            [path_to_MRI_library subject{winner} '/MEG/anatomy/' subject{winner} ...
            '.R.midthickness.4k_fs_LR.surf.gii']});
        
        mesh = ft_convert_units(mesh,'mm');
        
        % Transform 1 (MESH --> coreg via manual marking of fiducial points)
        mesh.pos = ft_warp_apply(rmatx,mesh.pos);
        mesh.pos = ft_warp_apply(initial_mri_realign{winner},mesh.pos);
        
        % Scale
        mesh.pos = ft_warp_apply([scaling_factor_all(winner) 0 0 0;0 scaling_factor_all(winner) 0 0; 0 0 scaling_factor_all(winner) 0; 0 0 0 1],mesh.pos);

        %transform 2 (MESH --> coreg via ICP adjustment)
        mesh.pos = ft_warp_apply(trans_matrix,mesh.pos);
        
        %ft_determine_coordsys(mri_realigned2,'interactive','no'); hold on;
        figure;ft_plot_sens(grad_trans);
        ft_plot_headshape(headshape_downsampled) %plot headshape
        ft_plot_mesh(mesh,'facealpha',0.8); camlight; hold on; view([100 4]);
        print('headmodel_3D_cortical_mesh_quality','-dpng');
        
        %% SAVE
        fprintf('\nSaving the necessary data\n');
        
        save headmodel headmodel
        save trans_matrix trans_matrix
        save grad_trans grad_trans
        save sourcemodel3d sourcemodel3d
        save mesh mesh
        
        
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

    function [headshape_downsampled] = downsample_headshape(path_to_headshape,numvertices)
        % Get headshape
        headshape = ft_read_headshape(path_to_headshape);
        % Convert to cm
        headshape = ft_convert_units(headshape,'cm');
        % Convert to BESA co-ordinates
%         headshape.pos = cat(2,fliplr(headshape.pos(:,1:2)),headshape.pos(:,3));
%         headshape.pos(:,2) = headshape.pos(:,2).*-1;
        
        % Get indices of facial points (up to 4cm above nasion)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Is 4cm the correct distance?
        % Possibly different for child system?
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        count_facialpoints = find(headshape.pos(:,3)<4);
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
        
        % Add the facial info back in
        headshape.pos = vertcat(headshape.pos,facialpoints);
        % Add in names of the fiducials from the sensor
        headshape.fid.label = {'NASION','LPA','RPA'};
        
        % Convert fiducial points to BESA
%         headshape.fid.pos = cat(2,fliplr(headshape.fid.pos(:,1:2)),headshape.fid.pos(:,3));
%         headshape.fid.pos(:,2) = headshape.fid.pos(:,2).*-1;
        
        % Plot for quality checking
        figure;ft_plot_headshape(headshape) %plot headshape
        view(0,0);
        print('headshape_quality2','-dpdf');
        
        % Export filename
        headshape_downsampled = headshape;
        
    end
    
    function [headshape_downsampled] = downsample_headshape_noface(path_to_headshape,numvertices,sensors)
        % Get headshape
        headshape = ft_read_headshape(path_to_headshape);
        % Convert to cm
        headshape = ft_convert_units(headshape,'cm');
        % Convert to BESA co-ordinates
%         headshape.pos = cat(2,fliplr(headshape.pos(:,1:2)),headshape.pos(:,3));
%         headshape.pos(:,2) = headshape.pos(:,2).*-1;
        
        % Get indices of facial points (up to 4cm above nasion)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Is 4cm the correct distance?
        % Possibly different for child system?
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        count_facialpoints = find(headshape.pos(:,3)<4);
        if isempty(count_facialpoints)
            disp('CANNOT FIND ANY FACIAL POINTS');
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
        
        % Add in names of the fiducials from the sensor
        headshape.fid.label = {'NASION','LPA','RPA'};
        
        % Convert fiducial points to BESA
%         headshape.fid.pos = cat(2,fliplr(headshape.fid.pos(:,1:2)),headshape.fid.pos(:,3));
%         headshape.fid.pos(:,2) = headshape.fid.pos(:,2).*-1;
%         
        % Plot for quality checking
        figure;ft_plot_sens(sensors) %plot channel position : between the 1st and 2nd coils
        ft_plot_headshape(headshape) %plot headshape
        view(0,0);
        print('headshape_quality2','-dpdf');
        
        % Export filename
        headshape_downsampled = headshape;
        
    end
    
end



