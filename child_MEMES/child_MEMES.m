function child_MEMES(dir_name,elpfile,hspfile,confile,...
    mrkfile,path_to_MRI_library,bad_coil,varargin)
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

% Deal with variable inputs
if isempty(varargin)
    transform_sensors = 'no';
    include_face        = 'yes';
    sens_coreg_method = 'rot3dfit';
else
    transform_sensors = varargin{1};
    include_face      = varargin{2};
    sens_coreg_method = varargin{3};
end

fprintf(['\nThis is MEMES for child MEG data v0.1\n\nMake sure you have ',...
    'asked Robert for the mesh, headmodel and sourcemodel library\n\n']);
% Check inputs
disp('Performing input check');
assert(length(bad_coil)<3,'You need at least 3 good coils for accurate alignment\n');
assert(path_to_MRI_library(end) == '/','path_to_MRI_library needs to end with /\n');

% CD to right place
cd(dir_name); fprintf('\nCDd to the right place\n');

% Get Polhemus Points
disp('Reading elp and hspfile');
[shape]  = parsePolhemus(elpfile,hspfile);
shape   = ft_convert_units(shape,'mm');
%
% Read the grads from the con file
disp('Reading Sensors');
grad_con = ft_read_sens(confile); %in cm, load grads
grad_con = ft_convert_units(grad_con,'mm'); %in mm

% Perform MEG sensor transform if specified in the
switch transform_sensors
    case 'no'
        % If the sensors have aready been transformed using MEG160 then don't do
        % the transform(!)
        disp('Assuming the data has already been transformed with MEG 160');
        grad_trans = grad_con;
        save grad_trans grad_trans
        
    case 'yes'
        % If user specified the sensors need to be transformed
        disp('Transforming the data based on data in the mrk file');
        
        % Read mrk_file
        disp('Reading the mrk file');
        mrk      = ft_read_headshape(mrkfile,'format','yokogawa_mrk');
        mrk      = ft_convert_units(mrk,'mm'); %in mm
        
        %% Perform Realighment Using Paul's Fancy Functions
        if strcmp(bad_coil,'')
            disp('NO BAD MARKERS');
            markers                     = mrk.fid.pos([2 3 1 4 5],:);%reorder mrk to match order in shape
            
            % If the user specifed to use the icp sensor coregistratio approach use
            % this...
            if strcmp(sens_coreg_method,'icp')
                fids_2_use = shape.fid.pnt(4:end,:);
                % For some reason this works better with only 3 points... check to
                % make sure this works for all?
                [R, T, err, dummy, info]    = icp(fids_2_use(1:5,:)',...
                    markers(1:5,:)',100,'Minimize', 'point');
                meg2head_transm             = [[R T]; 0 0 0 1];%reorganise and make 4*4 transformation matrix
                % Otherwise use the original rot3dfit method
            else
                [R,T,Yf,Err]                = rot3dfit(markers,shape.fid.pnt(4:end,:));%calc rotation transform
                meg2head_transm             = [[R;T]'; 0 0 0 1];%reorganise and make 4*4 transformation matrix
            end
            
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
                
                % Now take out the bad coil from the shape variable to prevent bad
                % plotting - needs FIXING for 2 markers (note: Dec 18)
                
                grad_trans.fid              = shape; %add in the head information
            else
                [R,T,Yf,Err]                = rot3dfit(markers,fids_2_use);%calc rotation transform
                meg2head_transm             = [[R;T]'; 0 0 0 1];%reorganise and make 4*4 transformation matrix
                grad_trans                  = ft_transform_geometry_PFS_hacked(meg2head_transm,grad_con); %Use my hacked version of the ft function - accuracy checking removed not sure if this is good or not
                
                % Now take out the bad coil from the shape variable to prevent bad
                % plotting
                shape.fid.pnt(badcoilpos+3,:) = [];
                shape.fid.label(badcoilpos+3,:) = [];
                grad_trans.fid              = shape; %add in the head information
            end
            
            
        end
        
        % Save grad_trans
        fprintf('Saving grad_trans\n');
        save grad_trans grad_trans
        
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

fprintf(['Downsampling headshape information to %d points whilst ',...
    'preserving facial information\n'],100);
headshape_downsampled = downsample_headshape_child(hspfile,100,include_face);
headshape_downsampled = ft_convert_units(headshape_downsampled,'mm'); %in mm
disp('Saving headshape downsampled');
save headshape_downsampled headshape_downsampled

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

    function [headshape_downsampled] = downsample_headshape_child(path_to_headshape,...
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
        % Save a version for later
        headshape_orig = headshape;
        
        % Get indices of facial points (up to 3cm above nasion)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Is 3cm the correct distance?
        % Possibly different for child system?
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        count_facialpoints = find(headshape.pos(:,3)<3 & headshape.pos(:,1)>1);
        if isempty(count_facialpoints)
            disp('CANNOT FIND ANY FACIAL POINTS - COREG BY ICP MAY BE INACCURATE');
        else
            facialpoints = headshape.pos(count_facialpoints,:,:);
            rrr = 1:4:length(facialpoints);
            facialpoints = facialpoints(rrr,:); clear rrr;
        end
        
        % Remove facial points for now
        headshape.pos(count_facialpoints,:) = [];
        
        % Plot the facial and head points in separate colours
        figure;
        if isempty(count_facialpoints)
            disp('Not plotting any facial points')
        else
            ft_plot_mesh(facialpoints,'vertexcolor','r','vertexsize',10); hold on;
        end
        ft_plot_mesh(headshape.pos,'vertexcolor','k','vertexsize',10); hold on;
        view([90 0]);
        
        % Create mesh out of headshape downsampled to x points specified in the
        % function call
        cfg = [];
        %cfg.numvertices = 1000;
        cfg.method = 'headshape';
        cfg.headshape = headshape.pos;
        mesh = ft_prepare_mesh(cfg, headshape);
        
        %
        [decimated_headshape] = decimate_headshape(headshape, 'gridaverage');
        
        
        % Create figure for quality checking
        figure; subplot(2,2,1);ft_plot_mesh(mesh,'facecolor','k',...
            'facealpha',0.1,'edgealpha',0); hold on;
        ft_plot_mesh(headshape_orig.pos,'vertexcolor','r','vertexsize',2); hold on;
        ft_plot_mesh(decimated_headshape,'vertexcolor','b','vertexsize',10); hold on;
        view(-180,0);
        subplot(2,2,2);ft_plot_mesh(mesh,'facecolor','k',...
            'facealpha',0.1,'edgealpha',0); hold on;
        ft_plot_mesh(headshape_orig.pos,'vertexcolor','r','vertexsize',2); hold on;
        ft_plot_mesh(decimated_headshape,'vertexcolor','b','vertexsize',10); hold on;
        view(0,0);
        subplot(2,2,3);ft_plot_mesh(mesh,'facecolor','k',...
            'facealpha',0.1,'edgealpha',0); hold on;
        ft_plot_mesh(headshape_orig.pos,'vertexcolor','r','vertexsize',2); hold on;
        ft_plot_mesh(decimated_headshape,'vertexcolor','b','vertexsize',10); hold on;
        view(90,0);
        subplot(2,2,4);ft_plot_mesh(mesh,'facecolor','k',...
            'facealpha',0.1,'edgealpha',0); hold on;
        ft_plot_mesh(headshape_orig.pos,'vertexcolor','r','vertexsize',2); hold on;
        ft_plot_mesh(decimated_headshape,'vertexcolor','b','vertexsize',10); hold on;
        view(-90,0);
        
        print('headshape_quality','-dpng');
        
        % Replace headshape.pos with decimated pos
        headshape.pos = decimated_headshape;
        
        % Only include points facial points 2cm below nasion
        %rrr  = find(facialpoints(:,3) > -2);
        
        
        % Add the facial points back in (default) or leave out if user specified
        % 'no' in function call
        if strcmp(include_facial_points,'yes')
            try
                % Add the facial info back in
                % Only include points facial points 2cm below nasion
                headshape.pos = vertcat(headshape.pos,...
                    facialpoints(find(facialpoints(:,3) > -2),:));
            catch
                disp('Cannot add facial info back into headshape');
            end
        else
            headshape.pos = headshape.pos;
            disp('Not adding facial points back into headshape');
        end
        
        %Add in names of the fiducials from the sensor
        headshape.fid.label = {'NASION','LPA','RPA'};
        
        % Plot for quality checking
        
        view_angle = [-180, 0]
        figure;
        
        for angle = 1:length(view_angle)
            
            subplot(1,2,angle)
            ft_plot_headshape(headshape,'vertexcolor','k','vertexsize',12) %plot headshape
            hold on;
            ft_plot_headshape(headshape_orig,'vertexcolor','r','vertexsize',2) %plot headshape
            view(view_angle(angle),10);
        end
        
        print('headshape_quality2','-dpng');
        
        
        % Export filename
        headshape_downsampled = headshape;
        
        function [decimated_headshape] = decimate_headshape(headshape, method)
            
            % Convert to a pointcloud in MATLAB
            headshape_pc = pointCloud(headshape.pos);
            
            switch method
                case 'gridaverage'
                    decimated_headshape = pcdownsample(headshape_pc,'gridAverage',3);
                    
                case 'nonuniform'
                    decimated_headshape = pcdownsample(headshape_pc,...
                        'nonuniformGridSample',20);
            end
            
            decimated_headshape = decimated_headshape.Location;
            
        end
        
        
    end


end

