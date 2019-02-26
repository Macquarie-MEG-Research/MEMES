``` matlab
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
% - bad_coil            = list of bad coils (up to length of 2). Enter as:
%                         {'LPAred','RPAyel','PFblue','LPFwh','RPFblack'}
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
% sens_coreg_method     = method used to realign MEG sensors based on 5
%                       marker coils. Use 'rot3dfit' or 'icp'. For some
%                       reason the usual rot3dfit method seems to fail
%                       sometimes. Try using 'icp' in this case...

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
```
