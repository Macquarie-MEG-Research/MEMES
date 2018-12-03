%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script computes source analysis 
%
%
%
%


% subject = {'2842','2851','2850','2877','2883','2890',...
%     '2891','3010','3011','3065','3076','3088','3089','3095','3096'...
%     ,'3128','3146','3157'};


subject = {'2660','2708','2735','2852','2857','2870','2880','2881'...
    ,'2925','3003','3004','3051','3072','3091'};

subject = {'3076'};

%Load template sourcemodel
load('/Users/44737483/Documents/fieldtrip-20170501/template/sourcemodel/standard_sourcemodel3d8mm.mat');
template_grid = sourcemodel;
template_grid = ft_convert_units(template_grid,'mm');
clear sourcemodel;

for i = 1:length(subject)

%% Load variables computed earlier

%cd(sprintf('/Users/44737483/Documents/scripts_mcq/MEMES_data/%s',subject{i}));

load(sprintf('/Users/44737483/Documents/alien_data/%s/visual/data_clean_noICA.mat',subject{i})); 

cd(['/Users/44737483/Documents/alien_data/' subject{i} '/MEMES_average']);

load('headmodel.mat');
load('sourcemodel3d.mat');
load('grad_trans.mat');
load('shape.mat');
%load(sprintf('/Users/44737483/Documents/scripts_mcq/MEMES_data/%s/mesh.mat',subject{i}));

%load(sprintf('/Users/44737483/Documents/scripts_mcq/MEMES_data/%s/trans_matrix',subject{i}));

%% Downsample at filter data at required frequency
cfg = [];
cfg.resamplefs = 250;
cfg.detrend = 'no';
data_clean_noICA_250 = ft_resampledata(cfg,data_clean_noICA);

%% BP Filter 
cfg = [];
cfg.bpfilter = 'yes'
cfg.bpfreq = [40 70];    %band-pass filter in the required range
data_filtered = ft_preprocessing(cfg,data_clean_noICA_250); clear data_clean_noICA_250

%% Here we redefine trials based on the time-points of interest.
% Make sure the timepoints are of equivalent length

cfg = [];
cfg.toilim = [-1.5 -0.3];
datapre = ft_redefinetrial(cfg, data_filtered);
cfg.toilim = [0.3 1.5];
datapost = ft_redefinetrial(cfg, data_filtered);

%% Create leadfields in subject{i}'s brain warped to MNI space

% 
% % create the subject{i} specific grid, using the template grid that has just been created
% cfg                = [];
% cfg.grid.warpmni   = 'yes';
% cfg.grid.template  = template_grid;
% cfg.grid.nonlinear = 'yes'; % use non-linear normalization
% cfg.mri            = mri_realigned;
% cfg.grid.unit      ='cm';
% cfg.inwardshift = '1.5';
% grid               = ft_prepare_sourcemodel(cfg);
% figure;ft_plot_mesh((grid.pos(grid.inside,:)));
% 
% %grid.pos(:,2) = grid.pos(:,2).*-1;
% 
% %grid.pos = ft_warp_apply(mri_realigned.transform,grid.pos);
% grid.pos = ft_warp_apply(trans_matrix,grid.pos);


%% Create leadfield
cfg = [];
cfg.headmodel = headmodel;
cfg.grid = sourcemodel3d
cfg.grad = grad_trans;
cfg.normalize = 'yes'; % May not need this
lf = ft_prepare_leadfield(cfg);

% make a figure of the single subject{i} headmodel, and grid positions
figure; hold on;
ft_plot_vol(headmodel,  'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.4; camlight;
ft_plot_mesh(lf.pos(lf.inside,:),'vertexsize',5);
ft_plot_sens(grad_trans, 'style', 'r*')

% Here we are keeping all parts of the trial to compute the 
% covariance matrix --> common filter
cfg = [];
cfg.covariance = 'yes';
cfg.covariancewindow = [-1.5 1.5]
avg = ft_timelockanalysis(cfg,data_filtered);

% Time lock analysis for datapre and datapost period
cfg = [];
cfg.covariance='yes';
avgpre = ft_timelockanalysis(cfg,datapre);
avgpst = ft_timelockanalysis(cfg,datapost);

% Create a figure to illustrate the averaged timecourse
figure
plot(avg.time,avg.avg)

%% Do ya beamforming
% Source reconstruction for the whole trial
cfg=[];
cfg.grad = grad_trans;
cfg.method='lcmv';
cfg.grid=lf;
cfg.headmodel=headmodel;
cfg.lcmv.keepfilter='yes';
sourceavg=ft_sourceanalysis(cfg, avg);

% Now do beamforming for the two time points separately using the same spatial
% filter computed from the whole trial
cfg = [];
cfg.grad = grad_trans;
cfg.method='lcmv';
cfg.grid=lf;
cfg.grid.filter=sourceavg.avg.filter; %uses the grid from the whole trial average
cfg.headmodel=headmodel;
%Pre-grating
sourcepreS1 = ft_sourceanalysis(cfg, avgpre);
%Post-grating
sourcepstS1=ft_sourceanalysis(cfg, avgpst);

%If mesh
% standard_mesh = ft_read_headshape({['/Users/44737483/Documents/scripts_mcq/HCP_restingstate/megconnectome-3.0 2/template/Conte69.R.inflated.4k_fs_LR.surf.gii'],...
%     ['/Users/44737483/Documents/scripts_mcq/HCP_restingstate/megconnectome-3.0 2/template/Conte69.L.inflated.4k_fs_LR.surf.gii']});
% 
% standard_mesh = ft_convert_units(standard_mesh,'mm');

% Make sure your field positions match the template grid

sourcepreS1.pos=template_grid.pos; % right(?)
sourcepstS1.pos=template_grid.pos; % right(?)

save sourcepreS1 sourcepreS1
save sourcepstS1 sourcepstS1

% save sourcepreS1 sourcepreS1
% save sourcepstS1 sourcepstS1

%Plot the difference - not necessary but useful for illustration purposes
cfg = [];
cfg.parameter = 'avg.pow';
cfg.operation = '((x1-x2)./x2)*100';
sourceR=ft_math(cfg,sourcepstS1,sourcepreS1);

mri = ft_read_mri('/Users/44737483/Documents/fieldtrip-20170501/template/anatomy/single_subj_T1.nii');

% Interpolate onto SPM brain
cfg              = [];
cfg.voxelcoord   = 'no';
cfg.parameter    = 'pow';
cfg.interpmethod = 'nearest';
sourceI  = ft_sourceinterpolate(cfg, sourceR, mri);

%% Plot
cfg = [];
cfg.funparameter = 'pow';
ft_sourceplot(cfg,sourceI);
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
saveas(gcf,'sourceI.png');

% cfg = [];
% cfg.method         = 'surface';
% cfg.funparameter   = 'avg.pow';
% %cfg.funcolorlim    = [0 20];
% %cfg.downsample     = 6;
% cfg.projmethod     = 'nearest';
% cfg.surfinflated   = 'surface_inflated_both.mat';
% %cfg.surfdownsample = 10
% cfg.projthresh     = 0.2;
% cfg.camlight       = 'no';
% ft_sourceplot(cfg, sourceI);
% ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
% colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
% %view ([70 0 50])
% light ('Position',[70 0 50])
% material dull;
% drawnow;

%close all

end


subject = {'2704','2721','2842','2851','2850','2877','2883','2890',...
    '2891','3010','3011','3076','3088','3089','3095','3096'...
    ,'3128','3146'};

subject = {'2660','2708','2735','2852','2857','2870','2880','2881'...
    ,'2925','3003','3004','3051','3072','3091'};

grandavgA = []; grandavgB = [];

for i = 1:length(subject)

cd(sprintf('/Users/44737483/Documents/alien_data/%s/MEMES',subject{i}));

load('sourcepreS1.mat'); load('sourcepstS1.mat'); 

% sourceR = sourcepreS1;
% sourceR.avg.pow = (((sourcepstS1.avg.pow - sourcepreS1.avg.pow)./...
%     sourcepreS1.avg.pow).*100);
% 
% cfg               = [];
% cfg.method        = 'ortho';
% cfg.funparameter  = 'pow';
% cfg.location = 'max';
% %cfg.maskparameter = 'mask';%turn on to show mask
% ft_sourceplot(cfg,sourceR);
% ft_hastoolbox('brewermap', 1);    
% colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
% title(subject{i}); drawnow;

grandavgA{i} = sourcepreS1; grandavgB{i} = sourcepstS1;

end




cfg=[];
cfg.dim         = grandavgA{1}.dim;
cfg.method      = 'montecarlo';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.parameter   = 'pow';
cfg.correctm    = 'cluster';
cfg.computecritval = 'yes'
cfg.numrandomization = 4000;
cfg.alpha       = 0.4; % Set alpha level
%cfg.clusteralpha = 0.025;
cfg.tail        = 1;    % Two sided testing

% Design Matrix
nsubj=numel(grandavgA);
cfg.design(1,:) = [1:nsubj 1:nsubj];
cfg.design(2,:) = [ones(1,nsubj) ones(1,nsubj)*2];
cfg.uvar        = 1; % row of design matrix that contains unit variable (in this case: subjects)
cfg.ivar        = 2; % row of design matrix that contains independent variable (the conditions)

stat = ft_sourcestatistics(cfg,grandavgB{:}, grandavgA{:});
%save stat stat

% Show raw source level statistics
cfg               = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'stat';
cfg.location = 'max';
%cfg.maskparameter = 'mask';%turn on to show mask
ft_sourceplot(cfg,stat);
ft_hastoolbox('brewermap', 1);    
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap


%% Interpolate onto SPM T1 brain
mri = ft_read_mri('/Users/44737483/Documents/fieldtrip-20170501/template/anatomy/single_subj_T1.nii');

cfg              = [];
cfg.voxelcoord   = 'no';
cfg.parameter    = 'stat';
cfg.interpmethod = 'nearest';
statint  = ft_sourceinterpolate(cfg, stat, mri); %your FT variable corresponding to the subject specific nifti
cfg.parameter    = 'mask';
maskint  = ft_sourceinterpolate(cfg, stat,mri);
statint.mask = maskint.mask;

%% Plot interpolated results
cfg               = [];
%cfg.method        = 'slice';
cfg.funparameter  = 'stat';
cfg.maskparameter = 'mask';
cfg.funcolorlim = 'maxabs';
cfg.location = 'max';
cfg.funcolormap = 'jet';
ft_sourceplot(cfg,statint);


cd('/Users/44737483/Dropbox/RS PhD Documents/Results/MEMES');

%% 
cfg = [];
cfg.method         = 'surface';
cfg.funparameter   = 'stat';
%cfg.funcolorlim    = [0 20];
%cfg.downsample     = 6;
cfg.maskparameter = 'mask';
cfg.projmethod     = 'nearest';
cfg.surfinflated   = 'surface_inflated_both.mat';
%cfg.surfdownsample = 10
%cfg.projthresh     = 0.3;
cfg.camlight       = 'no';
ft_sourceplot(cfg, statint);
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
%view ([70 0 50])
%delete(findall(gcf,'Type','light'))
light ('Position',[0 0 50])
light ('Position',[0 -50 0])
material dull;
drawnow;
view([0 0]);
print('vis_gamma_yoko_13','-dpng','-r300');
view([-60 0]);
print('vis_gamma_yoko_13_2','-dpng','-r300');
view([60 0]);
print('vis_gamma_yoko_13_3','-dpng','-r300');

%% Export to nifti
cfg = [];
cfg.filetype = 'nifti';
cfg.filename = 'vis_gamma_yoko_13';
cfg.parameter = 'stat';
ft_sourcewrite(cfg,statint);


%OptionZ.FrameRate=15;OptionZ.Duration=5.5;OptionZ.Periodic=true;
%CaptureFigVid([0,0; 360,0;], 'Gamma_40_70_with_MEMES',OptionZ)


