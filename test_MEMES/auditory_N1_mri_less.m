%% Source analysis

subject = {'2704','2721','2842','2851','2850','2877','2883','2890',...
    '2891','3010','3011','3065','3076','3088','3089','3095','3096'...
    ,'3128','3146','3157'};

% Time*subject array
grandavgA = [];
grandavgB = [];

% Load template_grid (5,8,10mm) - here I use 8
load('/Users/44737483/Documents/fieldtrip-20170501/template/sourcemodel/standard_sourcemodel3d8mm.mat');
template_grid = sourcemodel;
template_grid = ft_convert_units(template_grid,'mm');
clear sourcemodel;

for sub = 1:length(subject)
    
    fprintf('\nProcessing subject %s\n',subject{sub});
    
    %dir_name
    dir_name = ['/Users/44737483/Documents/MMN_data/' subject{sub} '/'];
    
    MEMES_dir = ['/Users/44737483/Documents/MMN_data/'...
        subject{sub} '/MEMES'];
    
    fprintf('\nLoading Data for Subject %s\n',subject{sub});
    
    cd(dir_name);
    load('data_clean_noICA_no_noise.mat');
    load('lay.mat');
    
    % Go to MEMES directory and load the 
    cd(MEMES_dir);
    load('headmodel');
    load('sourcemodel3d');
    load('grad_trans');
    load('mri_realigned_MEMES.mat');
    
    %% LP filter your data to beautify
    cfg =[];
    cfg.bpfilter = 'yes';
    cfg.bpfreq = [1 30];
    cfg.detrend     = 'yes';
    cfg.demean = 'yes';
    cfg.baselinewindow  = [-0.3 -0.1];
    data_clean_noICA_lpfilter = ft_preprocessing(cfg,data_clean_noICA);

    %% Create leadfields and common filter
    cfg = [];
    cfg.grid = sourcemodel3d;
    cfg.grid.unit = 'mm'
    cfg.headmodel = headmodel;
    cfg.grad = grad_trans;
    cfg.normalize = 'yes'; % May not need this
    lf = ft_prepare_leadfield(cfg);
    
    % make a figure of the single subject headmodel, and grid positions
    figure; hold on;
    ft_plot_vol(headmodel,  'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
    ft_plot_mesh(lf.pos(lf.inside,:));
    ft_plot_sens(grad_trans, 'style', 'r*')
    
    % Create
    cfg.covariance       = 'yes';
    %cfg.covariancewindow = 'all';
    cfg.covariancewindow = [-0.15 0.15];
    all_data_cov         = ft_timelockanalysis(cfg,data_clean_noICA_lpfilter);
    
    cfg = [];
    cfg.showlabels = 'yes';
    cfg.fontsize = 6;
    cfg.layout = lay;
    %cfg.xlim = [-0.3 1];
    cfg.xlim = [0.05 0.15];
    figure;
    cfg.baseline = [-0.15 -0.05];
    ft_topoplotER(cfg,all_data_cov);
    print('N1_sensor_level','-dpng','-r100');

    
    % Perform Source Analysis to get common filter across whole trial
    cfg=[];
    cfg.grad = grad_trans;
    cfg.method='lcmv';
    cfg.grid = lf;
    cfg.headmodel=headmodel;
    cfg.lcmv.keepfilter='yes';
    cfg.lcmv.fixedori = 'yes';
    sourceavg=ft_sourceanalysis(cfg, all_data_cov);
    
    %% Baseline Beamforming
    cfg = [];
    cfg.latency = [-0.15 -0.05];
    pre_tone = ft_selectdata(cfg,data_clean_noICA_lpfilter);

    % Time-Lock analysis
    cfg = [];
    cfg.covariance       = 'yes';
    cfg.covariancewindow = 'all';
    cfg.grad = grad_trans;
    erp_pre_tone = ft_timelockanalysis(cfg,pre_tone);
    
    % Source analysis
    cfg = [];
    cfg.method='lcmv';
    cfg.grid=lf;
    cfg.grid.filter=sourceavg.avg.filter; %uses the grid from the whole trial average
    cfg.headmodel=headmodel;
    cfg.grad = grad_trans;
    source_erp_pre_tone = ft_sourceanalysis(cfg, erp_pre_tone);
    
    % Make sure the .pos field is equal to the template_grid.pos
    source_erp_pre_tone.pos = template_grid.pos;

    % Add subject field
    source_erp_pre_tone.subject = subject{sub};
    
    % Now add these to a subject*time array
    grandavgA{sub} = source_erp_pre_tone;
        
    % Select data from the time window
    cfg = [];
    cfg.latency = [0.05 0.15];
    post_tone = ft_selectdata(cfg,data_clean_noICA_lpfilter);
    
    % Time-Lock analysis
    cfg = [];
    cfg.covariance       = 'yes';
    cfg.covariancewindow = 'all';
    erp_post_tone = ft_timelockanalysis(cfg,post_tone);

    % Now do source analysis with the common filter
    cfg = [];
    cfg.method='lcmv';
    cfg.grid=lf;
    cfg.grid.filter=sourceavg.avg.filter; %uses the grid from the whole trial average
    cfg.headmodel=headmodel;
    cfg.grad = grad_trans;
    source_erp_post_tone = ft_sourceanalysis(cfg, erp_post_tone);

    % Make sure the .pos field is equal to the template_grid.pos
    source_erp_post_tone.pos = template_grid.pos;

    % Add subject field
    source_erp_post_tone.subject = subject{sub};
    
    % Now add these to a subject*time array
    grandavgB{sub} = source_erp_post_tone;
    
    sourceR = source_erp_post_tone;
    sourceR.avg.pow = (((source_erp_post_tone.avg.pow - ...
        source_erp_pre_tone.avg.pow)./source_erp_pre_tone.avg.pow).*100);
    
    cfg               = [];
    cfg.funparameter  = 'pow';
    cfg.funcolorlim = 'maxabs';
    cfg.location = 'max';
    %cfg.funcolormap = 'jet';
    ft_sourceplot(cfg,sourceR);
    ft_hastoolbox('brewermap', 1);
    colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
    print('N1_sourcelevel','-dpng','-r100');
    
    close all
    
end

%% STATS

cfg=[];
%cfg.dim         = grandavgA{1}.dim;
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

stat = rmfield(stat,'cfg');

% Show raw source level statistics
cfg               = [];
cfg.maskparameter = 'mask';%turn on to show mask
cfg.method        = 'ortho';
cfg.funparameter  = 'stat';
cfg.location = 'max';
ft_sourceplot(cfg,stat);
drawnow;

    
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
cfg.zlim = [0 6];
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
delete(findall(gcf,'Type','light'))
light ('Position',[90 -10 90])
light ('Position',[-90 -10 90])
material dull;
drawnow;
view([110 0]);
print('aud_M100_yoko_18','-dpng','-r300');
view([-110 0]);
print('aud_M100_yoko_18_2','-dpng','-r300');

%% Export to nifti
cfg = [];
cfg.filetype = 'nifti';
cfg.filename = 'vis_gamma_yoko_18';
cfg.parameter = 'stat';
ft_sourcewrite(cfg,statint);