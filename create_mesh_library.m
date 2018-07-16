% Scripts to create a mesh library from the HCP-MEG database

path_to_MRI_library = '/Users/44737483/Documents/scripts_mcq/MRI_database/';
cd(path_to_MRI_library);

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

initial_mri_realign = cell(1,length(subject));
scalp_library = cell(1,length(subject));
mesh_library = cell(1,length(subject));

% You also need a headshape preferably decimated and with facial information
% preserved (variable = headshape_downsampled)

% 90deg rotation matrix
rmatx = [6.12323399573677e-17,1,0,0;-1,6.12323399573677e-17,...
    0,0;0,0,1,0;0,0,0,1];

for i = 1:length(subject)
    
    fprintf('Adding subject %d of %d to the MRI database \n',i,length(subject));
    
    mri_file = [path_to_MRI_library subject{i} '/MEG/anatomy/T1w_acpc_dc_restore.nii.gz'];
    
    % Load in MRI
    mri_orig            = ft_read_mri(mri_file); % in mm, read in mri from DICOM
    mri_orig            = ft_convert_units(mri_orig,'mm');
    mri_orig.coordsys   = 'neuromag';
    
    % Rotate 90 degrees
    mri_realigned = ft_transform_geometry(rmatx,mri_orig);
    
%     cfg = [];
%     cfg.fiducial.nas = vox_fids(1,:);
%     cfg.fiducial.lpa = vox_fids(2,:);
%     cfg.fiducial.rpa = vox_fids(3,:);
%     cfg.coordsys = 'bti';
%     mri_realigned = ft_volumerealign(cfg, mri_orig);
%     
    % Give rough estimate of fiducial points using ICP on one example subject
    cfg                         = [];
    cfg.method                  = 'headshape';
    %cfg.viewmode                = 'ortho';
    cfg.headshape               = headshape_downsampled;
    cfg.headshape.interactive   = 'no';
    cfg.headshape.icp           = 'yes';
    cfg.coordsys                = 'bti';
    [mri_realigned]             = ft_volumerealign(cfg, mri_realigned);
    %
    initial_mri_realign{i}    = mri_realigned.cfg.transform_icp;
    
    % % check that the MRI is consistent after realignment
    ft_determine_coordsys(mri_realigned, 'interactive', 'no');
    hold on; % add the subsequent objects to the figure
    drawnow; % workaround to prevent some MATLAB versions (2012b and 2014b) from crashing
    ft_plot_headshape(headshape_downsampled);
    drawnow;
    
    cfg = [];
    cfg.output    = 'scalp';
    cfg.scalpsmooth = 5;
    cfg.scalpthreshold = 0.08;
    scalp  = ft_volumesegment(cfg, mri_realigned);
    
    scalp_library{i} = scalp;
    
    %% Create mesh out of scalp surface
    cfg = [];
    cfg.method = 'isosurface';
    cfg.numvertices = 10000;
    mesh = ft_prepare_mesh(cfg,scalp);
    mesh = ft_convert_units(mesh,'mm');
    %
    
    cd('/Users/44737483/Documents/scripts_mcq/MEMES_data/quality_check');
    
    figure;ft_plot_mesh(mesh,'facecolor',[238,206,179]./255,'EdgeColor','none','facealpha',0.4); 
    hold on; 
    camlight; hold on; view([-180,-10]);
    ft_plot_headshape(headshape_downsampled);
    print(sprintf('mesh_with_headshape_%s',subject{i}),'-dpng');
    close all;
    
    mesh_library{i} = mesh;    
    
end

cd(path_to_MRI_library);

save('initial_mri_realign_july','initial_mri_realign');
save('mesh_library_july','mesh_library');

load('initial_mri_realign_july');
load('mesh_library_july');