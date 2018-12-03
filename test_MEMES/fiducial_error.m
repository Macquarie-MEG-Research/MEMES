subject = {'2660','2708','2735','2852','2880','2881'...
    ,'2925','3004','3051','3072','2704','2721','2842','2851',...
    '2850','2877','2883','2890','2891','3010','3011','3065',...
    '3076','3088','3089','3095','3096','3128','3146'};

% 3091 removed = no visual gamma
% 3003 ""
% 2870 removed = too young for database
% 2857 removed = no visual gamma due to bad coreg?

cond = {'MEMES','MEMES_average','MEMES_scaled','MEMES_noface'};

dist_MEMES_real_all = zeros(length(cond),...
    length(subject));

for sub = 1:length(subject)
    
    cd(['/Users/44737483/Documents/alien_data/' subject{sub} '/visual']);
    
    load('headshape_downsampled.mat');
    
    fids_real = headshape_downsampled.fid.pos;
    clear headshape_downsampled
    
    for con = 1:length(cond);
        fprintf('Computing Subject %s Cond %s\n',subject{sub},cond{con});

        cd(['/Users/44737483/Documents/alien_data/' subject{sub} '/' cond{con}]);
        
        load('headshape_downsampled.mat');
    
    fids_MEMES = headshape_downsampled.fid.pos;
    
    pdist(
    
    end
end



