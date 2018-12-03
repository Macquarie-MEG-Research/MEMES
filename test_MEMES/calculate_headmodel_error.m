%% Compute Headmodel Error

subject = {'2660','2708','2735','2852','2880','2881'...
    ,'2925','3004','3051','3072','2704','2721','2842','2851',...
    '2850','2877','2883','2890','2891','3010','3011','3065',...
    '3076','3088','3089','3095','3096','3128','3146'};

cond = {'MEMES','MEMES_average','MEMES_scaled','MEMES_noface'};

dist_MEMES_real_all = zeros(length(cond),11000,...
    length(subject));

for sub = 1:length(subject)
    for con = 1:length(cond);
        fprintf('Computing Subject %s\n',subject{sub});
        
        cd(['/Users/44737483/Documents/alien_data/' subject{sub} '/' cond{con}]);
        
        load('sourcemodel3d.mat');
        headmodel_MEMES = sourcemodel3d;
        clear headmodel;
        
        cd(['/Users/44737483/Documents/alien_data/' subject{sub} '/visual']);
        
        load('sourcemodel3d.mat');
        load('headshape_downsampled.mat');
        
        headmodel_real = sourcemodel3d;
        clear headmodel;
        
        dist_MEMES_real = zeros(length(headmodel_MEMES.pos),1);
        
        for pos = 1:length(headmodel_MEMES.pos)
            dist_MEMES_real(pos) = pdist2(headmodel_real.pos(pos,:),...
                headmodel_MEMES.pos(pos,:));
        end
        
        dist_MEMES_real_all(con,:,sub) = dist_MEMES_real;
        
        
        %     figure;
        %     ft_plot_vol(headmodel_real,  'facecolor', 'r', 'edgecolor', 'none');alpha 0.1; camlight;
        %     ft_plot_vol(headmodel_MEMES,  'facecolor', 'b', 'edgecolor', 'none');alpha 0.1; camlight;
        %     title(sprintf('Subject %s . Mean error = %.3f',subject{sub},...
        %         min(dist_MEMES_real))); drawnow;
        %     ft_plot_headshape(headshape_downsampled);
        clear dist_MEMES_real
    end
end


dist_MEMES_real_all_avg = (squeeze(mean(dist_MEMES_real_all,2)));
figure; boxplot(dist_MEMES_real_all_avg');

% Export peak data
MEMES = (dist_MEMES_real_all_avg(1,:)./10)';
MEMES_average = (dist_MEMES_real_all_avg(2,:)./10)';
MEMES_scaled = (dist_MEMES_real_all_avg(3,:)./10)';
MEMES_noface = (dist_MEMES_real_all_avg(4,:)./10)';
group = vertcat(ones(10,1),ones(19,1)*2,ones(14,1)*3);

MEMES_data = (vertcat(MEMES,MEMES_average,MEMES_scaled,...
    MEMES_noface));

MEMES_variation = vertcat(ones(length(MEMES),1),ones(length(MEMES),1)*2,...
    ones(length(MEMES),1)*3,ones(length(MEMES),1)*4);

group_all = vertcat(group,group,group,group);

sourcemodel_data = table(MEMES_data,MEMES_variation,group_all);

cd('/Users/44737483/Documents/scripts_mcq/MEMES/test_MEMES');

writetable(sourcemodel_data,'sourcemodel_data_whole_brain.csv','Delimiter',...
    ',','QuoteStrings',true);

figure; ft_plot_mesh(headmodel_MEMES.pos(headmodel_MEMES.inside,:)...
    ,'vertexcolor','b','vertexsize',30);
ft_plot_mesh(headmodel_real.pos(headmodel_MEMES.inside,:),'vertexcolor',...
    'r','vertexsize',30);
view([176,25]);
print('source_model_whole_brain','-dpng','-r300');


