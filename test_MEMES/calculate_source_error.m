%% Compute Sourcemodel Error
subject = {'2660','2708','2735','2852','2880','2881'...
    ,'2925','3004','3051','3072','2704','2721','2842','2851',...
    '2850','2877','2883','2890','2891','3010','3011','3065',...
    '3076','3088','3089','3095','3096','3128','3146'};

% 3091 removed = no visual gamma
% 3003 ""
% 2870 removed = too young for database
% 2857 removed = no visual gamma due to bad coreg?

cond = {'MEMES','MEMES_average','MEMES_scaled','MEMES_noface'};
%cond = {'MEMES'};

dist_MEMES_real_all = zeros(length(cond),...
    length(subject));

pow_MEMES_real_all = zeros(length(cond),...
    length(subject));

perc_loop_all = zeros(length(cond),...
    length(subject),197);

for sub = 1:length(subject)
    
    
    cd(['/Users/44737483/Documents/alien_data/' subject{sub} '/visual']);
    
    load('sourcepreS1.mat');
    load('sourcepstS1.mat');
    
    cfg = [];
    cfg.parameter = 'avg.pow';
    cfg.operation = '((x1-x2)./x2)*100';
    sourceR_real=ft_math(cfg,sourcepstS1,sourcepreS1);
        
    for con = 1:length(cond);
        fprintf('Computing Subject %s Cond %s\n',subject{sub},cond{con});

        cd(['/Users/44737483/Documents/alien_data/' subject{sub} '/' cond{con}]);
        
        load('sourcepreS1.mat');
        load('sourcepstS1.mat');
        
        cfg = [];
        cfg.parameter = 'avg.pow';
        cfg.operation = '((x1-x2)./x2)*100';
        sourceR_MEMES =ft_math(cfg,sourcepstS1,sourcepreS1);
        
        % Distance between peak
        ind_max_real = find(max(sourceR_real.pow)==sourceR_real.pow);
        ind_max_MEMES = find(max(sourceR_MEMES.pow)==sourceR_MEMES.pow);
        dist = pdist2(sourceR_real.pos(ind_max_real,:),...
            sourceR_real.pos(ind_max_MEMES,:));
        dist_MEMES_real_all(con,sub) = dist;
        
        % Difference in peak power
        pow_MEMES_real_all(con,sub) = max(sourceR_real.pow) ...
            - max(sourceR_MEMES.pow);
        
        % Percentage of overlap at various thresholds
       
       pos_inside = sourceR_real.pos(sourceR_real.inside,:);
       pow_inside_real = (sourceR_real.pow(sourceR_real.inside,:));
       pow_inside_MEMES = (sourceR_MEMES.pow(sourceR_MEMES.inside,:));

       V = sort(sourceR_real.pow(sourceR_real.inside,:), 'descend');
       W = sort(sourceR_MEMES.pow(sourceR_MEMES.inside,:), 'descend'); 
       
       perc_loop = zeros(1,99); count = 1;
       
       for perc = 1:0.5:99
       
        [pos_real] = find(pow_inside_real<(V(round((length(V)./100).*perc))));
        [pos_MEMES] = find(pow_inside_MEMES<(W(round((length(W)./100).*perc))));
        
        common_element = pos_real(ismember(pos_real,pos_MEMES));
        
        
        perc_common = (((length(pos_real)-length(common_element))...
            ./length(pos_real)) .*100);
        
        perc_loop(count) = 100-perc_common;
        count= count+1;
        
       end

       perc_loop_all(con,sub,:) = perc_loop;
    end
end



figure;boxplot((dist_MEMES_real_all./10)');
figure;boxplot((pow_MEMES_real_all)');

% Export peak data
MEMES = (dist_MEMES_real_all(1,:)./10)';
MEMES_average = (dist_MEMES_real_all(2,:)./10)';
MEMES_scaled = (dist_MEMES_real_all(3,:)./10)';
MEMES_noface = (dist_MEMES_real_all(4,:)./10)';
group = vertcat(ones(10,1),ones(19,1)*2,ones(14,1)*3);

MEMES_data = (vertcat(MEMES,MEMES_average,MEMES_scaled,...
    MEMES_noface));

MEMES_variation = vertcat(ones(length(MEMES),1),ones(length(MEMES),1)*2,...
    ones(length(MEMES),1)*3,ones(length(MEMES),1)*4);

group_all = vertcat(group,group,group,group);

sourcemodel_data = table(MEMES_data,MEMES_variation,group_all);

cd('/Users/44737483/Documents/scripts_mcq/MEMES/test_MEMES');

writetable(sourcemodel_data,'sourcemodel_data.csv','Delimiter',...
    ',','QuoteStrings',true);

% 
MEMES = (pow_MEMES_real_all(1,:)./10)';
MEMES_average = (pow_MEMES_real_all(2,:)./10)';
MEMES_scaled = (pow_MEMES_real_all(3,:)./10)';
MEMES_noface = (pow_MEMES_real_all(4,:)./10)';
group = vertcat(ones(10,1),ones(19,1)*2,ones(14,1)*3);

MEMES_data = (vertcat(MEMES,MEMES_average,MEMES_scaled,...
    MEMES_noface));

MEMES_variation = vertcat(ones(length(MEMES),1),ones(length(MEMES),1)*2,...
    ones(length(MEMES),1)*3,ones(length(MEMES),1)*4);

group_all = vertcat(group,group,group,group);

sourcemodel_data = table(MEMES_data,MEMES_variation,group_all);

cd('/Users/44737483/Documents/scripts_mcq/MEMES/test_MEMES');

writetable(sourcemodel_data,'sourcemodel_data_pow.csv','Delimiter',...
    ',','QuoteStrings',true);

%%%%%%%%%%%%%%%%%%%%

load('perc_loop_all_elekta.mat');

perc_loop_all = horzcat(perc_loop_all,perc_loop_all_elekta);

perc_loop_all_mean = squeeze(mean(perc_loop_all,2));

colors_for_lines = {[0.2823529411764706, 0.47058823529411764, 0.8156862745098039];
 [0.9333333333333333, 0.5215686274509804, 0.2901960784313726];
 [0.41568627450980394, 0.8, 0.39215686274509803];
 [0.8392156862745098, 0.37254901960784315, 0.37254901960784315]};

figure; 
for i = 1:4
    plot([1:0.5:99],perc_loop_all_mean(i,:),'LineWidth',4,'Color',...
        colors_for_lines{i}); hold on;
end
set(gca,'FontSize',18);
ylabel({'% Source Points', 'Overlapping'},'FontSize',30);
xlabel('Threshold %','FontSize',30);
lgd = legend({'MEMES','Average','Scaled','No Face'},'FontSize',14,...
    'Location','northeast');
print('overlap_cluster','-dpng','-r300');





figure;
ft_plot_headshape(headshape_downsampled,'vertexsize',20);
view([90,0]);
print('headshape_no_face','-dpng','-r200');




