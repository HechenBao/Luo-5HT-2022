function [] = indv_stress(path_data) % function for loading raw data and pre-processing
%% loading data for individual mice
%path_data = pwd;                        % get pathway
baseline = dir(fullfile(path_data, ['*_baseline_*.mat']));     
% find .mat file contains "baseline"
load([path_data,'/',baseline.name]);    % load to workspace
coG_base_raw = coef(1,:);               % get baseline coef for G signal
coTd_base_raw = coef(2,:);              % get baseline coef for Td signal

% find .mat file contains "stress"
stress = dir(fullfile(path_data,['*_stress_*.mat']));            % find .mat file contains "stress"
load([path_data,'/',stress.name]);         % load to workspace
coG_stress_raw = coef(1,:);                % get stress coef for G signal
coTd_stress_raw = coef(2,:);               % get stress coef for Td signal

% find .mat file contains "post-stress"
post = dir(fullfile(path_data,['*_post_*.mat']));            % find .mat file contains "stress"
load([path_data,'/',post.name]);         % load to workspace
coG_post_raw = coef(1,:);                % get stress coef for G signal
coTd_post_raw = coef(2,:);               % get stress coef for Td signal


a = size(path_data,2);                   % measure the string of the pathway
name = path_data((a-4):end);            % use the folder name for the summarzied data
save(fullfile(path_data,[name '_raw' '.mat']),'coG_base_raw','coTd_base_raw',...
    'coG_stress_raw','coTd_stress_raw','coG_post_raw','coTd_post_raw')       % save useable variables
export_path = path_data(1:end-5);       % export path


%% downsample to 1 Hz
coG_base = getround(coG_base_raw,10);               % get the round data length for downsampling    
coG_base = sepblockfun(coG_base,[1,10],'mean');     % downsampled to 1 Hz

coTd_base = getround(coTd_base_raw,10);               
coTd_base = sepblockfun(coTd_base,[1,10],'mean');

coG_stress = getround(coG_stress_raw,10);               % get the round data length for downsampling    
coG_stress = sepblockfun(coG_stress,[1,10],'mean');     % downsampled to 1 Hz

coTd_stress = getround(coTd_stress_raw,10);               
coTd_stress = sepblockfun(coTd_stress,[1,10],'mean');

coG_post = getround(coG_post_raw,10);               % get the round data length for downsampling    
coG_post = sepblockfun(coG_post,[1,10],'mean');     % downsampled to 1 Hz

coTd_post = getround(coTd_post_raw,10);               
coTd_post = sepblockfun(coTd_post,[1,10],'mean');
%% polyfit and calculate baseline ratio
%[p1,p2,coG_fit,coTd_fit,base_ratio]=fitsensor(coG_base,coTd_base,1); % get baseline fitting curve

    saveas(gcf,fullfile(export_path,'Output_Plots',[name '_ployfit.png']))
    exportgraphics(gcf,fullfile(export_path,'Output_Plots',[name '_polyfit.pdf']),'ContentType','vector')
    close 
    
stitch_coG = [coG_base,coG_stress,coG_post];        % stitch baseline and stress, post-stress data
stitch_coTd = [coTd_base,coTd_stress,coTd_post];     % stitch baseline and stress, post-stress data

t = 1:length(stitch_coG);            % time of total session

stitch_coG_p1 = polyval(p1,t);       % total session, coG fitted line
stitch_coTd_p2 = polyval(p2,t);      % total session, coTd fitted line

stitch_coG_fit = stitch_coG-stitch_coG_p1+mean(coG_base);    % subtract fitting lines from raw data,add back mean of raw data for the y shift 
stitch_coTd_fit = stitch_coTd-stitch_coTd_p2+mean(coTd_base); 

stitch_ratio = stitch_coG_fit./stitch_coTd_fit;     % get the new corrected ratio for the total session
diff = stitch_ratio./mean(stitch_ratio(1:length(base_ratio))).*100; % normalize the whole session to the baseline


%% plot everything for double check
time = 1:length(stitch_ratio);
time = time./3600; % convert x axis to hour
    
figure('Position', [100 100 900 600])    
    subplot(4,1,1)
    plot(time,stitch_coG,'k')
    hold on
    plot(time,stitch_coG_p1,'b')
    hold on
    plot(time,stitch_coG_fit,'r')
    legend('raw','polyfit','fitted data','Location','eastoutside')
    legend('boxoff')
    title(['coG' '    y=' num2str(p1(1)) '*x+' num2str(p1(2))])
    xlim([0 3.5])
    box off

    subplot(4,1,2)
    plot(time,stitch_coTd,'k')
    hold on
    plot(time,stitch_coTd_p2,'b')
    hold on
    plot(time,stitch_coTd_fit,'r')
    legend('raw','polyfit','fitted data','Location','eastoutside')
    legend('boxoff')
    title(['coTd' '    y=' num2str(p2(1)) '*x+' num2str(p2(2))])
    xlim([0 3.5])
    box off
    
    subplot(4,1,3)
    plot(time,stitch_coG./stitch_coTd,'k')
    hold on
    plot(time,stitch_ratio,'r')
    legend('raw','corrected','Location','eastoutside')
    legend('boxoff')
    title('coG coTd ratio')
    ylabel('coG/coTd')    
    xlim([0 3.5])
    box off
    
    subplot(4,1,4)
    phase1 = length(coG_base);
    phase2 = phase1+length(coG_stress);
    phase3 = phase2+length(coG_post);
    plot(time(1:phase1),diff(1:phase1),'k',...
    time(phase1:phase2),diff(phase1:phase2), 'r',...
    time(phase2:phase3),diff(phase2:phase3), 'm')
    legend('dF/F0%','Location','eastoutside')
    legend('boxoff')
    title('dF/F0')
    xlabel('Time(hr)')
    ylabel('dF/F0%')
    xlim([0 3.5])
    box off
    
    saveas(gcf,fullfile(export_path,'Output_Plots',[name '_Traces.png']))
    exportgraphics(gcf,fullfile(export_path,'Output_Plots',[name '_Traces.pdf']),'ContentType','vector')
    close 
    
%% get histograms of three phases
diff_base = diff(1:phase1);
diff_stress = diff(phase1:phase2);
diff_post = diff(phase2:phase3);

figure;
    x0=0;
    y0=0;
    width=700;
    height=500;
    set(gcf,'position',[x0,y0,width,height])
    bin_step = 1;
    bins = [min(diff):bin_step:max(diff)];

H_base = histogram(diff_base,'BinEdges',bins,'Normalization','probability',...
    'FaceColor','k','FaceAlpha',0.6,'EdgeColor','none');
hold on
H_stress = histogram(diff_stress,'BinEdges',bins,'Normalization','probability',...
    'FaceColor','r','FaceAlpha',0.5,'EdgeColor','none');
hold on
H_post = histogram(diff_post,'BinEdges',bins,'Normalization','probability',...
    'FaceColor','m','FaceAlpha',0.5,'EdgeColor','none');
    xlabel('dF/F0%')
    ylabel('Probability%')
    title('dF/F0 histogram')
    legend('baseline','stress','post stress')
    legend('boxoff')
    box off

    saveas(gcf,fullfile(export_path, 'Output_Plots',[name '_Histogram.png']))
    exportgraphics(gcf,fullfile(export_path, 'Output_Plots',[name '_Histogram.pdf']),'ContentType','vector')

H_base_val = H_base.Values;
H_stress_val = H_stress.Values;
H_post_val = H_post.Values;
    close 

%% get useful variables for statistics analysis
% binsize = 30.*60;    % 30 mins as bin
% get the segmentation of bin data, and average of each bin
% [output_binstress,output_binstress_mean] = bindata(diff_stress,binsize); 

%output = [mean(diff_base),




%% save data and variables to the summarized file
% save([name '_processed' '.mat'])
% output_dir = ('D:\5HT sensor\GRAB5HT photometry\Output\stress');
% save([output_dir,'stress_male.mat'],'-append')

save(fullfile(path_data,[name '_processed.mat']))

    
clear
clc
end
