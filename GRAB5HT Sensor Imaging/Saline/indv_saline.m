function [] = indv_saline(path_data)
%% loading data for individual mice
%path_data = pwd;                        % get pathway
baseline = dir(fullfile(path_data, ['*_baseline_*.mat']));     
% find .mat file contains "baseline"
load([path_data,'/',baseline.name]);    % load to workspace
coG_base_raw = coef(1,:);               % get baseline coef for G signal
coTd_base_raw = coef(2,:);              % get baseline coef for Td signal

% make sure they are in correct order
saline = dir(fullfile(path_data,['*_saline_*.mat']));            % find .mat file contains "saline"
load([path_data,'/',saline.name]);         % load to workspace
coG_saline_raw = coef(1,:);                % get stress coef for G signal
coTd_saline_raw = coef(2,:);               % get stress coef for Td signal

a = size(path_data,2);                   % measure the string of the pathway
name = path_data((a-4):end);            % use the folder name for the summarzied data
save(fullfile(path_data,[name '_raw' '.mat']),'coG_base_raw','coTd_base_raw',...
    'coG_saline_raw','coTd_saline_raw')       % save useable variables
export_path = path_data(1:end-5);       % export path

%% downsample to 1 Hz
coG_base = getround(coG_base_raw,10);               % get the round data length for downsampling    
coG_base = sepblockfun(coG_base,[1,10],'mean');     % downsampled to 1 Hz

coTd_base = getround(coTd_base_raw,10);               
coTd_base = sepblockfun(coTd_base,[1,10],'mean');

coG_saline = getround(coG_saline_raw,10);               % get the round data length for downsampling    
coG_saline = sepblockfun(coG_saline,[1,10],'mean');     % downsampled to 1 Hz

coTd_saline = getround(coTd_saline_raw,10);               
coTd_saline = sepblockfun(coTd_saline,[1,10],'mean');

%% polyfit and calculate baseline ratio
%[p1,p2,coG_fit,coTd_fit,base_ratio]=fitsensor(coG_base,coTd_base,1); % get baseline fitting curve

    saveas(gcf,fullfile(export_path,'Output_Plots',[name '_ployfit.png']))
    exportgraphics(gcf,fullfile(export_path,'Output_Plots',[name '_polyfit.pdf']),'ContentType','vector')
    close 
    
stitch_coG = [coG_base,coG_saline];        % stitch baseline and saline data
stitch_coTd = [coTd_base,coTd_saline];     % stitch baseline and saline data

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
    xlim([0 1])
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
    xlim([0 1])
    box off
    
    subplot(4,1,3)
    plot(time,stitch_coG./stitch_coTd,'k')
    hold on
    plot(time,stitch_ratio,'r')
    legend('raw','corrected','Location','eastoutside')
    legend('boxoff')
    title('coG coTd ratio')
    ylabel('coG/coTd')    
    xlim([0 1])
    box off
    
    subplot(4,1,4)
    phase1 = length(coG_base);
    phase2 = phase1+length(coG_saline);       
    plot(time(1:phase1),diff(1:phase1),'k',...
    time(phase1:phase2),diff(phase1:phase2), 'g')
    legend('dF/F0%','Location','eastoutside')
    legend('boxoff')
    title('dF/F0')
    xlabel('Time(hr)')
    ylabel('dF/F0%')
    xlim([0 1])
    box off
    
    saveas(gcf,fullfile(export_path,'Output_Plots',[name '_Traces.png']))
    exportgraphics(gcf,fullfile(export_path,'Output_Plots',[name '_Traces.pdf']),'ContentType','vector')
    close 
    
%% get histograms of two phases
diff_base = diff(1:phase1);
diff_saline = diff(phase1:phase2);

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
H_saline = histogram(diff_saline,'BinEdges',bins,'Normalization','probability',...
    'FaceColor','g','FaceAlpha',0.5,'EdgeColor','none');
    xlabel('dF/F0%')
    ylabel('Probability%')
    title('dF/F0 histogram')
    legend('baseline','saline')
    legend('boxoff')
    box off

    saveas(gcf,fullfile(export_path, 'Output_Plots',[name '_Histogram.png']))
    exportgraphics(gcf,fullfile(export_path, 'Output_Plots',[name '_Histogram.pdf']),'ContentType','vector')

H_base_val = H_base.Values;
H_saline_val = H_saline.Values;
    close 

%% get useful variables for statistics analysis
% binsize = 30.*60;    % 30 mins as bin
% get the segmentation of bin data, and average of each bin
% [output_binsaline,output_binsaline_mean] = bindata(diff_saline,binsize); 

%output = [mean(diff_base),




%% save data and variables to the summarized file
% save([name '_processed' '.mat'])
% output_dir = ('D:\5HT sensor\GRAB5HT photometry\Output\saline');
% save([output_dir,'saline_male.mat'],'-append')

save(fullfile(path_data,[name '_processed.mat']))

    
clear
clc
end
