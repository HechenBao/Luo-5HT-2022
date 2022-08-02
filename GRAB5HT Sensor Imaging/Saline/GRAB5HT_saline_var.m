function [] = GRAB5HT_saline_var(path_folder) % batch process function
% get subfolders for batch process
% get a list of all files and folders in this folder.
files = dir(path_folder);

% run scripts on each subfolders
files = files(~ismember({files.name},{'.','..'}));
dirFlags = [files.isdir];   % Get a logical vector that tells which is a directory.
subFolders = files(dirFlags);   % Extract only those that are directories.

for i = 1:numel(subFolders)-1
    fullpath = fullfile(subFolders(i).folder, subFolders(i).name);
    processed = dir(fullfile(fullpath, ['*_processed*.mat']));     
    % find .mat file contains "processed"
    load([fullpath,'/',processed.name]);    % load to workspace
    mean_diff_base = mean(diff_base);
    mean_diff_saline = mean(diff_saline);
    %diff_saline_mm = movmean(diff_saline,60);
    S = stepinfo(diff_saline);
    RiseTime = S.RiseTime./60;
    PeakTime = S.PeakTime./60;
    Peak = S.Peak;
    SettlingMin = S.SettlingMin;
    
    Threshold2 = 2.*std(diff_base)+mean(diff_base);
    Threshold3 = 3.*std(diff_base)+mean(diff_base);
    diff_base_th2 = sum(diff_base>Threshold2)./length(diff_base).*100;
    diff_base_th2Peak = sum(diff_base(diff_base>Threshold2))./(sum(diff_base>Threshold2));
    diff_base_th3 = sum(diff_base>Threshold3)./length(diff_base).*100;
    diff_base_th3Peak = sum(diff_base(diff_base>Threshold3))./(sum(diff_base>Threshold3));
       
    output_saline(i,:) = [mean_diff_base,mean_diff_saline,RiseTime,PeakTime,Peak,SettlingMin,...
        Threshold2,diff_base_th2,diff_base_th2Peak,...
        Threshold3,diff_base_th3,diff_base_th3Peak];
end

header = {'baseline','saline','Rise Time','Peak Time','Peak','SettingMin',...
    'Threshold2','base_th2','base_th2Peak',...
    'Threshold3','base_th3','base_th3Peak'};
output_saline = [header; num2cell(output_saline)];
save(fullfile(path_folder,'output_saline.mat'),'output_saline')

end
