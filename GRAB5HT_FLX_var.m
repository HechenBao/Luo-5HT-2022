function [] = GRAB5HT_FLX_var(path_folder) % batch process function
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
    mean_diff_FLX = mean(diff_FLX);
    %diff_FLX_mm = movmean(diff_FLX,60);
    S = stepinfo(diff_FLX);
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
    
    output_FLX(i,:) = [mean_diff_base,mean_diff_FLX,RiseTime,PeakTime,Peak,SettlingMin,...
        Threshold2,diff_base_th2,diff_base_th2Peak,...
        Threshold3,diff_base_th3,diff_base_th3Peak];
end

header = {'baseline','FLX','Rise Time','Peak Time','Peak','SettingMin',...
    'Threshold2','base_th2','base_th2Peak',...
    'Threshold3','base_th3','base_th3Peak'};
output_FLX = [header; num2cell(output_FLX)];
save(fullfile(path_folder,'output_FLX.mat'),'output_FLX')

end
