function [] = GRAB5HT_saline_bp(path_folder) % batch process function
% get subfolders for batch process
% get a list of all files and folders in this folder.
files = dir(path_folder);

% run scripts on each subfolders
files = files(~ismember({files.name},{'.','..'}));
dirFlags = [files.isdir];   % Get a logical vector that tells which is a directory.
subFolders = files(dirFlags);   % Extract only those that are directories.

for i = 1:numel(subFolders)-1
    fullpath = fullfile(subFolders(i).folder, subFolders(i).name);
    indv_saline(fullpath)        
end
