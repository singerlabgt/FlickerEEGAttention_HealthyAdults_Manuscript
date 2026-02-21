function [newestFilePath] = getNewestFile()
%UNTITLED2 Returns the name of the newest .mat file in Current MATLAB Folder
%   Matty 2/8/24
d = dir('*.mat'); % get all .mat files in folder
filtered_d = d(~ismember({d.name}, {'.', '..'}));  % remove '..mat' and '...mat' files from list
[~, index]   = max([filtered_d.datenum]); % get index of newest file
newestFilePath = fullfile(filtered_d(index).folder, filtered_d(index).name) %#ok<NOPRT>
% newestFolder = filtered_d(index).name %#ok<NOPRT>
end