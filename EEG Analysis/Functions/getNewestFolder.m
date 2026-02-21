function [newestFolder] = getNewestFolder()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
d = dir;
filtered_d = d(~ismember({d.name}, {'.', '..'}));
[~, index]   = max([filtered_d.datenum]);
% youngestFile = fullfile(filtered_d(index).folder, filtered_d(index).name)
newestFolder = filtered_d(index).name %#ok<NOPRT>
end