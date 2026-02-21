function [EEG] = mergingEEGSets(datasetsToMerge,mergedSetFilename,mergedFilepathName)
%Append multiple EEG files into one file 1/5/24
%   Detailed explanation goes here

% Edit this line for how many files you are merging ([1 2 ...])
% datasetsToMerge = [1 2];
EEG = pop_mergeset( ALLEEG, datasetsToMerge, 0);

% Edit the name you want the file to save as: this will save over other files if you do not change the name, set to save only a .set file right now but if you wanted a .fdt and .set remove the 'save mode' 'one file' part
% mergedSetFilename = 'S083_F2';
% mergedFilepathName = 'C:\Users\smettupalli8\Box\Project_FlickerHealthyYoungAdults\EEG Data Analysis\Preprocessing\Input Datasets\Raw_SetFiles_AllEEG';
EEG = pop_saveset( EEG, 'filename',mergedSetFilename,'filepath', mergedFilepathName,'savemode', 'onefile');
end