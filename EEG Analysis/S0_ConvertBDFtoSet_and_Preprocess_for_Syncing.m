% Preprocessing Pipeline for ERPLAB Analysis - a plugin of EEGLAB

% First, open EEGLAB, then load all your BDF/,set files you wish to preprocess into EEGLAB using 
% the Biosig3.8.1 plugin. Select the first dataset and have that loaded in
% before moving on.  The Current Folder in MATLAB should be the Script
% folder with this file.
% 
% Next, place the locations folder "32BioSemiChannelCoordinates.ced" in a 
% subfolder called "Functions".
% 
% Then, customize which preprocessing you wish to perform using the menu below.

%%%%%% ------- Functions to Perform ------------ %%%%%% Set to 0 or 1
%% Standard Intro Preprocessing Steps (Generally does not change)
tStart = tic;
% The two lines below are for parallel loop
% delete(gcp('nocreate'));            %% Run this first, esp after running parfor to remove previous errors.
% parpool(12);                          %% Sets the number of cores to use e.g. 8-12 for parrallel processing

dateAndTime = char(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss'));
dateAndTimeAbbr = char(datetime('now','TimeZone','local','Format','yyMMdd_HHmmss'));

script_name = mfilename;  % mfilename is the name of the file that is being ran now
folderName = ['0_OutputSets/' script_name '_' dateAndTimeAbbr];     %Folder name to save outputs into

saveUnprocessedDataset  = 1;                             % 1 to save unprocessed dataset to output folder

datasetType = 2;                                         % 1 if datasets are BDF files, 2 if datasets are SET files
addChLocs = 1;
applyEventLists = 1;                                     %1 to apply an event list from a file to the datasets, 0 to skip
convertTriggers = 1;                                     % 1 to convert triggers from larger numbers, 0 for not converting triggers
    saveTriggerConvertedData=0;
eventListFileName = 'AttentionTaskEventListStimONLY.txt';        %File name of the events to be analyzed
remove78 = 1;                                            % Remove Channels EXG7 EXG8

filterHighPass = 1;                                      %1 to perform a high-pass filter on the data, 0 to skip
    filterHighPassFrequency = 1;                             %Frequency to perform high-pass filter on in Hz.  1 Hz for ICA

notchFilter = 0;                                         %1 to remove line noise from the data, 0 to skip
    notchFilterFrequency = 60;                               %Line noise frequencies to remove in Hz (60 Hz in USA)
saveFilteredData = 1;

extractEpochsEarly = 0;                                  %1 to extract epochs early (BEFORE bad channel removal)
epochBinRange = [-4000.0 1200.0];
%% Non-standard post-intro preprocessing steps
removeBadChs = 0;
    saveRemoveChsDataset = 1;

trimNonTaskEEG = 0;                                      % 1 to remove EEG data that was recording when the task was NOT running (include only Task EEG)

reref = 4;                                               %0 to skip re-referencing the datasets
                                                         %1 to re-rereference to the average,
                                                         %2 to re-reference to a common channel,
                                                         %3 to re-reference to the average of mastoid channels
                                                         % 4: find channel labels for
                                                         % mastoid channels
                                                         % and rerefer to
                                                         % the average
                                                         % mastoid channels
    commonChannel = 'ch32';                                  %Channel number of the common channel if 2 is selected
    commonChannelName = 'Cz';                                %Channel label of the common channel if 2 is selected
    saveRRdata = 1;                                          % Save rereferenced data

%% ICA Steps
ComputeICAweights = 1;                                         %1 to perform ICA Analysis on datasets, 0 to skip
    preICAEpochExtract = 0;                                  % Get epochs BEFORE running ICA
    preICAEpochBinRange = epochBinRange;                      % [-4000.0 1200.0];  Universal epoch for both Lu and ERP analysis
    saveICADatasets = 1;                                     %1 to save datasets with ICA weights, 0 to skip
saveICAFigures = 0;                                      %1 to save MATLAB figure and JPEG of plot of components, 0 to skip

ICLabelFlag = 1;
    saveICLabelFlag = 1;

removeICs = 1;                                           % Remove ICs relating to eye and muscle artifacts
    removeEyeArtifacts = 1;                                  % Not implemented yet
    removeMuscleArtifacts = 1;                               % Not implemented yet
    saveSubCompData = 1;                                     % Save EEG dataset with artifact ICs removed as .set
    saveSubCompICfigs = 0;                                   % Not implemented yet

%% Post Subtract Artifact IC Compenent Steps
interpolateBadCh = 0;                                    % 1 to Intperolate Bad channels that were removed

filterLowPass = 0;
    filterLowPassFrequency = 30;

extractEpochs = 0;                                       %1 to extract epochs from the datasets AFTER ICA, 0 to skip 
    postICAEpochBinRange = preICAEpochBinRange;
% epochBinRange = [-500.0 1200.0];                       %Range that for the epoch is extracted at for analysis

saveDatasets = 0;                                        %1 to save the datasets, 0 to skip
%% Artifact Detection Steps
    artifactDetection = 0;                                   %1 to perform artifact detection on the epochs, 0 to skip
        artifactDetectionWindow = 200;                           %Moving window range for artifact detection in ms
        artifactDetectionVoltage = 100;                          %Artifact detection voltage in microVolts
%% ERP Steps
computeERP = 0;                                          % 1 to compute the ERP of the datasets, 0 to skip

plotEachERP = 0;                                         % 1 to plot each ERP and save the figures, 0 to skip
    ERPplotchannels = [13];                                         %Channels to plot in each figure
    saveFigs = 1;                                        % 1 to save the ERP MATLAB figure if plotting is performed, 0 to skip
    saveJPGs = 1;                                        % 1 to save the ERP JPG figure if plotting is performed, 0 to skip

saveERP = 0;                                             % 1 to save the ERP, 0 to skip
%%%%%% ----------------------------------------- %%%%%%

%% Get # of datasets,  Create Output Folder
nDatasets = size(ALLEEG,2);  % number of datasets loaded into EEGLAB
if nDatasets == 0
    error('No datasets loaded.  Please load datasets into EEGLAB and try again.')
elseif (nDatasets < 2)
    folderName = 'OutputDatasets/TestOutputs';  % These will prevent creating a new folder if a batch of files aren't being run
end
%Creates folder if it does not already exist
if ~exist(folderName, 'dir')
   mkdir(folderName);
end

%% Variables for IC Component removal / artifact detection
acces = {};
rejs = {};
pacces = {};
prejs = {};
totalEpochsList ={};
ICSfiles = {};       % files that went thru IC subtraction
ADfiles = {};        % List of Artifact Detection files for excel output
comp = {};
eyeCount = {};
muscleCount = {};

files = '';
pStepsString =''; %string containing all abbreviation to all preprocessing step done on 
%% Variables for Channel Removal
fileRemovalInfoSub = {};
chRemovalFiles = {};
totalChsRemoved = {};
ChannelindexesRemoved = {};
RemovedChannelLabels = {};
KeptChannelIndices = {};
KeptChannelLabels = {};
% Lu parfor
% EEGfolder='Y:\singer\LuZhang\Project6-EEG\Data\ICA_EyeBlink_Removal\';
% FileList=dir([EEGfolder '*.set']);    %%%%EEG data
% parfor iDataset=1:length(FileList)

%% Iterate thru all participants loaded in EEGLAB
lastDataset = nDatasets; % set as nDatasets to end at the last participant loaded in EEGLAB
firstDataset = 1; % set as 1 to start from the first participant loaded in EEGLAB
if firstDataset>nDatasets % in case first dataset is hardcoded to an invalid number
    firstDataset = 1;
end
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',firstDataset,'study',0);  

for iDataset = firstDataset:lastDataset
    %if iDataset ~= 1
        % EEG=pop_loadset(FileStruct(iDataset).filename, Folder,'all');
        % [EEG2, ~, ~] = pop_newset(EEG, EEG, iDatasetCURRENTSET,'retrieve',iDataset,'study',0);
        %  save_parfor
        %end
    preprocessedOutput = '';
    preprocessedAbbrSteps = '';
%% Get filename, Channel Locations, and Event List

    %Extract the filename of the file that was loaded into EEGLAB
    if (datasetType == 1)
        filename = EEG.comments(16:end);
        [filepath, name, ext] = fileparts(filename);
    elseif (datasetType == 2)
        [filepath, name, ext] = fileparts(EEG.filename);
    end
    if contains(name,'Flicker','IgnoreCase',true)
        name = strrep(name,'Flicker','F');  %shortens name: "F"= 1h Flicker session
    end
    subjectNum = name(1:4);

    % preprocessedOutput = '';
    % preprocessedAbbrSteps = '';
    files = sprintf('%s%s\n', files, name);
    if (saveUnprocessedDataset==1)
        if ~exist([folderName '/OriginalDatasets/'], 'dir')  % added this in to check if the folder exists and if it doesn't it makes it 
            mkdir([folderName '/OriginalDatasets/'])
        end
       % EEG = pop_saveset(EEG,'filename', [folderName '/OriginalDatasets/' name '_0_OriginalDataset.set'],'savemode', 'onefile');
        EEG = pop_saveset(EEG,'filename', [folderName name '_0_OriginalDataset.set'],'savemode', 'onefile');
        preprocessedOutput = sprintf('%sOriginal Data Set Saved.\n', preprocessedOutput);
    end

    %% Add channel locations
    if (addChLocs == 1)
        EEG = pop_chanedit(EEG, 'lookup', [pwd '/32BioSemiChannelCoordinates.ced']);  %'/Functions/32BioSemiChannelCoordinates.ced'
        
        % Set EEG channel type for EXG electrodes since it does not load automatically from .ced file for some reason. MKA 10/9/23
        EEG=pop_chanedit(EEG, 'settype',{'33:36','EOG'});  % Set EXG1-EXG4 as EOG channel (eye/optical channels),
        EEG=pop_chanedit(EEG, 'settype',{'37:38','MAST'}); % Set EXG5 & EXG6 as MAST (Mastoid) channels 

        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, iDataset, 'overwrite','on','gui','off', 'setname', name);
        preprocessedOutput = sprintf('%sChannel Locations added from file: %s \n', preprocessedOutput, '32BioSemiChannelCoordinates.ced');
    end
     %% Converts Triggers (from larger numbers to 2,22,etc.)
    if (convertTriggers == 1)
        NEvents = size(EEG.event,2);
        for iEvent = 1:NEvents
            if (EEG.event(iEvent).type == 61442) | (EEG.event(iEvent).type == 49154) | (strcmp(EEG.event(iEvent).type, '61442')) | (strcmp(EEG.event(iEvent).type, '49154'))
                EEG.event(iEvent).type = 2;
                EEG.urevent(iEvent).type = 2;
            elseif (EEG.event(iEvent).type == 61462) | (EEG.event(iEvent).type == 49174) | (strcmp(EEG.event(iEvent).type, '61462')) | (strcmp(EEG.event(iEvent).type, '49174'))
                EEG.event(iEvent).type = 22;
                EEG.urevent(iEvent).type = 22;
            elseif (EEG.event(iEvent).type == 61443) | (EEG.event(iEvent).type == 49155) | (strcmp(EEG.event(iEvent).type, '61443')) | (strcmp(EEG.event(iEvent).type, '49155'))
                EEG.event(iEvent).type = 3;
                EEG.urevent(iEvent).type = 3;
            elseif (EEG.event(iEvent).type == 61444) | (EEG.event(iEvent).type == 49156) | (strcmp(EEG.event(iEvent).type, '61444')) | (strcmp(EEG.event(iEvent).type, '49156'))
                EEG.event(iEvent).type = 4;
                EEG.urevent(iEvent).type = 4;
            elseif (EEG.event(iEvent).type == 61445) | (EEG.event(iEvent).type == 49157) | (strcmp(EEG.event(iEvent).type, '61445')) | (strcmp(EEG.event(iEvent).type, '49157'))
                EEG.event(iEvent).type = 5;
                EEG.urevent(iEvent).type = 5;
            end
        end
        [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        
        if (saveTriggerConvertedData==1)
            EEG = pop_saveset( EEG, 'filename', [folderName '\' name '_1_TrigsConverted.set'],'savemode', 'onefile');  % save to run folder, save as one file, addTrigsConverted to file name -10/8/23
        end
        [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        preprocessedOutput = sprintf('%sTriggers Converted \n', preprocessedOutput);
    end 
    % Apply eventlists
    if (applyEventLists == 1)
        EEG  = pop_editeventlist( EEG , 'AlphanumericCleaning', 'on', 'BoundaryNumeric', { -99}, ...
            'BoundaryString', { 'boundary' }, 'List', eventListFileName, 'SendEL2', 'EEG', ...
            'UpdateEEG', 'code', 'Warning', 'off' );
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, iDataset, 'overwrite','on','gui','off');
        preprocessedOutput = sprintf('%sApplied Events List: %s \n', preprocessedOutput, eventListFileName);
    end
    
    %% Removed Channels EXG7 & EXG8 (these chs are typically unused)
    if (remove78 == 1)
        EEG = pop_select( EEG, 'nochannel',{'EXG7','EXG8'});
        preprocessedOutput = sprintf('%sChannels EXG7 & EXG8 removed.\n', preprocessedOutput);
    end
    %% High Pass Filter, Line Noise Removal, SAVE (1) Filtered Dataset
    % High Pass Filter
    if (filterHighPass == 1)
        EEG  = pop_basicfilter( EEG,  1:EEG.nbchan , 'Boundary', 'boundary', 'Cutoff',  filterHighPassFrequency, ...
            'Design', 'fir', 'Filter', 'highpass', 'Order',  36, 'RemoveDC', 'on' );
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, iDataset, 'overwrite','on','gui','off');
        preprocessedOutput = sprintf('%sHigh Pass Filter: %d Hz\n', preprocessedOutput, filterHighPassFrequency);
        preprocessedAbbrSteps = sprintf('%sHPF%d', preprocessedAbbrSteps, filterHighPassFrequency);
    end

    % Notch filter at 60 Hz (to remove line noise)
    if (notchFilter == 1)
        EEG  = pop_basicfilter( EEG,  1:EEG.nbchan , 'Boundary', 'boundary', 'Cutoff',  notchFilterFrequency, ...
            'Design', 'notch', 'Filter', 'PMnotch', 'Order',  180 );
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, iDataset, 'overwrite','on','gui','off');
        preprocessedOutput = sprintf('%sNotch Filter: %d Hz\n', preprocessedOutput, notchFilterFrequency);
        preprocessedAbbrSteps = sprintf('%sNF%d', preprocessedAbbrSteps, notchFilterFrequency);
    end
    
    if (saveFilteredData == 1)
        EEG = pop_saveset( EEG, 'filename', [folderName '/' name '_2_' preprocessedAbbrSteps '.set'],'savemode', 'onefile');
        preprocessedOutput = sprintf('%sFiltered Data Saved.\n', preprocessedOutput);
        [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    end
    
    %% Extract epochs early (BEFORE ChRemoval & ICA)
    if (extractEpochsEarly == 1)
        EEG = pop_epochbin( EEG , epochBinRange,  'pre');
        EEG = pop_saveset( EEG, 'filename', [folderName '/' name '_3_EarlyEpoched.set'],'savemode', 'onefile');
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, iDataset, 'overwrite','on','gui','off');
        preprocessedOutput = sprintf('%sEpochs Extracted Early. Bin Range: [%d %d]\n', preprocessedOutput, epochBinRange(1), epochBinRange(2));
        % preprocessedOutput = sprintf('%sEpochs extracted: %d\n', preprocessedOutput);
    end
    %% Identify/Remove Bad Channels, SAVE (2) "ChsRemoved"
    if (removeBadChs == 1)
        %dEEG = EEG;
        %dEEG = pop_select(dEEG, 'channel', {'Fp1', 'AF3', 'F7', 'F3', 'FC1', 'FC5', 'T7', 'C3', 'CP1', 'CP5', 'P7', 'P3', 'Pz', 'PO3', 'O1', 'Oz', 'O2', 'PO4', 'P4', 'P8', 'CP6', 'CP2', 'C4', 'T8', 'FC6', 'FC2', 'F4', 'F8', 'AF4', 'Fp2', 'Fz', 'Cz'});
        scalpEEG = pop_select(EEG, 'nochannel', {'EXG1','EXG2','EXG3','EXG4','EXG5','EXG6'});
        mastoidEEG= pop_select(EEG, 'channel', {'EXG5', 'EXG6'});
        %cleanedEEG = pop_clean_rawdata(EEG, 'FlatlineCriterion',5,'ChannelCriterion',0.8,'LineNoiseCriterion',4,'Highpass','off','BurstCriterion','off','WindowCriterion','off','BurstRejection','off','Distance','Euclidian');
        % channels_ignore = {'EXG1','EXG2','EXG3','EXG4','EXG5','EXG6'}; %'EXG1','EXG2','EXG3','EXG4','EXG5','EXG6'
        FlatlineCriterion = 60;
        LineNoiseCriterion = 4;
        ChannelCriterion = 0.6;
        cleanedScalpEEG = pop_clean_rawdata(scalpEEG, ...
            'FlatlineCriterion',FlatlineCriterion,...
            'ChannelCriterion',ChannelCriterion,...
            'LineNoiseCriterion',LineNoiseCriterion,...
            'Highpass','off',...
            'BurstCriterion','off',...
            'WindowCriterion','off',...
            'BurstRejection','off',...
            'Distance','Euclidian'); % 'channels_ignore',channels_ignore
        removedScalpChs = {cleanedScalpEEG.chaninfo.removedchans.labels; cleanedScalpEEG.chaninfo.removedchans.type}'; % get list of all removed channels and their channel type
        removedScalpEEGchIndex = ismember(removedScalpChs(:,2), 'EEG'); % get index of the removed EEG channels only
        removedScalpEEGChannels = removedScalpChs(removedScalpEEGchIndex,1); % get names of removed Scalp EEG channels only
        % removedScalpEEGChannels2 = removedScalpChs(ismember(removedScalpChs(:,2), 'EEG'),1);        
        
        cleanedMastoidEEG = pop_clean_rawdata(mastoidEEG,...
            'FlatlineCriterion',FlatlineCriterion,...
            'ChannelCriterion','off',...
            'LineNoiseCriterion','off',...
            'Highpass','off',...
            'BurstCriterion','off',...
            'WindowCriterion','off',...
            'BurstRejection','off',...
            'Distance','Euclidian');
        removedMastoidChs = {cleanedMastoidEEG.chaninfo.removedchans.labels; cleanedMastoidEEG.chaninfo.removedchans.type}'; % get list of all removed channels and their channel type
        removedMastoidMASTchIndex = ismember(removedMastoidChs(:,2), 'MAST'); % get index of the removed Mastoid channels only
        removedMastoidMASTChannels = removedMastoidChs(removedMastoidMASTchIndex,1);  % get names of removed Mastoid channels only
        
        % combine list of removed Scalp and Mastoid channels
        removedChList = [removedScalpEEGChannels; removedMastoidMASTChannels];

        % remove bad removed channels from orignal EEG dataset
        cleanedEEG = pop_select(EEG, 'nochannel', removedChList);
       
        rmChanIdx = find(~ismember({EEG.chanlocs.labels}, {cleanedEEG.chanlocs.labels})); % get Indexes of removed channels
        % rci = join(string({EEG.chaninfo.removedchans(2:end).urchan}))   
        EEG.chaninfo.removedchans(2:end).labels;   
        keptChanIdx = find(ismember({EEG.chanlocs.labels}, {cleanedEEG.chanlocs.labels}));
        rmChanLabels = {EEG.chanlocs(rmChanIdx').labels}';  %get labels of removed Channels
        if isempty(rmChanIdx)
            cleanedEEG.removedChs = struct('RmChanIdx', num2cell(0), 'rmChanLabels', 'No Channels Removed');
        else
            cleanedEEG.removedChs = struct('RmChanIdx', num2cell(rmChanIdx'), 'rmChanLabels', rmChanLabels);
        end
        %chRemovalFiles
        % fileRemovalInfo(iDataset,1) = {name};
        % 
        % %totalChsRemoved
        % fileRemovalInfo(iDataset,2) = {size(rmChanIdx,2)};
        % 
        % %totalChsKept
        % fileRemovalInfo(iDataset,3) = {size({cleanedEEG.chanlocs.labels},2)};
        % 
        % %ChannelindexesRemoved
        % fileRemovalInfo(iDataset,4) = {rmChanIdx};
        % 
        % %RemovedChannelLabels
        % fileRemovalInfo(iDataset,5) = {join(string(rmChanLabels'))};
        % 
        % %KeptChannelIndices
        % fileRemovalInfo(iDataset,6) = {keptChanIdx};
        % 
        % %KeptChannelLabels
        % fileRemovalInfo(iDataset,7) = {join(string({cleanedEEG.chanlocs(:).labels}))};

        % Changed iDataset to end+1 and end to accomodate preprocessing
        % that does not start at first loaded file
        %chRemovalFiles
        fileRemovalInfoSub(end+1,1) = {name};
        
        %totalChsRemoved
        fileRemovalInfoSub(end,2) = {size(rmChanIdx,2)};
        
        %totalChsKept
        fileRemovalInfoSub(end,3) = {size({cleanedEEG.chanlocs.labels},2)};
        
        %ChannelindexesRemoved
        fileRemovalInfoSub(end,4) = {rmChanIdx};
        
        %RemovedChannelLabels
        fileRemovalInfoSub(end,5) = {join(string(rmChanLabels'))};
        
        %KeptChannelIndices
        fileRemovalInfoSub(end,6) = {keptChanIdx};
        
        %KeptChannelLabels
        fileRemovalInfoSub(end,7) = {join(string({cleanedEEG.chanlocs(:).labels}))};

        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, cleanedEEG, iDataset,'overwrite','off','gui','off');
        

        if (saveRemoveChsDataset == 1)
            EEG = pop_saveset( EEG, 'filename', [folderName '\' name '_4_ChsRemoved.set'],'savemode', 'onefile');
            [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        end
        preprocessedOutput = sprintf('%sBad Channels Removed\n', preprocessedOutput);
        preprocessedOutput = sprintf('%sBad Channel Removal parameters:\n', preprocessedOutput);
        preprocessedOutput = sprintf('%sFlatlineCriterion: %d \n', preprocessedOutput,FlatlineCriterion);
        if exist('ChannelCriterion','var')
            preprocessedOutput = sprintf('%sChannelCriterion: %d \n', preprocessedOutput,ChannelCriterion);
        end
        %preprocessedOutput = sprintf('%sLineNoiseCriterion: %d \n', preprocessedOutput,LineNoiseCriterion);
    end
    EEGChIndices = find(strcmp({EEG.chanlocs.type}, 'EEG')==1); %finds all remaining EEG channel indicies

    %% Trim Data to exclude non-task EEG data, SAVE(3)
    if (trimNonTaskEEG == 1)
        events = [];
        for i = 1:length(EEG.event)
            events(end+1) = EEG.event(i).type;
        end
        breaks = find(events == 9);
        count = 0;
        for i = length(breaks):-1:1
            if (breaks(i) ~= length(EEG.event))
                x = EEG.event(breaks(i)).latency;
                y = EEG.event(breaks(i) + 1).latency;
                range = (y - 4096) - x;
                EEG = eeg_eegrej(EEG, [x x+range]);
                lengths(end+1) = {range/512};
                starts(end+1) = {x/512};
            end
        end
    
        trimHeaders = {{'Start Time'} {'Length'}};
        trimTable = [starts' lengths'];
        trimTable = [trimHeaders; trimTable];
        trimTable = cell2table(trimTable);
        writetable(trimTable, [folderName '\' name '_Trimmed.xls']);
    
        EEG = pop_saveset( EEG, 'filename', [folderName '\' name '_' dateAndTimeAbbr '_5_Trimmed.set']);

        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    end

    %% Re-reference the data, SAVE(4)
    if (reref == 1)
        %Re-Reference the data to the average
        EEG = pop_eegchanoperator( EEG, {  'nch1 = ch1 - ( avgchan( 1:32) ) Label Fp1',  'nch2 = ch2 - ( avgchan( 1:32) ) Label AF3',...
            'nch3 = ch3 - ( avgchan( 1:32) ) Label F7',  'nch4 = ch4 - ( avgchan( 1:32) ) Label F3',  'nch5 = ch5 - ( avgchan( 1:32) ) Label FC1',...
            'nch6 = ch6 - ( avgchan( 1:32) ) Label FC5',  'nch7 = ch7 - ( avgchan( 1:32) ) Label T7',  'nch8 = ch8 - ( avgchan( 1:32) ) Label C3',...
            'nch9 = ch9 - ( avgchan( 1:32) ) Label CP1',  'nch10 = ch10 - ( avgchan( 1:32) ) Label CP5',  'nch11 = ch11 - ( avgchan( 1:32) ) Label P7',...
            'nch12 = ch12 - ( avgchan( 1:32) ) Label P3',  'nch13 = ch13 - ( avgchan( 1:32) ) Label Pz',  'nch14 = ch14 - ( avgchan( 1:32) ) Label PO3',...
            'nch15 = ch15 - ( avgchan( 1:32) ) Label O1',  'nch16 = ch16 - ( avgchan( 1:32) ) Label Oz',  'nch17 = ch17 - ( avgchan( 1:32) ) Label O2',...
            'nch18 = ch18 - ( avgchan( 1:32) ) Label PO4',  'nch19 = ch19 - ( avgchan( 1:32) ) Label P4',  'nch20 = ch20 - ( avgchan( 1:32) ) Label P8',...
            'nch21 = ch21 - ( avgchan( 1:32) ) Label CP6',  'nch22 = ch22 - ( avgchan( 1:32) ) Label CP2',  'nch23 = ch23 - ( avgchan( 1:32) ) Label C4',...
            'nch24 = ch24 - ( avgchan( 1:32) ) Label T8',  'nch25 = ch25 - ( avgchan( 1:32) ) Label FC6',  'nch26 = ch26 - ( avgchan( 1:32) ) Label FC2',...
            'nch27 = ch27 - ( avgchan( 1:32) ) Label F4',  'nch28 = ch28 - ( avgchan( 1:32) ) Label F8',  'nch29 = ch29 - ( avgchan( 1:32) ) Label AF4',...
            'nch30 = ch30 - ( avgchan( 1:32) ) Label Fp2',  'nch31 = ch31 - ( avgchan( 1:32) ) Label Fz',  'nch32 = ch32 - ( avgchan( 1:32) ) Label Cz',...
            'nch33 = ch33 - ( avgchan( 1:32) ) Label EXG1',  'nch34 = ch34 - ( avgchan( 1:32) ) Label EXG2',  'nch35 = ch35 - ( avgchan( 1:32) ) Label EXG3',...
            'nch36 = ch36 - ( avgchan( 1:32) ) Label EXG4',  'nch37 = ch37 - ( avgchan( 1:32) ) Label EXG5',  'nch38 = ch38 - ( avgchan( 1:32) ) Label EXG6',...
            'nch39 = ch39 - ( avgchan( 1:32) ) Label EXG7',  'nch40 = ch40 - ( avgchan( 1:32) ) Label EXG8'} , 'ErrorMsg', 'popup', 'KeepChLoc',...
            'on', 'Saveas', 'off');
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, iDataset, 'overwrite','on','gui','off');
        preprocessedOutput = sprintf('%sRe-Reference: Average\n', preprocessedOutput);
        % add rereference abbr here
    elseif (reref == 2)
        %Re-reference the data to a common channel
        EEG = pop_eegchanoperator( EEG, {  ['nch1 = ch1 - ( ' commonChannel ' ) Label Fp1'],  ['nch2 = ch2 - ( ' commonChannel ' ) Label AF3'],  ...
            ['nch3 = ch3 - ( ' commonChannel ' ) Label F7'],  ['nch4 = ch4 - ( ' commonChannel ' ) Label F3'],  ['nch5 = ch5 - ( ' commonChannel ' ) Label FC1'],  ...
            ['nch6 = ch6 - ( ' commonChannel ' ) Label FC5'],  ['nch7 = ch7 - ( ' commonChannel ' ) Label T7'],  ['nch8 = ch8 - ( ' commonChannel ' ) Label C3'],  ...
            ['nch9 = ch9 - ( ' commonChannel ' ) Label CP1'],  ['nch10 = ch10 - ( ' commonChannel ' ) Label CP5'],  ['nch11 = ch11 - ( ' commonChannel ' ) Label P7'], ...
            ['nch12 = ch12 - ( ' commonChannel ' ) Label P3'],  ['nch13 = ch13 - ( ' commonChannel ' ) Label Pz'],  ['nch14 = ch14 - ( ' commonChannel ' ) Label PO3'],  ...
            ['nch15 = ch15 - ( ' commonChannel ' ) Label O1'],  ['nch16 = ch16 - ( ' commonChannel ' ) Label Oz'],  ['nch17 = ch17 - ( ' commonChannel ' ) Label O2'],  ...
            ['nch18 = ch18 - ( ' commonChannel ' ) Label PO4'],  ['nch19 = ch19 - ( ' commonChannel ' ) Label P4'],  ['nch20 = ch20 - ( ' commonChannel ' ) Label P8'],  ...
            ['nch21 = ch21 - ( ' commonChannel ' ) Label CP6'],  ['nch22 = ch22 - ( ' commonChannel ' ) Label CP2'],  ['nch23 = ch23 - ( ' commonChannel ' ) Label C4'],  ...
            ['nch24 = ch24 - ( ' commonChannel ' ) Label T8'],  ['nch25 = ch25 - ( ' commonChannel ' ) Label FC6'],  ['nch26 = ch26 - ( ' commonChannel ' ) Label FC2'],  ...
            ['nch27 = ch27 - ( ' commonChannel ' ) Label F4'],  ['nch28 = ch28 - ( ' commonChannel ' ) Label F8'],  ['nch29 = ch29 - ( ' commonChannel ' ) Label AF4'],  ...
            ['nch30 = ch30 - ( ' commonChannel ' ) Label Fp2'],  ['nch31 = ch31 - ( ' commonChannel ' ) Label Fz'],  ['nch32 = ch32 - ( ' commonChannel ' ) Label Cz'],  ...
            ['nch33 = ch33 - ( ' commonChannel ' ) Label EXG1'],  ['nch34 = ch34 - ( ' commonChannel ' ) Label EXG2'],  ['nch35 = ch35 - ( ' commonChannel ' ) Label EXG3'],  ...
            ['nch36 = ch36 - ( ' commonChannel ' ) Label EXG4'],  ['nch37 = ch37 - ( ' commonChannel ' ) Label EXG5'],  ['nch38 = ch38 - ( ' commonChannel ' ) Label EXG6'],  ...
            ['nch39 = ch39 - ( ' commonChannel ' ) Label EXG7'],  ['nch40 = ch40 - ( ' commonChannel ' ) Label EXG8']} , 'Saveas', 'off');
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, iDataset, 'overwrite','on','gui','off');
        preprocessedOutput = sprintf('%sRe-Reference: Common Channel (%s)\n', preprocessedOutput, commonChannelName);
        % add rereference abbr here

    elseif (reref == 3)
        %Re-reference the data to the mastoid channels
        EEG = pop_eegchanoperator( EEG, {  'nch1 = ch1 - ( (ch37 + ch38)/2 ) Label Fp1',  'nch2 = ch2 - ( (ch37 + ch38)/2 ) Label AF3',...
            'nch3 = ch3 - ( (ch37 + ch38)/2 ) Label F7',  'nch4 = ch4 - ( (ch37 + ch38)/2 ) Label F3',  'nch5 = ch5 - ( (ch37 + ch38)/2 ) Label FC1',...
            'nch6 = ch6 - ( (ch37 + ch38)/2 ) Label FC5',  'nch7 = ch7 - ( (ch37 + ch38)/2 ) Label T7',  'nch8 = ch8 - ( (ch37 + ch38)/2 ) Label C3',...
            'nch9 = ch9 - ( (ch37 + ch38)/2 ) Label CP1',  'nch10 = ch10 - ( (ch37 + ch38)/2 ) Label CP5',  'nch11 = ch11 - ( (ch37 + ch38)/2 ) Label P7',...
            'nch12 = ch12 - ( (ch37 + ch38)/2 ) Label P3',  'nch13 = ch13 - ( (ch37 + ch38)/2 ) Label Pz',  'nch14 = ch14 - ( (ch37 + ch38)/2 ) Label PO3',...
            'nch15 = ch15 - ( (ch37 + ch38)/2 ) Label O1',  'nch16 = ch16 - ( (ch37 + ch38)/2 ) Label Oz',  'nch17 = ch17 - ( (ch37 + ch38)/2 ) Label O2',...
            'nch18 = ch18 - ( (ch37 + ch38)/2 ) Label PO4',  'nch19 = ch19 - ( (ch37 + ch38)/2 ) Label P4',  'nch20 = ch20 - ( (ch37 + ch38)/2 ) Label P8',...
            'nch21 = ch21 - ( (ch37 + ch38)/2 ) Label CP6',  'nch22 = ch22 - ( (ch37 + ch38)/2 ) Label CP2',  'nch23 = ch23 - ( (ch37 + ch38)/2 ) Label C4',...
            'nch24 = ch24 - ( (ch37 + ch38)/2 ) Label T8',  'nch25 = ch25 - ( (ch37 + ch38)/2 ) Label FC6',  'nch26 = ch26 - ( (ch37 + ch38)/2 ) Label FC2',...
            'nch27 = ch27 - ( (ch37 + ch38)/2 ) Label F4',  'nch28 = ch28 - ( (ch37 + ch38)/2 ) Label F8',  'nch29 = ch29 - ( (ch37 + ch38)/2 ) Label AF4',...
            'nch30 = ch30 - ( (ch37 + ch38)/2 ) Label Fp2',  'nch31 = ch31 - ( (ch37 + ch38)/2 ) Label Fz',...
            'nch32 = ch32 - ( (ch37 + ch38)/2 ) Label Cz',  'nch33 = ch33 - ( (ch37 + ch38)/2 ) Label EXG1',  'nch34 = ch34 - ( (ch37 + ch38)/2 ) Label EXG2',...
            'nch35 = ch35 - ( (ch37 + ch38)/2 ) Label EXG3',  'nch36 = ch36 - ( (ch37 + ch38)/2 ) Label EXG4',  'nch37 = ch37 - ( (ch37 + ch38)/2 ) Label EXG5',...
            'nch38 = ch38 - ( (ch37 + ch38)/2 ) Label EXG6',  'nch39 = ch39 - ( (ch37 + ch38)/2 ) Label EXG7',  'nch40 = ch40 - ( (ch37 + ch38)/2 ) Label EXG8'} ,...
            'ErrorMsg', 'popup', 'KeepChLoc', 'on', 'Saveas', 'off' );
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, iDataset, 'overwrite','on','gui','off');
        preprocessedOutput = sprintf('%sRe-Reference: Average Mastoid\n', preprocessedOutput);
    elseif (reref == 4)
        LMastIndex = find(strcmp({EEG.chanlocs.labels}, 'EXG5')==1); % finds index of left mastoid channel
        RMastIndex = find(strcmp({EEG.chanlocs.labels}, 'EXG6')==1); % finds index of right mastoid channel
        EEG = pop_reref(EEG, [LMastIndex RMastIndex],'keepref','on');  %rereference to average mastoid
        preprocessedOutput = sprintf('%sRe-Reference: Average Mastoid - LeftMastIndex(EXG5)=[%d], RightMastIndex(EXG6)=[%d]\n', preprocessedOutput, LMastIndex, RMastIndex);        
    end
    if reref >0 && saveRRdata == 1
        EEG = pop_saveset( EEG, 'filename', [folderName '\' name '_6_RR.set'],'savemode', 'onefile');
        [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    end
       
    
    %% Perform ICA Analysis on the datasets, SAVE(5) data and figures
    if sum(ismember({EEG.chanlocs.type},'EEG'))>1     % ICA should not run if there's less than 2 'EEG' scalp channels
        if (ComputeICAweights == 1) %length(EEG.chaninfo.removedchans)<40)  
            if (preICAEpochExtract == 1)
                EEG = pop_epochbin( EEG , preICAEpochBinRange,  'pre');
                [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, iDataset, 'overwrite','on','gui','off');
                EEG = pop_saveset( EEG, 'filename', [folderName '\' name '_7_Epoched.set'],'savemode', 'onefile');
                [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
                preprocessedOutput = sprintf('%sPreICA Epoch Bin Range: [%d %d]\n', preprocessedOutput, preICAEpochBinRange(1), preICAEpochBinRange(2));
            end
    
            EEG = pop_runica(EEG, 'icatype', 'runica', 'chanind', EEGChIndices, 'extended',1,'interrupt','on');
            [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, iDataset, 'overwrite','on','gui','off');
        
            %Save the dataset with ICA weights
            if (saveICADatasets == 1)
                EEG = pop_saveset( EEG, 'filename', [folderName '\' name  '_8_ICA.set'],'savemode', 'onefile');
                [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
            end
            preprocessedOutput = sprintf('%sICA Ran\n', preprocessedOutput);
        end
    
        %Save ICA figures
        if (saveICAFigures == 1) %length(EEG.chaninfo.removedchans)<40)
            pop_topoplot(EEG, 0, EEGChIndices , name, [6 6] ,0,'electrodes','on');
            savefig([folderName '\' name '_' dateAndTimeAbbr '_ICA.fig']);
            saveas(gcf,[folderName '\' name '_' dateAndTimeAbbr '_ICA.jpg']);
            close(gcf);
            preprocessedOutput = sprintf('%sICA Figures Saved.\n', preprocessedOutput);
        end
        
        %% ICLabel and Flag, SAVE(6)
        if (ICLabelFlag ==1) %length(EEG.chaninfo.removedchans)<40)
            EEG = pop_iclabel(EEG, 'default');  % make IC labels based on icweights
            [ALLEEG EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, iDataset, 'overwrite','on','gui','off', 'setname', name);
            
            eyeFlagThres = 0.9;
            muscleFlagThres = 0.9;
            EEG = pop_icflag(EEG, [NaN NaN;muscleFlagThres 1;eyeFlagThres 1;NaN NaN;NaN NaN;NaN NaN;NaN NaN]);  % flag IC for removal (eye >90%, muscle >90%)
            [ALLEEG EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, iDataset, 'overwrite','on','gui','off', 'setname', name);
        
            if (saveICLabelFlag ==1)
                EEG = pop_saveset( EEG, 'filename', [folderName '\' name '_9_ICLabelFlag.set'],'savemode', 'onefile');
                [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
            end
            preprocessedOutput = sprintf('%sIC Components Labeled and Flagged\n', preprocessedOutput);
        end
    
        %% Remove Artifact IC Components, SAVE(7)
        if (removeICs == 1) %length(EEG.chaninfo.removedchans)<40)
            % ADfolderName = 'ArtifactDetection';              %Name of the folder you wish to save files into for artifact detection
            %tableName = 'ICArtifactSummary.xls';                     %Name of excel file to store data into for further anaylsis
                                                   %Limitations of file extensions are based on the MATLAB function writetable
            
            comp(end+1) = {sum(EEG.reject.gcompreject)};
            
            eyeCompSubThresh = 0.8;  % Threshold which eye components greater than that value will be removed.  Changed from default 0.9 to 0.8 on 7/17/23 MKA
            eyes = {};
            eyeCount(end+1) = {sum(EEG.etc.ic_classification.ICLabel.classifications(:,3) > eyeCompSubThresh)};
            eyeIndices = find(EEG.etc.ic_classification.ICLabel.classifications(:,3) > eyeCompSubThresh);
            for i = 1:length(eyeIndices)
                eyes = [eyes; {{'Eye'} {EEG.etc.ic_classification.ICLabel.classifications(eyeIndices(i),3)*100}}];
            end
            preprocessedOutput = sprintf('%sEye Component Subtraction Threshold: %d\n', preprocessedOutput, eyeCompSubThresh);
            
            muscleCompSubThresh = 0.8;  % Threshold which eye components greater than that value will be removed
            muscleCount(end+1) = {sum(EEG.etc.ic_classification.ICLabel.classifications(:,2) > muscleCompSubThresh)};
            muscles = {};
            muscleIndices = find(EEG.etc.ic_classification.ICLabel.classifications(:,2) > muscleCompSubThresh);
            % muscleCount(end+1) = length(muscleIndices);
            for i = 1:length(muscleIndices)
                muscles = [muscles; {{'Muscle'} {EEG.etc.ic_classification.ICLabel.classifications(muscleIndices(i),2)*100}}];
            end
            preprocessedOutput = sprintf('%sMuscle Component Subtraction Threshold: %d\n', preprocessedOutput, muscleCompSubThresh);
            
            EEG = pop_subcomp( EEG, [], 0); % subtract flaged IC components from EEG data 
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, iDataset, 'overwrite','on','gui','off', 'setname', name);
            preprocessedOutput = sprintf('%sFlagged IC Components Removed\n', preprocessedOutput);
            if (saveSubCompData == 1)
                EEG = pop_saveset( EEG, 'filename', [folderName '\' name '_10_SubIC.set'],'savemode', 'onefile');
                [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
            end
            ICSfiles(end+1) = {name};
            
            ICLabelTableName = [folderName '\' name '_ComponentDetails_' dateAndTimeAbbr '.xls'];
            ICLabelHeaders = EEG.etc.ic_classification.ICLabel.classes; %["Brain","Muscle","Eye","Heart","Line Noise","Channel Noise","Other"];
            ICLabelPercents = EEG.etc.ic_classification.ICLabel.classifications;
            ICLabelTable = array2table(ICLabelPercents);
            ICLabelTable.Properties.VariableNames = ICLabelHeaders;
            writetable(ICLabelTable, ICLabelTableName);
            %writetable(dT, [folderName '\' name '_ComponentDetails.xls']);
        end
    else  % runs if there is not enough channels to do ICA - creates empty row for IC Summary Table so equal length as other tables
        ICSfiles(end+1) = {name}; %add file name of file that did not have ICA done
        comp(end+1) = {[]}; %add empty cell to list as placeholder
        eyeCount(end+1) = {[]}; %add empty cell to list as placeholder
        muscleCount(end+1) = {[]}; %add empty cell to list as placeholder
    end

    %% Interpolate Bad Channels, SAVE(8)
    if (interpolateBadCh == 1)
    end

    %% Post ICA: Lowpass filter, Extract Epochs, SAVE(9) datasets   
    %Low Pass Filter
    if (filterLowPass == 1)
        EEG  = pop_basicfilter( EEG,  1:EEG.nbchan , 'Boundary', 'boundary', 'Cutoff',  filterLowPassFrequency, ...
            'Design', 'fir', 'Filter', 'lowpass', 'Order',  36 );
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, iDataset, 'overwrite','on','gui','off');
        preprocessedOutput = sprintf('%sHigh Low Filter: %d Hz\n', preprocessedOutput, filterLowPassFrequency);
    end

    % Post ICA Extract epochs
    if (extractEpochs == 1)
        EEG = pop_epochbin( EEG , postICAEpochBinRange,  'pre');
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, iDataset, 'overwrite','on','gui','off');
        % preprocessedOutput = sprintf('%s Post-ICA Epochs Extracted. Bin Range: %d\n', preprocessedOutput);
        preprocessedOutput = sprintf('%sPost-ICA Epochs Extracted. Bin Range: [%d %d]\n', preprocessedOutput, postICAEpochBinRange(1), postICAEpochBinRange(2));
    end

    %Save preprocessed datasets
    if (saveDatasets == 1)
        EEG = pop_saveset( EEG, 'filename', [folderName '\' name '_' dateAndTimeAbbr '_11_PostICA.set'],'savemode', 'onefile');
        [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    end
    
    %% Artifact Detection
    if (artifactDetection == 1)
        % if (extractEpochs||preICAEpochExtract ||extractEpochsEarly  == 0 )
        %     epochBinRange = preICAEpochBinRange;
        %     EEG = pop_epochbin( EEG , epochBinRange,  'pre');
        %     [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, iDataset, 'overwrite','on','gui','off');
        % end
        lowPassFiltHz = 30;
        Windowstep = 100;
        EEG  = pop_artmwppth( EEG , 'Channel',  EEGChIndices, 'Flag', [ 1 2], 'LowPass',  lowPassFiltHz, 'Threshold',  ...
            artifactDetectionVoltage, 'Twindow', epochBinRange, 'Windowsize',  artifactDetectionWindow, 'Windowstep',  Windowstep );
        EEG = pop_saveset( EEG, 'filename', [folderName '\' name '_12_ArtD.set'],'savemode', 'onefile');
            [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
            [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, iDataset, 'overwrite','on','gui','off');
        preprocessedOutput = sprintf('%sArtifact Detection Threshold: %d V\n', preprocessedOutput, artifactDetectionVoltage);
        preprocessedOutput = sprintf('%sArtifact Detection Window: %d ms\n', preprocessedOutput, artifactDetectionWindow);
        preprocessedOutput = sprintf('%sArtifact Detection Windowstep: %d ms\n', preprocessedOutput, Windowstep);

        %Summary of Artifact Detection
        [EEG, tprej, acce, rej, histoflags ] = pop_summary_AR_eeg_detection(EEG,'none');
    
        %Storing the data into arrays
        ADfiles(end+1) = {name};
        acces(end+1) = {acce(1)};
        rejs(end+1) = {rej(1)};
        pacces(end+1) = {(acce(1)/(acce(1)+(rej(1))))*100};
        prejs(end+1) = {(rej(1)/(acce(1)+(rej(1))))*100};
        [~, currentEpochs, ~] = size(EEG.epoch);
        totalEpochsList = [totalEpochsList, currentEpochs];

     
    end

    %% ERP - Computing and Saving
    %Compute average ERPs
    if (computeERP == 1)
        ERP = pop_averager( EEG , 'Criterion', 'good', 'ExcludeBoundary',...
         'on', 'SEM', 'on' );
    end

    %Plot Each ERP
    if (plotEachERP == 1)
        ERP = pop_ploterps( ERP,  1,  ERPplotchannels, 'Axsize', [ 0.05 0.08], 'BinNum', 'on', 'Blc', 'pre', 'Box', [1 1], 'ChLabel', 'on', ...
            'FontSizeChan',  10, 'FontSizeLeg',  12, 'FontSizeTicks',  10, 'LegPos', 'bottom', 'Linespec', {'k-' , 'r-' , 'b-' , 'g-' , 'c-' }, ...
            'LineWidth',  1, 'Maximize', 'on', 'Position', [ 42 29.6429 106.857 31.9286], 'SEM', 'on', 'Style', 'Classic', 'Tag', 'ERP_figure', 'Transparency',  0.7, 'xscale', ...
            [ -500.0 1200.0   -500:300:1200 ], 'YDir', 'normal', 'yscale', [ -10.0 20.0   -10:5:20 ] );
        
        if saveFigs == 1
            savefig([folderName '\' name '_' dateAndTimeAbbr '_.fig']);
        end
        if saveJPGs == 1
           saveas(gcf,[folderName '\' name '_' dateAndTimeAbbr '_.jpg']);
        end
        close(gcf);
    end

    %Save the ERP
    if (saveERP == 1)
        ERP = pop_savemyerp(ERP, 'erpname', name, 'filename', [folderName '\' name '_' dateAndTimeAbbr '_.erp'], 'Warning', 'on');
    end
    ERP = [];  %clear ERP variable?
    
    %% Check is last file, ifnot load next set. Display completion message.
    if iDataset ~= lastDataset
        % Set next dataset as current dataset, if not at the end
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, iDataset,'retrieve',iDataset+1,'study',0);
    end
    disp(['Preprocessing completed for: ' name]);
end

%% Create Channel Removal Summary Table
if (removeBadChs ==1)
    ChRTableName = ['ChannelRemovalSummary_' dateAndTimeAbbr '.xls'];
    ChRheaders = ["File","Total Channels Removed","Total Channels Kept","Removed Ch Indexes","Remove Channel Labels","Kept Channel Indices","Kept Channel Labels"];
    ChRTable = cell2table(fileRemovalInfoSub);
    ChRTable.Properties.VariableNames = ChRheaders;
    writetable(ChRTable, [folderName '\' ChRTableName]);
end
%% Create IC removal summary excel sheet
if (removeICs ==1)
    ICTableName = ['ICCompRemovalSummary_' dateAndTimeAbbr '.xls'];
    ICheaders = ["File", "Total Components Removed", "Eye Components Removed", "Muscle Components Removed"];
    ICcell = [ICSfiles' comp' eyeCount' muscleCount'];
    % ICTable = [ICheaders; ICTable];
    ICTable = cell2table(ICcell);
    ICTable.Properties.VariableNames = ICheaders;
    writetable(ICTable, [folderName '\' ICTableName]);
end

%% Create Excel file of summary statistics of artifact detection for each subject
if (artifactDetection == 1)
    ArtDtableName = ['ArtifactDetectionSummary_' dateAndTimeAbbr '.xls'];
    artdetParam = {{[dateAndTimeAbbr ' Artifact Detection Parameters:']} {['LPF:' num2str(30) 'Hz']} { ['artifactDetectionVoltage:' num2str(artifactDetectionVoltage) 'uV']} {['epochBinRange: ' num2str(epochBinRange) 'ms' ]} {['artifactDetectionWind: ' num2str(artifactDetectionWindow) 'ms']} {['numberOfEpoch']}};
    headers = {{'Files'} {'% Accepted'} {'% Rejected'} {'# Accepted'} {'# Rejected'} {'# Epochs'}};
    ADTable = [ADfiles' pacces' prejs' acces' rejs' totalEpochsList'];
    ADTable = [artdetParam; headers; ADTable];
    ADTable = cell2table(ADTable);
    writetable(ADTable, [folderName '\' ArtDtableName]);
end
%% Create output text
%Generate output message
if isempty(preprocessedOutput)
    preprocessedOutput = 'None';
end
outputMessage = sprintf('Preprocessing performed on %s\nScript used: %s \n\nPreprocessing steps performed:\n%s \n%d out of %d files preprocessed:\n%s', dateAndTime, mfilename, preprocessedOutput, iDataset,    nDatasets, files);
endDateAndTime = char(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss'));
endDuration = toc(tStart);
outputMessage = sprintf('%s\nPreprocessing completed at: %s\nPreprocess duration: %g seconds.', outputMessage, endDateAndTime, endDuration);
% Add toc

%Display the output message in the command window
disp(outputMessage);

%Create a .txt file with the output message
%folderName = 'OutputDatasets/TestOutputs/';
filePath = fullfile(folderName, ['Preprocessing_' dateAndTimeAbbr '_.txt']);

if ~exist(folderName, 'dir')
    mkdir(folderName);
end

fileID = fopen(filePath, 'w');
if fileID == -1
    error('Failed to open or create the file: %s', filePath);
else
    fprintf(fileID, '%s', outputMessage);
    fclose(fileID);
    fprintf('File created and written: %s\n', filePath);
end

%% Create Summary Spreadsheet for all Procedures run

% Something to note: this section is created assuming AD is always on, if
% you don't have AD on, say only have IC or BCR, etc. on. And still want to
% create a summary spreadsheet using this section, you need to look at the
% lines where "combineTable = ???" lies, and those if statements in the
% combinging BCR and IC section

%Adding AD content
if (artifactDetection||removeICs||removeBadChs)
    %Add column of filenames first

    if (artifactDetection == 1)
        AllSumtableName = ['AllSummary_' dateAndTimeAbbr '.xls'];
        %artdetParam = {{[dateAndTimeAbbr ' Artifact Detection Parameters:']} {['LPF:' num2str(30) 'Hz']} { ['artifactDetectionVoltage:' num2str(artifactDetectionVoltage) 'uV']} {['epochBinRange: ' num2str(epochBinRange) 'ms' ]} {['artifactDetectionWind: ' num2str(artifactDetectionWindow) 'ms']} {['numberOfEpoch']}};
        %ADheaders = {{'Files'} {'% Accepted'} {'% Rejected'} {'# Accepted'} {'# Rejected'} {'# Epochs'}};
        ADheaders4Sum = ["Files", "% Accepted", "% Rejected", "# Accepted", "# Rejected", "# Epochs"];
        ADcells = [ADfiles' pacces' prejs' acces' rejs' totalEpochsList'];
        %ADTable = [artdetParam; headers; ADTable];
        ADTable4Sum = cell2table(ADcells);
        ADTable4Sum.Properties.VariableNames = ADheaders4Sum;
        %writetable(ADTable, [folderName '\' ArtDtableName]);
        combineTable = ADTable4Sum;
    end
    
    %Adding IC content
    if (removeICs ==1)
        %ICTableName = ['ICCompRemovalSummary_' dateAndTimeAbbr '.xls'];
        ICheaders4Sum = ["Total Components Removed", "Eye Components Removed", "Muscle Components Removed"];
        ICcell4Sum = [comp' eyeCount' muscleCount'];
        %ICTable = [ICheaders; ICTable];
        ICTable4Sum = cell2table(ICcell4Sum);
        ICTable4Sum.Properties.VariableNames = ICheaders4Sum;
        %combineTable = [ADTable4Sum, ICTable4Sum];
        %writetable(ICTable, [folderName '\' ICTableName]);
    end
    
    %Adding BCR content
    if (removeBadChs ==1)
        %ChRTableName = ['ChannelRemovalSummary_' dateAndTimeAbbr '.xls'];
        ChRheaders = ["File","Total Channels Removed","Total Channels Kept","Removed Ch Indexes","Remove Channel Labels","Kept Channel Indices","Kept Channel Labels"];
        ChRTable = cell2table(fileRemovalInfoSub);
        %EmptyChR = array2table(repmat("", 2, 7));
        
        %ChRTable = [ChRheaders; fileRemovalInfo];
        ChRTable.Properties.VariableNames = ChRheaders;
        %EmptyChR.Properties.VariableNames = ChRTable.Properties.VariableNames;
        %ChRTableWEmpty = [EmptyChR; ChRTable];
        %writetable(ChRTable, [folderName '\' ChRTableName]);
        if (removeICs ==1)
            combineTable = [ADTable4Sum, ICTable4Sum, ChRTable];
        elseif (removeICs == 0)
            combineTable = [ADTable4Sum, ChRTable];
        end
    end
    
    
    writetable(combineTable, [folderName '\' AllSumtableName]);
    
    
    %Adding .txt content
    
    % fileID = fopen('Preprocessing_' dateAndTimeAbbr '_.txt');
    % lineID = fgetl(fileID);
    % while(ischar(lineID))
    AllSumTXTtableName = ['AllSummaryTXT_' dateAndTimeAbbr '.xls'];
    % Specify the paths to your .txt and .xls files
    txtFilePath = filePath;
    
    % Open the .txt file
    fid = fopen(txtFilePath, 'r');
    
    % Initialize a cell array to store the text data
    txtData = {};
    
    % Read the .txt file line by line using fopen
    line = fgetl(fid);
    while ischar(line)
        txtData = [txtData;{line}];
        line = fgetl(fid); % Read the next line
    end
    
    [~, colcombine] = size(combineTable);
    [rowtxt, ~] = size(txtData);
    emptytxt = repmat({''},rowtxt,(colcombine-1));
    txtcombine = [txtData, emptytxt];
    txtTable = cell2table(txtcombine);
    if (artifactDetection == 1)
        combinedheaders = ADheaders4Sum;
        if (removeICs ==1)
            combinedheaders = [ADheaders4Sum, ICheaders4Sum];
            if(removeBadChs ==1)
                combinedheaders = [ADheaders4Sum, ICheaders4Sum, ChRheaders];
            end
        elseif(removeICs ==0)
            if(removeBadChs ==1)
                combinedheaders = [ADheaders4Sum, ChRheaders];
            end
        end
    end
    
    txtTable.Properties.VariableNames = combinedheaders;
    ADICBCRTable2cell = table2cell(combineTable);
    txtTable2cell = table2cell(txtTable);
    combinedTableCA = [ADICBCRTable2cell; txtTable2cell];
    
    combinedData = cell2table(combinedTableCA);
    combinedData.Properties.VariableNames = combinedheaders;
    
    
    writetable(combinedData, [folderName '\' AllSumTXTtableName]);
end

%% Redraw
eeglab redraw % Updates list of datasets in dropdown menu of GUI
disp('Step0: ConvertBDFtoSet For Syncing is complete')