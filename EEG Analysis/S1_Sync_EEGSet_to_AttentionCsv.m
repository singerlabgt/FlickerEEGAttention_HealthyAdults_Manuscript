%% Syncing EEG files to PsychoPy Attention Task output .csv's
%% Matty Attokaren 3/2/24
% Need EEGlab 2022.1, may need to open eeglab before running as well

function S1_Sync_EEGSet_to_AttentionCsv(dataset_folder)
    tStart = tic;
    SummaryTextOutput = '';
    %% Change to Folder that this script is in.  This folder should also contain
    % Data Subset folder with contains .set files and attention task .csv's
    if nargin < 1
        dataset_folder = '1_All_Data_For_Syncing';
    end
    % Old Options for data folders:
    % cd ..;  #Go up one folder
    % cd 'Data_Subset'; % Go into Data_Subset folder which should include .set files and attention task .csv's
    % cd '1_Test_Data_Subset_ForStep1';
    % cd '1_Old_Data_Subset';

    % cd '1_All_Data_For_Syncing';
    % cd 1_Input_Sets_and_CSVs_For_Step4_PSDs
    % cd '1_Issue_Sets'
    cd(dataset_folder)

    % cd '1_All40LightAccCutSets'
    dataset_folder_message = sprintf('Input dataset & .csv folder: %s \n', dataset_folder);
    disp(dataset_folder_message);
    path_data = pwd; % pwd always returns path of current folder

    %% Create list of files
    list_of_csv=dir('*.csv'); %creating a list of csv file names
    %list2=dir('*.bdf') %creating a list of bdf file names
    list_of_sets=dir('*.set'); %creating a list of set file names
    if isempty(list_of_sets)
        list_of_bdfs = dir('*.bdf'); % Creating a list of bdf file names
    end
    %% Set up output folder
    cd ..;
    synced_outputs_folder = '1_Synced_Outputs';  %Save locations
    if ~exist(synced_outputs_folder, 'dir')
        mkdir(synced_outputs_folder)
    end
    cd(synced_outputs_folder);

    currDate = strrep(datestr(datetime), ':', '_');
    dateAndTimeAbbr = char(datetime('now','TimeZone','local','Format','yyMMdd_HHmmss'));
    scriptName = mfilename;
    synced_outputs_Subfolder = [currDate '_' scriptName];
    mkdir(synced_outputs_Subfolder);
    cd(synced_outputs_Subfolder);
    path_output = pwd;

    %% Looping through each subject
    subjectsList = '';
    filenum = 1;
    filemins=struct('subjectID',{},'filecount',{},'minimum',{},'result',{},'leng',{},'sessionLength',{},'nSession',{},'triggerType',{});

    nCsv = height(list_of_csv); % number of CSVs in list of CSVs
    for iCsv = 1:nCsv  % iCsv cooreposnd to index of CSV in list, which should correspond to an indiviudal subject
        %if iCsv == 40 %each iCsv correspond to a subjectID; these subjects are excluded due to different trigger values or other errors
        %continue;
        %end

        %% This block seems like it should not be in the loop
        clearvars -except iCsv filemins filenum path_data path_output list_of_csv list_of_sets scriptName dataset_folder currDate synced_outputs_Subfolder synced_outputs_folder SummaryTextOutput nCsv tStart dateAndTimeAbbr subjectsList
        cd(path_data)
        % list_of_csv=dir('*.csv'); %creating a list of csv file names
        % list_of_sets=dir('*.set');
        session_time = 0;

        %% Variables for each event
        % dotreact=[]; % clear up dotreact when analyzing a new subject
        % dotnew=[]; % clear up dotnew when analyzing a new subject
        % Find subject ID regardless of filename format, as long as subject
        % number is in the filename preceeded by "S". [1/4/24] MKA
        subjectID_pattern = "S" + digitsPattern(3);  % Pattern looks for capital "S" followed by three digits
        csvFilename = list_of_csv(iCsv).name;
        subjectID = char(extract(csvFilename,subjectID_pattern)) %#ok<NOPTS>
        subjectsList = sprintf('%s%s\n', subjectsList, subjectID);

        % Old version of subjectID (requires ID to be in specific position (10-13)
        % subjectIDold=list_of_csv(iCsv).name(10:13) %subjectID = the id of the current subject

        % if isempty(list_of_sets)
        %     listtemp = dir([subjectID, '*', '.bdf']);
        % else
        iSubject_SetList=dir([subjectID, '*.set']); %listtemp contains names of all set files for a current subject


        if isempty(iSubject_SetList)
            continue;
        end
        % read csv
        csv_filepath = fullfile(path_data,list_of_csv(iCsv).name); %reading csv file for a specific subject
        csv = readtable(csv_filepath); % change the file name here
        reactionKey = csv.reactionKey_rt; % reactionKey = raw reactionKey.rt data
        dotraw = csv.reactionKey_started; % dotraw = raw reactionKey.started data

        %% extracting a consecutive segment of eeg data
        set_count = height(iSubject_SetList); %set_count is the number of set files availabe for each subject
        % count = 1; %initializing count
        for iSet = 1:set_count %looping through each set file
            cd(path_data);
            dotreact=[]; % clear up dotreact when analyzing a new subject
            dotnew=[]; % clear up dotnew when analyzing a new subject
            min_sums=[];
            set_file = fullfile(path_data,iSubject_SetList(iSet).name); % finding the file name of a specific set file
            set_file = char(set_file);
            EEG = pop_loadset(set_file);% open set file using EEGLAB
            if isempty(EEG.event)
                continue
            end
            % netindex = listtemp(iSet).name(12);
            EEGeventsTable=struct2table(EEG.event); %converting EEG.event to a table
            numberOfEEGEvents = height(EEGeventsTable);
            latency=EEGeventsTable.latency;
            if any(ismember(EEGeventsTable.Properties.VariableNames,'edftype'))
                trigger=EEGeventsTable.edftype;
                % if isa(trigger, 'cell')
                %     trigger=cell2mat(trigger);
                % end
            end
            if any(ismember(EEGeventsTable.Properties.VariableNames,'type'))
                trigger=EEGeventsTable.type;
            end

            iTrigger = 1;
            row_index_eegraw = 1;
            column_index_eegraw = 1;
            eegraw = []; %eegraw is a table that contains all latency values with trigger 2. Each column is a continuous sample.
            while iTrigger<numberOfEEGEvents+1
                if trigger(iTrigger,1) == 2 || trigger(iTrigger,1) == 61442 || trigger(iTrigger,1) == 49154 % if trigger is 2, record the corresponding latency value to eegraw
                    eegraw(row_index_eegraw,column_index_eegraw) = latency(iTrigger,1);
                    iTrigger=iTrigger+1;
                    row_index_eegraw=row_index_eegraw+1;
                elseif trigger(iTrigger,1) == 13 || trigger(iTrigger,1) == 61453 || trigger(iTrigger,1) == 49161 %if trigger is 13, start recording latency values in a new column in eegraw
                    column_index_eegraw=column_index_eegraw+1;
                    iTrigger=iTrigger+1;
                    row_index_eegraw=1;
                elseif trigger(iTrigger,1) == 9 || trigger(iTrigger,1) == 61451 % if trigger is 9, start recording latency values in a new column in eegraw
                    column_index_eegraw=column_index_eegraw+1;
                    iTrigger=iTrigger+1;
                    row_index_eegraw=1;
                elseif trigger(iTrigger,1) == 61449
                    column_index_eegraw=column_index_eegraw+1;
                    iTrigger=iTrigger+1;
                    row_index_eegraw=1;
                elseif latency(iTrigger,1) == 170497 %used for subject 57 only
                    column_index_eegraw=column_index_eegraw+1;
                    iTrigger=iTrigger+1;
                    row_index_eegraw=1;;
                else %skip if trigger !=2 or 13 or 9
                    iTrigger=iTrigger+1;
                end
            end
            dot_trial_durations_seconds=dotraw(2:end)-dotraw(1:end-1); % dotdiff is the duration of each trial in seconds
            dot_trial_durations_samples = dot_trial_durations_seconds*512; % dot_trial_duration is the duration of each trial in samples (512 samples per sec) recorded in csv file

            %% synching eeg and dot task data
            %dotnew=[] % clear up dotnew when analyzing a new subject
            nColumns_eegraw = width(eegraw); % nColumns_eegraw=number of continuous recordings in a specific set file
            min_indexes=[]; % the index where the minimum sum happens (the sync point!)
            trials_per_segment=[];
            for iColumns_eegraw = 1:nColumns_eegraw
                etemp=eegraw(:,iColumns_eegraw); %extracting a continuous recording in a specific set file
                etemp = etemp(all(etemp,2),:); %deleting rows of 0s and the end
                eeg_trial_durations=etemp(2:end,1)-etemp(1:end-1,1); %eeg_trial_duration is the duration of each trial in set file
                session_time = session_time + sum(eeg_trial_durations);

                num_of_dot_trials = height(dot_trial_durations_samples);
                num_of_eeg_trials = height(eeg_trial_durations);
                % nloops_subtraction= num_of_dot_trials - num_of_eeg_trials + 1; % nloops_subtraction=number of loops for finding differences between dot_trial_duration and eeg_trial_duration
                nloops_subtraction= abs(num_of_dot_trials - num_of_eeg_trials) + 1; % nloops_subtraction=number of loops for finding differences between dot_trial_duration and eeg_trial_duration

                diff_bw_EEG_and_dot_durations = zeros(length(eeg_trial_durations),1); %initializing - an array that contains the difference in duration between eeg trials and dot trials
                eegdottrial_durations_diffs_sums = zeros(nloops_subtraction,1); %initializing

                for iLoops_subtraction = 1:nloops_subtraction
                    if num_of_dot_trials > num_of_eeg_trials
                        diff_bw_EEG_and_dot_durations = dot_trial_durations_samples(iLoops_subtraction:iLoops_subtraction+length(eeg_trial_durations)-1) - eeg_trial_durations; %calculate difference between dot_trial_duration and eeg_trial_duration
                        eegdottrial_durations_diffs_sums(iLoops_subtraction,1)=sum(abs(diff_bw_EEG_and_dot_durations),'omitnan'); %caculate sum of the differences, save the sum to list_of_sum
                    else % If num og eeg_trials is greater than num_of_dot_trials
                        diff_bw_EEG_and_dot_durations = eeg_trial_durations(iLoops_subtraction:iLoops_subtraction+length(dot_trial_durations_samples)-1) - dot_trial_durations_samples; %calculate difference between dot_trial_duration and eeg_trial_duration
                        eegdottrial_durations_diffs_sums(iLoops_subtraction,1)=sum(abs(diff_bw_EEG_and_dot_durations),'omitnan'); %caculate sum of the differences, save the sum to list_of_sum
                    end
                end
                [min_sum_of_diffs,min_ind] = min(eegdottrial_durations_diffs_sums,[],'omitnan'); %finding the minimum of the eegdottial duration difference sums
                min_sums(iColumns_eegraw,1)=min_sum_of_diffs;
                min_indexes(iColumns_eegraw,1)=min_ind;
                trials_per_segment(iColumns_eegraw,1)=length(etemp);
                dotnew=[dotnew;etemp]; % dotnew = synched behaviral data
            end
            %[xa,ya] = alignsignals(dot_trial_duration,eeg_trial_duration,5)

            %dotreact=[] % clear up dotreact when analyzing a new subject

            %calculating for reaction time
            for iColumns_eegraw=1:nColumns_eegraw
                etemp=eegraw(:,iColumns_eegraw);
                etemp = etemp(all(etemp,2),:);
                reactemp=[]; %reactemp = trial start time + reaction time in samples (512 samples per sec)
                reactemp= etemp + reactionKey(min_indexes(iColumns_eegraw):min_indexes(iColumns_eegraw)+trials_per_segment(iColumns_eegraw)-1)*512;
                dotreact=[dotreact;reactemp];
            end
            dotsynch=[]; % 1st column contains synched trial start time; 2nd column contains synched reaction time
            dotsynch(:,2)=dotreact();
            dotsynch(:,1)=dotnew();

            % update the 3rd column on dotsynch
            for iHeight_dotsynch = 1:height(dotsynch)
                if ~isnan(dotsynch(iHeight_dotsynch,2))
                    dotsynch(iHeight_dotsynch,3) = 1;
                else
                    dotsynch(iHeight_dotsynch,3) = 0;
                end
            end
            tooearly = csv.tooEarly_keys;
            indexct = 1; % indexct = index in the tooearlyindex which recording of the index of the premature hit
            for iLength_tooearly=1:length(tooearly)
                if tooearly(iLength_tooearly,1) ~= "None" & tooearly(iLength_tooearly,1) ~= ""
                    tooearlyindex = iLength_tooearly; %records the index of the premature hit
                    prematureindex = 0;
                    iindex = 1;
                    while iindex <= length(min_indexes)
                        if min_indexes(1,1) > tooearlyindex %if eeg recording starts after a premature hit
                            break
                        elseif iindex == length(min_indexes) & tooearlyindex > min_indexes(iindex) & tooearlyindex < min_indexes(iindex) + trials_per_segment(iindex) %if we reach the end of index
                            prematureindex = prematureindex + tooearlyindex - min_indexes(iindex) + 1;
                            dotsynch(prematureindex,3) = -1; %-1 = premature, 0 = miss, 1 = hit
                            break
                        elseif tooearlyindex == min_indexes(1) %
                            prematureindex = 1;
                            dotsynch(prematureindex,3) = -1;
                            break
                        elseif iindex == length(min_indexes) & tooearlyindex > min_indexes(iindex) & tooearlyindex > min_indexes(iindex) + trials_per_segment(iindex) %if eeg recordings ends beore a premature hit
                            break
                        elseif tooearlyindex > min_indexes(iindex)
                            prematureindex = prematureindex + trials_per_segment(iindex,1);
                            iindex = iindex + 1;
                        else
                            prematureindex = prematureindex - trials_per_segment(iindex-1,1) + tooearlyindex - min_indexes(iindex-1) + 1;
                            dotsynch(prematureindex,3) = -1; %-1 = premature, 0 = miss, 1 = hit
                            break
                        end
                    end
                end
            end
            result = 'synched';

            indices = find(dotsynch(:,1)==0);
            dotsynch(indices,:) = [];

            %omiting trials that cannot be sync'd

            for iMinimum = 1:length(min_sums)
                if min_sums(iMinimum) > 200
                    etemp=eegraw(:,iMinimum); %extracting a continuous recording in a specific set file
                    etemp = etemp(all(etemp,2),:); %deleting rows of 0s and the end
                    eeg_trial_durations=etemp(2:end,1)-etemp(1:end-1,1); %eeg_trial_duration is the duration of each trial in set file
                    length_subtraction_omit=trials_per_segment(iMinimum)-1; %height(dot_trial_duration)-height(eeg_trial_duration)+1; % nloops_subtraction=number of loops for finding differences between dot_trial_duration and eeg_trial_duration
                    diff_bw_EEG_and_dot_durations = []; %initializing
                    eegdottrial_durations_diffs_sums = []; %initializing
                    while length_subtraction_omit > 0
                        diff_bw_EEG_and_dot_durations = dot_trial_durations_samples(min_indexes(iMinimum):min_indexes(iMinimum)+length_subtraction_omit-1)-eeg_trial_durations(1:length_subtraction_omit); %calculate difference between dot_trial_durations and eeg_trial_durations
                        if sum(abs(diff_bw_EEG_and_dot_durations),'omitnan') > 200 %caculate sum of the differences, save the sum to list_of_sum
                            length_subtraction_omit = length_subtraction_omit - 1;
                        else
                            index_Omit = length_subtraction_omit + 2 + sum(trials_per_segment(1:iMinimum-1));
                            index_Omit_end = index_Omit + trials_per_segment(iMinimum) - 2 - length_subtraction_omit;
                            break;
                        end
                    end
                    if exist('index_Omit','var')==1
                        dotsynch(index_Omit:index_Omit_end,3) = -2; % -2 indicates a excluding trial
                        min_sums(iMinimum,1)=sum(abs(diff_bw_EEG_and_dot_durations),'omitnan');
                    end
                end
            end

            for iMinimum = 1:length(min_sums)
                if min_sums(iMinimum) > 200
                    result = 'not synched';
                    break;
                end
            end

            hits = find(dotsynch(:,3)==1);
            misses = find(dotsynch(:,3)==0);
            premature = find(dotsynch(:,3)==-1);

            %determining trigger type
            if sum(trigger==2)
                triggerType = 2;
            end
            if sum(trigger==49154)
                triggerType = 49154;
            end
            if sum(trigger==61442)
                triggerType = 61442;
            end

            %save dotsynch to a location in box
            cd(path_output);
            % save([subjectID,'_',netindex,'_synched.mat',],"dotsynch","min_sums","trials_per_segment","hits","misses","premature","result") % saving dotsynch to a specific location
            save([subjectID,'_synched.mat',],"dotsynch","min_sums","trials_per_segment","hits","misses","premature","result") % saving dotsynch to a specific location

            filemins(filenum).subjectID = subjectID;
            filemins(filenum).minimum = min_sums;
            filemins(filenum).filecount = iSet;
            filemins(filenum).result = result;
            filemins(filenum).leng = height(dotsynch);
            filemins(filenum).sessionLength = session_time/512;
            filemins(filenum).nSession = length(trials_per_segment);
            filemins(filenum).triggerType = triggerType;
            filenum = filenum + 1;
        end
    end



    %% debugging
    etemp=eegraw(:,1);
    etemp = etemp(all(etemp,2),:);
    eeg_trial_durations=etemp(2:end,1)-etemp(1:end-1,1);

    %%
    hits = find(dotsynch(:,3)==1);
    misses = find(dotsynch(:,3)==0);
    premature = find(dotsynch(:,3)==-1);
    %%
    EEGeventsTable=struct2table(EEG.event);
    trigger = EEGeventsTable.type;
    triggerdiff = [];
    for deb=1:length(trigger)
        if trigger(deb) == 61442
            continue
        elseif trigger(deb) == 61462
            continue
        else
            triggerdiff = [triggerdiff;trigger(deb)];
        end
    end
    %%
    triggerind = find(trigger~=49152 & trigger~=49154);
    %%
    triggertest = trigger;
    triggertest(triggertest(:,1)==49155)=0;

    %%
    for iMinimum = 1:length(min_sums)
        if min_sums(iMinimum) > 200
            etemp=eegraw(:,iMinimum);
            etemp = etemp(all(etemp,2),:);
            eeg_trial_durations=etemp(2:end,1)-etemp(1:end-1,1);
            for iDiff = 1:trials_per_segment(iMinimum)-1 %index(iMinimum):index(iMinimum)+leng(iMinimum)-1
                diffOmit = dot_trial_durations_samples(min_indexes(iMinimum)+iDiff-1) - eeg_trial_durations(iDiff);
                indexOmit = min_indexes(iMinimum)+iDiff-1;
                if diffOmit > 1000
                    break;
                end
            end
        end
    end
    %%
    for iMinimum = 1:length(min_sums)
        if min_sums(iMinimum) > 200
            etemp=eegraw(:,iMinimum); %extracting a continuous recording in a specific set file
            etemp = etemp(all(etemp,2),:); %deleting rows of 0s and the end
            eeg_trial_durations=etemp(2:end,1)-etemp(1:end-1,1); %eeg_trial_duration is the duration of each trial in set file
            length_subtraction_omit=trials_per_segment(iMinimum)-1; %height(dot_trial_duration)-height(eeg_trial_duration)+1; % nloops_subtraction=number of loops for finding differences between dot_trial_duration and eeg_trial_duration
            diff_bw_EEG_and_dot_durations = []; %initializing
            eegdottrial_durations_diffs_sums = []; %initializing
            while length_subtraction_omit > 0
                diff_bw_EEG_and_dot_durations = dot_trial_durations_samples(min_indexes(iMinimum):min_indexes(iMinimum)+length_subtraction_omit-1)-eeg_trial_durations(1:length_subtraction_omit); %calculate difference between dot_trial_duration and eeg_trial_duration
                if sum(abs(diff_bw_EEG_and_dot_durations),'omitnan') > 200 %calculate sum of the differences, save the sum to list_of_sum
                    length_subtraction_omit = length_subtraction_omit - 1;
                else
                    index_Omit = length_subtraction_omit + 2 + sum(trials_per_segment(1:iMinimum-1));
                    index_Omit_end = index_Omit + trials_per_segment(iMinimum) - 2 - length_subtraction_omit;
                    break;
                end
            end
            if exist('index_Omit','var')==1
                dotsynch(index_Omit:index_Omit_end,3) = -2; % -2 indicates a excluding trial
                min_sums(iMinimum,1)=sum(abs(diff_bw_EEG_and_dot_durations),'omitnan');
            end
        end
    end
    %%
    indices = find(dotsynch(:,1)==0);
    dotsynch(indices,:) = [];

    cd ..;
    cd ..;

    %% Save relevant variables for a summary of what this script did
    save([path_output '\SummaryFile'], 'list_of_csv', 'list_of_sets', 'scriptName', 'filenum', 'dataset_folder','path_data', 'path_output', 'currDate', 'synced_outputs_Subfolder', 'synced_outputs_folder', 'nCsv', 'set_count', 'iCsv')

    %% Create output text
    %Generate output message

    if isempty(SummaryTextOutput)
        SummaryTextOutput = 'None';
    end
    SummaryTextOutput = sprintf('Step1 started on %s\nScript used: %s \n', currDate, scriptName);
    SummaryTextOutput =  sprintf('%s \nInput dataset & .csv folder: %s \n', SummaryTextOutput,  dataset_folder);
    SummaryTextOutput =  sprintf('%s \n%d out of %d files synced:\n%s', SummaryTextOutput, iCsv, nCsv, subjectsList);
    endDateAndTime = char(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss'));
    endDuration = toc(tStart);
    SummaryTextOutput = sprintf('%s\nStep1Syncing completed at: %s\nSyncing duration: %g seconds.', SummaryTextOutput, endDateAndTime, endDuration);
    % Add toc

    %Display the output message in the command window
    disp(SummaryTextOutput);

    %Create a .txt file with the output message
    filePath = fullfile(path_output, [scriptName, dateAndTimeAbbr '_Summary.txt']);

    fileID = fopen(filePath, 'w');
    if fileID == -1
        error('Failed to open or create the file: %s', filePath);
    else
        fprintf(fileID, '%s', SummaryTextOutput);
        fclose(fileID);
        fprintf('File created and written: %s\n', filePath);
    end
end