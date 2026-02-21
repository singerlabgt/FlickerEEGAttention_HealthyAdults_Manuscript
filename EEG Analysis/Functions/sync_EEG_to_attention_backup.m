%% Syncing EEG files to PsychoPy Attention Task output .csv's
%% Matty Attokaren 1/2/24 Based off of Cecilia Xie's 2022 Singer Lab project
% Need EEGlab 2022
%% Change to Folder that this script is in.  This folder should also contain
% Data Subset folder with contains .set files and attention task .csv's
% cd ..;  #Go up one folder
% cd 'Data_Subset'; % Go into Data_Subset folder which should include .set files and attention task .csv's
cd 'Test_Data_Subset';
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
cd 'Output';
currDate = strrep(datestr(datetime), ':', '_');
mkdir(currDate);
cd(currDate);
path_output = pwd;
%% Looping through each subject
filenum = 1;
filemins=struct('subjectID',{},'filecount',{},'minimum',{},'result',{},'leng',{},'sessionLength',{},'nSession',{},'triggerType',{});

for iCsv = 1:height(list_of_csv)  % iCsv cooreposnd to index of CSV in list, which should correspond to an indiviudal subject
    %if iCsv == 40 %each iCsv correspond to a subjectID; these subjects are excluded due to different trigger values or other errors
        %continue;
    %end
    
    %% This block seems like it should not be in the loop
    % clearvars -except iCsv filemins filenum path_data path_output
    cd(path_data)
    % list_of_csv=dir('*.csv'); %creating a list of csv file names
    % list_of_sets=dir('*.set');
    session_time = 0;

    %% Variables for each event
    dotreact=[]; % clear up dotreact when analyzing a new subject
    dotnew=[]; % clear up dotnew when analyzing a new subject
    % Find subject ID regardless of filename format, as long as subject
    % number is in the filename preceeded by "S". [1/4/24] MKA
    subjectID_pattern = "S" + digitsPattern(3);  % Pattern looks for capital "S" followed by three digits
    csvFilename = list_of_csv(iCsv).name;
    subjectID = char(extract(csvFilename,subjectID_pattern))
    
    % Old version of subjectID (requires ID to be in specific position (10-13) 
    % subjectIDold=list_of_csv(iCsv).name(10:13) %subjectID = the id of the current subject
    
    % if isempty(list_of_sets)
    %     listtemp = dir([subjectID, '*', '.bdf']);
    % else
   listtemp=dir([subjectID, '*.set']); %listtemp contains names of all set files for a current subject


    if isempty(listtemp)
        continue;
    end
    % read csv
    csv_filepath = fullfile(path_data,list_of_csv(iCsv).name); %reading csv file for a specific subject
    csv = readtable(csv_filepath); % change the file name here
    reactionKey = csv.reactionKey_rt; % reactionKey = raw reactionKey.rt data
    dotraw = csv.reactionKey_started; % dotraw = raw reactionKey.started data

% extracting a consecutive segment of eeg data
    set_count = height(listtemp); %set_count is the number of set files availabe for each subject
    % count = 1; %initializing count
    for iSet = 1:set_count %looping through each set file 
        cd(path_data);
        dotreact=[]; % clear up dotreact when analyzing a new subject
        dotnew=[]; % clear up dotnew when analyzing a new subject
        minimum=[];
        set_file = fullfile(path_data,listtemp(iSet).name); % finding the file name of a specific set file
        set_file = char(set_file); 
        EEG = pop_loadset(set_file);% open set file using EEGLAB
        if isempty(EEG.event)
            continue
        end    
        netindex = listtemp(iSet).name(12);
        table=struct2table(EEG.event); %converting EEG.event to a table
        max = height(table);
        latency=table.latency;
        if any(ismember(table.Properties.VariableNames,'edftype')) 
            trigger=table.edftype;
        else
            trigger=table.type;
        end
      
        iTrigger = 1;
        row_index_eegraw = 1;
        column_index_eegraw = 1;
        eegraw = []; %eegraw is a table that contains all latency values with trigger 2. Each column is a continuous sample.
        while iTrigger<max+1
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
        dotdiff=dotraw(2:end)-dotraw(1:end-1);
        dot_trial_duration = dotdiff*512; % dot_trial_duration is the duration of each trial in samples (512 samples per sec) recorded in csv file 
  
    % synching eeg and dot task data
        %dotnew=[] % clear up dotnew when analyzing a new subject
        nColumns_eegraw = width(eegraw); % nColumns_eegraw=number of continuous recordings in a specific set file
        index=[];
        leng=[];
        for iColumns_eegraw = 1:nColumns_eegraw
            etemp=eegraw(:,iColumns_eegraw); %extracting a continuous recording in a specific set file
            etemp = etemp(all(etemp,2),:); %deleting rows of 0s and the end
            eeg_trial_duration=etemp(2:end,1)-etemp(1:end-1,1); %eeg_trial_duration is the duration of each trial in set file
            session_time = session_time + sum(eeg_trial_duration);
            nloops_subtraction=height(dot_trial_duration)-height(eeg_trial_duration)+1; % nloops_subtraction=number of loops for finding differences between dot_trial_duration and eeg_trial_duration
            diff = zeros(length(eeg_trial_duration),1); %initializing
            list_of_sum = zeros(nloops_subtraction,1); %initializing
            for iLoops_subtraction = 1:nloops_subtraction
                diff = dot_trial_duration(iLoops_subtraction:iLoops_subtraction+length(eeg_trial_duration)-1)-eeg_trial_duration; %calculate difference between dot_trial_duration and eeg_trial_duration
                list_of_sum(iLoops_subtraction,1)=sum(abs(diff),'omitnan'); %caculate sum of the differences, save the sum to list_of_sum
            end
            [m,ind] = min(list_of_sum,[],'omitnan'); %finding the minimum of list_of_sum
            minimum(iColumns_eegraw,1)=m;
            index(iColumns_eegraw,1)=ind;
            leng(iColumns_eegraw,1)=length(etemp);
            dotnew=[dotnew;etemp]; % dotnew = synched behaviral data
        end
        %[xa,ya] = alignsignals(dot_trial_duration,eeg_trial_duration,5)
    
    %dotreact=[] % clear up dotreact when analyzing a new subject
    
    %calculating for reaction time
        for iColumns_eegraw=1:nColumns_eegraw
            etemp=eegraw(:,iColumns_eegraw);
            etemp = etemp(all(etemp,2),:);
            reactemp=[]; %reactemp = trial start time + reaction time in samples (512 samples per sec)
            reactemp=etemp+reactionKey(index(iColumns_eegraw):index(iColumns_eegraw)+leng(iColumns_eegraw)-1)*512;
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
                while iindex <= length(index)
                    if index(1,1) > tooearlyindex %if eeg recording starts after a premature hit
                        break
                    elseif iindex == length(index) & tooearlyindex > index(iindex) & tooearlyindex < index(iindex) + leng(iindex) %if we reach the end of index
                        prematureindex = prematureindex + tooearlyindex - index(iindex) + 1;
                        dotsynch(prematureindex,3) = -1; %-1 = premature, 0 = miss, 1 = hit
                        break
                    elseif tooearlyindex == index(1) %
                        prematureindex = 1;
                        dotsynch(prematureindex,3) = -1;
                        break
                    elseif iindex == length(index) & tooearlyindex > index(iindex) & tooearlyindex > index(iindex) + leng(iindex) %if eeg recordings ends beore a premature hit
                        break
                    elseif tooearlyindex > index(iindex)
                        prematureindex = prematureindex + leng(iindex,1);
                        iindex = iindex + 1;
                    else 
                        prematureindex = prematureindex - leng(iindex-1,1) + tooearlyindex - index(iindex-1) + 1;
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

        for iMinimum = 1:length(minimum)
            if minimum(iMinimum) > 200
                    etemp=eegraw(:,iMinimum); %extracting a continuous recording in a specific set file
                    etemp = etemp(all(etemp,2),:); %deleting rows of 0s and the end
                    eeg_trial_duration=etemp(2:end,1)-etemp(1:end-1,1); %eeg_trial_duration is the duration of each trial in set file
                    length_subtraction_omit=leng(iMinimum)-1; %height(dot_trial_duration)-height(eeg_trial_duration)+1; % nloops_subtraction=number of loops for finding differences between dot_trial_duration and eeg_trial_duration
                    diff = []; %initializing
                    list_of_sum = []; %initializing
                    while length_subtraction_omit > 0
                        diff = dot_trial_duration(index(iMinimum):index(iMinimum)+length_subtraction_omit-1)-eeg_trial_duration(1:length_subtraction_omit); %calculate difference between dot_trial_duration and eeg_trial_duration
                        if sum(abs(diff),'omitnan') > 200; %caculate sum of the differences, save the sum to list_of_sum
                            length_subtraction_omit = length_subtraction_omit - 1;
                        else
                            index_Omit = length_subtraction_omit + 2 + sum(leng(1:iMinimum-1));
                            index_Omit_end = index_Omit + leng(iMinimum) - 2 - length_subtraction_omit;
                            break;
                        end
                    end
                    dotsynch(index_Omit:index_Omit_end,3) = -2; % -2 indicates a excluding trial
                    minimum(iMinimum,1)=sum(abs(diff),'omitnan');
            end
        end

        for iMinimum = 1:length(minimum)
            if minimum(iMinimum) > 200
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
        save([subjectID,'_',netindex,'_synched.mat',],"dotsynch","minimum","leng","hits","misses","premature","result") % saving dotsynch to a specific location
        filemins(filenum).subjectID = subjectID;
        filemins(filenum).minimum = minimum;
        filemins(filenum).filecount = iSet;
        filemins(filenum).result = result;
        filemins(filenum).leng = height(dotsynch);
        filemins(filenum).sessionLength = session_time/512;
        filemins(filenum).nSession = length(leng);
        filemins(filenum).triggerType = triggerType;
        filenum = filenum + 1;
    end
end

  

%% debugging
etemp=eegraw(:,1);
etemp = etemp(all(etemp,2),:);
eeg_trial_duration=etemp(2:end,1)-etemp(1:end-1,1);

%%
hits = find(dotsynch(:,3)==1);
misses = find(dotsynch(:,3)==0);
premature = find(dotsynch(:,3)==-1);
%%
table=struct2table(EEG.event);
trigger = table.type;
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
for iMinimum = 1:length(minimum)
    if minimum(iMinimum) > 200
        etemp=eegraw(:,iMinimum);
        etemp = etemp(all(etemp,2),:);
        eeg_trial_duration=etemp(2:end,1)-etemp(1:end-1,1);
        for iDiff = 1:leng(iMinimum)-1 %index(iMinimum):index(iMinimum)+leng(iMinimum)-1
            diffOmit = dot_trial_duration(index(iMinimum)+iDiff-1) - eeg_trial_duration(iDiff);
            indexOmit = index(iMinimum)+iDiff-1;
            if diffOmit > 1000
                break;
            end
        end
    end
end
%%
for iMinimum = 1:length(minimum)
    if minimum(iMinimum) > 200
            etemp=eegraw(:,iMinimum); %extracting a continuous recording in a specific set file
            etemp = etemp(all(etemp,2),:); %deleting rows of 0s and the end
            eeg_trial_duration=etemp(2:end,1)-etemp(1:end-1,1); %eeg_trial_duration is the duration of each trial in set file
            length_subtraction_omit=leng(iMinimum)-1; %height(dot_trial_duration)-height(eeg_trial_duration)+1; % nloops_subtraction=number of loops for finding differences between dot_trial_duration and eeg_trial_duration
            diff = []; %initializing
            list_of_sum = []; %initializing
            while length_subtraction_omit > 0
                diff = dot_trial_duration(index(iMinimum):index(iMinimum)+length_subtraction_omit-1)-eeg_trial_duration(1:length_subtraction_omit); %calculate difference between dot_trial_duration and eeg_trial_duration
                if sum(abs(diff),'omitnan') > 200 %calculate sum of the differences, save the sum to list_of_sum
                    length_subtraction_omit = length_subtraction_omit - 1;
                else
                    index_Omit = length_subtraction_omit + 2 + sum(leng(1:iMinimum-1));
                    index_Omit_end = index_Omit + leng(iMinimum) - 2 - length_subtraction_omit;
                    break;
                end
            end
            dotsynch(index_Omit:index_Omit_end,3) = -2; % -2 indicates a excluding trial
            minimum(iMinimum,1)=sum(abs(diff),'omitnan');
    end
end
%%
    indices = find(dotsynch(:,1)==0);
    dotsynch(indices,:) = [];