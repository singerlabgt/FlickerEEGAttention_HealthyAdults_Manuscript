%%  Adapated from Lu's Code by Matty - Checks how well behavior .csv syncs to EEG.set
% Must run from this folder: SensoryFx_EEG_Analysis
% For example: Box\Project_FlickerHealthyYoungAdults\EEG Data Analysis\EEGDataAnalysis2024\MattyRepo_EEGDataAnalysis\SensoryFx_EEG_Analysis
% Add getNewestFolder() to path
function S2_LoadEEG_CheckSync(dataset_folder)
    S2Start = tic;
    eeglab

    %% Set up folders
    if nargin < 1
        dataset_folder = '1_All_Data_For_Syncing';
    end

    % EEGfolder='C:\Users\lzhang481\Box\Project_FlickerHealthyYoungAdults\Data\EEG\Preprocessing\All Subjects Preprocessing\Preprocessed Files 11_21_2022\'
    % EEGfolderName = '1_Old_Data_Subset'; % '1_Test_Data_Subset_ForStep1';
    % EEGfolderName = '1_Test_Data_Subset_ForStep1';
    % EEGfolderName = '1_All40LightAccCutSets'
    % EEGfolderName = '1_Issue_Sets';
    % EEGfolderName = '1_All_Data_For_Syncing';
    % EEGfolderName = '1_Input_Sets_and_CSVs_For_Step4_PSDs';
    cd(dataset_folder);
    EEGfolderPath= pwd;
    EEGSet_FileList=dir('*.set');    %%%%EEG data

    % Synfolder='C:\Users\lzhang481\Box\Project_FlickerHealthyYoungAdults\Data\EEG\Data Analysis\Synching\Output\06-Dec-2022 16_02_07\';
    % Synfolder='Y:\singer\LuZhang\Project6-EEG\Data\Allsubj\';
    cd ..;
    Synfolder = '1_Synced_Outputs';
    % SynfolderName = '2_Synced_Outputs\18-Jan-2024 15_45_32';
    cd(Synfolder);
    newestFolder = getNewestFolder();
    cd(newestFolder);
    SyncFolderPath = pwd;
    SyncList=dir('*synched*'); %%%%Trial Info Data;

    % SaveResult='Y:\singer\LuZhang\Project6-EEG\Results\Step0-PreparingData\';
    cd ..;
    cd ..;
    SaveFolderName = '2_CheckSync_Outputs';
    if ~exist(SaveFolderName, 'dir')
        mkdir(SaveFolderName)
    end
    cd(SaveFolderName);

    currDate = strrep(datestr(datetime), ':', '_');
    scriptName = mfilename;
    synced_outputs_Subfolder = [currDate '_' scriptName '\'];
    mkdir(synced_outputs_Subfolder);
    cd(synced_outputs_Subfolder);
    save_output_path = pwd;

    cd ..;
    cd ..;
    % %% Change to Folder that this script is in.  This folder should also contain
    % % Data Subset folder with contains .set files and attention task .csv's
    % % cd ..;  #Go up one folder
    % % cd 'Data_Subset'; % Go into Data_Subset folder which should include .set files and attention task .csv's
    % cd '1_Test_Data_Subset_ForStep1';
    % path_data = pwd; % pwd always returns path of current folder
    %
    % %% Create list of files
    % list_of_csv=dir('*.csv'); %creating a list of csv file names
    % %list2=dir('*.bdf') %creating a list of bdf file names
    % list_of_sets=dir('*.set'); %creating a list of set file names
    % if isempty(list_of_sets)
    %     list_of_bdfs = dir('*.bdf'); % Creating a list of bdf file names
    % end
    % %% Set up output folder
    % cd ..;
    % synced_outputs_folder = '2_Synced_Outputs';
    % if ~exist(synced_outputs_folder, 'dir')
    %    mkdir(synced_outputs_folder)
    % end
    % cd(synced_outputs_folder);
    %
    % currDate = strrep(datestr(datetime), ':', '_');
    % mkdir(currDate);
    % cd(currDate);
    % path_output = pwd;
    %%
    %% %%Load trial information from Synfolder
    FileStruct={};
    nFiles = length(EEGSet_FileList);
    for iFile=1:nFiles
        EEGSet_FileList(iFile).Subj=EEGSet_FileList(iFile).name(1:4);  % Gets Subject number from file name
        % FileList(iFile).Session=FileList(iFile).name(12);
        iFilename = EEGSet_FileList(iFile).name;
        expression = '(?<=[Ff][a-zA-Z]*)\d_';  % finds a numeric which follows and 'F' or 'f' and is followed by and '_' (underscsore)
        session_number_index = regexp(iFilename,expression);
        EEGSet_FileList(iFile).Session = iFilename(session_number_index);
        % tempF=dir([SynfolderName '*' FileList(iFile).Subj '_' FileList(iFile).Session '_synched.mat']);
        tempF=dir([SyncFolderPath '\' EEGSet_FileList(iFile).Subj '_synched.mat']);
        if isempty(tempF)
            continue
        else
            SynTemp=load([tempF.folder,'\',tempF.name]);
            TempN1=fieldnames(EEGSet_FileList);
            TempN2=fieldnames(SynTemp);
            TempN3=[TempN1;TempN2];
            if iscell(FileStruct)
                FileStruct = FieldName2Struct(TempN3);
            end

            FileStruct(iFile)= MapFields1to2(SynTemp,EEGSet_FileList(iFile));
        end
    end

    Invalid=[];
    %% %%Load trial information from Synfolder

    for iFile=1:length(FileStruct)
        if isempty(FileStruct(iFile).name)
            Invalid=[Invalid;iFile];
        end
    end

    FileStruct(Invalid)=[];

    ValidList=[];
    FlickerList=[];
    clear TrialNmat EEGList

    %% Events
    for iFile=1:length(FileStruct)
        % for iFile=1:1

        tic
        clear EEG OUTEEG TrialType;
        %     iFile=33

        EEG=pop_loadset(FileStruct(iFile).name,FileStruct(iFile).folder,'all');
        if isempty(EEG)
            continue;
        end


        % eventCell=struct2cell(EEG.event);
        % eventID=table2array(cell2table(squeeze(eventCell(2,1,:)))); % type
        EEG_event_table = struct2table(EEG.event);
        eventID = EEG_event_table.type;

        % eventInd=table2array(cell2table(squeeze(eventCell(1,1,:)))); % latency - timestamp of each event     MKA - this gets bepoch colum from EEG.event.  What is bepoch?
        eventInd = EEG_event_table.latency;
        % eventInd=table2array(cell2table(squeeze(EEG.event.type));

        % EEGeventTable = struct2table(EEG.event)
        eventTs=eventInd/EEG.srate;
        IEI=diff(eventTs);  % MKA- IEI = inter-epoch interval?
        % %     Key1=find(eventID==2);    %%%%%Two types of labels, either 2 or 61442 were used for different recording days to represent event "color change".
        % %     Key2=find(eventID==61442);

        eventTemp=unique(eventID);
        ColorChangeLabel={};

        for iEventTemp=1:length(eventTemp)
            Key=find(eventID==eventTemp(iEventTemp));  % find the index of events that equal ievent
            Key(Key>length(IEI))=[];
            tempLag(iEventTemp)=nanmedian(IEI(Key));
            if tempLag(iEventTemp)>0.19&&tempLag(iEventTemp)<0.21
                ColorChangeLabel=[ColorChangeLabel num2str(eventTemp(iEventTemp))];  % MKA -  Why is color change label is empty?
            end
        end

        %     ColorChangeLabel={};
        %     if ~isempty(Key1)      %%%%%%%check the interval is almost 0.2 seconds between color changing time and color changing back
        %        if median(IEI(Key1))>0.19&&median(IEI(Key1))<0.21
        %           ColorChangeLabel=[ColorChangeLabel '2'];
        %        end
        %     end
        %
        %     if ~isempty(Key2)
        %        if median(IEI(Key2))>0.19&&median(IEI(Key2))<0.21n
        %           ColorChangeLabel=[ColorChangeLabel '61442'];
        %
        %        end
        %     end

        %%%%%61442

        eventOriginalCount=0;
        for iEvent=1:length(EEG.event)
            for iColorChangeLabel=1:length(ColorChangeLabel)
                if EEG.event(iEvent).type==str2num(ColorChangeLabel{iColorChangeLabel})
                    eventOriginalCount=eventOriginalCount+1;
                end
            end
        end
        [OUTEEG, IndInlucdEpoch{iFile}] = pop_epoch(EEG, ColorChangeLabel, [-4 1],'epochinfo','yes');

        %     if iFile<=31
        %     [OUTEEG, IndInlucdEpoch{iFile}] = pop_epoch(EEG, {'2'}, [-2 0],'epochinfo','yes');
        %     else
        %     [OUTEEG, IndInlucdEpoch{iFile}] = pop_epoch(EEG, {'61442'}, [-2 0],'epochinfo','yes');
        %
        %     end
        numAllEpoch(iFile)=eventOriginalCount;
        numIncludeEpoch(iFile)=length(IndInlucdEpoch{iFile});
        [ChN(iFile),SampleN(iFile),TrialN(iFile)]=size(OUTEEG.data); % from EEG - TrialN should = numIncludeEpoch
        TrialNmat(iFile)=length(FileStruct(iFile).hits)+length(FileStruct(iFile).misses)+length(FileStruct(iFile).premature);

        TrialNmatMax(iFile)=max(union(union(FileStruct(iFile).hits,FileStruct(iFile).misses),FileStruct(iFile).premature));

        %     TrialNmat(iFile)=max(FileStruct(iFile).hits)+max(FileStruct(iFile).misses)+max(FileStruct(iFile).premature);

        TrialType(FileStruct(iFile).hits)=1;
        TrialType(FileStruct(iFile).misses)=0;
        TrialType(FileStruct(iFile).premature)=-1;
        EEGList{iFile}=OUTEEG;
        EEGList{iFile}.TrialType=TrialType;
        toc
        %     if length(TrialType)==TrialR
        %        ValidList=[ValidList;iFile];
        %     end

    end

    % TrialNVector=[TrialN(:) TrialNmat(:)]
    %% Counts
    CountEvent=[];
    CountName={};
    ReacT={};
    for iFile=1:length(FileStruct)
        if ~isempty(FileStruct(iFile).dotsynch)
            CountEvent(iFile,1)=length(FileStruct(iFile).dotsynch(:,1));  % should be from behavior
            CountEvent(iFile,4:6)=[length(FileStruct(iFile).hits) length(FileStruct(iFile).misses) length(FileStruct(iFile).premature)];

            CountEvent(iFile,3)=sum(CountEvent(iFile,4:6));
            ReacT{iFile}=(FileStruct(iFile).dotsynch(:,2)-FileStruct(iFile).dotsynch(:,1))/512;  %%%Reaction time in seconds.
        end
        if ~isempty(IndInlucdEpoch{iFile})
            CountEvent(iFile,10)= max(IndInlucdEpoch{iFile});  % col 10 max index of include trial should be = or smaller than num epoch, from EEG
        end
    end
    CountName{1}='DIMdotsynch';
    CountName{4}='Numhits';
    CountName{5}='Nummiss';
    CountName{6}='Numpremature';

    CountName{2}='MaxInd456';
    CountName{3}='NumInd456';

    CountName{7}='numEEGeventAll';
    CountName{8}='numEEGeventInc';
    CountName{9}='EEGDataDim';
    CountName{10}='MAXEEGeventInc';

    CountEvent(:,2)=TrialNmatMax(:);

    CountEvent(:,7)=numAllEpoch(:); % from EEG
    CountEvent(:,8)=numIncludeEpoch(:); % from EEG
    CountEvent(:,9)=TrialN(:); % from EEG

    % CountEvent(:,10)=numAllEpoch;

    TestTable=array2table(CountEvent,'VariableNames',CountName);
    TestTable(:,11) = {EEGSet_FileList.Subj}';
    A=TestTable(:,[1:2 7:10]);

    % important comparisons
    numIncludeEpoch_TrialN_diff = sum(abs(CountEvent(:,8)-CountEvent(:,9)))   %#ok<NOPRT> %%%Expected being 0, where dimension of eeg data match length of inlcuded index
    numAllEEGEpoch_vs_maxIndexofIncludeTrial_diff = sum(abs(CountEvent(:,7)-CountEvent(:,10)))   %#ok<NASGU,NOPRT> %%%Expected being 0, where dimension of eeg data match length of inlcuded index
    DIMdotsynch_vs_numAllEEGEpoch_diff = sum(abs(CountEvent(:,1)-CountEvent(:,7)))   %#ok<NOPRT> % %%%Expected being 0, where dimension of eeg data match length of inlcuded index
    % 1 from behavior, 7 from EEG, maybe


    % I1=find(CountEvent(:,1)-CountEvent(:,2)>0)

    % FileStruct(I1).name
    %
    % C:\Users\lzhang481\Box\Project_FlickerHealthyYoungAdults\Data\EEG\Data Analysis\Synching\Output\30-Nov-2022 15_00_48

    channels = {'Fp1', 'AF3', 'F7', 'F3', 'FC1', 'FC5', 'T7', 'C3', 'CP1', 'CP5', 'P7', 'P3', 'Pz', 'PO3', 'O1', 'Oz', 'O2', 'PO4', 'P4', 'P8', 'CP6', 'CP2', 'C4', 'T8', 'FC6', 'FC2', 'F4', 'F8', 'AF4', 'Fp2', 'Fz', 'Cz', 'EXG1', 'EXG2', 'EXG3', 'EXG4', 'EXG5', 'EXG6', 'EXG7', 'EXG8'};
    EEGchInd=1:32;
    EEGch=channels(EEGchInd);

    ChNTotal=length(channels);
    NeedFields={'labels','theta','radius'};
    tempN=fieldnames(EEGList{1}.chanlocs);
    tempN1=tempN;
    tempN2=tempN;
    tempN1([1 2 11 12])=[];
    NeedI=setdiff(1:length(tempN),[1 2 12]);
    tempN2([2 12])=[];

    ChanPos=FieldName2Struct(tempN1);
    TempTable=cell2table(squeeze(struct2cell(EEGList{iFile}.chanlocs)));
    ChInd=table2array(cell2table(table2array(TempTable(11,:))));
    IndEx=find(ChInd>=33);

    ChIndWritten=[];
    DataTemp=zeros(length(tempN1),length(EEGchInd))+nan;

    %%
    for iFile=1:length(EEGList)
        TempTable=cell2table(squeeze(struct2cell(EEGList{iFile}.chanlocs)));
        ChI=table2array(cell2table(table2array(TempTable(11,:))));
        IndEx=find(ChI>=33);

        Invalid=find(ChI>=33);
        ChIndV=ChI;
        ChIndV(Invalid)=[];
        TempTable(:,Invalid)=[];

        InforAdd=TempTable;
        InforAdd([1 2 11 12],:)=[];

        [ChAdd,I1]=setdiff(ChIndV,ChIndWritten);
        temp=table2array(cell2table(table2array(InforAdd)));
        DataTemp(:,ChAdd)=temp(:,I1);
        ChIndWritten=find(isnan(DataTemp(1,:))==0);
        if length(ChIndWritten)==length(EEGchInd)
            break
        end
    end


    for ifield=1:length(tempN1)
        ChanPos=setfield(ChanPos,{1,1},tempN1{ifield},DataTemp(ifield,:));
    end
    ChanPos.labels=EEGch;

    ChanEEGLab=rmfield(ChanPos,'labels');
    tempName=fieldnames(ChanEEGLab);
    for iCh=1:length(EEGchInd)
        for iN=1:length(tempName)
            dataTemp=getfield(ChanPos,{1},tempName{iN});
            ChanEEGLab=setfield(ChanEEGLab,{1,iCh},tempName{iN},dataTemp(iCh));
        end
        ChanEEGLab(iCh).labels=EEGch{iCh};
        ChanEEGLab(iCh).urchan=iCh;
    end

    %% Plots channel locations
    % figure;
    % plot(ChanPos.Y,ChanPos.X,'r.');
    % for iCh=1:length(ChanPos.X)
    %     text(ChanPos.Y(iCh),ChanPos.X(iCh),EEGch{iCh});
    % end

    %% Channel Count (?)
    ChCount=[];
    ChIndList={};
    for iFile=1:length(EEGList)
        if ~isempty(EEGList{iFile})
            ChCount(iFile)=length(EEGList{iFile}.chanlocs); % generally 38 for new semiprocessed EEGs
            tempData=EEGList{iFile}.data;
            AllChData=zeros(ChNTotal,size(tempData,2),size(tempData,3))+nan;

            for iCh=1:ChCount(iFile)
                [~,ChIndList{iFile}(iCh)]=ismember(EEGList{iFile}.chanlocs(iCh).labels,channels);
            end
            AllChData(ChIndList{iFile},:,:)=tempData;
            EEGList{iFile}.AllChData=AllChData;

            TrialType=zeros(length(FileStruct(iFile).dotsynch(:,1)),1)+nan;
            TrialType(FileStruct(iFile).hits)=1;
            TrialType(FileStruct(iFile).misses)=0;
            TrialType(FileStruct(iFile).premature)=-1;
            ReacTime=FileStruct(iFile).dotsynch(:,2)-FileStruct(iFile).dotsynch(:,1);

            TrialTypeEEG=TrialType(IndInlucdEpoch{iFile});
            ReacTime=ReacTime(IndInlucdEpoch{iFile});

            if size(tempData,3)~=length(TrialTypeEEG)
                disp('Dimention non matching for EEG and Behavior trials');
                iFile
                TrialTypeEEG=[];
            else
                EEGList{iFile}.TrialType=TrialTypeEEG;
                EEGList{iFile}.ReactionTime=ReacTime;
            end
            clear tempData AllChData TrialType TrialTypeEEG;
        end
    end

    for iFile=1:length(EEGList)
        test(iFile)=length(EEGList{iFile}.TrialType)-length(EEGList{iFile}.ReactionTime);
    end
    %% Save
    saveStart = tic;

    save([save_output_path '\Summary.mat'],'EEGfolderPath','EEGSet_FileList', 'SyncFolderPath', 'SyncList', 'save_output_path', 'scriptName', 'nFiles', 'newestFolder')
    save([save_output_path '\S2Data.mat'],'FileStruct','EEGList','TestTable','ChCount','ChIndList','ChanPos','ChanEEGLab','IndInlucdEpoch','-v7.3')
    
    saveDuration = toc(saveStart);
    
    endDateAndTime = char(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss'));
    fprintf('Step2 completed at: %s\nSyncing duration: %g seconds.', endDateAndTime, saveDuration);

    %% Create output text
    %Generate output message
    SummaryTextOutput='';
    SummaryTextOutput = sprintf('Step2 started on %s\nScript used: %s \n', currDate, scriptName);
    SummaryTextOutput =  sprintf('%s \nInput dataset: %s \nFrom folder: %s', SummaryTextOutput,  dataset_folder,EEGfolderPath);
    % SummaryTextOutput =  sprintf('%s \n%d out of %d files synced:\n%s', SummaryTextOutput, iCsv, nCsv, subjectsList);

    endDuration = toc(S2Start);
    SummaryTextOutput = sprintf('%s\nStep2 CheckSync completed at: %s\nRun duration: %g seconds.', SummaryTextOutput, endDateAndTime, endDuration);

    %Display the output message in the command window
    disp(SummaryTextOutput);

    %Create a .txt file with the output message
    filePath = fullfile(save_output_path, [scriptName, currDate '_Summary.txt']);

    fileID = fopen(filePath, 'w');
    if fileID == -1
        error('Failed to open or create the file: %s', filePath);
    else
        fprintf(fileID, '%s', SummaryTextOutput);
        fclose(fileID);
        fprintf('File created and written: %s\n', filePath);
    end
end