%% Step 4 - Does PSD.  (derived from Step2_WPLI by Lu Z)TrialType
% Must add GenMatCode-main and all subfolders to path before running
S4Start = tic;
currDate = strrep(datestr(datetime), ':', '_');
currDate = datestr(datetime, 'yy-mm-dd_HHMMSSFFF');
scriptName = mfilename;
%% Load All EEGs
% LoadPath='Y:\singer\LuZhang\Project6-EEG\Results\Step0-PreparingData\';
cd("2_CheckSync_Outputs\")

%% Get newest file
newestFolder = getNewestFolder();
%cd(newestFolder);
%newestEEGData = getNewestFile();  % Automatically gets the newest .mat file in folder
%newestEEGData=
%% Select file to load
inputDatafile = '01-Apr-2024 13_38_57EEGDataStep2CheckSyncOutput.mat' %#ok<NOPTS> % can replace variable with hardcoded EEG file
tic
load(inputDatafile)
%load('01-Apr-2024 13_38_57EEGDataStep2CheckSyncOutput.mat')

%% Create folder
disp(['EEG data loaded:', inputDatafile])
toc
% load('02-Feb-2024 17_46_26EEGData4EyeBlinkRemoval.mat')
cd ..
%cd ..

SaveFolder=['4_COH\' currDate '_' scriptName];
mkdir(SaveFolder);

save([SaveFolder '\Step4DataInputSummary.mat'],'inputDatafile','currDate','scriptName')
%% Create output text
%Generate output message
SummaryTextOutput='';
% if isempty(SummaryTextOutput)
%     SummaryTextOutput = 'None';
% end

SummaryTextOutput = sprintf('Step4 started on %s\nScript used: %s \n', currDate, scriptName);
SummaryTextOutput =  sprintf('%s \nInput dataset: %s \n', SummaryTextOutput,  inputDatafile);
% SummaryTextOutput =  sprintf('%s \n%d out of %d files synced:\n%s', SummaryTextOutput, iCsv, nCsv, subjectsList);
% endDateAndTime = char(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss'));
% endDuration = toc(tStart);
% SummaryTextOutput = sprintf('%s\nStep1Syncing completed at: %s\nSyncing duration: %g seconds.', SummaryTextOutput, endDateAndTime, endDuration);
% Add toc

%Display the output message in the command window
disp(SummaryTextOutput);

%Create a .txt file with the output message
filePath = fullfile(SaveFolder, [scriptName, currDate '_Summary.txt']);

fileID = fopen(filePath, 'w');
if fileID == -1
    error('Failed to open or create the file: %s', filePath);
else
    fprintf(fileID, '%s', SummaryTextOutput);
    fclose(fileID);
    fprintf('File created and written: %s\n', filePath);
end
%% Set Channels, Locations for Heatmap
% (Necessary if already done in Step 0?)
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
TempTable=cell2table(squeeze(struct2cell(EEGList{1}.chanlocs)));
ChInd=table2array(cell2table(table2array(TempTable(11,:))));
IndEx=find(ChInd>=33);

ChIndWritten=[];
DataTemp=zeros(length(tempN1),length(EEGchInd))+nan;
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

%% Assign Condition Groups  *** Run to set GroupName

FlickerSubj{4}=[20:27 30:33 44:51 53 55:58 62:64 66:69 71:73 75 77:81 83:87];   % All control subjects (Random and Light) no cuts
GroupName{4}='BothControls'; % 'BothControls' = Random and Light together.  previously 'Random' or 'Light'

FlickerSubj{5}=[20:27 30:33 44:51 53 55:58 62:63 66 68 69 71:73 75 77:81 83:87];   % All control subjects minus 64 and 67 (for RT)
GroupName{5}='BothControlsRT'; % 'BothControls' = Random and Light together.  previously 'Random' or 'Light'

FlickerSubj{1}=[10:19 34:42 52 54 59:61 65 70 74 76 82];      % SubjID of 40Hz group.
GroupName{1}='40Hz';

FlickerSubj{2} = [62:64 66:69 71:73 75 77:81 83:87];  % All SubjID for Light Group (no cuts)
GroupName{2}='Light';

FlickerSubj{3}=[20:27 30:33 44:51 53 55:58];   % SubjID of Random flicker group
GroupName{3} = 'Random';

FlickerSubj{6} = [62:63 66 68 69 71:73 75 77:81 83:87];  % All SubjID for Light Group for RT (removed s064 and s067 due to incorrect average RTs
GroupName{6}='LightRT'; % No 64 67 for WPLI vs RT
% Updated Groups (includes 2023 EEGs)
% FlickerSubj{1}=[20:28 30:33 45 48:50 56:58];   % (20 total) SubjID of Random flicker group (Accuracy cuts)
% FlickerSubj{2}=[10:19 39 40 52 54 59 60 74 76 82];      %  (20 total) % SubjID of 40Hz group. Newest 40 Hz EEGs added (2023) (Accuracy Cut
% FlickerSubj{3}=[];  % check if ID is still matching

%% Get Accuracies for all subjects
Acc=zeros(length(FileStruct),1)+nan;
SubjsAvgRT = zeros(length(FileStruct),1)+nan;
for iFile=1:length(FileStruct)
    if ~isempty(FileStruct(iFile).Subj)
        %% May need to uncomment out the below
        SubjID(iFile) = str2num(FileStruct(iFile).Subj(end-1:end));
        Acc(iFile)=length(FileStruct(iFile).hits)/(length(FileStruct(iFile).hits)+length(FileStruct(iFile).misses));

        %% Calculate Avg RT per subject using data found in FileStruct->dotsynch
        isHit = (FileStruct(iFile).dotsynch(:,3)==1);  % get logical index for all trials that are Hits (misses and premature hits will = 0)
        colorchangeTimes = FileStruct(iFile).dotsynch(isHit,1); % get time of color change for Hit trials only
        subjRTtimes = FileStruct(iFile).dotsynch(isHit,2); % get RT time for Hit trials (this should ignore premature hits)
        subjRTduration =  subjRTtimes - colorchangeTimes; % subtract RT time from color change time to get duration of RT (AKA the reaction time)
        subjRTdurInSecs = subjRTduration/512; % sample rate is generally 512 samples per second
        SubjsAvgRT(iFile) = mean(subjRTdurInSecs);
    end
end

%% SubjG assignment for WPLI calculation
SubjG{1}=[];
SubjG{2}=[];
SubjG{3}=[];
SubjG{4}=[];
SubjG{5}=[];
SubjG{6}=[];

[~,SubjG{1},~]=intersect(SubjID,FlickerSubj{1}); % 40Hz
[~,SubjG{2},~]=intersect(SubjID,FlickerSubj{2}); % Light
[~,SubjG{3},~]=intersect(SubjID,FlickerSubj{3}); % Random
[~,SubjG{4},~]=intersect(SubjID,FlickerSubj{4}); % Both Controls
[~,SubjG{5},~]=intersect(SubjID,FlickerSubj{5}); % BothControlsRT
[~,SubjG{6},~]=intersect(SubjID,FlickerSubj{6}); % LightRT
%% Set Trial Groups (Hit, Miss, Etc.)
ChCount=[];
ChIndList={};
TrialType{1}=1; %%%%%Hit trial
TrialType{2}=0; %%%%%Miss trial
TrialType{3}=[0 1]; %%%%%All trial
% TrialType{4}=-1; %%%%%Premature trial

TrialTypeName{1}='Hit'; %%%%%Hit trial
TrialTypeName{2}='Miss'; %%%%%Miss trial
TrialTypeName{3}='HitAndMiss'; %%%%%All trial
% TrialTypeName{4}='Premature'; %%%%%Premature trial


%% PSD parameters
psdParameter.Fs=512;
psdParameter.window=1024;  % can increase window to 1024 from 512 8/18/24 - less smoth, but higher res
psdParameter.noverlap = psdParameter.window/2;        % can change overlap to 256 or half window size
psdParameter.nfft=1024;     % decrease nfft from 1024 to 512 - less smooth but decrease resolution
nFre=psdParameter.nfft/2+1; % modified pwelch see Lu genmat code on git

%% Epoch, Samp Rate, Make Save Trial Folder
% (When is data epoched? What if data is already epoched from preprocessing)
DataTimeRange=[-4 1];     %%%4 seconds before color-change and 1s after.
AnaRange=[-4 0]; %%4s before color-change
% AnaRange=[0 1]; %%1s after color change

SampRate=512;
SampI=(AnaRange-DataTimeRange(1))*SampRate;
SampI=SampI(1)+1:SampI(2);


% clear CohGroup
parfor iFile=1:length(EEGList)
    ChTempN(iFile)=size(EEGList{iFile}.AllChData,1);
end

SavePath = ['4_COH\' currDate '_' scriptName '\'];
if ~exist(SavePath, 'dir')
    mkdir(SavePath)
end

%% Get coherence for all pairs of channels PER SUBJECT!
% This takes a looong time! (Data from here goes in TrialCrossSpec)

% SaveTrialSubj='Y:\singer\LuZhang\Project6-EEG\Results\Step2-COH\TrialCrossSpec\';
SaveTrialSubj=['4_COH\' currDate '_' scriptName '\TrialCrossSpec\'];
mkdir(SaveTrialSubj);
parpool(12)

tic
for iFile=1:length(EEGList)
    SaveTemp=[SaveTrialSubj EEGList{iFile}.filename(1:4) '\'];
    mkdir(SaveTemp);
    if ~isempty(EEGList{iFile})
        ChTempN=size(EEGList{iFile}.AllChData,1);
        TrialTypeTemp=EEGList{iFile}.TrialType;
        TrialI=[];
        %%  Calculate for pair of channels
        tic
        for iCh=1:ChTempN
            for jCh=iCh:ChTempN
                clear TempTrial1 tempSig1 TempTrial2 tempSig2 TrialSpec
                tempSig1=squeeze(EEGList{iFile}.AllChData(iCh,SampI,:));
                tempSig2=squeeze(EEGList{iFile}.AllChData(jCh,SampI,:));

                if sum(sum(isnan(tempSig1)))>1||sum(sum(isnan(tempSig2)))>1
                    continue
                end
                % tic
                for iTrial=1:length(TrialTypeTemp)  % loop by trial
                    if length(TrialTypeTemp)>1
                        TempTrial1(iTrial).Data=tempSig1(:,iTrial);
                        TempTrial2(iTrial).Data=tempSig2(:,iTrial);

                    else
                        TempTrial1(iTrial).Data=tempSig1;
                        TempTrial2(iTrial).Data=tempSig2;

                    end
                    TempTrial1(iTrial).Time=([1:length(TempTrial1(iTrial).Data)]-1)/512;
                    TempTrial2(iTrial).Time=([1:length(TempTrial2(iTrial).Data)]-1)/512;


                end
                %                  clear TrialSpec
                %%%Old version to calculate CrossSpec,tested equal to
                %%%new version
                %                  [TrialSpec1.Sxy,TrialSpec1.Sxx,TrialSpec1.Syy,TrialSpec1.w,TrialSpec1.options,ValidIndex]=coh_TrialData(TempTrial1,TempTrial2,psdParameter);
                %%%Old version to calculate CrossSpec,tested equal to
                %                    %%%new version
                %                   psdParameter.noverlap=500;
                %                   psdParameter.nfft=512;
                %                   psdParameter.window=512;
                %                   nFre=psdParameter.nfft/2+1;

                [TrialSpec.Sxy,TrialSpec.Sxx,TrialSpec.Syy,TrialSpec.w,TrialSpec.options,ValidIndex]=crossspec_EqualTriL(TempTrial1,TempTrial2,psdParameter);  %% find in genmat code
                save([SaveTemp 'Ch' num2str(iCh) 'Ch' num2str(jCh) '.mat'],'TrialSpec','ValidIndex','psdParameter');

                %                  a=crossspec_Trial(TrialSpec);
                %                  figure;
                %                  plot(a.Fre,abs(((a.wpli))))
                %                  figure;
                %                  plot(a.Fre,abs(mean((a.wpli(1:30,:)))))
                %                  figure;
                %                  plot(a.Fre,abs(mean((a.wpli(1:30,:)))))

                %                  figure;
                %                  plot(TempTrial1(4).Data);hold on;plot(TempTrial2(4).Data,'r.')
                for iTrialType=1:length(TrialType)  % group: hit, miss, etc.
                    TrialI=[];
                    parfor j=1:length(TrialType{iTrialType})
                        TrialI=union(TrialI,find(TrialTypeTemp==TrialType{iTrialType}(j)));
                    end
                    TrialI=intersect(TrialI,ValidIndex);
                    if ~isempty(TrialI)
                        %                    CohGroup{iGroup,iFile}{iCh,jCh}=Coh_TrialIndex(TrialSpec1,TrialI);
                        CohGroup{iTrialType,iFile}{iCh,jCh}=crossspec_TrialIndex(TrialSpec,TrialI);

                        %% %confirmed Old and New version of Cross-Spectrum results in same coherence results.
                        %                    D1=Coh_TrialIndex(TrialSpec1,TrialI);
                        %                    D2=crossspec_TrialIndex(TrialSpec,TrialI);
                        %                    figure;
                        %                    plot(D1.Fre,(D1.Cxy));hold on;
                        %                    plot(D2.Fre,(D2.Cxy),'r.');hold on;c
                        %%%confirmed Old and New version of Cross-Spectrum results in same coherence results.

                    end
                end
                % toc
                %%                  figure;
                %                  Temp=CohGroup{1,iFile}{iCh,jCh};
                %                  subplot(2,1,1)
                %                  plot(Temp.Fre,(Temp.Cxy));
                %                  subplot(2,1,2)
                %
                %                  plot(Temp.Fre,log(abs(Temp.Pxx)));
                %                  hold on;
                %                  plot(Temp.Fre,log(abs(Temp.Pyy)));
            end
        end
        toc
    end
end
toc
% Check if everything above works ***


%% Save workspace to be used for Step 5: WPLITrial Group.  This step takes a long time, creates 20gb file!
COHSaveFileName =['COHdata_forWPLITrialType_' currDate '_' scriptName '.mat'];
COH_Save_Path = [SavePath COHSaveFileName];
save(COH_Save_Path,'-v7.3')  % Check this
% load(COH_Save_Path)  % loading should not be necessary as all variables in workspace should be in that save file

%% For visualization - Set-up - Must run before plotting anything below
% (of what? COH,WPLI and PSD?)
load('chanPosColin27');

%
Fre=CohGroup{1,1}{1,2}.Fre;  % may need to import CohGroup from a COHdata file
FBand=[2 100]; % consider changing [1 100] to [2 100] due to normalization
FreInd=find(Fre>=FBand(1)&Fre<=FBand(2));
Fplot=Fre(FreInd);

PSDall=zeros(length(FileStruct),length(FreInd),ChNTotal,length(TrialType))+nan;
COHall=zeros(length(FileStruct),length(FreInd),ChNTotal,ChNTotal-1,length(TrialType))+nan;
WPLIall=zeros(length(FileStruct),length(FreInd),ChNTotal,ChNTotal-1,length(TrialType))+nan;

%
for iTrialType=1:length(TrialType)
    for iFile=1:length(FileStruct)
        for iCh=1:size(CohGroup{iTrialType,iFile},1)
            for jCh=iCh+1:size(CohGroup{iTrialType,iFile},2)
                if isempty(CohGroup{iTrialType,iFile}{iCh,jCh})
                    continue;
                end
                if iCh==1
                    PSDall(iFile,:,iCh,iTrialType)=CohGroup{iTrialType,iFile}{iCh,jCh}.Pxx(FreInd);
                    PSDall(iFile,:,jCh,iTrialType)=CohGroup{iTrialType,iFile}{iCh,jCh}.Pyy(FreInd);
                end
                COHall(iFile,:,iCh,jCh,iTrialType)=CohGroup{iTrialType,iFile}{iCh,jCh}.Cxy(FreInd);
                WPLIall(iFile,:,iCh,jCh,iTrialType)=CohGroup{iTrialType,iFile}{iCh,jCh}.wpli(FreInd);
            end
        end
    end
end

%% Peak Alpha WPLI Distribution Histogram
% Get all alpha values
alphaLowerLimitFreqHz = 8;
alphaUpperLimitFreqHz = 13;

% Find indices of values in Fre between 8 and 13 (inclusive)
allAlphaFreIndices = find(Fre >= alphaLowerLimitFreqHz & Fre <= alphaUpperLimitFreqHz);
HitMissTrialType = 3;
WPLIallAlpha = squeeze(WPLIall(:,allAlphaFreIndices,1:32,1:32,HitMissTrialType)); % size(WPLIall) ans = 67 197 40 39 3- 67subs x allFreqs x iCh x jCh x TrialType

% Find the peak value within alpha (dimension 2)
[peakAlphaValues, peakAlphaIndices] = max(WPLIallAlpha, [], 2);

% Reshape the result to 3D
PeakWPLIallAlpha = squeeze(peakAlphaValues);

%% Plot distribution of PeakWPLIallAlpha
% Flatten PeakWPLIallAlpha to 1D for distribution analysis
PeakWPLIallAlphaFlat = PeakWPLIallAlpha(:);

% Remove NaN values
PeakWPLIallAlphaFlat = PeakWPLIallAlphaFlat(~isnan(PeakWPLIallAlphaFlat));

% Calculate the total number of data points (channel pairs)
nDataPoints = numel(PeakWPLIallAlphaFlat);

% Compute top percentiles
top25Percent = prctile(PeakWPLIallAlphaFlat, 75);
top10Percent = prctile(PeakWPLIallAlphaFlat, 90);
top5Percent = prctile(PeakWPLIallAlphaFlat, 95);
top1Percent = prctile(PeakWPLIallAlphaFlat, 99);
top0_1Percent = prctile(PeakWPLIallAlphaFlat, 99.9);
top0_01Percent = prctile(PeakWPLIallAlphaFlat, 99.99);

% Count data points greater than each percentile
countAbove25Percent = sum(PeakWPLIallAlphaFlat > top25Percent);
countAbove10Percent = sum(PeakWPLIallAlphaFlat > top10Percent);
countAbove5Percent = sum(PeakWPLIallAlphaFlat > top5Percent);
countAbove1Percent = sum(PeakWPLIallAlphaFlat > top1Percent);
countAbove0_1Percent = sum(PeakWPLIallAlphaFlat > top0_1Percent);
countAbove0_01Percent = sum(PeakWPLIallAlphaFlat > top0_01Percent);

% Display the results
fprintf('Top 25%% WPLI value: %.4f, Data points above: %d\n', top25Percent, countAbove25Percent);
fprintf('Top 10%% WPLI value: %.4f, Data points above: %d\n', top10Percent, countAbove10Percent);
fprintf('Top 5%% WPLI value: %.4f, Data points above: %d\n', top5Percent, countAbove5Percent);
fprintf('Top 1%% WPLI value: %.4f, Data points above: %d\n', top1Percent, countAbove1Percent);
fprintf('Top 0.1%% WPLI value: %.4f, Data points above: %d\n', top0_1Percent, countAbove0_1Percent);
fprintf('Top 0.01%% WPLI value: %.4f, Data points above: %d\n', top0_01Percent, countAbove0_01Percent);

% Plot the histogram of PeakWPLIallAlpha
figure;
histogram(PeakWPLIallAlphaFlat, 'Normalization', 'probability', 'BinWidth', 0.02);
hold on;

% Add vertical lines for top percentiles
xline(top25Percent, '--k', 'Top 25%', 'LineWidth', 1.5);
xline(top10Percent, '--c', 'Top 10%', 'LineWidth', 1.5);
xline(top5Percent, '--r', 'Top 5%', 'LineWidth', 1.5);
xline(top1Percent, '--g', 'Top 1%', 'LineWidth', 1.5);
xline(top0_1Percent, '--b', 'Top 0.1%', 'LineWidth', 1.5);
xline(top0_01Percent, '--m', 'Top 0.01%', 'LineWidth', 1.5);

% Label the axes
xlabel('WPLI');
ylabel('Fraction');

% Add a title including the total number of data points
title(['Distribution of Peak Alpha WPLI (Total data points: ', num2str(nDataPoints), ')']);

% Improve plot appearance
grid on;


% % Display the size of the resulting array
% disp('Size of the resulting 3D array:');
% disp(size(PeakWPLIallAlpha));
%% Plot the histogram of PeakWPLIallAlpha (SuppFig4A MS version)
% May need to run the previous section first! MKA 2025-03-18
figure;
histogram(PeakWPLIallAlphaFlat, 'Normalization', 'probability', 'BinWidth', 0.02);
hold on;

% Add vertical lines for top percentiles
xline(top25Percent, '--k', 'Top 25%', 'LineWidth', 1.5);

% Label the axes
xlabel('WPLI');
ylabel('Fraction');

% Add a title including the total number of data points
title(['Distribution of Peak Alpha WPLI (Total data points: ', num2str(nDataPoints), ')']);

% Improve plot appearance
grid on;

saveas(gcf, fullfile('Fig4Panels', 'SuppFig4A_PeakWPLI_Alpha.svg'));

%% Peak Alpha WPLI 40Hz vs Light ALL CHANNELS (not used) - search "stats preceding fig4D" for signif channels only
% Find the peak WPLI value within alpha (dimension 2)
[peakAlphaValues, peakAlphaIndices] = max(WPLIallAlpha, [], 2, "includemissing");

alphaFreqHz = 8:.5:13;

% Initialize a copy of peakAlphaIndices
peakAlphaIndicesNaN = peakAlphaIndices;

% If all alpha WPLI values are NaN, then set the max index to NaN
peakAlphaIndicesNaN(all(isnan(WPLIallAlpha),2)) = NaN;

% Reshape the result to 3D
PeakAlphaIndicesNaN3D = squeeze(peakAlphaIndicesNaN);

% % Check with single participant
% singleparticiantAllChPairPeakAlpha = squeeze(PeakAlphaIndicesNaN3D(1,:,:))
% 
% % Convert Indices to Correct Corresponding Frequnecy (Hz)
% % Find valid indices (values between 1 and 11)
% validIdx = singleparticiantAllChPairPeakAlpha >= 1 & singleparticiantAllChPairPeakAlpha <= 11;
% 
% % Initialize the output array with NaN, preserving original shape
% singleallchpfreqs = NaN(size(singleparticiantAllChPairPeakAlpha));
% 
% % Perform mapping only for valid indices
% singleallchpfreqs(validIdx) = alphaFreqHz(singleparticiantAllChPairPeakAlpha(validIdx));

% Convert Indices (3D) to Correct Corresponding Frequncy (Hz) 
validIdx = PeakAlphaIndicesNaN3D >= 1 & PeakAlphaIndicesNaN3D <= 11;

% Initialize the output array with NaN, preserving original shape
PeakAlphaFreqs = NaN(size(PeakAlphaIndicesNaN3D));

% Perform mapping only for valid indices
PeakAlphaFreqs(validIdx) = alphaFreqHz(PeakAlphaIndicesNaN3D(validIdx));
% singleallchpfreqs = alphaFreqHz(singleparticiantAllChPairPeakAlpha)

GroupSubjs_40Hz = 1; % 40Hz group
GroupSubjs_Light = 6; % LightRT group

% % Extract Peak Alpha WPLI into groups: 40 Hz & Light
% % PeakWPLIallAlpha: 67 subs x 32ch x 32ch
% PeakAlphaWPLI_40HzGroup = PeakWPLIallAlpha(SubjG{GroupSubjs_40Hz},:,:);
% PeakAlphaWPLI_LightGroup = PeakWPLIallAlpha(SubjG{GroupSubjs_Light},:,:);
% % single particiapn
% singledudeWPLIallpair = squeeze(PeakWPLIallAlpha(1,:,:))


% Extract Peak Alpha WPLI freq into groups: 40 Hz & Light
% PeakWPLIallAlpha: 67 subs x 32ch x 32ch
PeakAlphaFreqs_40HzGroup = PeakAlphaFreqs(SubjG{GroupSubjs_40Hz},:,:);
PeakAlphaFreqs_LightGroup = PeakAlphaFreqs(SubjG{GroupSubjs_Light},:,:);

% Get average peak alpha WPLI frequncy for each group for all channel pairs
MeanPeakAlphaFreq_40Hz = squeeze(mean(PeakAlphaFreqs_40HzGroup, 1));
MeanPeakAlphaFreq_Light = squeeze(mean(PeakAlphaFreqs_LightGroup, 1));

% Display average peak alpha frequency tables (32x32)
ChannelLabels = {'Fp1', 'AF3', 'F7', 'F3', 'FC1', 'FC5', 'T7', 'C3', 'CP1', 'CP5', ...
                 'P7', 'P3', 'Pz', 'PO3', 'O1', 'Oz', 'O2', 'PO4', 'P4', 'P8', ...
                 'CP6', 'CP2', 'C4', 'T8', 'FC6', 'FC2', 'F4', 'F8', 'AF4', 'Fp2', 'Fz', 'Cz'};

fprintf('\nMean Peak Alpha Frequency (40Hz Group):\n');
disp(array2table(MeanPeakAlphaFreq_40Hz, 'VariableNames', ChannelLabels, 'RowNames', ChannelLabels));

fprintf('\nMean Peak Alpha Frequency (Light Group):\n');
disp(array2table(MeanPeakAlphaFreq_Light, 'VariableNames', ChannelLabels, 'RowNames', ChannelLabels));


% Flatten both 32 x 32 arrays into two 1D-arrays (492 channel pairs each)
UpperTriIdx = find(triu(ones(32, 32), 1)); % Indices of upper triangular elements
FlattenedAlphaFreq_40Hz = MeanPeakAlphaFreq_40Hz(UpperTriIdx);
FlattenedAlphaFreq_Light = MeanPeakAlphaFreq_Light(UpperTriIdx);

% Compute average peak alpha frequency across all channel pairs for each group
AvgPeakAlphaFreq_40Hz = mean(FlattenedAlphaFreq_40Hz);
AvgPeakAlphaFreq_Light = mean(FlattenedAlphaFreq_Light);

fprintf('\nAverage Peak Alpha Frequency Across All Channel Pairs:\n');
fprintf('40Hz Group: %.4f Hz\n', AvgPeakAlphaFreq_40Hz);
fprintf('Light Group: %.4f Hz\n', AvgPeakAlphaFreq_Light);


% Test for normality (to decide if t-test or ranksum)
[H_40Hz, p_40Hz] = kstest(FlattenedAlphaFreq_40Hz);
[H_Light, p_Light] = kstest(FlattenedAlphaFreq_Light);

fprintf('\nNormality Test Results:\n');
fprintf('40Hz Group: H = %d, p = %.16f\n', H_40Hz, p_40Hz);
fprintf('Light Group: H = %d, p = %.16f\n', H_Light, p_Light);

% Decide on statistical test
if H_40Hz == 0 && H_Light == 0
    % Normally distributed: Use independent t-test
    [h_ttest, p_ttest] = ttest2(FlattenedAlphaFreq_40Hz, FlattenedAlphaFreq_Light);
    test_used = 't-test';
    p_value = p_ttest;
else
    % Non-normally distributed: Use Wilcoxon rank-sum test
    [p_ranksum, h_ranksum] = ranksum(FlattenedAlphaFreq_40Hz, FlattenedAlphaFreq_Light);
    test_used = 'Wilcoxon rank-sum test';
    p_value = p_ranksum;
end

% Display results
fprintf('Statistical Test Used: %s\n', test_used);
fprintf('p-value: %.16f\n', p_value);

% Perform one-sided Wilcoxon rank-sum test
[p_ranksum_right, h_ranksum_right] = ranksum(FlattenedAlphaFreq_40Hz, FlattenedAlphaFreq_Light, 'tail', 'right');
[p_ranksum_left, h_ranksum_left] = ranksum(FlattenedAlphaFreq_40Hz, FlattenedAlphaFreq_Light, 'tail', 'left');

% Display results
fprintf('\nOne-Sided Wilcoxon Rank-Sum Test Results:\n');
fprintf('H0: 40Hz <= Light | p-value (right-tailed, 40Hz > Light): %.16f\n', p_ranksum_right);
fprintf('H0: 40Hz >= Light | p-value (left-tailed, 40Hz < Light): %.16f\n', p_ranksum_left);


% Create Violin plots showing distribution
% Create figure
figure;
hold on;

% Combine data for violin plot
groupLabels = [repmat({'40Hz'}, length(FlattenedAlphaFreq_40Hz), 1); ...
               repmat({'Light'}, length(FlattenedAlphaFreq_Light), 1)];
data = [FlattenedAlphaFreq_40Hz; FlattenedAlphaFreq_Light];

% Create violin plot
violinplot(data, groupLabels);

% Format plot
title('Violin Plot of Peak Alpha Frequency');
ylabel('Peak Alpha Frequency (Hz)');
xlabel('Group');
ylim([8 13]); % Set y-axis range from 8 Hz to 13 Hz
grid on;
hold off;




%% PSD related parameters
LogPSDraw=log(abs(PSDall));
NoiseInd=find(Fplot>=58&Fplot<=62);
NormBandI=setdiff(1:length(Fplot),NoiseInd);
PSDall=PSDall./repmat(nansum(PSDall(:,NormBandI,:,:),2),1,length(Fplot),1,1);

LogPSD=log(abs(PSDall));

PlotColor2=[1 0 0;0 0 1];
% ParamPSD.ANOVAstats='Anova';
ParamPSD.PlotType=3;
ParamPSD.SigPlot='Anova';
ParamPSD.SigPlot='Ttest';
ParamPSD.CorrName='fdr';   %%%methold for multi-compairson
ParamPSD.Q=0.1;
ParamPSD.Ytick=[0 0.002 0.004];
ParamPSD.LegendShow=0;
ParamPSD.Legend=[];
ParamPSD.TimeRepeatAnova=1;
ParamPSD.GroupRepeatAnova=0;
ParamPSD.RepeatAnova=0;
ParamPSD.TimeCol=Fplot;
ParamPSD.Paired=1;
ParamPSD.BinName='Fre';
ParamPSD.Bin=Fplot;
ParamPSD.TimeComparison=0;
ParamPSD.statisP=1; % 1 to do stats and plot. 0 will do stats, but not plot, will be faster.  Uses R.  If error, set as 0
ParamPSD.Ytick=[-8:4:0];
ParamPSD.Crit_p=0.05;


% One color for ea of the 3 groups.  Blue for 40, Gold/yellow/orange for Light, Red for Random
FlickerColor=[31 125 184; 219 129 50; 150 27 27]/255;
FlickerColor=[31 125 184; 150 27 27; 219 129 50]/255; %40, Random, Light
% FlickerColor=[0.5 0.5 0.5;0.9 0.1 0.3];

load('Functions\GenMatCode-main\Plotfun\Color\colorMapPN.mat')
load('Functions\GenMatCode-main\Plotfun\Color\colorMapPNraw.mat')

%% WPLI - Weight Phase Lag Index - Parameters
ParamWPLI=ParamPSD;
ParamWPLI.Ytick= [0:0.1:0.2]; %#ok<NBRAK2>
SubSaveWPLI=[SavePath 'WPLI\'];
ParamWPLI.SigPlot='Anova';
mkdir(SubSaveWPLI)
ParamWPLI.statisP=1;
P.xLeft=0.01;        %%%%%%Left Margin
P.xRight=0.01;       %%%%%%Right Margin
P.yTop=0.01;         %%%%%%Top Margin
P.yBottom=0.01;      %%%%%%Bottom Margin
P.xInt=0.005;        %%%%%%Width-interval between subplots
P.yInt=0.005;        %%%%%%Height-interval between subplots

%% WPLI - Weight Phase Lag Index - Calculation *** Fig4b  - this takes a long time
WPLIStimGroupIndices = [1 3 6]; % 1=40, 3=Random, 6=LightRT
SubjGWPLI = SubjG(WPLIStimGroupIndices);  % Which three groups to include
GroupNameWPLI = GroupName(WPLIStimGroupIndices);
todayDate = datestr(now, 'yymmdd');
for iTrialType=3%1:length(TrialType)
    SaveTemp=[SubSaveWPLI todayDate '\' TrialTypeName{iTrialType} '\'];
    mkdir(SaveTemp)
    SubSaveFig=[SaveTemp 'Chan\'];
    mkdir(SubSaveFig)

    %     CH-Ch WPLI plot.tif figure;
    iPlot=0;
    alphaPeakAmplitudeList = zeros(length(EEGchInd),length(EEGchInd),length(SubjGWPLI));
    alphaPeakFrequencyList = zeros(length(EEGchInd),length(EEGchInd),length(SubjGWPLI));
    % alphaPeakAmplitudeListEmpty = double.empty(length(EEGchInd),length(EEGchInd),length(SubjG),0)
    for iCh=1:length(EEGchInd)
        for jCh=iCh+1:length(EEGchInd)
            clear DataPlot
            for iStimGroup=1:length(SubjGWPLI)
                DataPlot{iStimGroup}= squeeze(WPLIall(SubjGWPLI{iStimGroup},:,EEGchInd(iCh),EEGchInd(jCh),iTrialType));
                Invalid=isnan(DataPlot{iStimGroup}(:,1));
                DataPlot{iStimGroup}(Invalid,:)=[];
            end
            if isempty(DataPlot{1})||isempty(DataPlot{2})
                continue;
            end
            iPlot=iPlot+1;
            subplotLU(length(EEGchInd),length(EEGchInd),iCh,jCh,P);
            ParamWPLI.PathSave=[SaveTemp 'Light40Rand' EEGch{iCh} '-' EEGch{jCh}];

            % [~,COHComStatis{iCom,iCh,jCh}]=RateHist_GroupPlot(Fplot,DataPlot,FlickerColor,ParamCOH);
            tic
            RateHist_GroupPlot(Fplot,DataPlot,FlickerColor,ParamWPLI);
            % toc
            text(50,ParamWPLI.Ytick(end),[EEGch{iCh} '-' EEGch{jCh}]);
            set(gca,'xlim',FBand,'xtick',[0:20:120],'xticklabel',[],'yticklabel',[]);
        end
    end
    %% Permutation to set WPLI threshold ***fig4bc
    % Create "significant" WPLI threhold curve using permutation of current channel pair data.
    % -MKA 2024-12-11

    %Set up save folder
    saveDate = datestr(datetime, 'yy-mm-dd_HHMMSSFFF');
    SaveTemp=[SaveTrialSubj 'WPLIPermutation_' saveDate '\'];
    mkdir(SaveTemp);

    % Define parameters
    nPermutations = 10000; % Number of permutations, 10k or 1mil
    nSubjs = length(SubjID);  % total number of subjects in analysis
    % Define the frequency range of interest: 2-55Hz and every half frequency inbetween
    % frequencies = [2:0.5:55]; % start at 1 or 2 hz? Cut off before 60 Hz

    % Pre-allocate storage for maximum WPLI curves across frequencies
    % WPLIperms = zeros(nPermutations, length(frequencies));
    % WPLIallPerm=zeros(length(FileStruct),length(FreInd),ChNTotal,ChNTotal-1,length(TrialType))+nan;

    %% Preallocate CohGroupPerm as a cell array
    CohGroupPerm = cell(3, nPermutations);

    % Begin permutations
    tStartPerm = tic;
    for iPermutation = 1:nPermutations
        % Step 1: Randomly select two subjects and a channel pair
        randSubjs = randperm(nSubjs, 2); %randomly choose a random Subject A and Subject B from all three groups.

        SubAInd = randSubjs(1);
        SubAData = EEGList{1,SubAInd};

        SubBInd = randSubjs(2);
        SubBData = EEGList{1,SubBInd};

        randChs = randperm(32,2); %randomly choose a chan Chi from subject i,  Chj for subject j
        Chi = randChs(1);
        Chj = randChs(2);

        % Step 2: Find the lower number of trials between the two selected subjects
        nTrials_SubA = SubAData.trials;
        nTrials_SubB = SubBData.trials;
        minTrials = min(nTrials_SubA, nTrials_SubB);
        minTrialsIndex = 1:minTrials;

        % Step 3: Get `minTrials` from each subject
        clear TempTrial1 DataA TempTrial2 DataB TrialSpec
        DataA=squeeze(SubAData.AllChData(Chi,SampI,minTrialsIndex)); % analog to "tempSigA"
        DataB=squeeze(SubBData.AllChData(Chj,SampI,minTrialsIndex)); % analog to "tempSigB"

        %% Step 4: Calculate WPLI for selected trials
        % WPLIperms(iPermuation, :) = calculateWPLI(DataA, DataB);

        % Preallocate TempTrial1 and TempTrial2 as structure arrays "FOR SPEED"
        TempTrial1(minTrials).Data = []; % Preallocate Data field
        TempTrial1(minTrials).Time = []; % Preallocate Time field
        TempTrial2(minTrials).Data = []; % Preallocate Data field
        TempTrial2(minTrials).Time = []; % Preallocate Time field

        for iTrial=1:minTrials  % loop by trial
            if minTrials>1
                TempTrial1(iTrial).Data=DataA(:,iTrial);
                TempTrial2(iTrial).Data=DataB(:,iTrial);
            else
                TempTrial1(iTrial).Data=DataA;
                TempTrial2(iTrial).Data=DataB;
            end
            TempTrial1(iTrial).Time=([1:length(TempTrial1(iTrial).Data)]-1)/512;
            TempTrial2(iTrial).Time=([1:length(TempTrial2(iTrial).Data)]-1)/512;
        end

        [TrialSpec.Sxy,TrialSpec.Sxx,TrialSpec.Syy,TrialSpec.w,TrialSpec.options,ValidIndex]=crossspec_EqualTriL(TempTrial1,TempTrial2,psdParameter);  %% find in genmat code
        % save([SaveTemp 'Ch' num2str(Chi) 'Ch' num2str(Chj) '.mat'],'TrialSpec','ValidIndex','psdParameter');

        jTrialType = iTrialType;
        for jTrialType=3%1:length(TrialType)  % group: hit, miss, etc.
            %iTrialType hardcoded to 3; 3=Hit&Miss (all trials)
            TrialI=[];

            parfor j=1:length(TrialType{jTrialType})  % starting pool takes time; does this need to be parfor?
                TrialI=union(TrialI,find(TrialTypeTemp==TrialType{jTrialType}(j)));
            end

            TrialI=intersect(TrialI,ValidIndex);

            if ~isempty(TrialI)
                %  CohGroup{iGroup,iFile}{iCh,jCh}=Coh_TrialIndex(TrialSpec1,TrialI);
                CohGroupPerm{iPermutation}=crossspec_TrialIndex(TrialSpec,TrialI);
                % CohGroup -> {3x67 cell} -> {3 trialTypes x 67 subjects}
            end
        end
    end
    elapsedTime = toc(tStartPerm);
    disp(['Elapsed time for script: ', num2str(elapsedTime), ' seconds']);

    %Save Permuation variable
    filename = [nPermutations 'CohGroupPerm_', datestr(now, 'yyyy-mm-dd'), '.mat'];
    save(filename, 'CohGroupPerm');

    %% Get WPLI from coherence
    Fre=CohGroupPerm{1}.Fre;  % may need to import CohGroup from a COHdata file
    FBand=[2 100]; % consider changing [1 100] to [2 100] due to normalization
    FreInd=find(Fre>=FBand(1)&Fre<=FBand(2));
    Fplot=Fre(FreInd);

    PSDallPerm=zeros(length(FileStruct),length(FreInd),ChNTotal,length(TrialType))+nan;
    COHallPerm=zeros(length(FileStruct),length(FreInd),ChNTotal,ChNTotal-1,length(TrialType))+nan;
    WPLIallPerm=zeros(length(FileStruct),length(FreInd),ChNTotal,ChNTotal-1,length(TrialType))+nan;

    for iPermutation = 1:nPermutations
        if isempty(CohGroupPerm{iPermutation})
            continue;
        end
        if iCh==1
            PSDallPerm(iPermutation,:,iCh)=CohGroupPerm{iPermutation}.Pxx(FreInd);
            PSDallPerm(iPermutation,:,jCh)=CohGroupPerm{iPermutation}.Pyy(FreInd);
        end
        COHallPerm(iPermutation,:,iCh,jCh)=CohGroupPerm{iPermutation}.Cxy(FreInd);
        WPLIallPerm(iPermutation,:,iCh,jCh)=CohGroupPerm{iPermutation}.wpli(FreInd);
    end
end

%% Step 5: Determine the significance threshold for each frequency

% EXTRACT relevant data for wpli threshold
HitMissDataset = CohGroupPerm(3,:);
% Initialize an output cell array of the same size.
extractedStructs = cell(1, numel(HitMissDataset));
% Extract the (31, 32) cell data from each cell in CohGroupPerm.
extractedStructs = cellfun(@(x) x{31, 32}, HitMissDataset, 'UniformOutput', false);
% Convert the extracted 'wpli' data to a 10,000x513 double matrix.
wpliMatrix = cell2mat(cellfun(@(x) x.wpli, extractedStructs, 'UniformOutput', false)');
wpliFreqs = extractedStructs{1,1}.Fre;

% FIND WPLI threshold (value_of_topk)
% Sort each column of wpliMatrix in descending order
sortedMatrix = sort(wpliMatrix, 1, 'descend');

% Extract the top_kth greatest value from each column
p_value_threshold = 0.01; % p-val Significance threshold analogus to alpha value of 0.01
top_0_01 = ceil(nPermutations * p_value_threshold); % Top K for p-value cutoff
value_of_top0_01 = sortedMatrix(top_0_01, :);

p_value_threshold2 = 0.001; % p-val Significance threshold analogus to alpha value of 0.01
top_0_001 = ceil(nPermutations * p_value_threshold2); % Top K for p-value cutoff
value_of_top0_001 = sortedMatrix(top_0_001, :);

p_value_threshold3 = 0.0001; % p-val Significance threshold analogus to alpha value of 0.01
top_0_0001 = ceil(nPermutations * p_value_threshold3); % Top K for p-value cutoff
value_of_top0_0001 = sortedMatrix(top_0_0001, :);

%% Step 6: Plot permutation WPLI threshold curve
% Plot wpliFreqs (x-axis) against value100 (y-axis)
figure;
hold on;
upperFreqLimit = 55;  % Upper frequency limit of plot in (hz).  Otherwise plot will go to >250 Hz
upperFreqIndex = find(wpliFreqs >= upperFreqLimit, 1, 'first');
% Plot the first line
plot(wpliFreqs, value_of_top0_01, 'LineWidth', 2, 'DisplayName', ...
    ['p = ' num2str(p_value_threshold) ', Top K = ' num2str(top_0_01)]);

% Plot the second line
plot(wpliFreqs, value_of_top0_001, 'LineWidth', 2, 'DisplayName', ...
    ['p = ' num2str(p_value_threshold2) ', Top K = ' num2str(top_0_001)]);

% Plot the third line
plot(wpliFreqs, value_of_top0_0001, 'LineWidth', 2, 'DisplayName', ...
    ['p = ' num2str(p_value_threshold3) ', Top K = ' num2str(top_0_0001)]);

% Add a legend
legend('show', 'Location', 'best');

% Add axis labels and title
xlabel('Frequency (Hz)', 'FontSize', 12);
ylabel('WPLI', 'FontSize', 12);
title([num2str(size(wpliMatrix, 1)) ' Permutations. Threshold Curves for Different p-Values'], 'FontSize', 14);

% Set y-axis to start at zero
xlim([0 55])
% ylim([0, max(max(value_of_top0_01(1:upperFreqIndex)), max(value_of_top0_001(1:upperFreqIndex)), max(value_of_top0_0001(1:upperFreqIndex)))]);
ylim([0, max([max(value_of_top0_01(1:upperFreqIndex)), ...
    max(value_of_top0_001(1:upperFreqIndex)), ...
    max(value_of_top0_0001(1:upperFreqIndex))])]);


% plot(wpliFreqs(1:index), value_of_topk(1:index), 'LineWidth', 2);
% plot(wpliFreqs(1:index), value_of_topk2(1:index), 'LineWidth', 2);
%
% legend()
%
% % Add axis labels and title
% xlabel('Frequency (Hz)', 'FontSize', 12);
% ylabel('WPLI', 'FontSize', 12);
% title([num2str(nPermutations) ' permutations.  Top ' num2str(top_k) ' WPLI value. pvalue threshold of ' num2str(p_value_threshold)], 'FontSize', 14);
%
% % Set y-axis to start at zero
% ylim([0, max(value_of_topk(1:index))]);

grid on;  % improve readability of plot
hold off;
% wpli_threshold_curve = zeros(32, 32, length(frequencies));
% for freq_idx = 1:length(frequencies)
%     % Sort WPLI values for the current frequency across all permutations
%     sorted_values = sort(perm_idx_WPLICurve(:, freq_idx), 'descend');
%
%     % Find the WPLI value corresponding to the top `top_k` value
%     wpli_threshold_curve(iCh, jCh, freq_idx) = sorted_values(top_k);
% end

%% Real > Permutated WPLI
% get real WPLI - from
load('10kCohGroupPerm_2024-12-17.mat')
upperFreqLimit = 55;  % Upper frequency limit of plot in (hz).  Otherwise plot will go to >250 Hz
upperFreqIndex = find(wpliFreqs >= upperFreqLimit, 1, 'first');
lowerFreqLimit = 2;
lowerFreqIndex = find(wpliFreqs >= lowerFreqLimit, 1, 'first');

averagePermutatedWPLI_0to55_top0_01 = mean(value_of_top0_01(1:upperFreqIndex), 'omitnan');
averagePermutatedWPLI_0to55_top0_001 = mean(value_of_top0_001(1:upperFreqIndex), 'omitnan');
averagePermutatedWPLI_0to55_top0_0001 = mean(value_of_top0_0001(1:upperFreqIndex), 'omitnan');

averagePermutatedWPLI_2to55_top0_01 = mean(value_of_top0_01(lowerFreqIndex:upperFreqIndex), 'omitnan');
averagePermutatedWPLI_2to55_top0_001 = mean(value_of_top0_001(lowerFreqIndex:upperFreqIndex), 'omitnan');
averagePermutatedWPLI_2to55_top0_0001 = mean(value_of_top0_0001(lowerFreqIndex:upperFreqIndex), 'omitnan');

permutatedWPLI = value_of_top0_001;

%% Initialize the average WPLI storage
AvgRealWPLI = cell(length(EEGchInd), length(EEGchInd)); % Cell array for all channel pairs
% AvgRealWPLI is a 32x32 cell containing 3x197 doubles (average WPLI per stim group (3) per freq (197)

% Loop through all channel pairs 32x32
for iCh = 1:length(EEGchInd)
    for jCh = iCh+1:length(EEGchInd) % Avoid duplicates and diagonal (upper triangle)

        % Initialize a 2D array to store averages for this channel pair
        AvgRealWPLI{iCh, jCh} = zeros(length(SubjGWPLI), size(WPLIall, 2)); % Rows: StimGroups, Cols: Frequencies

        % Loop through stimulus groups
        for iStimGroup = 1:length(SubjGWPLI)
            % Extract WPLI data for this stimulus group and channel pair
            DataPlot = squeeze(WPLIall(SubjGWPLI{iStimGroup}, :, EEGchInd(iCh), EEGchInd(jCh), iTrialType));

            % Remove invalid rows containing NaN
            Invalid = isnan(DataPlot(:, 1));
            DataPlot(Invalid, :) = [];

            % Compute the average across all valid rows for this stim group
            if ~isempty(DataPlot) % Ensure there is data after removing NaN
                AvgRealWPLI{iCh, jCh}(iStimGroup, :) = mean(DataPlot, 1); % Average across rows
            end
        end
    end
end
% outcome: 
% AvgRealWPLI is a 32x32 cell containing 3x197 doubles (average WPLI per stim group (3) per freq (197)
%% Peak WPLI Extraction
%%% set up
BOIAlphas=[1 4 8 8 10 13 30 39.5 43;4 8 13 10 13 30 37 41.5 100];
BOI2HzAlphas=[2 4 8 8 10 13 30 39.5 43 55;4 8 13 10 13 30 37 41.5 55 100];
BOI2Hz=[2 4 8 13 30 39.5 43 55;4 8 13 30 37 41.5 55 100]; % No LOWER or UPPER ALPHA
BandName={'Delta','Theta','Alpha','LowAlpha','HighAlpha','Beta','Gamma-1','Gamma-E', 'Gamma-55','Gamma-2'};
BandHzNameAlphas={'1-4 Hz','4-8 Hz','8-13Hz','8-10Hz','10-13Hz','13-30 Hz','30-37 Hz','39-41 Hz','43-55 Hz', '55-100 Hz'};
BandHzName2HzAlphas={'2-4 Hz','4-8 Hz','8-13Hz','8-10Hz','10-13Hz','13-30 Hz','30-37 Hz','39-41 Hz','43-55 Hz', '55-100 Hz'};
BandHzName2HzAlphas={'2-4 Hz','4-8 Hz','8-13Hz','13-30 Hz','30-37 Hz','39-41 Hz','43-55 Hz', '55-100 Hz'};


lowerboundHz = 2;
upperboundHZ = 39; % can change to 55 Hz

%%%
% $% Filter frequency bands up to upperbound Hz
BOI_filtered = BOI2Hz(:, BOI2Hz(2, :) <= upperboundHZ);
BandName_filtered = BandName(BOI2Hz(2, :) <= upperboundHZ);
BandHzName_filtered = BandHzName2Hz(BOI2Hz(2, :) <= upperboundHZ);
% 
% Preallocate results (cell array for flexibility)
numBands = size(BOI_filtered, 2);
[numChannels, ~] = size(AvgRealWPLI);
peakWPLI = cell(numChannels, numChannels);

% %% TESTING ONLY - Create BOIsFreqIndices - comment when fixed
% % Example: Fre is the frequency vector (197x1 double), BOI_filtered contains the bands up to 55 Hz
% BOIsFreqIndices = cell(1, numBands);
% BOIsFreqValues = cell(1, numBands); % To store frequencies for each band
% 
% % Calculate the indices and corresponding frequencies for each band
% for b = 1:numBands
%     BOIsFreqIndices{b} = find(Fplot >= BOI_filtered(1, b) & Fplot <= BOI_filtered(2, b));
%     BOIsFreqValues{b} = Fplot(BOIsFreqIndices{b}); % Map indices to frequencies
% end
% 
% % Display the indices and their corresponding frequencies for each band
% for b = 1:numBands
%     fprintf('Band: %s (%s)\n', BandName_filtered{b}, BandHzName_filtered{b});
%     fprintf('Indices: %s\n', mat2str(BOIsFreqIndices{b}));
%     fprintf('Frequencies (Hz): %s\n', mat2str(BOIsFreqValues{b}));
%     fprintf('\n');
% end

%%% Loop over channel pairs
for iCh = 1:length(EEGchInd)
    for jCh = 1:length(EEGchInd)
        if isempty(AvgRealWPLI{iCh, jCh})
            continue; % Skip empty cells
        end
        data = AvgRealWPLI{iCh, jCh}; % data: 3x197 double (stimulation groups x frequencies)

        % Preallocate storage for this channel pair
        peakWPLI{iCh, jCh} = zeros(size(data, 1), numBands); % stim groups x frequency bands

        % Loop over bands of interest
        for iBand = 1:numBands
            % Find indices corresponding to the frequency band
            freqIndices = Fplot >= BOI_filtered(1, iBand) & Fplot <= BOI_filtered(2, iBand);

            % Get peak WPLI for each stimulation group
            for group = 1:size(data, 1)
                peakWPLI{iCh, jCh}(group, iBand) = max(data(group, freqIndices));
                % peakWPLI: 32x32 cell of 3x5 double (3 stim groups, 5 BOIs)
            end
        end
    end
end

% Result is stored in peakWPLI{i, j}(group, iBand), where:
% i, j = channel indices
% group = stimulation group
% iBand = frequency band index

% %% TESTING ONLY - Max WPLI Frequency per band
% % Example variables
% maxWPLIFreqHz = cell(length(EEGchInd), numChannels); % To store the frequencies of max WPLI
% 
% % Loop through all channel pairs
% for i = 1:numChannels
%     for j = 1:numChannels
%         if isempty(AvgRealWPLI{i, j})
%             continue; % Skip empty cells
%         end
% 
%         data = AvgRealWPLI{i, j}; % 3x197 double (stimulation groups x frequencies)
%         maxWPLIFreqHz{i, j} = zeros(size(data, 1), numBands); % Preallocate for groups x bands
% 
%         % Loop through stimulation groups
%         for group = 1:size(data, 1)
%             % Loop through each band
%             for iBand = 1:numBands
%                 freqIndices = BOIsFreqIndices{iBand}; % Get indices for the band
%                 [~, maxIndex] = max(data(group, freqIndices)); % Find index of max value
%                 maxWPLIFreqHz{i, j}(group, iBand) = Fplot(freqIndices(maxIndex)); % Map to frequency in Hz
%             end
%         end
%     end
% end
% 
% % Display example output for a specific channel pair
% exampleChannelPair = [1, 2]; % Change to any valid pair
% if ~isempty(maxWPLIFreqHz{exampleChannelPair(1), exampleChannelPair(2)})
%     fprintf('Max WPLI frequencies for channel pair (%d, %d):\n', exampleChannelPair(1), exampleChannelPair(2));
%     for group = 1:numGroups
%         fprintf('  %s:\n', GroupNameWPLI{group}); % Use the group name
%         for iBand = 1:numBands
%             fprintf('    Band: %s (%s) - Max WPLI at %.2f Hz\n', ...
%                 BandName_filtered{iBand}, BandHzName_filtered{iBand}, ...
%                 maxWPLIFreqHz{exampleChannelPair(1), exampleChannelPair(2)}(group, iBand));
%         end
%     end
% end

%%% Plot fraction of channel pairs with WPLI greater than p=0.01 permutated WPLI

% Initialize fraction of channel pairs exceeding the threshold
numGroups = 3; % Number of stimulation groups
numBands = size(BOI_filtered, 2); % Number of frequency bands
fractionExceed = zeros(numGroups, numBands);

averagePermutatedWPLI_2to55_top0_0001 = 0.0894;
% averagePermutatedWPLI_2to55_top0_001 = 0.04
% 0.1202; % top quartile threshold
averagePermutatedWPLIvalue = averagePermutatedWPLI_2to55_top0_0001; % top quartile threshold % averagePermutatedWPLI_2to55_top0_001; 

% Count total non-empty cells
numChannels = size(AvgRealWPLI, 1);
totalPairs = 0;

for i = 1:numChannels
    for j = 1:numChannels
        if ~isempty(AvgRealWPLI{i, j})
            totalPairs = totalPairs + 1;
        end
    end
end

% Loop through stimulation groups and frequency bands
for iGroup = 1:numGroups
    for iBand = 1:numBands
        exceedCount = 0;

        % Loop over all channel pairs
        for i = 1:numChannels
            for j = 1:numChannels
                if isempty(AvgRealWPLI{i, j})
                    continue; % Skip empty cells
                end

                % Check if the value for the group and band exceeds the threshold
                if peakWPLI{i, j}(iGroup, iBand) > averagePermutatedWPLIvalue %
                    exceedCount = exceedCount + 1;
                end
            end
        end

        % Calculate the fraction for this group and band
        fractionExceed(iGroup, iBand) = exceedCount / totalPairs;
    end
end
%% Peak Alpha WPLI 40Hz vs Light SIGNIFICANT CH ONLY ( stats preceding fig4D stats)
%%% Groups
GroupSubjs_40Hz = 1; % 40Hz group
GroupSubjs_Light = 6; % LightRT group

%%% Get all alpha values
alphaLowerLimitFreqHz = 8;
alphaUpperLimitFreqHz = 13;

% Find indices of values in Fre between 8 and 13 (inclusive)
allAlphaFreIndices = find(Fre >= alphaLowerLimitFreqHz & Fre <= alphaUpperLimitFreqHz);
HitMissTrialType = 3;
WPLIallAlpha = squeeze(WPLIall(:,allAlphaFreIndices,1:32,1:32,HitMissTrialType)); % size(WPLIall) ans = 67 197 40 39 3- 67subs x allFreqs x iCh x jCh x TrialType
%WPLIallAlpha: (67subs, x 11 freqs x 32 x 32chs)

%%% Find the peak WPLI value within alpha (dimension 2) and the corresponding Index (correspond to
% Freq Hz) of the peak alpha WPLI value
[peakAlphaValues4D, peakAlphaIndices4D] = max(WPLIallAlpha, [], 2);
peakAlphaValues = squeeze(peakAlphaValues4D);
peakAlphaIndices = squeeze(peakAlphaIndices4D);

%%%%%%%%%%%%%% pre-process peak alpha WLPI values
% peakAlphaValues = squeeze(peakAlphaValues);  % 4d to 3d

% test with single
peakAlphaValueSingleSub = squeeze(peakAlphaValues(1,:,:));

% Extract Peak Alpha WPLI values for each sub into groups: 40 Hz & Light
PeakAlphaWPLIs_40HzGroup = peakAlphaValues(SubjG{GroupSubjs_40Hz},:,:);
PeakAlphaWPLIs_LightGroup = peakAlphaValues(SubjG{GroupSubjs_Light},:,:);

% Get average peak alpha WPLI  for each group for all channel pairs
MeanAlphaPeakWPLI_40Hz = squeeze(mean(PeakAlphaWPLIs_40HzGroup, 1));
MeanAlphaPeakWPLI_Light = squeeze(mean(PeakAlphaWPLIs_LightGroup, 1));

%%% Get alpha peak WPLI frequency from peakWPLI_Freq (see: %% Plot bar graph of # channels exceeding
%%% p=0.0001 WPLI threshold (New Fig4C MS) MKA 2025-02-06)

% Preallocate with NaN
peakAlphaFreq_40HzGroup = NaN(numChannels, numChannels);
peakAlphaFreq_LightGroup = NaN(numChannels, numChannels);

% Apply cellfun with error handling
validCells = ~cellfun(@isempty, peakWPLI_Freq); % Logical mask for non-empty cells
peakAlphaFreq_40HzGroup(validCells) = cellfun(@(x) x(1,3), peakWPLI_Freq(validCells));
peakAlphaFreq_LightGroup(validCells) = cellfun(@(x) x(1,3), peakWPLI_Freq(validCells));



% %%%%%%%%%%%% Get frequency value of the peak alpha WPLI
% % Initialize a copy of peakAlphaIndices
% peakAlphaIndicesNaN = peakAlphaIndices;
% 
% % test with single participant
% singleSubPeakAlphaIndexAllChPairs = squeeze(peakAlphaIndices(1,:,:));  % should contain values between 1 and 11
% 
% % Convert Indices to Correct Corresponding Frequnecy (Hz)
% % Find valid indices (values between 1 and 11)
% validIdx = singleSubPeakAlphaIndexAllChPairs >= 1 & singleSubPeakAlphaIndexAllChPairs <= 11;
% 
% % Initialize the output array with NaN, preserving original shape
% singleallchpfreqs = NaN(size(singleSubPeakAlphaIndexAllChPairs));
% 
% % Perform mapping only for valid indices
% singleallchpfreqs(validIdx) = alphaFreqHz(singleSubPeakAlphaIndexAllChPairs(validIdx));
% 
% % If all alpha WPLI values are NaN, then set the max index to NaN
% peakAlphaIndicesNaN(all(isnan(WPLIallAlpha),2)) = NaN;
% 
% % % Reshape the result to 3D: 67x32x32 double
% % PeakAlphaIndicesNaN3D = squeeze(peakAlphaIndicesNaN);
% 
% % Convert Indices (3D) to Correct Corresponding Frequency (Hz) 
% validIdx = peakAlphaIndicesNaN >= 1 & peakAlphaIndicesNaN <= 11;
% 
% % Initialize the output array with NaN, preserving original shape
% PeakAlphaFreqs = NaN(size(peakAlphaIndicesNaN));
% 
% % Perform mapping only for valid indices
% PeakAlphaFreqs(validIdx) = alphaFreqHz(peakAlphaIndicesNaN(validIdx));
% 
% 
% % Extract Peak Alpha WPLI frequency into groups: 40 Hz & Light
% % PeakWPLIallAlpha: 67 subs x 32ch x 32ch
% PeakAlphaFreqs_40HzGroup = PeakAlphaFreqs(SubjG{GroupSubjs_40Hz},:,:);
% PeakAlphaFreqs_LightGroup = PeakAlphaFreqs(SubjG{GroupSubjs_Light},:,:);
% 
% % Get average peak alpha WPLI frequncy for each group for all channel pairs
% MeanPeakAlphaFreq_40Hz = squeeze(mean(PeakAlphaFreqs_40HzGroup, 1));
% MeanPeakAlphaFreq_Light = squeeze(mean(PeakAlphaFreqs_LightGroup, 1));

%%% Find channel pairs that exceed WPLI threshold
% WPLI threshold:
averagePermutatedWPLI_2to55_top0_0001 = 0.0894;
averagePermutatedWPLIvalue = averagePermutatedWPLI_2to55_top0_0001;

% display the averagePermutatedWPLIvalue
fprintf("averagePermutatedWPLI_2to55_top0_0001: %.4f\n", averagePermutatedWPLIvalue);


% create logical mask of all channel pairs that have a peak band WPLI that exceeds the WPLI threshold (e.g. 0.0894) for 40 Hz and Light
% Groups. peakWPLI: 32x32 cell of 3x5 doubles
% SignificantWPLI_40Hz = MeanAlphaPeakWPLI_40Hz > averagePermutatedWPLIvalue;
% SignificantWPLI_Light = MeanAlphaPeakWPLI_Light > averagePermutatedWPLIvalue;
SignificantWPLI_40Hz = exceedMask(:,:,1,3);
SignificantWPLI_Light = exceedMask(:,:,3,3);


% count the total number of channel pairs that exceed the WPLI threshold, display the result
numSignificantPairs_40Hz = sum(sum(SignificantWPLI_40Hz));
numSignificantPairs_Light = sum(sum(SignificantWPLI_Light));

fprintf('\nNumber of Significant Channel Pairs (exceeding WPLI threshold) in 40Hz Group: %d\n', numSignificantPairs_40Hz);
fprintf('Number of Significant Channel Pairs (exceeding WPLI threshold) in Light Group: %d\n', numSignificantPairs_Light);


%%% use logical mask to extract peak alpha frequencies of signficant WPLI (exceeding threshold) channel pairs for 40 Hz and Light from
% , display the result
% SignificantAlphaFreqs_40Hz = MeanPeakAlphaFreq_40Hz;
SignificantAlphaFreqs_40Hz = peakAlphaFreq_40HzGroup .* SignificantWPLI_40Hz;
SignificantAlphaFreqs_40Hz(~SignificantWPLI_40Hz) = NaN;

SignificantAlphaFreqs_Light = peakAlphaFreq_LightGroup .* SignificantWPLI_Light;
SignificantAlphaFreqs_Light(~SignificantWPLI_Light) = NaN;



% Compute average peak alpha frequency across all channel pairs that exceed WPLI thresold for each
% group, display the result
MeanSignificantAlphaFreq_40Hz = mean(SignificantAlphaFreqs_40Hz(SignificantAlphaFreqs_40Hz > 0));
MeanSignificantAlphaFreq_Light = mean(SignificantAlphaFreqs_Light(SignificantAlphaFreqs_Light > 0));

fprintf('\nMean Peak Alpha Frequency (Significant 40Hz Group): %.2f Hz\n', MeanSignificantAlphaFreq_40Hz);
fprintf('Mean Peak Alpha Frequency (Significant Light Group): %.2f Hz\n', MeanSignificantAlphaFreq_Light);

% Compute median values
MedianSignificantAlphaFreq_40Hz = median(SignificantAlphaFreqs_40Hz(SignificantAlphaFreqs_40Hz > 0));
MedianSignificantAlphaFreq_Light = median(SignificantAlphaFreqs_Light(SignificantAlphaFreqs_Light > 0));

fprintf('\nMedian Peak Alpha Frequency (Significant 40Hz Group): %.2f Hz\n', MedianSignificantAlphaFreq_40Hz);
fprintf('Median Peak Alpha Frequency (Significant Light Group): %.2f Hz\n', MedianSignificantAlphaFreq_Light);
% Compute 25th and 75th percentiles
Q1_40Hz = prctile(SignificantAlphaFreqs_40Hz(SignificantAlphaFreqs_40Hz > 0), 25);
Q3_40Hz = prctile(SignificantAlphaFreqs_40Hz(SignificantAlphaFreqs_40Hz > 0), 75);
Q1_Light = prctile(SignificantAlphaFreqs_Light(SignificantAlphaFreqs_Light > 0), 25);
Q3_Light = prctile(SignificantAlphaFreqs_Light(SignificantAlphaFreqs_Light > 0), 75);

fprintf('\n25th Percentile (Significant 40Hz Group): %.2f Hz\n', Q1_40Hz);
fprintf('75th Percentile (Significant 40Hz Group): %.2f Hz\n', Q3_40Hz);
fprintf('25th Percentile (Significant Light Group): %.2f Hz\n', Q1_Light);
fprintf('75th Percentile (Significant Light Group): %.2f Hz\n', Q3_Light);

% Compute STE (Standard Error of the Mean)
std_40Hz = nanstd(SignificantAlphaFreqs_40Hz(:)); % Standard deviation
std_Light = nanstd(SignificantAlphaFreqs_Light(:));

N_40Hz = sum(~isnan(SignificantAlphaFreqs_40Hz(:))); % Sample size
N_Light = sum(~isnan(SignificantAlphaFreqs_Light(:)));

STE_40Hz = std_40Hz / sqrt(N_40Hz);
STE_Light = std_Light / sqrt(N_Light);

fprintf('\nStandard Error of the Mean (STE) - 40Hz Group: %.4f Hz\n', STE_40Hz);
fprintf('Standard Error of the Mean (STE) - Light Group: %.4f Hz\n', STE_Light);

% Test the distribution of alpha peak frequencies for normality, display results
[h_40Hz, p_40Hz] = kstest(SignificantAlphaFreqs_40Hz(:));
[h_Light, p_Light] = kstest(SignificantAlphaFreqs_Light(:));
fprintf('\nKolmogorov–Smirnov Normality Test (40Hz Group): p = %.8f\n', p_40Hz);
fprintf('Kolmogorov–Smirnov Normality Test (Light Group): p = %.8f\n', p_Light);

% Create violin plot of 40 Hz and Light alpha peak frequency distributions
figure;
v = violinplot([SignificantAlphaFreqs_40Hz(:), SignificantAlphaFreqs_Light(:)], {'40Hz', 'Light'});
ylabel('Peak Alpha Frequency (Hz)');
title(sprintf('Violin Plot of Peak Alpha Frequency Distributions\nWPLI threshold: %.4f', averagePermutatedWPLIvalue));

% Compute means and number of non-NaN data points
mean_40Hz = nanmean(SignificantAlphaFreqs_40Hz(:));
mean_Light = nanmean(SignificantAlphaFreqs_Light(:));
N_40Hz = sum(~isnan(SignificantAlphaFreqs_40Hz(:)));
N_Light = sum(~isnan(SignificantAlphaFreqs_Light(:)));

% Annotate mean values on the plot
hold on;
plot(1, mean_40Hz, 'kd', 'MarkerFaceColor', 'k', 'MarkerSize', 8); % Mean for 40Hz
plot(2, mean_Light, 'kd', 'MarkerFaceColor', 'k', 'MarkerSize', 8); % Mean for Light

% Annotate median values on the plot
plot(1, MedianSignificantAlphaFreq_40Hz, 'bs', 'MarkerFaceColor', 'b', 'MarkerSize', 8); % Median for 40Hz
plot(2, MedianSignificantAlphaFreq_Light, 'bs', 'MarkerFaceColor', 'b', 'MarkerSize', 8); % Median for Light

% Annotate 25th and 75th percentiles on the plot
plot(1, Q1_40Hz, 'm^', 'MarkerFaceColor', 'm', 'MarkerSize', 6); % 25th Percentile 40Hz
plot(1, Q3_40Hz, 'm^', 'MarkerFaceColor', 'm', 'MarkerSize', 6); % 75th Percentile 40Hz
plot(2, Q1_Light, 'm^', 'MarkerFaceColor', 'm', 'MarkerSize', 6); % 25th Percentile Light
plot(2, Q3_Light, 'm^', 'MarkerFaceColor', 'm', 'MarkerSize', 6); % 75th Percentile Light

% Display text for mean, median, and percentiles
text(1, mean_40Hz + 0.2, sprintf('Mean: %.2f Hz', mean_40Hz), 'HorizontalAlignment', 'center');
text(2, mean_Light + 0.2, sprintf('Mean: %.2f Hz', mean_Light), 'HorizontalAlignment', 'center');
text(1, MedianSignificantAlphaFreq_40Hz - 0.2, sprintf('Median: %.2f Hz', MedianSignificantAlphaFreq_40Hz), 'HorizontalAlignment', 'center', 'Color', 'b');
text(2, MedianSignificantAlphaFreq_Light - 0.2, sprintf('Median: %.2f Hz', MedianSignificantAlphaFreq_Light), 'HorizontalAlignment', 'center', 'Color', 'b');
text(1, Q1_40Hz - 0.3, sprintf('Q1: %.2f Hz', Q1_40Hz), 'HorizontalAlignment', 'center', 'Color', 'm');
text(1, Q3_40Hz + 0.3, sprintf('Q3: %.2f Hz', Q3_40Hz), 'HorizontalAlignment', 'center', 'Color', 'm');
text(2, Q1_Light - 0.3, sprintf('Q1: %.2f Hz', Q1_Light), 'HorizontalAlignment', 'center', 'Color', 'm');
text(2, Q3_Light + 0.3, sprintf('Q3: %.2f Hz', Q3_Light), 'HorizontalAlignment', 'center', 'Color', 'm');

% Display sample sizes
text(1, min(ylim) + 1, sprintf('N = %d', N_40Hz), 'HorizontalAlignment', 'center');
text(2, min(ylim) + 1, sprintf('N = %d', N_Light), 'HorizontalAlignment', 'center');

hold off;

% Rank sum test: testing if Light has significantly greater alpha peak frequency than 40 Hz
[p_ranksum, h_ranksum] = ranksum(SignificantAlphaFreqs_40Hz(:), SignificantAlphaFreqs_Light(:), 'tail', 'left');

fprintf('\nRank Sum Test p-value: %.16f\n', p_ranksum);

% Display peak alpha frequency of channel pairs exceeding WPLI threshold in tables (32x32)
ChannelLabels = {'Fp1', 'AF3', 'F7', 'F3', 'FC1', 'FC5', 'T7', 'C3', 'CP1', 'CP5', ...
                 'P7', 'P3', 'Pz', 'PO3', 'O1', 'Oz', 'O2', 'PO4', 'P4', 'P8', ...
                 'CP6', 'CP2', 'C4', 'T8', 'FC6', 'FC2', 'F4', 'F8', 'AF4', 'Fp2', 'Fz', 'Cz'};

fprintf('\nMean Peak Alpha Frequency (Significant channel pairs only - in 40Hz Group):\n');
disp(array2table(SignificantAlphaFreqs_40Hz, 'VariableNames', ChannelLabels, 'RowNames', ChannelLabels));

fprintf('\nMean Peak Alpha Frequency (Significant channel pairs only - in Light Group):\n');
disp(array2table(SignificantAlphaFreqs_Light, 'VariableNames', ChannelLabels, 'RowNames', ChannelLabels));

fprintf('\nMean Peak Alpha WPLI (All channel pairs in 40Hz Group):\n');
disp(array2table(MeanAlphaPeakWPLI_40Hz, 'VariableNames', ChannelLabels, 'RowNames', ChannelLabels));

fprintf('\nMean Peak Alpha WPLI (All channel pairs in Light Group):\n');
disp(array2table(MeanAlphaPeakWPLI_Light, 'VariableNames', ChannelLabels, 'RowNames', ChannelLabels));

%% Peak Alpha WPLI 40Hz vs Light SIGNIFICANT CH ONLY (abandoned approached)
% Initialize logical masks for valid channel pairs
numChannels = 32;
validChannels_40Hz = false(numChannels, numChannels);
validChannels_Light = false(numChannels, numChannels);
group40HzIndex = 1;
groupLightIndex = 3;

% Identify valid channels where WPLI exceeds threshold in each group
for i = 1:numChannels
    for j = 1:numChannels
        if isempty(peakWPLI{i, j})
            continue; % Skip empty cells
        end

        % Check if peak WPLI is greater than threshold for 40Hz and Light groups
        if peakWPLI{i, j}(group40HzIndex, 1) > averagePermutatedWPLIvalue
            validChannels_40Hz(i, j) = true;
        end
        if peakWPLI{i, j}(groupLightIndex, 1) > averagePermutatedWPLIvalue
            validChannels_Light(i, j) = true;
        end
    end
end

% Find common channel pairs where both groups exceed the threshold
validChannels = validChannels_40Hz & validChannels_Light;

% Collect flattened peak alpha frequencies for valid channels only
selectedFreqs_40Hz = [];
selectedFreqs_Light = [];

for i = 1:numChannels
    for j = 1:numChannels
        if validChannels(i, j)
            selectedFreqs_40Hz = [selectedFreqs_40Hz; PeakAlphaFreqs_40HzGroup(:, i, j)];
            selectedFreqs_Light = [selectedFreqs_Light; PeakAlphaFreqs_LightGroup(:, i, j)];
        end
    end
end
%% Plot the results (no counts)
figure;
for iGroup = 1:numGroups
    subplot(1, numGroups, iGroup);
    bar(fractionExceed(iGroup, :));
    title(GroupNameWPLI{iGroup}); % Use group name as title
    xlabel('Frequency Band');
    ylabel('Fraction of Channel Pairs');

    % Customize x-axis labels with both BandName and BandHzName
    xticks(1:numBands);
    xticklabels(arrayfun(@(iBand) [BandName_filtered{iBand}, ' (', BandHzName_filtered{iBand}, ')'], ...
        1:numBands, 'UniformOutput', false));
    xtickangle(45);

    ylim([0 1]); % Fractions are between 0 and 1
end

% Add super title
sgtitle('Fraction of Channel Pairs Exceeding average permutated p=0.01 WPLI value');

%% Same as above but with counts above each bar
figure;
for iGroup = 1:numGroups
    subplot(1, numGroups, iGroup);
    barHandle = bar(fractionExceed(iGroup, :));
    title(GroupNameWPLI{iGroup}); % Use group name as title
    xlabel('Frequency Band');
    ylabel('Fraction of Channel Pairs');

    % Customize x-axis labels with both BandName and BandHzName
    xticks(1:numBands);
    xticklabels(arrayfun(@(iBand) [BandName_filtered{iBand}, ' (', BandHzName_filtered{iBand}, ')'], ...
        1:numBands, 'UniformOutput', false));
    xtickangle(45);

    ylim([0 1]); % Fractions are between 0 and 1

    %%% Add exceed count above each bar
    barHeights = fractionExceed(group, :);
    for iBand = 1:numBands
        % Position the text slightly above the bar height
        text(iBand, barHeights(iBand) + 0.02, num2str(round(barHeights(iBand) * totalPairs)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
end

% Add super title
sgtitle('Fraction of Channel Pairs Exceeding average permutated p=0.0001 WPLI value');

%% Bar graph with counts below
figure;
for group = 1:numGroups
    subplot(1, numGroups, group);
    barHandle = bar(fractionExceed(group, :));
    title(GroupNameWPLI{group}); % Use group name as title
    xlabel('Frequency Band');
    ylabel('Fraction of Channel Pairs');

    % Customize x-axis labels with both BandName and BandHzName
    xticks(1:numBands);
    xticklabels(arrayfun(@(iBand) [BandName_filtered{iBand}, ' (', BandHzName_filtered{iBand}, ')'], ...
        1:numBands, 'UniformOutput', false));
    xtickangle(45);

    ylim([0 1]); % Fractions are between 0 and 1

    % Add exceed count below the top of each bar
    barHeights = fractionExceed(group, :);
    for iBand = 1:numBands
        % Position the text slightly below the bar height
        text(iBand, barHeights(iBand) - 0.02, num2str(round(barHeights(iBand) * totalPairs)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10, 'VerticalAlignment', 'top');
    end
end

% Add super title
sgtitle('Fraction of Channel Pairs Exceeding average permutated p=0.01 WPLI value');

%% Plot the results 2025-01-08 MKA fig4D supp? - top quartile
% BOI=[1 4 8 13 30 39.5;4 8 13 30 37 41.5];
% BOI2HzFiveBands=[2 4 8 13 39.5;4 8 13 30 41.5];
% BandNameFiveBands={'Delta','Theta','Alpha','Beta','Gamma-E'};
% BandHzName={'1-4 Hz','4-8 Hz','8-13Hz','13-30 Hz','39-41 Hz'};
% BandHzName2HzFiveBands={'2-4 Hz','4-8 Hz','8-13Hz','13-30 Hz','39-41 Hz'};

BOI2Hz=[2 4 8 8 10 13 30 39.5 43 55;4 8 13 10 13 30 37 41.5 55 100];
BandName={'Delta','Theta','Alpha','LowAlpha','HighAlpha','Beta','Gamma-1','Gamma-E', 'Gamma-55','Gamma-2'};
BandHzName2Hz={'2-4 Hz','4-8 Hz','8-13Hz','8-10Hz','10-13Hz','13-30 Hz','30-37 Hz','39-41 Hz','43-55 Hz', '55-100 Hz'};

numBands = size(BOI2Hz, 2);

%%% Get WPLI for each stim group for each BOI and ch Pair
%%%% Loop thru all ch pairs
for i = 1:numChannels
    for j = 1:numChannels
        if isempty(AvgRealWPLI{i, j})
            continue; % Skip empty cells
        end
        data = AvgRealWPLI{i, j}; % 3x197 double (stimulation groups x frequencies)

        % Preallocate storage for this channel pair
        peakWPLI{i, j} = zeros(size(data, 1), numBands); % stim groups x frequency bands

        % Loop over bands of interest
        for iBand = 1:numBands
            % Find indices corresponding to the frequency band
            freqIndices = Fplot >= BOI2Hz(1, iBand) & Fplot <= BOI2Hz(2, iBand);

            % Get peak WPLI for each stimulation group
            for group = 1:size(data, 1)
                peakWPLI{i, j}(group, iBand) = max(data(group, freqIndices));
            end
        end
    end
end

% Result is stored in peakWPLI{i, j}(group, iBand), where:
% i, j = channel indices
% group = stimulation group
% iBand = frequency band index

%%% Plot fraction of channel pairs with WPLI greater than p=0.01 permutated WPLI
% Example variables (replace these with actual data)
% averagePermutatedWPLItop0_01 = 0.5; % Replace with actual value
% GroupNameWPLI = {'Group 1', 'Group 2', 'Group 3'}; % Replace with actual group names

% Initialize fraction of channel pairs exceeding the threshold
numGroups = 3; % Number of stimulation groups
numBands = size(BOI2Hz, 2)-1; % Number of frequency bands
fractionExceed = zeros(numGroups, numBands);

% averagePermutatedWPLI_2to55_top0_0001 = 0.0894
% averagePermutatedWPLI_2to55_top0_001 = 0.04
averagePermutatedWPLIvalue =  0.1202; % top quartile threshold %

% Count total non-empty cells
numChannels = size(AvgRealWPLI, 1);
totalPairs = 0;

for i = 1:numChannels
    for j = 1:numChannels
        if ~isempty(AvgRealWPLI{i, j})
            totalPairs = totalPairs + 1;
        end
    end
end

% Loop through stimulation groups and frequency bands
for group = 1:numGroups
    for iBand = 1:numBands
        exceedCount = 0;

        % Loop over all channel pairs
        for i = 1:numChannels
            for j = 1:numChannels
                if isempty(AvgRealWPLI{i, j})
                    continue; % Skip empty cells
                end

                % Check if the value for the group and band exceeds the threshold
                if peakWPLI{i, j}(group, iBand) > averagePermutatedWPLIvalue %averagePermutatedWPLItop0_01
                    exceedCount = exceedCount + 1;
                end
            end
        end

        % Calculate the fraction for this group and band
        fractionExceed(group, iBand) = exceedCount / totalPairs;
    end
end
figure;
for iGroup = 1:numGroups
    subplot(1, numGroups, iGroup);
    barHandle = bar(fractionExceed(iGroup, :));
    title(GroupNameWPLI{iGroup}); % Use iGroup name as title
    xlabel('Frequency Band');
    ylabel('Fraction of Channel Pairs');

    % Customize x-axis labels with both BandName and BandHzName
    xticks(1:numBands);
    xticklabels(arrayfun(@(iBand) [BandName{iBand}, ' (', BandHzName2Hz{iBand}, ')'], ...
        1:numBands, 'UniformOutput', false));
    xtickangle(45);

    ylim([0 1]); % Fractions are between 0 and 1

    % Add exceed count conditionally above or below the bar
    barHeights = fractionExceed(iGroup, :);
    for iBand = 1:numBands
        exceedCount = round(barHeights(iBand) * totalPairs); % Calculate exceed count
        if barHeights(iBand) < 0.2
            % Place count above the bar
            text(iBand, barHeights(iBand) + 0.04, num2str(exceedCount), ...
                'HorizontalAlignment', 'center', 'FontSize', 10);
        else
            % Place count below the top of the bar
            text(iBand, barHeights(iBand) - 0.01, num2str(exceedCount), ...
                'HorizontalAlignment', 'center', 'FontSize', 10, 'VerticalAlignment', 'top');
        end
    end

    % Store percentage of Alpha band exceedance for this iGroup
    alphaExceedPercent(iGroup) = barHeights(3) * 100; % Alpha band is iBand = 4
end

% Add super title
sgtitle(['Fraction of Channel Pairs Exceeding top quartile WPLI value: ' num2str(averagePermutatedWPLIvalue)]);

% Calculate and display the average percentage of Alpha band exceedance
averageAlphaPercent = mean(alphaExceedPercent);
disp(['Average Percentage of Channel Pairs Exceeding Threshold for Alpha Band: ', num2str(averageAlphaPercent), '%']);

% Display percentages for each group
for group = 1:numGroups
    disp(['Group ', GroupNameWPLI{group}, ': ', num2str(alphaExceedPercent(group)), '%']);
end
%% Plot bar graph of # channels exceeding p=0.0001 WPLI threshold (New Fig4C MS) MKA 2025-02-06
BOI2HzFiveBands=[2 4 8 13 30; 4 8 13 30 37];
BandNameFiveBands={'Delta','Theta','Alpha','Beta','Slow Gamma'};
BandHzName2HzFiveBands={'2-4 Hz','4-8 Hz','8-13Hz','13-30 Hz','30-37 Hz'};
numBands = size(BOI2HzFiveBands, 2);

%%% Initialize storage for peak frequencies
peakWPLI_Freq = cell(numChannels, numChannels);

%%% Loop over channel pairs
for i = 1:numChannels
    for j = 1:numChannels
        if isempty(AvgRealWPLI{i, j})
            continue; % Skip empty cells
        end
        AvgRealWPLI_ijChPair = AvgRealWPLI{i, j}; % AvgRealWPLI_ijChPair = 3x197 double (stimulation groups x frequencies)

        % Preallocate storage for this channel pair
        peakWPLI{i, j} = zeros(size(AvgRealWPLI_ijChPair, 1), numBands); % stim groups x frequency bands
        peakWPLI_Freq{i, j} = zeros(size(AvgRealWPLI_ijChPair, 1), numBands); % Store peak frequency
        

        % Loop over bands of interest
        for iBand = 1:numBands
            % Find indices corresponding to the frequency band
            iBandFreqIndices = Fplot >= BOI2HzFiveBands(1, iBand) & Fplot <= BOI2HzFiveBands(2, iBand);
            Fplot_iBand = Fplot(iBandFreqIndices); % Extract frequencies in this band

            % Get peak WPLI for each stimulation group
            for iGroup = 1:size(AvgRealWPLI_ijChPair, 1)
                AvgRealWPLI_ijChPair_iBand_iGroup = AvgRealWPLI_ijChPair(iGroup, iBandFreqIndices); % 

                % Find peak WPLI value
                [peakWPLI{i, j}(iGroup, iBand), maxIdx] = max(AvgRealWPLI_ijChPair_iBand_iGroup);

                % Store corresponding frequency
                peakWPLI_Freq{i, j}(iGroup, iBand) = Fplot_iBand(maxIdx);
            end
        end
    end
end

% Result is stored in peakWPLI{i, j}(group, iBand), where:
% i, j = channel indices
% group = stimulation group
% iBand = frequency band index

%%% Plot fraction of channel pairs with WPLI greater than p=0.01 permutated WPLI
% Example variables (replace these with actual data)
% averagePermutatedWPLItop0_01 = 0.5; % Replace with actual value
% GroupNameWPLI = {'Group 1', 'Group 2', 'Group 3'}; % Replace with actual group names

% Initialize fraction of channel pairs exceeding the threshold
numGroups = 3; % Number of stimulation groups
numBands = size(BOI2HzFiveBands, 2); % Number of frequency bands
fractionExceed = zeros(numGroups, numBands);

% averagePermutatedWPLI_2to55_top0_0001 = 0.0894
% averagePermutatedWPLI_2to55_top0_001 = 0.04
averagePermutatedWPLIvalue = averagePermutatedWPLI_2to55_top0_0001; % top quartile threshold % averagePermutatedWPLI_2to55_top0_001;

% Count total non-empty cells
numChannels = size(AvgRealWPLI, 1);
totalPairs = 0;

%%% count the number of total channel pairs based on number of AvgRealWPLI values (use mask)
for i = 1:numChannels
    for j = 1:numChannels
        if ~isempty(AvgRealWPLI{i, j})
            totalPairs = totalPairs + 1;
        end
    end
end

totalPairsMaskSum = sum(~cellfun(@isempty, AvgRealWPLI), 'all');


% Preallocate arrays
peakWPLIarray = nan(numChannels, numChannels, numGroups, numBands);
exceedMask = false(numChannels, numChannels, numGroups, numBands);
fractionExceed = zeros(numGroups, numBands);

% Create a logical mask for non-empty channel pairs
validPairs = ~cellfun(@isempty, AvgRealWPLI);


% Loop through stimulation groups and frequency bands
for iGroup = 1:numGroups
    for iBand = 1:numBands
        % % Extract peakWPLI values into an array (set NaN for empty pairs)
        % peakWPLIarray(validPairs, iGroup, iBand) = cellfun(@(x) x(iGroup, iBand), peakWPLI(validPairs));
        % peakWPLIarray(~validPairs) = NaN; % Ensure empty cells remain NaN

        % Extract peakWPLI values into an array (set NaN for empty pairs)
        tempValues = nan(numChannels, numChannels); % Temporary storage
        tempValues(validPairs) = cellfun(@(x) x(iGroup, iBand), peakWPLI(validPairs));

        % Store in preallocated array
        peakWPLIarray(:, :, iGroup, iBand) = tempValues;

        % Create a logical mask for values exceeding the threshold
        exceedMask(:, :, iGroup, iBand) = peakWPLIarray(:, :, iGroup, iBand) > averagePermutatedWPLIvalue;

        % Count the number of exceeding pairs
        exceedCount = sum(exceedMask(:, :, iGroup, iBand), 'all');

        % Calculate the fraction
        fractionExceed(iGroup, iBand) = exceedCount / totalPairs;
    end
end

% Loop through stimulation groups and frequency bands
for group = 1:numGroups
    for iBand = 1:numBands
        exceedCount = 0;

        % Loop over all channel pairs
        for i = 1:numChannels
            for j = 1:numChannels
                if isempty(AvgRealWPLI{i, j})
                    continue; % Skip empty cells
                end

                % Check if the value for the group and band exceeds the threshold
                if peakWPLI{i, j}(group, iBand) > averagePermutatedWPLIvalue %averagePermutatedWPLItop0_01
                    exceedCount = exceedCount + 1;
                end
            end
        end

        % Calculate the fraction for this group and band
        fractionExceed(group, iBand) = exceedCount / totalPairs;
    end
end

%%% Plot High Functional Connectivity Bar Plot *** fig4c
figure;
for iGroup = 1:numGroups
    subplot(1, numGroups, iGroup);
    barHandle = bar(fractionExceed(iGroup, :));
    if iGroup ~= 3
        title(GroupNameWPLI{iGroup}); % Use iGroup name as title
    else
        title("Light"); % Use "Light" instead of "LightRT"
    end
    if iGroup ==2
        xlabel('Frequency Band');
    end
    if iGroup == 1
        ylabel('Fraction of Channel Pairs');
    end
    
    % Customize x-axis labels with both BandName and BandHzName
    iBand = 1;
    xticks(1:numBands);
    xticklabels(arrayfun(@(iBand) [BandNameFiveBands{iBand}], ...
        1:numBands, 'UniformOutput', false));
    xtickangle(45);

    ylim([0 0.55]); % Fractions are between 0 and 1

    % Add exceed count conditionally above or below the bar
    barHeights = fractionExceed(iGroup, :);
    for iBand = 1:numBands
        exceedCount = round(barHeights(iBand) * totalPairs); % Calculate exceed count
        if barHeights(iBand) < 0.2
            % Place count above the bar
            text(iBand, barHeights(iBand) + 0.06500000000000000000000001, num2str(exceedCount), ...
                'HorizontalAlignment', 'center', 'FontSize', 10, 'VerticalAlignment', 'top'); %
        else
            % Place count below the top of the bar
            text(iBand, barHeights(iBand) + 0.065, num2str(exceedCount), ...
                'HorizontalAlignment', 'center', 'FontSize', 10, 'VerticalAlignment', 'top');
        end
    end

    % Store percentage of Alpha band exceedance for this iGroup
    alphaExceedPercent(iGroup) = barHeights(3) * 100; % Alpha band is iBand = 4
end

% Add super title
sgtitle(['Fraction of Channel Pairs Out of 492 Exceeding permutated p=0.0001 WPLI value: ' num2str(averagePermutatedWPLIvalue)], 'FontSize', 6);

% Calculate and display the average percentage of Alpha band exceedance
averageAlphaPercent = mean(alphaExceedPercent);
disp(['Average Percentage of Channel Pairs Exceeding Threshold for Alpha Band: ', num2str(averageAlphaPercent), '%']);

% Display percentages for each group
for group = 1:numGroups
    disp(['Group ', GroupNameWPLI{group}, ': ', num2str(alphaExceedPercent(group)), '%']);
end

% Set figure size
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 4 2.5]); % 6x4 inches

% Save as SVG
saveDateFig4C = datestr(datetime, 'yy-mm-dd_HHMMSSFFF');
saveas(gcf, ['Fig4C_WPLI_exceed' saveDateFig4C '.svg']);
saveas(gcf, ['Fig4C_WPLI_exceed' saveDateFig4C '.png']);
%% Get Alpha Peak Amplitude and Alpha Peak Frequency *** fig4b
iWPLICh = 1; % 1 = FP1
jWPLICh = 5; % 5 = FC1

for iStimGroup = 1:length(SubjGWPLI)
    groupMean = mean(DataPlot{iStimGroup});
    alphaMean = groupMean(13:23);
    [alphaPeakAmplitude, alphaPeakIndex] = max(groupMean(13:23));
    alphaFreqs = Fplot(13:23); % 13 to 23 should be 8 to 13 Hz
    alphaPeakFrequency = alphaFreqs(alphaPeakIndex);
    alphaPeakAmplitudeList(iWPLICh,jWPLICh,iStimGroup) = alphaPeakAmplitude;
    alphaPeakFrequencyList(iWPLICh,jWPLICh,iStimGroup) = alphaPeakFrequency;
end

% Create the actual plot fig4b (commented for getting alpha measurements) UNCOMMENT TO PLOT!
figure;
xrange = [0 37]; % x axis range (frequency)
subplot('Position',[0.1 0.1 0.88 0.88])
% Cut data down to 55Hz
% Fplot55=Fplot(1:107);
% DataPlot55 = cellfun(@(x) x(1:107), DataPlot, 'UniformOutput', false); % Works like this-> DataPlot55=DataPlot(1:107);
RateHist_GroupPlot(Fplot,DataPlot,FlickerColor,ParamWPLI);  %% Lu's function

% Add a dotted line at averagePermutatedWPLIvalue
yline(averagePermutatedWPLIvalue, '--', 'Color', 'k', 'LineWidth', 1.2);

text(range(xrange)/2-5,ParamWPLI.Ytick(end),[EEGch{iWPLICh} '-' EEGch{jWPLICh}],'FontSize', 10);
% set(gca,'xlim',[0 100],'xtick',[0:20:120],'ylim',[0 0.3],'ytick',ParamWPLI.Ytick);
set(gca,'xlim',xrange,'xtick',[1 4 8 13 30 37 40 50 60 80 100],'ylim',[0 0.15],'ytick',[0 0.05 0.1 0.15]);
xlabel('Frequency Hz')
% xlim(xrange);

ylabel('WPLI')
title(['EEG Channel Pair: ' EEGch{iWPLICh} '-' EEGch{jWPLICh}]);
ax=gca;
ax.XGrid = 'on';
ax.YGrid = 'off';
ax.FontSize = 7;  % Set the desired font size for tick marks
LuFontStandard
papersizePX=[0 0 8 8];
papersizePX=1.3*[0 0 5 3.2]; % 09/09/24
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
saveas(gcf,[SubSaveFig '40LightRandBands_' EEGch{iWPLICh} '-' EEGch{jWPLICh}],'svg');
% saveas(gcf,[SubSaveFig '40LightRandBands_' EEGch{iWPLICh} '-' EEGch{jWPLICh}],'tiff');

close all

nonzerosGroup1= nonzeros(alphaPeakAmplitudeList(:,:,1));
nonzerosGroup2= nonzeros(alphaPeakAmplitudeList(:,:,2));
nonzerosGroup3= nonzeros(alphaPeakAmplitudeList(:,:,3));
AlphaPeakSaveFilename = fullfile(SaveFolder, 'AlphaPeakAmpAndFreq') %#ok<NOPTS>
save(AlphaPeakSaveFilename,"alphaPeakAmplitudeList","alphaPeakFrequencyList");
%     papersizePX=[0 0 6*length(EEGchInd) 6*length(EEGchInd)];
%     set(gcf, 'PaperUnits', 'centimeters');
%     set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
%     saveas(gcf,[SaveTemp num2str(FBand(1)) '-' num2str(FBand(2)) 'HzAllCh40Random' TrialTypeName{iCom}],'pdf');
%     saveas(gcf,[SaveTemp num2str(FBand(1)) '-' num2str(FBand(2)) 'HzAllCh40Random' TrialTypeName{iCom}],'png');
%     saveas(gcf,[SaveTemp num2str(FBand(1)) '-' num2str(FBand(2)) 'HzAllCh40Random' TrialTypeName{iCom} '.eps'],'epsc');

close all

close all

%% Frequency Band Definition, for maps - May need to start running from here for Fig4
BOI=[1 4 8 13 30 39.5 43;4 8 13 30 37 41.5 100];
BandName={'Delta','Theta','Alpha','Beta','Gamma-1','Gamma-E','Gamma-2'};
% Gamma-E is  40
BandHzName={'1-4 Hz','4-8 Hz','8-13Hz','13-30 Hz','30-37 Hz','39-41 Hz','43-100 Hz'};

FreqFunc{1}=@nanmean;
FreqFunc{2}=@nanmedian;
FreqFunc{3}=@nanmax;
FreqFuncNames={'mean','median','peak'}; % Check three different spots in band - previously called 'FunGroupName'

%comparison between groups w/ freq data
diffTPmap=zeros(length(ChanEEGLab),length(ChanEEGLab),size(BOI,2),length(TrialType)); % for T test
diffRPmap=diffTPmap;
diffTmap=diffTPmap;

rSpear=diffTPmap;
pSpear=diffTPmap;
% % TNodeTh=10;

% FC_BrainEEGLu(ChanPos,AdjWeight,NodeWeight,Param)

%% Preallocate ChannelPairName Cell Array
signifChPairNameTopo = cell(3,3,3,3,2);
signifChPairNameTopo{3,3,3,3,2} = [];
% WPLIBResults = struct();
%% FC parameters for plotting; (may contain p-value variable)
clear FCpara
FCpara.ColorMap=colorMapPN; %%%Color map for correlation link
%% Create orange indigo colormap
% Number of colors in the colormap
n = 64;

% Define orange and indigo RGB values
orange = [1, 0.75, 0];
indigo = [0.2, 0.1, 1];

% Create a colormap by interpolating between orange and indigo
custom_cmap = [linspace(indigo(1), orange(1), n)', ...
    linspace(indigo(2), orange(2),  n)', ...
    linspace(indigo(3), orange(3), n)'];

% Apply the custom colormap
% colormap(custom_cmap);

FCpara.ColorMap=custom_cmap; %%%Color map for correlation link

%% Other parameters
FCpara.NodeColor=[0.8 0.8 0.8]; %%%Node Color of Nodes for FC, not important, it is actually defined in ChanPos
FCpara.Clim=[-1 1];    %%%Color Limit FC
% FCpara.MarkerSize=8; %%%MarkerSize of scatter
% FCpara.EdgeColor=[1 0 0]; %%% This is not needed as ColorMap field and Clim field would determine the edge color
FCpara.EdgeTh=0.1; %%%
FCpara.NodeTh=0.05; %%%

FCparaT=FCpara;  %%%%T test parameters
FCparaT.EdgeTh=3; %%%
FCparaT.NodeTh=0.001; %%%

FCspear=FCpara;  %%%%Spearman r parameters
FCspear.EdgeTh=0.1; %%%  This is the plotting threshold
FCspear.NodeTh=0.001; %%%

WPLIEdgeTh=0.6;
WPLINodeTh=0.05;
TEdgeTh=1;
TNodeTh=0.001;
pTEdgeTh=0.05;
pSpearEdgeTh=0.05; %8/8/24 0.05 to 0.1
rSpearEdgeTh=0.1;
rSpearNodeTh=0.001;

ScaleWPLI=0.5;
ScaleT=1;
ScaleSpear=0.2;

% BOI=[5;15];
BOI=[1 4 8 13 30 39.5 43;4 8 13 30 37 41.5 100];
BandName={'Delta','Theta','Alpha','Beta','Gamma-1','Gamma-E','Gamma-2'};
BandHzName={'1-4 Hz','4-8 Hz','8-13Hz','13-30 Hz','30-37 Hz','39-41 Hz','43-100 Hz'};

%% FC calculation - WPLI-Behavior - Spearman Topo Plots *** (old fig4c) fig4d fig4e fig4f
FCpara.EdgeTh=0.1; %%%
BOI=[2 4 8 8 10 13 30 39 43;4 8 13 10 13 30 37 41 100]; %[1 4 8 8 10 13 30 39.5 43;4 8 13 10 13 30 37 41.5 100];
BandName={'Delta','Theta','Alpha','LowAlpha','HighAlpha','Beta','Gamma-1','Gamma-E','Gamma-2'};
BandHzName={'2-4 Hz','4-8 Hz','8-13Hz','8-10Hz','10-13Hz','13-30 Hz','30-37 Hz','39-41 Hz','43-100 Hz'};
for iFreqFunc=[3] %1:length(FreqFuncNames) % mean = 1, median =2, peak = 3
    todayDate = datestr(now, 'yymmdd');
    for iTrialType=3 %1:length(TrialType) % 1=Hit, 2=Miss, 3=HitANDMiss
        %% Create save folder
        SaveTemp=[SubSaveWPLI TrialTypeName{iTrialType} '\'];
        SaveTemp=[SaveTemp todayDate '_' FreqFuncNames{iFreqFunc} '_p' num2str(pSpearEdgeTh) '_colored\' ];
        mkdir(SaveTemp)

        %% WPLI figure -  ALL groups (6 groups) & ALL BOI - *** (old fig4c) fig4d Complete version - see MS version below
        % Fig 4 C and D are composed of panels created by this section, cut and pasted together in
        % illustrator

        % FCpara.EdgeTh=0.06;
        % FCpara.EdgeTh=0.15; % previous arbitray threshold
        FCpara.EdgeTh=0.1202; % top quartile threshold
        % FCpara.EdgeTh=averagePermutatedWPLI_2to55_top0_0001; % = 0.0894 - top 0.0001 permutation threshold
        % FCpara.EdgeTh=averagePermutatedWPLI_2to55_top0_001; % = 0.0484 -  top 0.001 permutation threshold

        FCgroups = [1,3,6]; % [1,3,6] = [40hz , Random, LightRT]

        figure;
        nGroups = length(FCgroups); % (1:3) Just first 3 groups
        for iBOI=1:size(BOI,2)-1
            NeedI=find(Fplot>=BOI(1,iBOI)&Fplot<BOI(2,iBOI));
            for iStimGroup=1:length(FCgroups)
                if iFreqFunc==3
                    Tdata{iStimGroup,iBOI,iTrialType}=squeeze(FreqFunc{iFreqFunc}(WPLIall(SubjG{FCgroups(iStimGroup)},NeedI,EEGchInd,EEGchInd,iTrialType),[],2));
                    MapGroup{iStimGroup,iBOI,iTrialType}=squeeze(nanmean(FreqFunc{iFreqFunc}(WPLIall(SubjG{FCgroups(iStimGroup)},NeedI,EEGchInd,EEGchInd,iTrialType),[],2),1));
                else
                    Tdata{iStimGroup,iBOI,iTrialType}=squeeze(FreqFunc{iFreqFunc}(WPLIall(SubjG{FCgroups(iStimGroup)},NeedI,EEGchInd,EEGchInd,iTrialType),2));
                    MapGroup{iStimGroup,iBOI,iTrialType}=squeeze(nanmean(FreqFunc{iFreqFunc}(WPLIall(SubjG{FCgroups(iStimGroup)},NeedI,EEGchInd,EEGchInd,iTrialType),2),1));
                end

                EdgeColor=[0.8 0.8 0.8];
                %              FC_BrainEEGLu(ChanPosColin27,MapGroup{iG,iFF,iCom},[],EdgeTh,NodeTh,EdgeColor,[])
                %              axis off

                subplotLU(nGroups,size(BOI,2),iStimGroup,iBOI);

                % WPLIsForPlot(iStimGroup,iBOI,iTrialType) = MapGroup{iStimGroup,iBOI,iTrialType};
                % nChPairAboveThreshold(iStimGroup,iBOI,iTrialType) = sum(WPLIsForPlot>FCpara.EdgeTh);

                % Calculate the number of channel pairs with WPLI > EdgeTh
                currentMap = MapGroup{iStimGroup, iBOI, iTrialType};
                nChPairAboveThreshold = sum(currentMap(:) > FCpara.EdgeTh);
                AboveThreholdChannelPairAllGroups(iStimGroup,iBOI,iTrialType) = nChPairAboveThreshold;

                %%% The plotting function
                FC_BrainEEGLu(ChanPosColin27,MapGroup{iStimGroup,iBOI,iTrialType},[],FCpara)
                axis off

                % Add text below the plot
                text(0.5, -0.15, sprintf('Pairs > Th: %d', nChPairAboveThreshold), ...
                    'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 10);

                if iStimGroup==nGroups
                    xlabel(BandName{iBOI});
                    text(-0.1,0,0.1,[BandName{iBOI} ' (' BandHzName{iBOI} ')'],'horizontalalignment','center','verticalalignment','top','fontsize',10)
                end

                if iBOI==1
                    ylabel(GroupName{FCgroups(iStimGroup)})
                    yt=text(0,0.1,0.1,GroupName{FCgroups(iStimGroup)},'horizontalalignment','center','verticalalignment','bottom','fontsize',10,'rotation',90);
                end
            end
        end


        % Define the edge threshold title
        edgeThresholdTitle = sprintf('Edge Threshold: %.4f', FCpara.EdgeTh);
        % Add the title displaying the edge threshold
        sgtitle(edgeThresholdTitle, 'FontSize', 6, 'FontWeight', 'bold', 'Interpreter', 'none');

        % % Add the title displaying the edge threshold and position it higher
        % titleHandle = sgtitle(edgeThresholdTitle, 'FontSize', 12);
        % titleHandle.Position = [0.5, 0.98, 0]; % [x, y, z] position in normalized figure units


        %          subplot('position',[0.5 0.51 0.3 0.01]);
        %      b=colorbar('southoutside');
        %      set(gca,'xtick',[],'ytick',[])
        %      set(b,'position',[0.5 0.5 0.3 0.03],'Limits',[0 1],'Ticks',[0 1],'Ticklabels',PowerLab);
        %      xlabel(b,'Log Normalized Power')
        LuFontStandard;
        papersizePX=[0 0 6*size(BOI,2) 6*nGroups+3];
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));

        % Adjust title spacing relative to the subplots
        t = sgtitle(edgeThresholdTitle, 'FontSize', 6);
        % t.Position(2) = t.Position(2) + 0.03; % Raise the title slightly

        saveWPLIFigName = ['WPLI40_' num2str(nGroups) 'Groups_' num2str(100*FCpara.EdgeTh) 'E-2EdgeThrs620.png'];
        % saveas(gcf,[SaveTemp saveWPLIFigName],'pdf');
        saveas(gcf,[SaveTemp saveWPLIFigName],'png');
        saveas(gcf,[SaveTemp saveWPLIFigName],'svg');
        % saveas(gcf,[SaveTemp saveWPLIFigName],'epsc');

        FCpara.EdgeTh=0.1; % reset
        %% WPLI figure ***  fig4d MS version
        % Fig 4 C and D are composed of panels created by this section, cut and pasted together in
        % illustrator
        % Lower Alpha and Upper Alpha only
        BOI=[8 10; 10 13]; %[1 4 8 8 10 13 30 39.5 43;4 8 13 10 13 30 37 41.5 100];
        BandName={'LowerAlpha','UpperAlpha'};
        BandHzName={'8-10Hz','10-13Hz'};
        % FCpara.EdgeTh=0.06;
        % FCpara.EdgeTh=0.15; % previous arbitray threshold
        FCpara.EdgeTh=0.1202; % top quartile threshold
        % FCpara.EdgeTh=averagePermutatedWPLI_2to55_top0_0001; % = 0.0894 - top 0.0001 permutation threshold
        % FCpara.EdgeTh=averagePermutatedWPLI_2to55_top0_001; % = 0.0484 -  top 0.001 permutation threshold

        FCgroups = [1,3,6]; % [1,3,6] = [40hz , Random, LightRT]

        figure;
        nGroups = length(FCgroups); % (1:3) Just first 3 groups
        for iBOI=1:size(BOI,2) % just lower alpha and upper alpha
            NeedI=find(Fplot>=BOI(1,iBOI)&Fplot<BOI(2,iBOI));
            for iStimGroup=1:length(FCgroups)
                if iFreqFunc==3
                    Tdata{iStimGroup,iBOI,iTrialType}=squeeze(FreqFunc{iFreqFunc}(WPLIall(SubjG{FCgroups(iStimGroup)},NeedI,EEGchInd,EEGchInd,iTrialType),[],2));
                    MapGroup{iStimGroup,iBOI,iTrialType}=squeeze(nanmean(FreqFunc{iFreqFunc}(WPLIall(SubjG{FCgroups(iStimGroup)},NeedI,EEGchInd,EEGchInd,iTrialType),[],2),1));
                else
                    Tdata{iStimGroup,iBOI,iTrialType}=squeeze(FreqFunc{iFreqFunc}(WPLIall(SubjG{FCgroups(iStimGroup)},NeedI,EEGchInd,EEGchInd,iTrialType),2));
                    MapGroup{iStimGroup,iBOI,iTrialType}=squeeze(nanmean(FreqFunc{iFreqFunc}(WPLIall(SubjG{FCgroups(iStimGroup)},NeedI,EEGchInd,EEGchInd,iTrialType),2),1));
                end

                EdgeColor=[0.8 0.8 0.8];
                %              FC_BrainEEGLu(ChanPosColin27,MapGroup{iG,iFF,iCom},[],EdgeTh,NodeTh,EdgeColor,[])
                %              axis off

                % subplotLU(nGroups,size(BOI,2),iStimGroup,iBOI); % old 3x2
                % subplotLU(1, nGroups*size(BOI,2), 1, iStimGroup+3*(iBOI-1)); % 1x6 grid
                subplotLU(size(BOI,2), nGroups, iBOI, iStimGroup); % 2x3 grid

                % WPLIsForPlot(iStimGroup,iBOI,iTrialType) = MapGroup{iStimGroup,iBOI,iTrialType};
                % nChPairAboveThreshold(iStimGroup,iBOI,iTrialType) = sum(WPLIsForPlot>FCpara.EdgeTh);

                % Calculate the number of channel pairs with WPLI > EdgeTh
                currentMap = MapGroup{iStimGroup, iBOI, iTrialType};
                nChPairAboveThreshold = sum(currentMap(:) > FCpara.EdgeTh);
                AboveThreholdChannelPairAllGroups(iStimGroup,iBOI,iTrialType) = nChPairAboveThreshold;

                %%% The plotting function
                FC_BrainEEGLu(ChanPosColin27,MapGroup{iStimGroup,iBOI,iTrialType},[],FCpara)
                axis off

                % Add text below the plot
                % text(0.5, -0.15, sprintf('Pairs > Th: %d', nChPairAboveThreshold), ...
                %     'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 10);

                if iStimGroup==2
                    xlabel(BandName{iBOI});
                    text(0.125,0,0.1,[BandName{iBOI} ' (' BandHzName{iBOI} ')'],'horizontalalignment','center','verticalalignment','top','fontsize',18)
                end

                    xlabel(GroupName{FCgroups(iStimGroup)})
                    yt=text(-0.125,0,0.1,GroupName{FCgroups(iStimGroup)},'horizontalalignment','center','verticalalignment','bottom','fontsize',14);

            end
        end


        % Define the edge threshold title
        % edgeThresholdTitle = sprintf('Edge Threshold: %.4f', FCpara.EdgeTh);
        % % Add the title displaying the edge threshold
        % sgtitle(edgeThresholdTitle, 'FontSize', 6, 'FontWeight', 'bold', 'Interpreter', 'none');
        % t.Position(2) = t.Position(2) + 0.05; % Move title slightly up

        % % Add the title displaying the edge threshold and position it higher
        % titleHandle = sgtitle(edgeThresholdTitle, 'FontSize', 12);
        % titleHandle.Position = [0.5, 0.98, 0]; % [x, y, z] position in normalized figure units


        %          subplot('position',[0.5 0.51 0.3 0.01]);
        %      b=colorbar('southoutside');
        %      set(gca,'xtick',[],'ytick',[])
        %      set(b,'position',[0.5 0.5 0.3 0.03],'Limits',[0 1],'Ticks',[0 1],'Ticklabels',PowerLab);
        %      xlabel(b,'Log Normalized Power')
        LuFontStandard;
        % papersizePX=[0 0 6*size(BOI,2)*4 6/3*nGroups+3]; % 1x6
        % 3x2 grid: papersizePX=[0 0 6*size(BOI,2) 6*nGroups+3];
        papersizePX=[0 0 nGroups*6 size(BOI,2)*6+2.5]; % 2x3 : [0 0 width{x} heigth{y}]
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));

        % Adjust title spacing relative to the subplots
        % t = sgtitle(edgeThresholdTitle, 'FontSize', 6);
        % t.Position(2) = t.Position(2) + 0.03; % Raise the title slightly

        
        % saveas(gcf,[SaveTemp saveWPLIFigName],'pdf');
        saveDateFig4D = datestr(datetime, 'yy-mm-dd_HHMMSSFFF');
        saveWPLIFigName = ['WPLI40_' num2str(nGroups) 'Groups_' saveDateFig4D '_' num2str(100*FCpara.EdgeTh) 'E-2EdgeThrs620.png'];
        saveas(gcf,['Fig4Panels/' saveWPLIFigName ],'png');
        % saveas(gcf,['Fig4Panels/' saveWPLIFigName],'svg');
        % saveas(gcf,[SaveTemp saveWPLIFigName],'epsc');

        FCpara.EdgeTh=0.1; % reset
        %% Perform Chi-Squared on proportions on elevated channels: 40vL & 40VRand
        % Data
        total_pairs = 496;
        elevated_40Hz = 93;
        elevated_Light = 32;

        % Define contingency table
        table_40Hz_Light = [elevated_40Hz, total_pairs - elevated_40Hz;
            elevated_Light, total_pairs - elevated_Light];

        % % Perform chi-squared test and compute effect size
        % [chi2_40Hz_Light, p_40Hz_Light, V_40Hz_Light] = analyzeChiSquared(table_40Hz_Light, total_pairs);
        % 
        % % Display results
        % fprintf('Results for 40Hz vs Light:\n');
        % fprintf('Chi-squared (X^2): %.2f\n', chi2_40Hz_Light);
        % fprintf('p-value: %.4f\n', p_40Hz_Light);
        % fprintf('Cramér''s V (Effect size): %.4f\n', V_40Hz_Light);

        %% Chi Squared Tests: 40vLight, 40vRandom for Lower and Upper Alpha Fig4D stats - OLD
        % Contingency tables for the tests

        % Lower Alpha:
        % 40Hz vs Light
        observed_LowerAlpha_40vLight = [93, 403;
            32, 464];

        % Lower Alpha: 40Hz vs Random
        observed_LowerAlpha_40vRandom = [93, 403;
            7, 489];

        % Upper Alpha: 
        % 40Hz vs Light
        observed_UpperAlpha_40vLight = [63, 433;
            7, 489];

        % Upper Alpha: 40Hz vs Random
        observed_UpperAlpha_40vRandom = [63, 433;
            267, 229];

        % Perform chi-squared tests
        fprintf('Lower Alpha (40Hz vs Light):\n');
        [chi2_LA_40vLight, p_LA_40vLight, V_LA_40vLight, dof_LA_40vLight, total_LA_40vLight] = chi_squared_test(observed_LowerAlpha_40vLight);

        fprintf('\nLower Alpha (40Hz vs Random):\n');
        [chi2_LA_40vRandom, p_LA_40vRandom, V_LA_40vRandom, dof_LA_40vRandom, total_LA_40vRandom] = chi_squared_test(observed_LowerAlpha_40vRandom);

        fprintf('\nUpper Alpha (40Hz vs Light):\n');
        [chi2_UA_40vLight, p_UA_40vLight, V_UA_40vLight, dof_UA_40vLight, total_UA_40vLight] = chi_squared_test(observed_UpperAlpha_40vLight);

        fprintf('\nUpper Alpha (40Hz vs Random):\n');
        [chi2_UA_40vRandom, p_UA_40vRandom, V_UA_40vRandom, dof_UA_40vRandom, total_UA_40vRandom] = chi_squared_test(observed_UpperAlpha_40vRandom);

        % FDR correction:
        % Lower Alpha p-values
        p_values_LowerAlpha = [p_LA_40vLight, p_LA_40vRandom];

        % Upper Alpha p-values
        p_values_UpperAlpha = [p_UA_40vLight, p_UA_40vRandom];

        % Apply FDR correction to Lower Alpha and Upper Alpha
        adjusted_p_LowerAlpha = fdr_correction(p_values_LowerAlpha);
        adjusted_p_UpperAlpha = fdr_correction(p_values_UpperAlpha);

        % Display results
        disp('FDR-corrected p-values for Lower Alpha:');
        disp(adjusted_p_LowerAlpha);

        disp('FDR-corrected p-values for Upper Alpha:');
        disp(adjusted_p_UpperAlpha);

        % Display FDR-corrected p-values for Lower Alpha in scientific notation
        disp('FDR-corrected p-values for Lower Alpha (scientific notation):');
        fprintf('%.15e\n', adjusted_p_LowerAlpha);

        % Display FDR-corrected p-values for Upper Alpha in scientific notation
        disp('FDR-corrected p-values for Upper Alpha (scientific notation):');
        fprintf('%.15e\n', adjusted_p_UpperAlpha);
        %% Chi-squared test of top quartile ch pairs & permutated channel pairs - no longer needed 1/16/24
        % % Assuming AboveThresholdChannelPairAllGroupsTopQuart and
        % % AboveThresholdChannelPairAllGroupsPermutated are already loaded.
        % 
        % % Set denominator for proportions
        % denominator = 496;
        % 
        % % Extract dimensions
        % dims = size(AboveThreholdChannelPairAllGroupsTopQuart);
        % 
        % % Initialize matrices for storing results
        % chi2_stat = zeros(dims); % Chi-squared statistic
        % p_value = zeros(dims);   % P-value
        % h_test = zeros(dims);    % Hypothesis test result (1: reject null, 0: fail to reject)
        % 
        % % Loop through each element
        % for i = 1:dims(1)
        %     for j = 1:dims(2)
        %         for k = 1:dims(3)
        %             % Observed data
        %             obs1 = AboveThreholdChannelPairAllGroupsTopQuart(i,j,k);
        %             obs2 = AboveThreholdChannelPairAllGroupsPermutated(i,j,k);
        % 
        %             % Proportions
        %             prop1 = obs1 / denominator;
        %             prop2 = obs2 / denominator;
        % 
        %             % Pooled proportion under null hypothesis
        %             pooled_p = (obs1 + obs2) / (2 * denominator);
        % 
        %             % Expected counts under null hypothesis
        %             exp1 = pooled_p * denominator;
        %             exp2 = pooled_p * denominator;
        % 
        %             % Chi-squared statistic for this pair
        %             chi2_stat(i,j,k) = ((obs1 - exp1)^2 / exp1) + ((obs2 - exp2)^2 / exp2);
        % 
        %             % Degrees of freedom
        %             df = 1;
        % 
        %             % Compute p-value
        %             p_value(i,j,k) = 1 - chi2cdf(chi2_stat(i,j,k), df);
        % 
        %             % Hypothesis test: reject null if p < 0.05
        %             h_test(i,j,k) = p_value(i,j,k) < 0.05;
        %         end
        %     end
        % end
        % 
        % % Display results for inspection
        % disp('Chi-squared statistics:');
        % disp(chi2_stat);
        % 
        % disp('P-values (3D matrix):');
        % disp(p_value);
        % 
        % disp('Hypothesis test results (3D matrix, 1: reject null, 0: fail to reject):');
        % disp(h_test);

        %% WPLI figure (Both Controls)
        % figure;
        % for iBOI=1:size(BOI,2)
        %     NeedI=find(Fplot>=BOI(1,iBOI)&Fplot<BOI(2,iBOI));
        %     for iStimGroup=1:length(SubjG)
        %         if iFreqFunc==3
        %             Tdata{iStimGroup,iBOI,iTrialType}=squeeze(FreqFunc{iFreqFunc}(WPLIall(SubjG{iStimGroup},NeedI,EEGchInd,EEGchInd,iTrialType),[],2));
        %             MapGroup{iStimGroup,iBOI,iTrialType}=squeeze(nanmean(FreqFunc{iFreqFunc}(WPLIall(SubjG{iStimGroup},NeedI,EEGchInd,EEGchInd,iTrialType),[],2),1));
        %         else
        %             Tdata{iStimGroup,iBOI,iTrialType}=squeeze(FreqFunc{iFreqFunc}(WPLIall(SubjG{iStimGroup},NeedI,EEGchInd,EEGchInd,iTrialType),2));
        %             MapGroup{iStimGroup,iBOI,iTrialType}=squeeze(nanmean(FreqFunc{iFreqFunc}(WPLIall(SubjG{iStimGroup},NeedI,EEGchInd,EEGchInd,iTrialType),2),1));
        %         end
        %
        %         EdgeColor=[0.8 0.8 0.8];
        %         %              FC_BrainEEGLu(ChanPosColin27,MapGroup{iG,iFF,iCom},[],EdgeTh,NodeTh,EdgeColor,[])
        %         %              axis off
        %
        %         subplotLU(2,size(BOI,2),iStimGroup,iBOI);
        %         FC_BrainEEGLu(ChanPosColin27,MapGroup{iStimGroup,iBOI,iTrialType},[],FCpara)
        %         axis off
        %         if iStimGroup==2
        %             xlabel(BandName{iBOI});
        %             text(-0.1,0,0.1,[BandName{iBOI} ' (' BandHzName{iBOI} ')'],'horizontalalignment','center','verticalalignment','top','fontsize',10)
        %         end
        %
        %         if iBOI==1
        %             ylabel(GroupName{iStimGroup})
        %             yt=text(0,0.1,0.1,GroupName{iStimGroup},'horizontalalignment','center','verticalalignment','bottom','fontsize',10,'rotation',90);
        %         end
        %     end
        % end
        % %          subplot('position',[0.5 0.51 0.3 0.01]);
        % %      b=colorbar('southoutside');
        % %      set(gca,'xtick',[],'ytick',[])
        % %      set(b,'position',[0.5 0.5 0.3 0.03],'Limits',[0 1],'Ticks',[0 1],'Ticklabels',PowerLab);
        % %      xlabel(b,'Log Normalized Power')
        % LuFontStandard;
        % papersizePX=[0 0 6*size(BOI,2) 6*2+3];
        % set(gcf, 'PaperUnits', 'centimeters');
        % set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
        %
        % WPLI40BothControlsFigName = 'WPLI40BothControls620';
        % % saveas(gcf,[SaveTemp WPLI40BothControlsFigName],'pdf');
        % saveas(gcf,[SaveTemp WPLI40BothControlsFigName],'png');
        % % saveas(gcf,[SaveTemp WPLI40BothControlsFigName],'epsc');

        %% WPLI Difference Figure (VARIABLES NEED TO BE RENAMED IN THIS SECTION)
        % figure;
        % pTEdgeTh=0.1;
        % %comparison between groups w/ freq data
        % diffTPmap=zeros(length(ChanEEGLab),length(ChanEEGLab),size(BOI,2),length(TrialType)); % for T test
        % diffRPmap=diffTPmap;
        % diffTmap=diffTPmap;
        %
        % % Group selection: 1=40Hz, 2=Light, 3=Random, 4=LightRT
        % Group1 = 1;
        % Group2 = 2;
        %
        % for iBOI=1:size(BOI,2)-1 % minus 1 to remove BOI with 60 Hz
        %     for iCh=1:length(ChanEEGLab)
        %         for jCh=iCh+1:length(ChanEEGLab)
        %             [~,diffTPmap(iCh,jCh,iBOI,iTrialType),~,stats]=ttest2(Tdata{Group1,iBOI,iTrialType}(:,iCh,jCh),Tdata{Group2,iBOI,iTrialType}(:,iCh,jCh));
        %             [diffRPmap(iCh,jCh,iBOI,iTrialType),~,~]=ranksum(Tdata{Group1,iBOI,iTrialType}(:,iCh,jCh),Tdata{Group2,iBOI,iTrialType}(:,iCh,jCh));
        %             diffTmap(iCh,jCh,iBOI,iTrialType)=stats.tstat;
        %         end
        %     end
        %     subplotLU(1,size(BOI,2),1,iBOI);
        %     Adj=diffTmap(:,:,iBOI,iTrialType);
        %     AdjP=diffTPmap(:,:,iBOI,iTrialType);
        %     Adj(AdjP>pTEdgeTh)=0;
        %     FC_BrainEEGLu(ChanPosColin27,Adj,[],FCparaT)
        %     axis off
        %
        %     if iBOI==1
        %         yt=text(-0.0,0.1,0.1,['T, Sig-Diff FC, ' GroupName{Group1} '-' GroupName{Group2} ],'horizontalalignment','center','verticalalignment','bottom','fontsize',10,'rotation',90);
        %     end
        %
        %     % %          subplotLU(2,size(BOI,2),2,iFF);
        %     % %           xlabel(BName{iFF});
        %     % %           text(0,-0.55,[BName{iFF} ' (' BName2{iFF} ')'],'horizontalalignment','center','verticalalignment','top','fontsize',10)
        %     % %              if iFF==1
        %     % %                 a=ylabel('40Hz-Random');
        %     % % %                 a.Position=[0.01 0.5 0.03 0.4];
        %     % % %                 a.verticalalignment='middle';
        %     % % %                 set(a,'Position',[0.01 0.5 0.03 0.4],'Verticalalignment','middle')
        %     % %                 set(a,'Verticalalignment','middle')
        %     % %                 yt=text(-120,0,'P, 40-Rand.','horizontalalignment','center','verticalalignment','bottom','fontsize',10,'rotation',90);
        %     % %              end
        %
        %     %                 xlabel(BName{iFF});
        %     text(-0.1,0,0.1,[BandName{iBOI} ' (' BandHzName{iBOI} ')'],'horizontalalignment','center','verticalalignment','top','fontsize',10)
        % end
        %
        % LuFontStandard;
        % papersizePX=[0 0 6*size(BOI,2) 6+2];
        % set(gcf, 'PaperUnits', 'centimeters');
        % set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
        % WPLIDiffFigName = ['WPLIDiff ' GroupName{Group1} '-' GroupName{Group2} ' p' num2str(pTEdgeTh)];
        % sgtitle(WPLIDiffFigName)
        %
        %
        % % saveas(gcf,[SaveTemp WPLIDiffBothControlsfigName],'pdf');
        % saveas(gcf,[SaveTemp WPLIDiffFigName '.png'],'png');
        % % saveas(gcf,[SaveTemp WPLIDiffBothControlsfigName],'epsc');

        %% Define Group Set Names - beginning of fig4e fig4f
        GroupSetsName{1} ='40andL'; % 40 and Light
        GroupSetsName{2}='40andR'; % 40 and Random
        GroupSetsName{3}='All3Groups'; % All three groups
        %% Loop thru all groups *** fig4e fig4f
        pSpearEdgeTh = 0.1;
        for iGroupSet = [1 2]% 1:length(GroupSetsName)
            if iGroupSet == 1
                %% 40&LightRT
                % included subs

                Group1 = 1; % 1= 40Hz flicker group
                Group2 = 6; % 2 = LightRT group
                dataName = [TrialTypeName{iTrialType} ' ' FreqFuncNames{iFreqFunc} ' p' num2str(pSpearEdgeTh) ' ' GroupName{Group1} GroupName{Group2}];

                IncludedSubj=union(SubjG{Group1},SubjG{Group2}); %40 + LightRT

            elseif iGroupSet == 2
                %% 40andR
                Group1 = 1; % 1= 40Hz flicker group
                Group2 = 3; % 3 = Random group
                dataName = [TrialTypeName{iTrialType} ' ' FreqFuncNames{iFreqFunc} ' p' num2str(pSpearEdgeTh) ' ' GroupName{Group1} GroupName{Group2}];

                IncludedSubj=union(SubjG{Group1},SubjG{Group2}); %40 + Random

            elseif iGroupSet == 3
                %% All3Groups
                Group1 = 1; % 1= 40Hz flicker group
                Group2 = 6; % 6 = LightRT
                Group3 = 3; % 3 = Random group
                dataName = [TrialTypeName{iTrialType} FreqFuncNames{iFreqFunc} GroupName{Group1} GroupName{Group2} GroupName{Group3}];

                IncludedSubj=union(SubjG{Group1},SubjG{Group2},SubjG{Group3}); %40 + LightRT + Random

            end
            %% Topo plot ALPHA GENERAL
            %% Define Bands of Interest
            % BOI=[8 8 10; 13 10 13];
            % BandName={'Alpha','LowAlpha','HighAlpha'};
            % BandHzName={'8-13Hz','8-10Hz','10-13Hz'};
            BOI=[ 8 10; 10 13];
            BandName={'LowAlpha','HighAlpha'};
            BandHzName={'8-10Hz','10-13Hz'};
            %% Spearman Acc-FC
            SpearmanAccFctopo = figure;
            for iBOI=1:size(BOI,2)
                NeedI=find(Fplot>=BOI(1,iBOI)&Fplot<=BOI(2,iBOI));
                if iFreqFunc==3
                    WPLItemp=squeeze(FreqFunc{iFreqFunc}(WPLIall(IncludedSubj,NeedI,EEGchInd,EEGchInd,iTrialType),[],2));
                else
                    WPLItemp=squeeze(FreqFunc{iFreqFunc}(WPLIall(IncludedSubj,NeedI,EEGchInd,EEGchInd,iTrialType),2));
                end
                Acctemp=Acc(IncludedSubj);

                %% Include only subjects with high accuracy (>80%) - Accuracy cuts
                highAccThreshold = 0.8;
                highAccSubjectsIndex = Acctemp>highAccThreshold;
                WPLItemp = WPLItemp(highAccSubjectsIndex,:,:);
                Acctemp80 = Acctemp(highAccSubjectsIndex);

                %% loop thru each channel - calculate Spearman correlation R and p-value
                for iCh=1:length(ChanEEGLab)-1
                    %              for jCh=2:length(ChanEEGLab)
                    [rSpear(iCh,iCh+1:end,iBOI,iTrialType),pSpear(iCh,iCh+1:end,iBOI,iTrialType)]=corr(squeeze(WPLItemp(:,iCh,iCh+1:end)),Acctemp80,'type','spearman','rows','pairwise');
                    %              end
                end
                subplotLU(1,size(BOI,2),1,iBOI); %can change to 2 for second row

                Adj=rSpear(:,:,iBOI,iTrialType); %32x32 of r values?
                AdjP=pSpear(:,:,iBOI,iTrialType); %32x32 of p values?
                Adj(AdjP>pSpearEdgeTh)=0;  % deletes all r values of ch-pairs with p-value > threshold
                signifchPairRAccTopo = Adj~=0; % creates 32x32 logical of significant channel pairs
                signifChPairNameTopo{iFreqFunc,iTrialType,iGroupSet,iBOI,1} = getListofSignifChPairs(signifchPairRAccTopo,EEGch); % get list of channel names not equal to zero

                %% Now, plot only for significant positive WPLI-accuracy correlations
                % posWPLIAccChPairFolder = ['PositiveAccCorr\' GroupSetsName{iGroupSet} '\' BName{iBOI} '\'];
                % mkdir([SaveTemp posWPLIAccChPairFolder])
                % for iCh = 1:length(EEGchInd)
                %     for jCh = iCh+1:length(EEGchInd)
                %         % Check if the channel pair has a positive significant correlation
                %         if Adj(iCh, jCh) > 0
                %             clear DataPlot;
                %
                %             for iStimGroup = 1:length(SubjGWPLI)
                %                 DataPlot{iStimGroup} = squeeze(WPLIall(SubjGWPLI{iStimGroup}, :, EEGchInd(iCh), EEGchInd(jCh), iTrialType));
                %                 Invalid = isnan(DataPlot{iStimGroup}(:, 1));
                %                 DataPlot{iStimGroup}(Invalid, :) = [];
                %             end
                %
                %             if isempty(DataPlot{1}) || isempty(DataPlot{2})
                %                 continue;
                %             end
                %
                %             iPlot = iPlot + 1;
                %
                %             % Plot WPLI data for the channel pair
                %             figure;
                %             xrange = [0 50]; % Frequency range for x-axis
                %             subplot('Position', [0.1 0.1 0.88 0.88]);
                %             RateHist_GroupPlot(Fplot, DataPlot, FlickerColor, ParamWPLI); % Custom function for plotting
                %
                %             text(range(xrange)/2 - 5, ParamWPLI.Ytick(end), [EEGch{iCh} '-' EEGch{jCh}], 'FontSize', 10); % Add channel pair label
                %             set(gca, 'xlim', xrange, 'xtick', [1 4 8 13 30 40 50 60 80 100], 'ylim', [0 0.2], 'ytick', ParamWPLI.Ytick);
                %             xlabel('Frequency (Hz)');
                %             ylabel('WPLI');
                %
                %             % Formatting for grid and font size
                %             ax = gca;
                %             ax.XGrid = 'on';
                %             ax.YGrid = 'off';
                %             ax.FontSize = 7;
                %
                %             % Set paper size and save the figure
                %             LuFontStandard; % Custom function for standard fonts
                %             papersizePX = 1.3 * [0 0 5 3.2]; % Custom figure size
                %             set(gcf, 'PaperUnits', 'centimeters');
                %             set(gcf, 'PaperPosition', papersizePX, 'PaperSize', papersizePX(3:4));
                %
                %             % Save the figure in SVG format
                %             % saveas(gcf, [SaveTemp 'PositiveAccCorr\' '40LightRandBands_' EEGch{iCh} '-' EEGch{jCh}], 'svg');
                %             saveas(gcf, [SaveTemp posWPLIAccChPairFolder 'WPLI_40LR_' EEGch{iCh} '-' EEGch{jCh}], 'png');
                %         end
                %     end
                % end
                % figure(SpearmanAccFctopo);

                %% PLotting and axises  ****fig4e
                FC_BrainEEGLu(ChanPosColin27,Adj,[sum(Adj,1)/2],FCspear) % plotting function! (3rd parameter is node weight)

                nSignifChPairs = sum(AdjP<pSpearEdgeTh & AdjP>0, 'all');

                axis off
                if iBOI==1
                    % yt=text(0,0.1,0.1,'Sig-Corr. FC-Acc','horizontalalignment','center','verticalalignment','bottom','fontsize',10,'rotation',90);
                end
                % text(-0.1,0,0.1,[BandName{iBOI} ' (' BandHzName{iBOI} ')'],'horizontalalignment','center','verticalalignment','top','fontsize',10)
                text(-0.11,0,-0.1,['nChpairs=' num2str(nSignifChPairs)],'horizontalalignment','center','verticalalignment','top','fontsize',10)


                if iBOI==2
                    % title(dataName,  'Units', 'normalized', 'Position', [0.5, 0.9, 0]) % MKA 9/23
                    % title(dataName,  'Units', 'normalized', 'Position', [0, 0.9, 0])
                end


                %% Display the list of significant channel pairs below the subplot
                % signifChPairNames = getListofSignifChPairs(signifchPairRAccTopo, EEGch);
                % % text(0.5, -0.2, strjoin(signifChPairNames, ', '), 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 8);
                % % verticalSignifChPairs = strjoin(signifChPairNames, '\n');
                % % text(-0.11, -0.2, verticalSignifChPairs, 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 8);
                %
                % % New subplot for channel pair names (as three columns)
                % subplotLU(2, size(BOI, 2), 2, iBOI);  % New row (2nd row) for the names
                %
                % % Divide the list into three columns
                % numNames = length(signifChPairNames);
                % numPerCol = ceil(numNames / 3);
                %
                % % Split the list of significant channel pairs into three columns
                % col1 = signifChPairNames(1:numPerCol);
                % col2 = signifChPairNames(numPerCol+1:min(2*numPerCol, numNames));
                % col3 = signifChPairNames(2*numPerCol+1:end);
                %
                % % Prepare the text to display in columns
                % colText = sprintf('%s\n', col1{:});
                % colText2 = sprintf('%s\n', col2{:});
                % colText3 = sprintf('%s\n', col3{:});
                %
                % % Display the three columns of significant channel pairs
                % text(0.2, 0.5, colText, 'Units', 'normalized', 'HorizontalAlignment', 'left', 'FontSize', 6);
                % text(0.5, 0.5, colText2, 'Units', 'normalized', 'HorizontalAlignment', 'left', 'FontSize', 6);
                % text(0.8, 0.5, colText3, 'Units', 'normalized', 'HorizontalAlignment', 'left', 'FontSize', 6);
                % axis off
            end

            %% Set up and save figure
            LuFontStandard;
            papersizePX=[0 0 6*size(BOI,2) 6+2]; % size of paper - height and width of brain plot 9/10. gain increase the papersize to try to increase the resolution
            set(gcf, 'PaperUnits', 'centimeters');
            set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));


            SpearmanFCAccName1 = ['SpearmanFCAcc_ALPHA' '_' dataName];
            % saveas(gcf,[SaveTemp SpearmanFCAccName1],'pdf');
            saveas(gcf,[SaveTemp SpearmanFCAccName1 '.png'],'png');
            print(gcf, '-dsvg', [SaveTemp SpearmanFCAccName1 '.svg'] , '-r0');  % -r0 ensures full vector output.
            % saveas(gcf,[SaveTemp SpearmanFCAccName1],'svg');  % -r0 ensures full vector output.);
            % saveas(gcf,[SaveTemp SpearmanFCAccName1 '.eps'],'epsc');
            close all

            %% Spearman RT-FC - Topo plot GENERAL ALPHA   ***fig4f
            SpearmanRTFCtopo = figure;

            for iBOI=1:size(BOI,2)
                NeedI=find(Fplot>=BOI(1,iBOI)&Fplot<=BOI(2,iBOI));
                if iFreqFunc==3
                    WPLItemp=squeeze(FreqFunc{iFreqFunc}(WPLIall(IncludedSubj,NeedI,EEGchInd,EEGchInd,iTrialType),[],2));
                else
                    WPLItemp=squeeze(FreqFunc{iFreqFunc}(WPLIall(IncludedSubj,NeedI,EEGchInd,EEGchInd,iTrialType),2));
                end
                RTtemp=SubjsAvgRT(IncludedSubj);

                %% Include only subjects with high accuracy (>80%) - Accuracy cuts
                highAccThreshold = 0.8;
                highAccSubjectsIndex = Acctemp>highAccThreshold;
                WPLItemp = WPLItemp(highAccSubjectsIndex,:,:);
                RTtemp = RTtemp(highAccSubjectsIndex);

                %%
                for iCh=1:length(ChanEEGLab)-1
                    %              for jCh=2:length(ChanEEGLab)
                    [rSpear(iCh,iCh+1:end,iBOI,iTrialType),pSpear(iCh,iCh+1:end,iBOI,iTrialType)]=corr(squeeze(WPLItemp(:,iCh,iCh+1:end)),RTtemp,'type','spearman','rows','pairwise');
                    %              end
                end
                subplotLU(1,size(BOI,2),1,iBOI); %can change to 2 for second row

                Adj=rSpear(:,:,iBOI,iTrialType); %32x32 of r values?
                AdjP=pSpear(:,:,iBOI,iTrialType); %32x32 of p values?
                Adj(AdjP>pSpearEdgeTh)=0;
                signifchPairRRTTopo = Adj~=0; % creates 32x32 logical of significant channel pairs
                signifChPairNameTopo{iFreqFunc,iTrialType,iGroupSet,iBOI,2} = getListofSignifChPairs(signifchPairRRTTopo,EEGch); % get list of channel names not equal to zero

                %% Now, plot only for significant negative WPLI-accuracy correlations
                % negWPLIRTChPairFolder = ['NegativeRTCorr\' GroupSetsName{iGroupSet} '\' BName{iBOI} '\'];
                % mkdir([SaveTemp negWPLIRTChPairFolder])
                % for iCh = 1:length(EEGchInd)
                %     for jCh = iCh+1:length(EEGchInd)
                %         % Check if the channel pair has a positive significant correlation
                %         if Adj(iCh, jCh) < 0
                %             clear DataPlot;
                %
                %             for iStimGroup = 1:length(SubjGWPLI)
                %                 DataPlot{iStimGroup} = squeeze(WPLIall(SubjGWPLI{iStimGroup}, :, EEGchInd(iCh), EEGchInd(jCh), iTrialType));
                %                 Invalid = isnan(DataPlot{iStimGroup}(:, 1));
                %                 DataPlot{iStimGroup}(Invalid, :) = [];
                %             end
                %
                %             if isempty(DataPlot{1}) || isempty(DataPlot{2})
                %                 continue;
                %             end
                %
                %             iPlot = iPlot + 1;
                %
                %             % Plot WPLI data for the channel pair
                %             figure;
                %             xrange = [0 50]; % Frequency range for x-axis
                %             subplot('Position', [0.1 0.1 0.88 0.88]);
                %             RateHist_GroupPlot(Fplot, DataPlot, FlickerColor, ParamWPLI); % Custom function for plotting
                %
                %             text(range(xrange)/2 - 5, ParamWPLI.Ytick(end), [EEGch{iCh} '-' EEGch{jCh}], 'FontSize', 10); % Add channel pair label
                %             set(gca, 'xlim', xrange, 'xtick', [1 4 8 13 30 40 50 60 80 100], 'ylim', [0 0.2], 'ytick', ParamWPLI.Ytick);
                %             xlabel('Frequency (Hz)');
                %             ylabel('WPLI');
                %
                %             % Formatting for grid and font size
                %             ax = gca;
                %             ax.XGrid = 'on';
                %             ax.YGrid = 'off';
                %             ax.FontSize = 7;
                %
                %             % Set paper size and save the figure
                %             LuFontStandard; % Custom function for standard fonts
                %             papersizePX = 1.3 * [0 0 5 3.2]; % Custom figure size
                %             set(gcf, 'PaperUnits', 'centimeters');
                %             set(gcf, 'PaperPosition', papersizePX, 'PaperSize', papersizePX(3:4));
                %
                %             % Save the figure in SVG format
                %             % saveas(gcf, [SaveTemp 'PositiveAccCorr\' '40LightRandBands_' EEGch{iCh} '-' EEGch{jCh}], 'svg');
                %             saveas(gcf, [SaveTemp negWPLIRTChPairFolder 'WPLI_40LR_' EEGch{iCh} '-' EEGch{jCh}], 'png');
                %         end
                %     end
                % end
                % figure(SpearmanRTFCtopo);

                %% Plot & axises  ***fig4f
                FC_BrainEEGLu(ChanPosColin27,Adj,[sum(Adj,1)/2],FCspear)

                nSignifChPairs = sum(AdjP<pSpearEdgeTh & AdjP>0, 'all');

                axis off
                if iBOI==1
                    % yt=text(0,0.1,0.1,'Sig-Corr. FC-RT','horizontalalignment','center','verticalalignment','bottom','fontsize',10,'rotation',90);
                end
                % text(-0.1,0,0.1,[BandName{iBOI} ' (' BandHzName{iBOI} ')'],'horizontalalignment','center','verticalalignment','top','fontsize',10)
                text(-0.11,0,-0.1,['nChpairs=' num2str(nSignifChPairs)],'horizontalalignment','center','verticalalignment','top','fontsize',10)

                if iBOI==2
                    % title(dataName,  'Units', 'normalized', 'Position', [0.5, 0.9, 0]) % MKA 9/23
                    % title(dataName,  'Units', 'normalized', 'Position', [0, 0.9, 0])
                end

                %% Display the list of significant channel pairs below the subplot
                % signifChPairNames = getListofSignifChPairs(signifchPairRRTTopo, EEGch);
                %
                % % New subplot for channel pair names (as three columns)
                % subplotLU(2, size(BOI, 2), 2, iBOI);  % New row (2nd row) for the names
                %
                % % Divide the list into three columns
                % numNames = length(signifChPairNames);
                % numPerCol = ceil(numNames / 3);
                %
                % % Split the list of significant channel pairs into three columns
                % col1 = signifChPairNames(1:numPerCol);
                % col2 = signifChPairNames(numPerCol+1:min(2*numPerCol, numNames));
                % col3 = signifChPairNames(2*numPerCol+1:end);
                %
                % % Prepare the text to display in columns
                % colText = sprintf('%s\n', col1{:});
                % colText2 = sprintf('%s\n', col2{:});
                % colText3 = sprintf('%s\n', col3{:});
                %
                % % Display the three columns of significant channel pairs
                % text(0.2, 0.5, colText, 'Units', 'normalized', 'HorizontalAlignment', 'left', 'FontSize', 6);
                % text(0.5, 0.5, colText2, 'Units', 'normalized', 'HorizontalAlignment', 'left', 'FontSize', 6);
                % text(0.8, 0.5, colText3, 'Units', 'normalized', 'HorizontalAlignment', 'left', 'FontSize', 6);
                % axis off
            end

            LuFontStandard;
            papersizePX=[0 0 6*size(BOI,2) 6+2];
            set(gcf, 'PaperUnits', 'centimeters');
            set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));

            SpearmanFCRTName2 = ['SpearmanFCRT_ALPHA'  '_' dataName];
            % saveas(gcf,[SaveTemp SpearmanFCRTName2],'pdf');
            saveas(gcf,[SaveTemp SpearmanFCRTName2 '.png'],'png');
            saveas(gcf,[SaveTemp SpearmanFCRTName2 '.svg'],'svg');
            % saveas(gcf,[SaveTemp SpearmanFCRTName2 '.eps'],'epsc');
            close all
        end
        %% Spearman Acc-FC - Topo plot (pre 6/18/24) 40 vs Light (All BANDS)
        % figure;
        %
        % Group1 = 1; % 1= 40Hz flicker group
        % Group2 = 2; % 2 = Light group
        % dataName = [TrialTypeName{iTrialType} FreqFuncNames{iFreqFunc} GroupName{Group1} GroupName{Group2}];
        %
        % IncludedSubj=union(SubjG{Group1},SubjG{Group2}); %40 vs Light
        % for iBOI=1:size(BOI,2)
        %     NeedI=find(Fplot>=BOI(1,iBOI)&Fplot<=BOI(2,iBOI));
        %     if iFreqFunc==3
        %         WPLItemp=squeeze(FreqFunc{iFreqFunc}(WPLIall(IncludedSubj,NeedI,EEGchInd,EEGchInd,iTrialType),[],2));
        %     else
        %         WPLItemp=squeeze(FreqFunc{iFreqFunc}(WPLIall(IncludedSubj,NeedI,EEGchInd,EEGchInd,iTrialType),2));
        %     end
        %     Acctemp=Acc(IncludedSubj);
        %
        %     for iCh=1:length(ChanEEGLab)-1
        %         %              for jCh=2:length(ChanEEGLab)
        %         WPLIijChPair = squeeze(WPLItemp(:,iCh,iCh+1:end));
        %         [rSpear(iCh,iCh+1:end,iBOI,iTrialType),pSpear(iCh,iCh+1:end,iBOI,iTrialType)]=corr(WPLIijChPair,Acctemp,'type','spearman','rows','pairwise');
        %         %              end
        %     end
        %     subplotLU(1,size(BOI,2),1,iBOI);
        %
        %     Adj=rSpear(:,:,iBOI,iTrialType);
        %     AdjP=pSpear(:,:,iBOI,iTrialType);
        %     Adj(AdjP>pSpearEdgeTh)=0;
        %     FC_BrainEEGLu(ChanPosColin27,Adj,[],FCspear)
        %     axis off
        %     if iBOI==1
        %         yt=text(0,0.1,0.1,'Sig-Corr. FC-Acc','horizontalalignment','center','verticalalignment','bottom','fontsize',10,'rotation',90);
        %     end
        %     text(-0.1,0,0.1,[BandName{iBOI} ' (' BandHzName{iBOI} ')'],'horizontalalignment','center','verticalalignment','top','fontsize',10)
        %
        %     if iBOI==4
        %         title(dataName,  'Units', 'normalized', 'Position', [0.5, 0.9, 0])
        %     end
        % end
        %
        % LuFontStandard;
        % papersizePX=[0 0 6*size(BOI,2) 6+2];
        % set(gcf, 'PaperUnits', 'centimeters');
        % set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
        %
        %
        % SpearmanFCAccName1 = ['SpearmanFCAcc' dataName];
        % saveas(gcf,[SaveTemp SpearmanFCAccName1],'pdf');
        % saveas(gcf,[SaveTemp SpearmanFCAccName1],'png');
        % saveas(gcf,[SaveTemp SpearmanFCAccName1 '.eps'],'epsc');
        % close all
        %% Spearman Acc-FC - Topo plot (pre 6/18/24) (40 vs Light) ALPHA
        % BOI=[8 8 10; 13 10 13];
        % BandName={'Alpha','LowAlpha','HighAlpha'};
        % % Gamma-E is  40
        % BandHzName={'8-13Hz','8-10Hz','10-13Hz',};
        % figure;
        %
        % Group1 = 1; % 1= 40Hz flicker group
        % Group2 = 2; % 2 = Light group
        % dataName = [TrialTypeName{iTrialType} FreqFuncNames{iFreqFunc} GroupName{Group1} GroupName{Group2}];
        %
        % IncludedSubj=union(SubjG{Group1},SubjG{Group2}); %40 + Light
        % for iBOI=1:size(BOI,2)
        %     NeedI=find(Fplot>=BOI(1,iBOI)&Fplot<=BOI(2,iBOI));
        %     if iFreqFunc==3
        %         WPLItemp=squeeze(FreqFunc{iFreqFunc}(WPLIall(IncludedSubj,NeedI,EEGchInd,EEGchInd,iTrialType),[],2));
        %     else
        %         WPLItemp=squeeze(FreqFunc{iFreqFunc}(WPLIall(IncludedSubj,NeedI,EEGchInd,EEGchInd,iTrialType),2));
        %     end
        %     Acctemp=Acc(IncludedSubj);
        %
        %     for iCh=1:length(ChanEEGLab)-1
        %         %              for jCh=2:length(ChanEEGLab)
        %         [rSpear(iCh,iCh+1:end,iBOI,iTrialType),pSpear(iCh,iCh+1:end,iBOI,iTrialType)]=corr(squeeze(WPLItemp(:,iCh,iCh+1:end)),Acctemp,'type','spearman','rows','pairwise');
        %         %              end
        %     end
        %     subplotLU(1,size(BOI,2),1,iBOI);
        %
        %     Adj=rSpear(:,:,iBOI,iTrialType);
        %     AdjP=pSpear(:,:,iBOI,iTrialType);
        %     Adj(AdjP>pSpearEdgeTh)=0;
        %     FC_BrainEEGLu(ChanPosColin27,Adj,[],FCspear)
        %
        %     axis off
        %     if iBOI==1
        %         yt=text(0,0.1,0.1,'Sig-Corr. FC-Acc','horizontalalignment','center','verticalalignment','bottom','fontsize',10,'rotation',90);
        %     end
        %     text(-0.1,0,0.1,[BandName{iBOI} ' (' BandHzName{iBOI} ')'],'horizontalalignment','center','verticalalignment','top','fontsize',10)
        %
        %     if iBOI==2
        %         title(dataName,  'Units', 'normalized', 'Position', [0.5, 0.9, 0])
        %     end
        % end
        %
        % LuFontStandard;
        % papersizePX=[0 0 6*size(BOI,2) 6+2];
        % set(gcf, 'PaperUnits', 'centimeters');
        % set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
        %
        %
        % SpearmanFCAccName1 = ['SpearmanFCAccALPHA_p' num2str(pSpearEdgeTh*100) '_' dataName];
        % % saveas(gcf,[SaveTemp SpearmanFCAccName1],'pdf');
        % saveas(gcf,[SaveTemp SpearmanFCAccName1],'png');
        % % saveas(gcf,[SaveTemp SpearmanFCAccName1 '.eps'],'epsc');
        % close all
        %% Spearman Acc-FC - Topo plot (40 vs Random) All BANDS
        % BOI=[1 4 8 13 30 39.5 43;4 8 13 30 37 41.5 100];
        % BandName={'Delta','Theta','Alpha','Beta','Gamma-1','Gamma-E','Gamma-2'};
        % % Gamma-E is  40
        % BandHzName={'1-4 Hz','4-8 Hz','8-13Hz','13-30 Hz','30-37 Hz','39-41 Hz','43-100 Hz'};
        % figure;
        %
        % Group1 = 1; % 1= 40Hz flicker group
        % Group2 = 3; % 3 = Random group
        % dataName = [TrialTypeName{iTrialType} FreqFuncNames{iFreqFunc} GroupName{Group1} GroupName{Group2}];
        %
        % IncludedSubj=union(SubjG{1},SubjG{3}); % 40 vs Random
        % for iBOI=1:size(BOI,2)
        %     NeedI=find(Fplot>=BOI(1,iBOI)&Fplot<=BOI(2,iBOI));
        %     if iFreqFunc==3
        %         WPLItemp=squeeze(FreqFunc{iFreqFunc}(WPLIall(IncludedSubj,NeedI,EEGchInd,EEGchInd,iTrialType),[],2));
        %
        %     else
        %         WPLItemp=squeeze(FreqFunc{iFreqFunc}(WPLIall(IncludedSubj,NeedI,EEGchInd,EEGchInd,iTrialType),2));
        %     end
        %     Acctemp=Acc(IncludedSubj);
        %
        %     for iCh=1:length(ChanEEGLab)-1
        %         %              for jCh=2:length(ChanEEGLab)
        %         [rSpear(iCh,iCh+1:end,iBOI,iTrialType),pSpear(iCh,iCh+1:end,iBOI,iTrialType)]=corr(squeeze(WPLItemp(:,iCh,iCh+1:end)),Acctemp,'type','spearman','rows','pairwise');
        %         %              end
        %     end
        %     subplotLU(1,size(BOI,2),1,iBOI);
        %
        %     Adj=rSpear(:,:,iBOI,iTrialType);
        %     AdjP=pSpear(:,:,iBOI,iTrialType);
        %     Adj(AdjP>pSpearEdgeTh)=0;
        %     FC_BrainEEGLu(ChanPosColin27,Adj,[],FCspear)
        %     axis off
        %     if iBOI==1
        %         yt=text(0,0.1,0.1,'Sig-Corr. FC-Acc','horizontalalignment','center','verticalalignment','bottom','fontsize',10,'rotation',90);
        %     end
        %     text(-0.1,0,0.1,[BandName{iBOI} ' (' BandHzName{iBOI} ')'],'horizontalalignment','center','verticalalignment','top','fontsize',10)
        %     if iBOI==4
        %         title(dataName,  'Units', 'normalized', 'Position', [0.5, 0.9, 0])
        %     end
        % end
        %
        % LuFontStandard;
        % papersizePX=[0 0 6*size(BOI,2) 6+2];
        % set(gcf, 'PaperUnits', 'centimeters');
        % set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
        %
        % SpearmanFCAccfigName = ['SpearmanFCAcc' dataName];
        % saveas(gcf,[SaveTemp SpearmanFCAccfigName],'pdf');
        % saveas(gcf,[SaveTemp SpearmanFCAccfigName],'png');
        % saveas(gcf,[SaveTemp SpearmanFCAccfigName 'eps'],'epsc');
        % close all
        %% Spearman Acc-FC - Topo plot  (40+RANDOM) ALPHA
        % BOI=[8 8 10; 13 10 13];
        % BandName={'Alpha','LowAlpha','HighAlpha'};
        % % Gamma-E is  40
        % BandHzName={'8-13Hz','8-10Hz','10-13Hz',};
        % figure;
        %
        % Group1 = 1; % 1= 40Hz flicker group
        % Group2 = 3; % 2 = Random group
        % dataName = [TrialTypeName{iTrialType} FreqFuncNames{iFreqFunc} GroupName{Group1} GroupName{Group2}];
        %
        % IncludedSubj=union(SubjG{Group1},SubjG{Group2}); %40 vs Random
        % for iBOI=1:size(BOI,2)
        %     NeedI=find(Fplot>=BOI(1,iBOI)&Fplot<=BOI(2,iBOI));
        %     if iFreqFunc==3
        %         WPLItemp=squeeze(FreqFunc{iFreqFunc}(WPLIall(IncludedSubj,NeedI,EEGchInd,EEGchInd,iTrialType),[],2));
        %
        %     else
        %         WPLItemp=squeeze(FreqFunc{iFreqFunc}(WPLIall(IncludedSubj,NeedI,EEGchInd,EEGchInd,iTrialType),2));
        %     end
        %     Acctemp=Acc(IncludedSubj);
        %
        %     for iCh=1:length(ChanEEGLab)-1
        %         %              for jCh=2:length(ChanEEGLab)
        %         [rSpear(iCh,iCh+1:end,iBOI,iTrialType),pSpear(iCh,iCh+1:end,iBOI,iTrialType)]=corr(squeeze(WPLItemp(:,iCh,iCh+1:end)),Acctemp,'type','spearman','rows','pairwise');
        %         %              end
        %     end
        %     subplotLU(1,size(BOI,2),1,iBOI);
        %
        %     Adj=rSpear(:,:,iBOI,iTrialType);
        %     AdjP=pSpear(:,:,iBOI,iTrialType);
        %     Adj(AdjP>pSpearEdgeTh)=0;
        %     FC_BrainEEGLu(ChanPosColin27,Adj,[],FCspear)
        %     axis off
        %     if iBOI==1
        %         yt=text(0,0.1,0.1,'Sig-Corr. FC-Acc','horizontalalignment','center','verticalalignment','bottom','fontsize',10,'rotation',90);
        %     end
        %     text(-0.1,0,0.1,[BandName{iBOI} ' (' BandHzName{iBOI} ')'],'horizontalalignment','center','verticalalignment','top','fontsize',10)
        %
        %     if iBOI==2
        %         title(dataName,  'Units', 'normalized', 'Position', [0.5, 0.9, 0])
        %     end
        % end
        %
        % LuFontStandard;
        % papersizePX=[0 0 6*size(BOI,2) 6+2];
        % set(gcf, 'PaperUnits', 'centimeters');
        % set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
        %
        %
        % SpearmanFCAccName1 = ['SpearmanFCAccALPHA' dataName];
        % % saveas(gcf,[SaveTemp SpearmanFCAccName1],'pdf');
        % saveas(gcf,[SaveTemp SpearmanFCAccName1],'png');
        % % saveas(gcf,[SaveTemp SpearmanFCAccName1 '.eps'],'epsc');
        % close all

        %% Spearman RT-FC - Topo plot (40+Light)  ALL BANDS
        % BOI=[1 4 8 13 30 39.5 43;4 8 13 30 37 41.5 100];
        % BandName={'Delta','Theta','Alpha','Beta','Gamma-1','Gamma-E','Gamma-2'};
        % % Gamma-E is  40
        % BandHzName={'1-4 Hz','4-8 Hz','8-13Hz','13-30 Hz','30-37 Hz','39-41 Hz','43-100 Hz'};
        %
        % figure;
        %
        % Group1 = 1; % 1= 40Hz flicker group
        % Group2 = 2; % 2 = Light group
        % dataName = [TrialTypeName{iTrialType} FreqFuncNames{iFreqFunc} GroupName{Group1} GroupName{Group2}];
        %
        % IncludedSubj=union(SubjG{Group1},SubjG{Group2}); %40 vs Light
        % for iBOI=1:size(BOI,2)
        %     NeedI=find(Fplot>=BOI(1,iBOI)&Fplot<=BOI(2,iBOI));
        %     if iFreqFunc==3
        %         WPLItemp=squeeze(FreqFunc{iFreqFunc}(WPLIall(IncludedSubj,NeedI,EEGchInd,EEGchInd,iTrialType),[],2));
        %
        %     else
        %         WPLItemp=squeeze(FreqFunc{iFreqFunc}(WPLIall(IncludedSubj,NeedI,EEGchInd,EEGchInd,iTrialType),2));
        %     end
        %     RTtemp=SubjsAvgRT(IncludedSubj);
        %     %Acctemp=Acc(IncludedSubj);
        %
        %     for iCh=1:length(ChanEEGLab)-1
        %         %              for jCh=2:length(ChanEEGLab)
        %         [rSpear(iCh,iCh+1:end,iBOI,iTrialType),pSpear(iCh,iCh+1:end,iBOI,iTrialType)]=corr(squeeze(WPLItemp(:,iCh,iCh+1:end)),RTtemp,'type','spearman','rows','pairwise');
        %         %              end
        %     end
        %     subplotLU(1,size(BOI,2),1,iBOI);
        %
        %     Adj=rSpear(:,:,iBOI,iTrialType);
        %     AdjP=pSpear(:,:,iBOI,iTrialType);
        %     Adj(AdjP>pSpearEdgeTh)=0;
        %     FC_BrainEEGLu(ChanPosColin27,Adj,[],FCspear)
        %     axis off
        %     if iBOI==1
        %         yt=text(0,0.1,0.1,'Sig-Corr. FC-RT','horizontalalignment','center','verticalalignment','bottom','fontsize',10,'rotation',90);
        %     end
        %     text(-0.1,0,0.1,[BandName{iBOI} ' (' BandHzName{iBOI} ')'],'horizontalalignment','center','verticalalignment','top','fontsize',10)
        %
        %     if iBOI==4
        %         title(dataName,  'Units', 'normalized', 'Position', [0.5, 0.9, 0])
        %     end
        % end
        %
        % LuFontStandard;
        % papersizePX=[0 0 6*size(BOI,2) 6+2];
        % set(gcf, 'PaperUnits', 'centimeters');
        % set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
        %
        %
        % SpearmanFCRTName1 = ['SpearmanFCRT' dataName];
        % saveas(gcf,[SaveTemp SpearmanFCRTName1],'pdf');
        % saveas(gcf,[SaveTemp SpearmanFCRTName1],'png');
        % saveas(gcf,[SaveTemp SpearmanFCRTName1 '.eps'],'epsc');
        % close all

        %% Spearman RT-FC - Topo plot (40 vs Random)  ALL BANDS
        % BOI=[1 4 8 13 30 39.5 43;4 8 13 30 37 41.5 100];
        % BandName={'Delta','Theta','Alpha','Beta','Gamma-1','Gamma-E','Gamma-2'};
        % % Gamma-E is  40
        % BandHzName={'1-4 Hz','4-8 Hz','8-13Hz','13-30 Hz','30-37 Hz','39-41 Hz','43-100 Hz'};
        % figure;
        %
        % Group1 = 1; % 1= 40Hz flicker group
        % Group2 = 3; % 3 = Random group
        % dataName = [TrialTypeName{iTrialType} FreqFuncNames{iFreqFunc} GroupName{Group1} GroupName{Group2}];
        %
        % IncludedSubj=union(SubjG{Group1},SubjG{Group2}); %40 vs Light
        % for iBOI=1:size(BOI,2)
        %     NeedI=find(Fplot>=BOI(1,iBOI)&Fplot<=BOI(2,iBOI));
        %     if iFreqFunc==3
        %         WPLItemp=squeeze(FreqFunc{iFreqFunc}(WPLIall(IncludedSubj,NeedI,EEGchInd,EEGchInd,iTrialType),[],2));
        %
        %     else
        %         WPLItemp=squeeze(FreqFunc{iFreqFunc}(WPLIall(IncludedSubj,NeedI,EEGchInd,EEGchInd,iTrialType),2));
        %     end
        %     RTtemp=SubjsAvgRT(IncludedSubj);
        %     %Acctemp=Acc(IncludedSubj);
        %
        %     for iCh=1:length(ChanEEGLab)-1
        %         %              for jCh=2:length(ChanEEGLab)
        %         [rSpear(iCh,iCh+1:end,iBOI,iTrialType),pSpear(iCh,iCh+1:end,iBOI,iTrialType)]=corr(squeeze(WPLItemp(:,iCh,iCh+1:end)),RTtemp,'type','spearman','rows','pairwise');
        %         %              end
        %     end
        %     subplotLU(1,size(BOI,2),1,iBOI);
        %
        %     Adj=rSpear(:,:,iBOI,iTrialType);
        %     AdjP=pSpear(:,:,iBOI,iTrialType);
        %     Adj(AdjP>pSpearEdgeTh)=0;
        %     FC_BrainEEGLu(ChanPosColin27,Adj,[],FCspear)
        %     axis off
        %     if iBOI==1
        %         yt=text(0,0.1,0.1,'Sig-Corr. FC-RT','horizontalalignment','center','verticalalignment','bottom','fontsize',10,'rotation',90);
        %     end
        %     text(-0.1,0,0.1,[BandName{iBOI} ' (' BandHzName{iBOI} ')'],'horizontalalignment','center','verticalalignment','top','fontsize',10)
        %
        %     if iBOI==4
        %         title(dataName,  'Units', 'normalized', 'Position', [0.5, 0.9, 0])
        %     end
        % end
        %
        % LuFontStandard;
        % papersizePX=[0 0 6*size(BOI,2) 6+2];
        % set(gcf, 'PaperUnits', 'centimeters');
        % set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
        %
        %
        % SpearmanFCRTName2 = ['SpearmanFCRT' dataName];
        % saveas(gcf,[SaveTemp SpearmanFCRTName2],'pdf');
        % saveas(gcf,[SaveTemp SpearmanFCRTName2],'png');
        % saveas(gcf,[SaveTemp SpearmanFCRTName2 '.eps'],'epsc');
        % close all

        %% Spearman RT-FC - Topo plot (40+Light) ALPHA ONLY
        % BOI=[8 8 10; 13 10 13];
        % BandName={'Alpha','LowAlpha','HighAlpha'};
        % % Gamma-E is  40
        % BandHzName={'8-13Hz','8-10Hz','10-13Hz',};
        %
        % figure;
        %
        % Group1 = 1; % 1= 40Hz flicker group
        % Group2 = 2; % 2 = Light group
        % dataName = [TrialTypeName{iTrialType} FreqFuncNames{iFreqFunc} GroupName{Group1} GroupName{Group2}];
        %
        % IncludedSubj=union(SubjG{Group1},SubjG{Group2}); %40 vs Light
        % for iBOI=1:size(BOI,2)
        %     NeedI=find(Fplot>=BOI(1,iBOI)&Fplot<=BOI(2,iBOI));
        %     if iFreqFunc==3
        %         WPLItemp=squeeze(FreqFunc{iFreqFunc}(WPLIall(IncludedSubj,NeedI,EEGchInd,EEGchInd,iTrialType),[],2));
        %
        %     else
        %         WPLItemp=squeeze(FreqFunc{iFreqFunc}(WPLIall(IncludedSubj,NeedI,EEGchInd,EEGchInd,iTrialType),2));
        %     end
        %     RTtemp=SubjsAvgRT(IncludedSubj);
        %     %Acctemp=Acc(IncludedSubj);
        %
        %     for iCh=1:length(ChanEEGLab)-1
        %         %              for jCh=2:length(ChanEEGLab)
        %         [rSpear(iCh,iCh+1:end,iBOI,iTrialType),pSpear(iCh,iCh+1:end,iBOI,iTrialType)]=corr(squeeze(WPLItemp(:,iCh,iCh+1:end)),RTtemp,'type','spearman','rows','pairwise');
        %         %              end
        %     end
        %     subplotLU(1,size(BOI,2),1,iBOI);
        %
        %     Adj=rSpear(:,:,iBOI,iTrialType);
        %     AdjP=pSpear(:,:,iBOI,iTrialType);
        %     Adj(AdjP>pSpearEdgeTh)=0;
        %     FC_BrainEEGLu(ChanPosColin27,Adj,[],FCspear)
        %     axis off
        %     if iBOI==1
        %         yt=text(0,0.1,0.1,'Sig-Corr. FC-RT','horizontalalignment','center','verticalalignment','bottom','fontsize',10,'rotation',90);
        %     end
        %     text(-0.1,0,0.1,[BandName{iBOI} ' (' BandHzName{iBOI} ')'],'horizontalalignment','center','verticalalignment','top','fontsize',10)
        %
        %     if iBOI==4
        %         title(dataName,  'Units', 'normalized', 'Position', [0.5, 0.9, 0])
        %     end
        % end
        %
        % LuFontStandard;
        % papersizePX=[0 0 6*size(BOI,2) 6+2];
        % set(gcf, 'PaperUnits', 'centimeters');
        % set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
        %
        %
        % SpearmanFCRTName1 = ['SpearmanFCRTALPHA' dataName];
        % % saveas(gcf,[SaveTemp SpearmanFCRTName1],'pdf');
        % saveas(gcf,[SaveTemp SpearmanFCRTName1],'png');
        % saveas(gcf,[SaveTemp SpearmanFCRTName1 '.eps'],'epsc');
        % close all
        %% Spearman RT-FC - Topo plot (40 vs Random) ALPHA ONLY
        % BOI=[8 8 10; 13 10 13];
        % BandName={'Alpha','LowAlpha','HighAlpha'};
        % % Gamma-E is  40
        % BandHzName={'8-13Hz','8-10Hz','10-13Hz',};
        %
        % figure;
        %
        % Group1 = 1; % 1= 40Hz flicker group
        % Group2 = 3; % 3 = Random group
        % dataName = [TrialTypeName{iTrialType} FreqFuncNames{iFreqFunc} GroupName{Group1} GroupName{Group2}];
        %
        % IncludedSubj=union(SubjG{Group1},SubjG{Group2}); %40 vs Random
        % for iBOI=1:size(BOI,2)
        %     NeedI=find(Fplot>=BOI(1,iBOI)&Fplot<=BOI(2,iBOI));
        %     if iFreqFunc==3
        %         WPLItemp=squeeze(FreqFunc{iFreqFunc}(WPLIall(IncludedSubj,NeedI,EEGchInd,EEGchInd,iTrialType),[],2));
        %
        %     else
        %         WPLItemp=squeeze(FreqFunc{iFreqFunc}(WPLIall(IncludedSubj,NeedI,EEGchInd,EEGchInd,iTrialType),2));
        %     end
        %     RTtemp=SubjsAvgRT(IncludedSubj);
        %     %Acctemp=Acc(IncludedSubj);
        %
        %     for iCh=1:length(ChanEEGLab)-1
        %         %              for jCh=2:length(ChanEEGLab)
        %         [rSpear(iCh,iCh+1:end,iBOI,iTrialType),pSpear(iCh,iCh+1:end,iBOI,iTrialType)]=corr(squeeze(WPLItemp(:,iCh,iCh+1:end)),RTtemp,'type','spearman','rows','pairwise');
        %         %              end
        %     end
        %     subplotLU(1,size(BOI,2),1,iBOI);
        %
        %     Adj=rSpear(:,:,iBOI,iTrialType);
        %     AdjP=pSpear(:,:,iBOI,iTrialType);
        %     Adj(AdjP>pSpearEdgeTh)=0;
        %     FC_BrainEEGLu(ChanPosColin27,Adj,[],FCspear)
        %     axis off
        %     if iBOI==1
        %         yt=text(0,0.1,0.1,'Sig-Corr. FC-RT','horizontalalignment','center','verticalalignment','bottom','fontsize',10,'rotation',90);
        %     end
        %     text(-0.1,0,0.1,[BandName{iBOI} ' (' BandHzName{iBOI} ')'],'horizontalalignment','center','verticalalignment','top','fontsize',10)
        %
        %     if iBOI==2
        %         title(dataName,  'Units', 'normalized', 'Position', [0.5, 0.9, 0])
        %     end
        % end
        %
        % LuFontStandard;
        % papersizePX=[0 0 6*size(BOI,2) 6+2];
        % set(gcf, 'PaperUnits', 'centimeters');
        % set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
        %
        %
        % SpearmanFCRTName2 = ['SpearmanFCRTALPHA' dataName];
        % % saveas(gcf,[SaveTemp SpearmanFCRTName2],'pdf');
        % saveas(gcf,[SaveTemp SpearmanFCRTName2],'png');
        % % saveas(gcf,[SaveTemp SpearmanFCRTName2 '.eps'],'epsc');
        % close all
    end
end
close all

%% PSD fig2b
SubSavePSD=[SavePath 'PSD\'];
PowerLim=[-7 -2];
PowerLab={'-7' '-2'};
TLim=[-5 5];
TLab={'-5' '5'};
PLim=[-4 4];
PLab={'10e-4' '10e-0'};

for iFreqFunc=1:length(FreqFuncNames)
    clear MapGroup diffTmap diffMap Tdata;
    clear diffTPmap diffRPmap diffTmap

    for iTrialType=1:length(TrialType)
        % iCom=1;
        SaveTemp=[SubSavePSD TrialTypeName{iTrialType} '\'];
        SaveTemp=[SaveTemp FreqFuncNames{iFreqFunc} '\' ];
        mkdir(SaveTemp)

        % Band of interest
        BOI=[5;15];
        BOI=[1 4 8 13 30 39.5 43;4 8 13 30 37 41.5 100];
        BandName={'Delta','Theta','Alpha','Beta','Gamma-1','Gamma-E','Gamma-2'};
        BandHzName={'1-4 Hz','4-8 Hz','8-13Hz','13-30 Hz','30-37 Hz','39-41 Hz','43-100 Hz'};
        %

        figure;
        IncludedSubj=union(SubjG{1},SubjG{2});
        for iBOI=1:size(BOI,2)
            NeedI=find(Fplot>=BOI(1,iBOI)&Fplot<=BOI(2,iBOI));

            %          for iCh=1:length(ChanEEGLab)
            if iFreqFunc==3
                PSDtemp=squeeze(FreqFunc{iFreqFunc}(LogPSD(IncludedSubj,NeedI,EEGchInd,iTrialType),[],2));
            else
                PSDtemp=squeeze(FreqFunc{iFreqFunc}(LogPSD(IncludedSubj,NeedI,EEGchInd,iTrialType),2));
            end
            Acctemp=Acc(IncludedSubj); % Behavior accuracy

            % Power and acuracy correlation
            [rSpear(:,iBOI,iTrialType),pSpear(:,iBOI,iTrialType)]=corr(PSDtemp,Acctemp,'type','spearman','rows','pairwise');

            %          end
            %% Visualize the data with topoplot in eeglab
            subplotLU(2,size(BOI,2),1,iBOI);
            topoplot(rSpear(:,iBOI,iTrialType), ChanEEGLab,'colormap',colorMapPN,'maplimits',[-1 1]);
            if iBOI==1
                yt=text(-0.55,0,'EEG-Behavior R','horizontalalignment','center','verticalalignment','bottom','fontsize',10,'rotation',90);
            end
            if iBOI==size(BOI,2)
                b=colorbar('southoutside');set(b,'position',[0.52 0.93 0.2 0.01],'xtick',[-1 1],'xticklabel',{'-1' '1'},'xlim',[-1 1]);
                xlabel(b,'PSD-Acc Correlation','verticalalignment','top')
            end


            subplotLU(2,size(BOI,2),2,iBOI);
            topoplot(log10(pSpear(:,iBOI,iTrialType)), ChanEEGLab,'colormap',colorMapPN,'maplimits',[-4 4]); % p-value
            xlabel(BandName{iBOI});
            text(0,-0.55,[BandName{iBOI} ' (' BandHzName{iBOI} ')'],'horizontalalignment','center','verticalalignment','top','fontsize',10)
            if iBOI==1
                a=ylabel('40Hz-BothControls');
                %                 set(a,'position',[0.01 0.5 0.03 0.4],'verticalalignment','middle')
                set(a,'verticalalignment','middle')
                yt=text(-0.55,0,'P EEG-Behavior','horizontalalignment','center','verticalalignment','bottom','fontsize',10,'rotation',90);
            end

            if iBOI==size(BOI,2)
                %  c=colorbar('southoutside');set(c,'position',[0.52 0.46 0.2 0.03],'xtick',[-4 0],'xticklabel',{'10e-4' '10e-0'},'xlim',[-4 0]);
                %  xlabel(c,'P values','verticalalignment','top')
                c=colorbar('southoutside');
                set(gca,'xtick',[],'ytick',[])
                set(c,'position',[0.52 0.51 0.2 0.01],'Limits',[PLim(1) 0],'ticks',[PLim(1) 0],'ticklabels',PLab);
                xlabel(c,'P values','verticalalignment','top')

            end

        end
        LuFontStandard;
        papersizePX=[0 0 6*size(BOI,2) 6*2+3];
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));

        saveas(gcf,[SaveTemp 'SpearmanPSDAcc'],'pdf');
        saveas(gcf,[SaveTemp 'SpearmanPSDAcc'],'png');
        saveas(gcf,[SaveTemp 'SpearmanPSDAcc.eps'],'epsc');
        close all

        %% Scatter plot with behavior - PSDAcc fig
        IncludedSubj=[SubjG{1}(:);SubjG{2}(:)];
        FlickerID=[zeros(size(SubjG{1}(:)))+1;zeros(size(SubjG{2}(:)))+2];
        Param.Corr='Spearman';  %%%Type of correlation, see Matlab function corr for more details
        Param.Pth=0.05; %%%threshold of Pvalue
        Param.ColorMap=colorMapPN; %%%Color map for correlation link
        Param.NodeColor=repmat([0.8 0.8 0.8],6,1); %%%Color of Nodes for correlation link plot
        Param.Clim=[-0.6 0.6];    %%%Color Limit for Correlation
        Param.Title='Pool All Sample'; %%Any title for label the figure
        Param.MarkerSize=8; %%%MarkerSize of scatter
        Param.SubjIDColor=FlickerColor;
        Param.SubjID=FlickerID;


        for iBOI=1:size(BOI,2)
            NeedI=find(Fplot>=BOI(1,iBOI)&Fplot<=BOI(2,iBOI));
            figure;

            %          for iCh=1:length(ChanEEGLab)
            if iFreqFunc==3
                PSDtemp=squeeze(FreqFunc{iFreqFunc}(LogPSD(IncludedSubj,NeedI,EEGchInd,iTrialType),[],2));
            else
                PSDtemp=squeeze(FreqFunc{iFreqFunc}(LogPSD(IncludedSubj,NeedI,EEGchInd,iTrialType),2));
            end

            Acctemp=Acc(IncludedSubj);

            tempName=EEGch;
            tempName{end+1}='Acc';
            multiCorr2GroupSubplot(6,6,[PSDtemp Acctemp],size(PSDtemp,2)+1,tempName,Param)
            LuFontStandard;

            papersizePX=[0 0 6*6 6*6];
            set(gcf, 'PaperUnits', 'centimeters');
            set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));

            saveas(gcf,[SaveTemp Param.Corr BandName{iBOI} BandHzName{iBOI} 'PSDAcc'],'pdf');
            saveas(gcf,[SaveTemp Param.Corr BandName{iBOI} BandHzName{iBOI} 'PSDAcc'],'png');
            saveas(gcf,[SaveTemp Param.Corr BandName{iBOI} BandHzName{iBOI} 'PSDAcc.eps'],'epsc');
            close all

        end

    end
    %
    %     close all
    %

end

FreqFunc{1}=@nanmean;
FreqFunc{2}=@nanmedian;
FreqFunc{3}=@nanmax;
FreqFuncNames={'mean','median','peak'};

diffTPmap=zeros(length(ChanEEGLab),length(ChanEEGLab),size(BOI,2),length(TrialType));
diffRPmap=diffTPmap;
diffTmap=diffTPmap;

rSpear=diffTPmap;
pSpear=diffTPmap;
% % TNodeTh=10;

COHEdgeTh=0.2;
COHNodeTh=0.01;
TEdgeTh=3;
TNodeTh=0.001;
pTEdgeTh=0.05;
pSpearEdgeTh=0.05;
rSpearEdgeTh=0.1;
rSpearNodeTh=0.001;

ScaleCOH=0.5;
ScaleT=1;
ScaleSpear=0.2;

close all

%% PSD all three groups together plot- Ranktest  fig2b
SubSavePSD=[SavePath 'PSD\'];
% SubSaveCOH=[SavePath 'COH\'];
ParamPSD.Ytick = [-8.1:3:-2.1];
ParamPSD.SigPlot='Ranktest';

saveDate = datestr(datetime, 'yy-mm-dd_HHMMSSFFF');

SubSavePSD=[SubSavePSD 'Ranktest_' saveDate '\'];
mkdir(SubSavePSD)


for iTrialType=3 %1:length(TrialType)
    SaveTemp=[SubSavePSD TrialTypeName{iTrialType} '\'];
    mkdir(SaveTemp)
    SubSaveFig=[SaveTemp 'Chan\'];
    mkdir(SubSaveFig)
    FBand=[1 55];

    %% HzAllChPSD40Random
    figure;
    for iCh=1:length(EEGchInd)
        clear DataPlot
        for iStimGroup=1:length(SubjG)
            DataPlot{iStimGroup}= squeeze(LogPSD(SubjG{iStimGroup},:,EEGchInd(iCh),iTrialType));
            Invalid=isnan(DataPlot{iStimGroup}(:,1));
            DataPlot{iStimGroup}(Invalid,:)=[];
        end
        if isempty(DataPlot{1})||isempty(DataPlot{2})
            continue;
        end
        %          subplotLUpage(6,6,iCh);
        ParamPSD.PathSave=[SaveTemp '40LightRandCh' EEGch{iCh}];

        % figure;
        % subplot('Position',[0.1 0.1 0.88 0.88])
        [~,FlickComStatis{iTrialType,iCh}]=RateHist_GroupPlot(Fplot,DataPlot,FlickerColor,ParamPSD); %plot and stats,
        text(27.5,-2,EEGch{iCh}); % text(27.5,0.3,EEGch{iCh}); %pre 6/14/24
        set(gca,'xlim',FBand,'xtick',[1 4 8 13 30 37 39 41 43 55],'ylim',[-7 -2],'ytick',[-7 -6 -5 -4 -3 -2]);
        xlabel('Frequency Hz')
        ylabel('Normalized Power (Log)')
        ax=gca;
        ax.XGrid = 'on';
        ax.YGrid = 'off';
        LuFontStandard
        papersizePX=[0 0 12 12];
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
        saveas(gcf,[SubSaveFig 'ThreeGroup_' EEGch{iCh}],'tiff');
        saveas(gcf,[SubSaveFig 'ThreeGroup_' EEGch{iCh}],'png');
        saveas(gcf,[SubSaveFig 'ThreeGroup_' EEGch{iCh}],'svg');
        saveas(gcf,[SubSaveFig 'ThreeGroup_' EEGch{iCh} '.eps'],'epsc');

        close all
    end
    %
    papersizePX=[0 0 6*6 6*6];
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
    saveas(gcf,[SaveTemp num2str(FBand(1)) '-' num2str(FBand(2)) 'HzAllChPSD40LightRand' TrialTypeName{iTrialType}],'pdf');
    saveas(gcf,[SaveTemp num2str(FBand(1)) '-' num2str(FBand(2)) 'HzAllChPSD40LightRand' TrialTypeName{iTrialType}],'png');
    saveas(gcf,[SaveTemp num2str(FBand(1)) '-' num2str(FBand(2)) 'HzAllChPSD40LightRand' TrialTypeName{iTrialType}],'svg');
    saveas(gcf,[SaveTemp num2str(FBand(1)) '-' num2str(FBand(2)) 'HzAllChPSD40LightRand' TrialTypeName{iTrialType} '.eps'],'epsc');


end
close all
ParamPSD.Ytick=[-8:4:0]; % reset this param to original value to not mess up other functions: ParamPSD.Ytick=[-8:4:0];

%% Set-up FOOOF save folder
fooof_starttime = datestr(datetime, 'yy-mm-dd_HHMMSSFFF');
fooof_save_folder=['FOOOF Results' '\' fooof_starttime '\'];
mkdir(fooof_save_folder)
%% Plot single FOOOF PSD
% set up fooof input variables:
allFreqs = Fplot; % = row vector of frequency values
iStimGroup = 3; % 1=40; 2=Light, 3=Random;
trialtype = 3; % 3=HIT&MISS
singleFOOOFEEGCh = 1;
subjectNumInStimGroup = 1;
input_power_spectrum = PSDall(SubjG{iStimGroup}(subjectNumInStimGroup),:,singleFOOOFEEGCh,trialtype); % size(PSDall) = 67 197 40 3
% PSDall(participants, frequencies, channels, trialtypes)

f_range = [2 55]; % f_range = fitting range - !!different upper ranges yield different results
% set default settings values  - % settings = fooof model settings, in a struct
settings = struct(...
    'peak_width_limits', [0.5, 12], ...
    'max_n_peaks', Inf, ...
    'min_peak_height', 0.0, ...
    'peak_threshold', 2.0, ...
    'aperiodic_mode', 'fixed', ...
    'verbose', true);
return_model = 1; % return_model = boolean of whether to return the FOOOF model fit, optional

% run fooof:
firstsubjectfirstCh = fooof(Fplot, input_power_spectrum, f_range, settings, return_model);

% Plot FOOOF results:
% Extract data from the results structure
fooof_freqs = firstsubjectfirstCh.freqs; % Frequency values
power_spectrum = firstsubjectfirstCh.power_spectrum; % Original power spectrum
% fooofed_spectrum = firstsubjectfirstCh.fooofed_spectrum; % Full FOOOF fit
ap_fit = firstsubjectfirstCh.ap_fit; % Aperiodic fit
difference_spectrum = power_spectrum - ap_fit;
% fooofed_difference_spectrum = fooofed_spectrum - ap_fit;


% Plot settings
figure;
hold on;

% Plot the original power spectrum
plot(fooof_freqs, power_spectrum, 'k', 'LineWidth', 1.5, 'DisplayName', 'Power Spectrum');

% Plot the FOOOFed spectrum (from fooof results)
% plot(fooof_freqs, fooofed_spectrum, 'r', 'LineWidth', 1.5, 'DisplayName', 'FOOOFed Spectrum');

% Plot the aperiodic fit
plot(fooof_freqs, ap_fit, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Aperiodic Fit');

% Plot the difference spectrum (subtracted in MATLAB)
plot(fooof_freqs, difference_spectrum, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Adjusted Spectrum (PS-ApFit)');

% % Plot the FOOOFed difference spectrum (subtracted in MATLAB)
% plot(fooof_freqs, fooofed_difference_spectrum, '--', 'LineWidth', 1.5, 'DisplayName', 'fooofedSpect minus apfit');

% Set logarithmic scale for frequency
% set(gca, 'XScale', 'log'); % Logarithmic x-axis
% set(gca, 'YScale', 'log'); % Logarithmic y-axis

% Add labels, legend, and grid
xlabel('Frequency (Hz)');
ylabel('Power');
legend_handle = legend('show');
% legend('show');
grid on;

% Position the legend in the middle right of the plot
set(legend_handle, 'Location', 'east');

fooof_figname = ['S' num2str(SubjG{iStimGroup}(subjectNumInStimGroup)) '_FOOOF_' EEGch{singleFOOOFEEGCh} '_' num2str(f_range(1)) '-' num2str(f_range(2)) ' Hz'];
title(['S' num2str(SubjG{iStimGroup}(subjectNumInStimGroup)) ' (' GroupName{iStimGroup} ') - Ch:' EEGch{singleFOOOFEEGCh} ' - '  num2str(f_range(1)) '-' num2str(f_range(2)) ' Hz']);
hold off;

% Save the figure as a single file
saveas(gcf, [fooof_save_folder fooof_figname '.png']);

%% FOOOF PSD - prep for fig2c
% set up fooof input variables:
allFreqs = Fplot; % = row vector of frequency values
f_range = [2 45]; % f_range = fitting range
% set default settings values - % settings = fooof model settings, in a struct
settings = struct(...
    'peak_width_limits', [0.5, 12], ...
    'max_n_peaks', Inf, ...
    'min_peak_height', 0.0, ...
    'peak_threshold', 2.0, ...
    'aperiodic_mode', 'fixed', ...
    'verbose', true);
return_model = 1; % return_model = boolean of whether to return the FOOOF model fit, optional
trialType = 3;
input_power_spectrum = PSDall(SubjG{iStimGroup}(1),:,1,trialType);
dummy_fooof = fooof(Fplot, input_power_spectrum, f_range, settings, return_model);
dummy_fooof.difference_spectrum = dummy_fooof.power_spectrum - dummy_fooof.ap_fit;

nSubjPSDs = length(PSDall(:,1,1,3));  % effectively gets the total number of EEG subjects that have a PSD
nCh = length(ChanEEGLab);
emptyStruct = struct(); % Create an empty structure
clear allFooofResults allFoooFDiffPSD
allFooofResults = repmat(dummy_fooof, nSubjPSDs, nCh); % Replicate the empty structure 32 times


for iSub = 1:nSubjPSDs
    for iCh=1:length(ChanEEGLab)
        input_power_spectrum = PSDall(iSub,:,iCh,3); % power_spectrum = row vector of power values
        % PSDall(participants, frequencies, channels, trialtypes)

        % run fooof:
        iSubiChFoofResults = fooof(Fplot, input_power_spectrum, f_range, settings, return_model);
        % Extract data from the results structure
        fooof_freqs = iSubiChFoofResults.freqs; % Frequency values
        power_spectrum = iSubiChFoofResults.power_spectrum; % Original power spectrum
        fooofed_spectrum = iSubiChFoofResults.fooofed_spectrum; % Full FOOOF fit
        ap_fit = iSubiChFoofResults.ap_fit; % Aperiodic fit
        iSubiChFoofResults.difference_spectrum = power_spectrum - ap_fit;
        allFooofResults(iCh) = iSubiChFoofResults;
        allFooofDiffPSD(iSub,:,iCh,trialType) = iSubiChFoofResults.difference_spectrum;
    end
end
%% Plot average FOOOF results:
% % Extract data from the results structure
% fooof_freqs = iChFooofResults.freqs; % Frequency values
% power_spectrum = iChFooofResults.power_spectrum; % Original power spectrum
% fooofed_spectrum = iChFooofResults.fooofed_spectrum; % Full FOOOF fit
% ap_fit = iChFooofResults.ap_fit; % Aperiodic fit
% difference_spectrum = power_spectrum - ap_fit;
% % fooofed_difference_spectrum = fooofed_spectrum - ap_fit;
%
% % Plot settings
% figure;
% hold on;
%
% % Plot the original power spectrum
% plot(fooof_freqs, power_spectrum, 'k', 'LineWidth', 1.5, 'DisplayName', 'Power Spectrum');
%
% % Plot the FOOOFed spectrum (from fooof results)
% plot(fooof_freqs, fooofed_spectrum, 'r', 'LineWidth', 1.5, 'DisplayName', 'FOOOFed Spectrum');
%
% % Plot the aperiodic fit
% plot(fooof_freqs, ap_fit, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Aperiodic Fit');
%
% % Plot the difference spectrum (subtracted in MATLAB)
% plot(fooof_freqs, difference_spectrum, 'g--', 'LineWidth', 1.5, 'DisplayName', 'powerSpect minus apfit');
%
% % % Plot the FOOOFed difference spectrum (subtracted in MATLAB)
% % plot(fooof_freqs, fooofed_difference_spectrum, '--', 'LineWidth', 1.5, 'DisplayName', 'fooofedSpect minus apfit');
%
% % Set logarithmic scale for frequency
% % set(gca, 'XScale', 'log'); % Logarithmic x-axis
% % set(gca, 'YScale', 'log'); % Logarithmic y-axis
%
% % Add labels, legend, and grid
% xlabel('Frequency (Hz)');
% ylabel('Power');
% legend('show');
% grid on;
% title(['FOOOF Analysis Results: ' num2str(f_range(1)) '-' num2str(f_range(2)) ' Hz']);
% hold off;

%%
% Predefine the number of stimulation groups, subjects, and channels
nStimGroups = 3;
nFreqs = length(fooof_freqs);
nCh = 32;

% Preallocate for average difference spectra per channel and stim group
avgDiffSpectraChannels = zeros(nStimGroups, nFreqs, nCh);

% Compute the average difference spectrum for each channel and stim group
for iStimGroup = 1:nStimGroups
    for iCh = 1:nCh
        % Extract difference spectra for the current stim group and channel
        diffSpectraGroup = allFooofDiffPSD(SubjG{iStimGroup},:,iCh,trialtype); % Subj x Freq
        % Average across subjects
        avgDiffSpectraChannels(iStimGroup, :, iCh) = mean(diffSpectraGroup, 1, 'omitnan');
    end
end

% Plot the average difference spectra for each channel, one plot per stim group
for iStimGroup = 1:nStimGroups
    figure;
    hold on;
    colors = lines(nCh); % Generate distinct colors for each channel
    for iCh = 1:nCh
        plot(fooof_freqs, avgDiffSpectraChannels(iStimGroup, :, iCh), 'LineWidth', 1.5, ...
            'DisplayName', EEGch{iCh}, 'Color', colors(iCh, :));
    end

    % Customize the plot
    xlabel('Frequency (Hz)');
    ylabel('Difference Spectrum Power');
    title(['Average Difference Spectrum by Channel - Stim Group ' GroupName{iStimGroup}]);
    legend('show', 'Location', 'eastoutside');
    grid on;
    hold off;

    % Save the figure
    saveas(gcf, [fooof_save_folder 'Avg_Diff_Spectrum_by_Channel_StimGroup' GroupName{iStimGroup} '.png']);
end

%% FOOOFed Spectrum Plot + settings
figure;
hold on;

% Plot the original power spectrum
plot(fooof_freqs, power_spectrum, 'k', 'LineWidth', 1.5, 'DisplayName', 'Power Spectrum');

% Plot the FOOOFed spectrum (from fooof results)
% plot(fooof_freqs, fooofed_spectrum, 'r', 'LineWidth', 1.5, 'DisplayName', 'FOOOFed Spectrum');

% Plot the aperiodic fit
plot(fooof_freqs, ap_fit, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Aperiodic Fit');

% Plot the difference spectrum (subtracted in MATLAB)
plot(fooof_freqs, difference_spectrum, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Adjusted Spectrum (PS-ApFit)');

% % Plot the FOOOFed difference spectrum (subtracted in MATLAB)
% plot(fooof_freqs, fooofed_difference_spectrum, '--', 'LineWidth', 1.5, 'DisplayName', 'fooofedSpect minus apfit');

% Set logarithmic scale for frequency
% set(gca, 'XScale', 'log'); % Logarithmic x-axis
% set(gca, 'YScale', 'log'); % Logarithmic y-axis

% Add labels, legend, and grid
xlabel('Frequency (Hz)');
ylabel('Power');
legend_handle = legend('show');
% legend('show');
grid on;

% Position the legend in the middle right of the plot
set(legend_handle, 'Location', 'east');

fooof_figname = ['S' num2str(SubjG{iStimGroup}(subjectNumInStimGroup)) '_FOOOF_' EEGch{singleFOOOFEEGCh} '_' num2str(f_range(1)) '-' num2str(f_range(2)) ' Hz'];
title(['S' num2str(SubjG{iStimGroup}(subjectNumInStimGroup)) ' (' GroupName{iStimGroup} ') - Ch:' EEGch{singleFOOOFEEGCh} ' - '  num2str(f_range(1)) '-' num2str(f_range(2)) ' Hz']);
hold off;

% Save the figure as a single file
saveas(gcf, [fooof_save_folder fooof_figname '.png']);
%% PSD topographical Plots on Brain (heatmaps) - group difference, PSD-Acc correlation brain *** fig2c fig2d fig3
PowerLim=[-7 -2];
PowerLab={'-7' '-2'};
TLim=[-5 5];
TLab={'-5' '5'};
PLim=[-4 4];
PLab={'10e-4' '10e-0'};
fooof_freqs = firstsubjectfirstCh.freqs; % Frequency values

for iFreqFunc=3 %1:length(FreqFuncNames)
    clear MapGroup diffTmap diffMap Tdata;
    clear diffTPmap diffRPmap diffTmap

    for iTrialType=3%1:length(TrialType)
        % iCom=1;
        PSDsaveDate = datestr(datetime, 'yy-mm-dd_HHMMSSFFF');

        SaveTemp=['Fig2Results\' PSDsaveDate '\'];
        SaveTemp=[SaveTemp FreqFuncNames{iFreqFunc} '\' ];
        mkdir(SaveTemp)

        BOI=[5;15];
        BOI=[1 4 8 13 30 39.5 43;4 8 13 30 37 41.5 100];
        BandName={'Delta','Theta','Alpha','Beta','Gamma-1','Gamma-E','Gamma-2'};
        BandHzName={'1-4 Hz','4-8 Hz','8-13Hz','13-30 Hz','30-37 Hz','39-41 Hz','43-100 Hz'};
        %
        %% PSD on brain ***old Fig2C
        % figure;
        for iBOI=1:size(BOI,2)
            NeedI=find(Fplot>=BOI(1,iBOI)&Fplot<=BOI(2,iBOI));

            for iStimGroup=1:length(SubjG)
                if iFreqFunc==3
                    Tdata{iStimGroup,iBOI,iTrialType}=squeeze(FreqFunc{iFreqFunc}(LogPSD(SubjG{iStimGroup},NeedI,EEGchInd,iTrialType),[],2));
                    MapGroup{iStimGroup,iBOI,iTrialType}=squeeze(nanmean(FreqFunc{iFreqFunc}(LogPSD(SubjG{iStimGroup},NeedI,EEGchInd,iTrialType),[],2),1));

                else
                    Tdata{iStimGroup,iBOI,iTrialType}=squeeze(FreqFunc{iFreqFunc}(LogPSD(SubjG{iStimGroup},NeedI,EEGchInd,iTrialType),2));
                    MapGroup{iStimGroup,iBOI,iTrialType}=squeeze(nanmean(FreqFunc{iFreqFunc}(LogPSD(SubjG{iStimGroup},NeedI,EEGchInd,iTrialType),2),1));
                end
                subplotLU(2,size(BOI,2),iStimGroup,iBOI); %
                [~,~,~,xmesh,ymesh]=topoplot(MapGroup{iStimGroup,iBOI,iTrialType}, ChanEEGLab,'colormap',jet,'maplimits',PowerLim);
                if iStimGroup==2
                    xlabel(BandName{iBOI});
                    text(0,-0.55,[BandName{iBOI} ' (' BandHzName{iBOI} ')'],'horizontalalignment','center','verticalalignment','top','fontsize',10)

                end
                if iBOI==1
                    ylabel(GroupName{iStimGroup})
                    yt=text(-0.55,0,GroupName{iStimGroup},'horizontalalignment','center','verticalalignment','bottom','fontsize',10,'rotation',90);
                end
            end
        end
        % subplot('position',[0.5 0.51 0.3 0.01]);
        % b=colorbar('southoutside');
        % set(gca,'xtick',[],'ytick',[])
        % set(b,'position',[0.5 0.5 0.3 0.03],'Limits',[0 1],'Ticks',[0 1],'Ticklabels',PowerLab);
        % xlabel(b,'Log Normalized Power')
        % LuFontStandard;
        % papersizePX=[0 0 6*size(BOI,2) 6*2+3];
        % set(gcf, 'PaperUnits', 'centimeters');
        % set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
        % saveas(gcf,[SaveTemp 'GroupPSDonBrain'],'pdf');
        % saveas(gcf,[SaveTemp 'GroupPSDonBrain'],'png');
        % saveas(gcf,[SaveTemp 'GroupPSDonBrain.eps'],'epsc');
        %% FOOOF PSD on brain *** old Fig2C
        % % Ensure the desired colormap is loaded
        % % batlow is part of the cmocean or scientific colormaps package
        % if exist('batlow', 'file') == 2
        %     cmap = batlow;  % Load batlow if available
        % else
        %     cmap = parula;  % Fallback to parula if batlow isn't available
        % end
        % 
        % figure;
        % BOI=[1 4 8 13 30 39.5 43;4 8 13 30 37 41.5 100];
        % BandName={'Delta','Theta','Alpha','Beta','Gamma-1','Gamma-E','Gamma-2'};
        % BandHzName={'1-4 Hz','4-8 Hz','8-13Hz','13-30 Hz','30-37 Hz','39-41 Hz','43-100 Hz'};
        % PowerLim = [-0.5 1];
        % for iBOI=1:size(BOI,2)
        %     % NeedI=find(Fplot>=BOI(1,iBOI)&Fplot<=BOI(2,iBOI));
        %     NeedI=find(fooof_freqs>=BOI(1,iBOI)&fooof_freqs<=BOI(2,iBOI))';
        % 
        %     for iStimGroup=1:length(SubjG)
        %         if iFreqFunc==3
        %             FooofTdata{iStimGroup,iBOI,iTrialType}=squeeze(FreqFunc{iFreqFunc}(allFooofDiffPSD(SubjG{iStimGroup},NeedI,EEGchInd,iTrialType),[],2));
        %             FooofMapGroup{iStimGroup,iBOI,iTrialType}=squeeze(nanmean(FreqFunc{iFreqFunc}(allFooofDiffPSD(SubjG{iStimGroup},NeedI,EEGchInd,iTrialType),[],2),1));
        % 
        %         else
        %             FooofTdata{iStimGroup,iBOI,iTrialType}=squeeze(FreqFunc{iFreqFunc}(allFooofDiffPSD(SubjG{iStimGroup},NeedI,EEGchInd,iTrialType),2));
        %             FooofMapGroup{iStimGroup,iBOI,iTrialType}=squeeze(nanmean(FreqFunc{iFreqFunc}(allFooofDiffPSD(SubjG{iStimGroup},NeedI,EEGchInd,iTrialType),2),1));
        %         end
        %         subplotLU(3,size(BOI,2),iStimGroup,iBOI); %
        %         [~,~,~,xmesh,ymesh]=topoplot(FooofMapGroup{iStimGroup,iBOI,iTrialType}, ChanEEGLab,'colormap',parula,'maplimits',PowerLim);
        %         if iStimGroup==3
        %             xlabel(BandName{iBOI});
        %             text(0,-0.55,[BandName{iBOI} ' (' BandHzName{iBOI} ')'],'horizontalalignment','center','verticalalignment','top','fontsize',10)
        % 
        %         end
        %         if iBOI==1
        %             ylabel(GroupName{iStimGroup})
        %             yt=text(-0.55,0,GroupName{iStimGroup},'horizontalalignment','center','verticalalignment','bottom','fontsize',10,'rotation',90);
        %         end
        %     end
        % end
        % % subplot('position',[0.5 0.51 0.3 0.01]);
        % b=colorbar('southoutside');
        % set(gca,'xtick',[],'ytick',[])
        % % set(b,'position',[0.5 0.5 0.3 0.03],'Limits',[0 1],'Ticks',[0 1],'Ticklabels',PowerLab);
        % xlabel(b,'(not Log Normalized) Power')
        % LuFontStandard;
        % papersizePX=[0 0 6*size(BOI,2) 6*2+3];
        % set(gcf, 'PaperUnits', 'centimeters');
        % set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
        % % saveas(gcf,[SaveTemp 'FOOOFPSD'],'pdf');
        % saveas(gcf,[SaveTemp 'FOOOFPSD_oldColorbarRange'],'png');
        % % saveas(gcf,[SaveTemp 'PSD40BothControls.eps'],'epsc');
        %% Updated FOOOF color map 1/22/25 *** Fig2C
        % % Compute global 10% minimum and 90% maximum across all data
        % allData = [];  % Initialize an empty array to collect all data values
        % for iBOI = 1:size(BOI, 2)
        %     for iStimGroup = 1:length(SubjG)
        %         % Collect all FooofMapGroup data into a single array
        %         allData = [allData; FooofMapGroup{iStimGroup, iBOI, iTrialType}(:)];
        %     end
        % end
        % 
        % % Compute 10% and 90% percentiles
        % cbarMin = prctile(allData, 10);
        % cbarMax = prctile(allData, 90);
        % 
        % % Update PowerLim based on the computed values
        % FOOOFPowerLim = [cbarMin, cbarMax];
        % 
        % figure;
        % for iBOI = 1:size(BOI, 2)
        %     NeedI = find(fooof_freqs >= BOI(1, iBOI) & fooof_freqs <= BOI(2, iBOI))';
        % 
        %     for iStimGroup = 1:3 %length(SubjG)
        %         subplotLU(3, size(BOI, 2), iStimGroup, iBOI);
        %         [~, ~, ~, xmesh, ymesh] = topoplot(FooofMapGroup{iStimGroup, iBOI, iTrialType}, ChanEEGLab, ...
        %             'colormap', parula, 'maplimits', FOOOFPowerLim);
        % 
        %         if iStimGroup == 3
        %             xlabel(BandName{iBOI});
        %             text(0, -0.55, [BandName{iBOI} ' (' BandHzName{iBOI} ')'], ...
        %                 'horizontalalignment', 'center', 'verticalalignment', 'top', 'fontsize', 10);
        %         end
        %         if iBOI == 1
        %             ylabel(GroupName{iStimGroup});
        %             text(-0.55, 0, GroupName{iStimGroup}, ...
        %                 'horizontalalignment', 'center', 'verticalalignment', 'bottom', ...
        %                 'fontsize', 10, 'rotation', 90);
        %         end
        %     end
        % end
        % 
        % % b = colorbar('southoutside');
        % % caxis(FOOOFPowerLim);  % Set colorbar limits to match PowerLim
        % % xlabel(b, '(not Log Normalized) Power');
        % 
        % LuFontStandard;
        % papersizePX = [0 0 6 * size(BOI, 2) 6 * 2 + 3];
        % set(gcf, 'PaperUnits', 'centimeters');
        % set(gcf, 'PaperPosition', papersizePX, 'PaperSize', papersizePX(3:4));
        % saveas(gcf, [SaveTemp 'FOOOFPSD_updatedColorbarRange'], 'png');
        % 
        %% Plot and save FOOOF vertical colorbar separately *** Fig2C
        % figure;
        % 
        % % Create a dummy image to generate a colorbar
        % imagesc([0 1]);  % Placeholder data
        % colormap(parula);  % Use the same colormap
        % caxis(FOOOFPowerLim);   % Apply the same color limits
        % 
        % % Customize colorbar
        % b = colorbar('eastoutside');  % Set colorbar orientation to vertical
        % 
        % % Set ticks and format tick labels with 2 significant figures
        % tickValues = linspace(FOOOFPowerLim(1), FOOOFPowerLim(2), 5);  % Generate 5 evenly spaced tick values
        % set(b, 'Ticks', tickValues);  % Set tick positions
        % set(b, 'TickLabels', arrayfun(@(x) sprintf('%.2g', x), tickValues, 'UniformOutput', false));  % Format tick labels
        % 
        % % Add label to the colorbar
        % ylabel(b, '(not Log Normalized) Power', 'fontsize', 12, ...
        %     'rotation', 270, 'VerticalAlignment', 'bottom', ...
        %     'HorizontalAlignment', 'center');
        % 
        % % Adjust figure layout to fit the colorbar and label
        % set(gca, 'Visible', 'off');  % Hide axes
        % set(gcf, 'PaperUnits', 'centimeters');
        % set(gcf, 'PaperPosition', [0 0 2 10]);  % Adjust to fit vertical colorbar
        % set(gcf, 'PaperSize', [2 10]);
        % 
        % % Save colorbar as a PNG file
        % saveas(gcf, [SaveTemp 'FOOOF_Colorbar'], 'png');

        %% Group differences on brain *** old fig2d ? suppfig2
        figure;
        for iBOI=1:size(BOI,2)

            for iCh=1:length(ChanEEGLab)
                [~,diffTPmap(iCh,iBOI,iTrialType),~,stats]=ttest2(Tdata{2,iBOI,iTrialType}(:,iCh),Tdata{1,iBOI,iTrialType}(:,iCh));
                [diffRPmap(iCh,iBOI,iTrialType),~,~]=ranksum(Tdata{2,iBOI,iTrialType}(:,iCh),Tdata{1,iBOI,iTrialType}(:,iCh));
                diffTmap(iCh,iBOI,iTrialType)=stats.tstat;
            end
            subplotLU(2,size(BOI,2),1,iBOI);
            topoplot(diffTmap(:,iBOI,iTrialType), ChanEEGLab,'colormap',colorMapPN,'maplimits',TLim);
            if iBOI==1
                yt=text(-0.55,0,['T, ' GroupName{Group1} '-' GroupName{Group2}],'horizontalalignment','center','verticalalignment','bottom','fontsize',10,'rotation',90);
            end
            if iBOI==size(BOI,2)
                b=colorbar('southoutside');
                %                  set(b,'position',[0.52 0.93 0.2 0.03],'xtick',[-6 6],'xticklabel',{'-6' '6'},'xlim',[-6 6]);
                set(gca,'xtick',[],'ytick',[])
                set(b,'position',[0.52 0.93 0.2 0.01],'ticks',TLim,'ticklabels',TLab);

                xlabel(b,'T statistics','verticalalignment','top')
            end


            subplotLU(2,size(BOI,2),2,iBOI);
            topoplot(log10(diffTPmap(:,iBOI,iTrialType)), ChanEEGLab,'colormap',colorMapPN,'maplimits',PLim);
            xlabel(BandName{iBOI});
            text(0,-0.55,[BandName{iBOI} ' (' BandHzName{iBOI} ')'],'horizontalalignment','center','verticalalignment','top','fontsize',10)
            if iBOI==1
                a=ylabel('40Hz-BothControls');
                %                 a.Position=[0.01 0.5 0.03 0.4];
                %                 a.verticalalignment='middle';
                %                 set(a,'Position',[0.01 0.5 0.03 0.4],'Verticalalignment','middle')
                set(a,'Verticalalignment','middle')
                yt=text(-0.55,0,['P, ' GroupName{Group1} '-' GroupName{Group2}],'horizontalalignment','center','verticalalignment','bottom','fontsize',10,'rotation',90);
            end

            if iBOI==size(BOI,2)
                c=colorbar('southoutside');
                %                     set(c,'position',[0.52 0.46 0.2 0.03],'xtick',[-4 0],'xticklabel',{'10e-4' '10e-0'},'xlim',[-4 0]);
                set(gca,'xtick',[],'ytick',[])
                set(c,'position',[0.52 0.51 0.2 0.01],'Limits',[PLim(1) 0],'ticks',[PLim(1) 0],'ticklabels',PLab);
                xlabel(c,'P values','verticalalignment','top')
            end
        end

        LuFontStandard;
        papersizePX=[0 0 6*size(BOI,2) 6*2+3];
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));

        saveas(gcf,[SaveTemp 'PSDDiff40BothControls'],'pdf');
        saveas(gcf,[SaveTemp 'PSDDiff40BothControls'],'png');
        saveas(gcf,[SaveTemp 'PSDDiff40BothControls.eps'],'epsc');

        %% Group differences on brain - Calculate t-test stat and plot *** fig2d - suppfig2
        % Define group pairs to compare % see: open: GroupName
        groupPairs = [
            1, 2;  % First pair: Group1 = 1 (40 Hz), Group2 = 2 (Light)
            1, 3   % Second pair: Group1 = 1, Group2 = 3 (Random)
            ];

        % Channel subset Fp1, Cz, Oz
        % Channels of interest (indices)
        PSDttestChSubset = [1, 32, 16]; % [1, 32, 16] = [Fp1, Cz, Oz]
        
        % Initialize a matrix to store p-values for FDR correction
        allPValues = [];

        % Loop through each group pair
        for iPair = 1:size(groupPairs, 1)
            Group1 = groupPairs(iPair, 1);
            Group2 = groupPairs(iPair, 2);
            dataName = [TrialTypeName{iTrialType} FreqFuncNames{iFreqFunc} GroupName{Group1} GroupName{Group2}];

            figure;
            
            for iBOI=1:size(BOI,2)
                %% calculate stats test
                for iCh=1:length(ChanEEGLab)
                    [~,diffTPmap(iCh,iBOI,iTrialType),~,stats]=ttest2(Tdata{Group1,iBOI,iTrialType}(:,iCh),Tdata{Group2,iBOI,iTrialType}(:,iCh)); % two-sided ttest
                    [diffRPmap(iCh,iBOI,iTrialType),~,~]=ranksum(Tdata{Group1,iBOI,iTrialType}(:,iCh),Tdata{Group2,iBOI,iTrialType}(:,iCh));
                    diffTmap(iCh,iBOI,iTrialType)=stats.tstat;
                end

                % Collect p-values for channels of interest
                for chIdx = 1:length(PSDttestChSubset)
                    iCh = PSDttestChSubset(chIdx);
                    % allPValues = [p-value, t-stat, groupPair, BOI, Ch]
                    allPValues(end+1, :) = [diffTPmap(iCh, iBOI, iTrialType), diffTmap(iCh, iBOI, iTrialType), iPair, iBOI, iCh]; % Store p-value with metadata
                end

                %% Plot T-stat
                subplotLU(2,size(BOI,2),1,iBOI);
                topoplot(diffTmap(:,iBOI,iTrialType), ChanEEGLab,'colormap',colorMapPN,'maplimits',TLim);
                if iBOI==1
                    yt=text(-0.55,0,['T, ' GroupName{Group1} '-' GroupName{Group2}],'horizontalalignment','center','verticalalignment','bottom','fontsize',10,'rotation',90);
                end
                if iBOI==size(BOI,2)
                    b=colorbar('southoutside');
                    %                  set(b,'position',[0.52 0.93 0.2 0.03],'xtick',[-6 6],'xticklabel',{'-6' '6'},'xlim',[-6 6]);
                    set(gca,'xtick',[],'ytick',[])
                    set(b,'position',[0.52 0.93 0.2 0.01],'ticks',TLim,'ticklabels',TLab);

                    xlabel(b,'T statistics','verticalalignment','top')
                end

                %% Plot p-value
                subplotLU(2,size(BOI,2),2,iBOI);
                topoplot(log10(diffTPmap(:,iBOI,iTrialType)), ChanEEGLab,'colormap',colorMapPN,'maplimits',PLim);
                xlabel(BandName{iBOI});
                text(0,-0.55,[BandName{iBOI} ' (' BandHzName{iBOI} ')'],'horizontalalignment','center','verticalalignment','top','fontsize',10)
                if iBOI==1
                    a=ylabel([GroupName{Group1} '-' GroupName{Group2}]);
                    %                 a.Position=[0.01 0.5 0.03 0.4];
                    %                 a.verticalalignment='middle';
                    %                 set(a,'Position',[0.01 0.5 0.03 0.4],'Verticalalignment','middle')
                    set(a,'Verticalalignment','middle')

                    yt=text(-0.55,0,['P, ' GroupName{Group1} '-' GroupName{Group2}],'horizontalalignment','center','verticalalignment','bottom','fontsize',10,'rotation',90);
                end

                if iBOI==size(BOI,2)
                    c=colorbar('southoutside');
                    %                     set(c,'position',[0.52 0.46 0.2 0.03],'xtick',[-4 0],'xticklabel',{'10e-4' '10e-0'},'xlim',[-4 0]);
                    set(gca,'xtick',[],'ytick',[])
                    set(c,'position',[0.52 0.51 0.2 0.01],'Limits',[PLim(1) 0],'ticks',[PLim(1) 0],'ticklabels',PLab);
                    xlabel(c,'P values','verticalalignment','top')
                end
            end

            LuFontStandard;
            papersizePX=[0 0 6*size(BOI,2) 6*2+3];
            set(gcf, 'PaperUnits', 'centimeters');
            set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));

            saveas(gcf,[SaveTemp 'PSDDiff' GroupName{Group1} '-' GroupName{Group2}],'pdf');
            saveas(gcf,[SaveTemp 'PSDDiff' GroupName{Group1} '-' GroupName{Group2}],'png');
            saveas(gcf,[SaveTemp 'PSDDiff' GroupName{Group1} '-' GroupName{Group2} '.eps'],'epsc');
        end
        %% PSD ttest FDR Correction fig2d stats
        % Filter p-values for the first 5 bands
        numBandsToInclude = 5; % Only include the first 5 bands

        % allPValues = [p-value, t-stat, groupPiar, BOI, Ch]
        filteredPValues = allPValues(allPValues(:, 4) <= numBandsToInclude, :);

        % Extract only the raw p-values for FDR correction
        pValuesForFDR = filteredPValues(:, 1);

        % % Perform FDR Correction (three different methods below)
        BHFDR = mafdr(pValuesForFDR, 'BHFDR', true); % Benjamini-Hochberg FDR correction
        [FDRST,Q,aPrioriProb] = mafdr(pValuesForFDR); % Storey-Tibshirani method (2002)
        
        % fdh_bh, parameters - gets the critical p value
        q=0.05;
        method='pdep';
        report='yes';
        [fdrbh_h, fdrbh_crit_p, fdrbh_adj_p]=fdr_bh(pValuesForFDR,q,method,report);

        % Add FDR corrected p-values back to the filtered results
        filteredPValuesPlusFDR = [filteredPValues, fdrAdjustedPValues];

        % Display Results
        fprintf('Channel\tGroup Pair\tBand\tT-stat\tRaw P-value\tFDR Corrected P-value\n');
        for i = 1:size(filteredPValuesPlusFDR, 1)
            % Get channel name from the struct
            chName = ChanEEGLab(filteredPValuesPlusFDR(i, 5)).labels; % Access the 'labels' field of the struct

            % Extract group indices
            group1Idx = groupPairs(filteredPValuesPlusFDR(i, 3), 1);
            group2Idx = groupPairs(filteredPValuesPlusFDR(i, 3), 2);

            % Get group pair names
            groupPair = sprintf('%s-%s', GroupName{group1Idx}, GroupName{group2Idx});

            % Get band name
            band = BandName{filteredPValuesPlusFDR(i, 4)}; % Assuming BandName is a cell array
            
            % Get t-test stat
            tstat =  filteredPValuesPlusFDR(i, 2);               

            % Get raw and FDR-corrected p-values
            rawP = filteredPValuesPlusFDR(i, 1);
            fdrP = filteredPValuesPlusFDR(i, 6);

            % Print result
            fprintf('%s\t%s\t%s\t%.4f\t%.4f\t%.4f\n', chName, groupPair, band, tstat, rawP, fdrP);
        end

        % Display total number of comparisons
        numComparisons = size(filteredPValuesPlusFDR, 1);
        fprintf('Total number of comparisons used for FDR correction: %d\n', numComparisons);


        %% NEW PSD-ACC and PSD_RT 2025-02-04 MKA
        %% % Set-up for PSD-Beh - Spearman PSD-Acc and PSD-RT ***fig3b & fig3c
        SaveTemp = ['Fig3Results\' PSDsaveDate '\'];
        mkdir(SaveTemp)

        % PSD groups: 1=40Hz, 3=Random, 6=LightRT
        PSDBehaviorGroup1 = 1;
        PSDBehaviorGroup2 = 3;
        PSDBehaviorGroup3 = 6;
        allGroupNames = [GroupName{PSDBehaviorGroup1} GroupName{PSDBehaviorGroup2} GroupName{PSDBehaviorGroup3}];
        dataName = [TrialTypeName{iTrialType} FreqFuncNames{iFreqFunc} allGroupNames];
        IncludedSubj = union(SubjG{PSDBehaviorGroup1}, union(SubjG{PSDBehaviorGroup2}, SubjG{PSDBehaviorGroup3}));

        figstatsChSubset = PSDttestChSubset;

        % Create separate figures for PSD-Acc and PSD-RT
        figureAcc = figure;
        figureRT = figure;

        Acctemp = Acc(IncludedSubj);
        RTtemp = SubjsAvgRT(IncludedSubj);

        %%% Variables for scatterplot:
        FlickerID=[zeros(size(SubjG{PSDBehaviorGroup1}(:)))+1;zeros(size(SubjG{PSDBehaviorGroup2}(:)))+2;zeros(size(SubjG{PSDBehaviorGroup3}(:)))+3];
        Param.Corr='Spearman';  %%%Type of correlation, see Matlab function corr for more details
        Param.Pth=0.05; %%%threshold of Pvalue
        Param.EdgeColor=colorMapPN; %%%Color map for correlation link
        Param.NodeColor=repmat([0.8 0.8 0.8],6,1); %%%Color of Nodes for correlation link plot
        Param.Clim=[-1 1];    %%%Color Limit for Correlation
        Param.Title='Pool All Sample'; %%Any title for label the figure
        Param.MarkerSize=8; %%%MarkerSize of scatter
        Param.SubjIDColor=FlickerColor;
        Param.SubjID=FlickerID;

        %%% Loop thru BOI for Fig3b and Fig3c heatmaps & plot
        for iBOI = 1:size(BOI, 2)
            NeedI = find(Fplot >= BOI(1, iBOI) & Fplot <= BOI(2, iBOI));
            
            if iFreqFunc == 3
                PSDtemp = squeeze(FreqFunc{iFreqFunc}(LogPSD(IncludedSubj, NeedI, EEGchInd, iTrialType), [], 2));
            else
                PSDtemp = squeeze(FreqFunc{iFreqFunc}(LogPSD(IncludedSubj, NeedI, EEGchInd, iTrialType), 2));
            end

            %% Scatterplot of correlation - indiv behavior and PSD fig3a+
            %% % PSD_Acc scatter
            PSDAcc_scatterplot = figure;
            tempName=EEGch;
            tempName{end+1}='Acc';
            multiCorr2GroupSubplot(6,6,[PSDtemp Acctemp],size(PSDtemp,2)+1,tempName,Param)
            LuFontStandard;

            papersizePX=[0 0 6*6 6*6];
            set(gcf, 'PaperUnits', 'centimeters');
            set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));

            % saveas(gcf,[SaveTemp Param.Corr BName{iBOI} BName2{iBOI} 'PSDAcc_' allGroupNames],'pdf');
            saveas(gcf,[SaveTemp Param.Corr BName{iBOI} BName2{iBOI} 'PSDAcc_' allGroupNames],'png');
            saveas(gcf,[SaveTemp Param.Corr BName{iBOI} BName2{iBOI} 'PSDAcc_' allGroupNames '.eps'],'epsc');

           
            %% % PSD_RT scatter
            PSDRT_scatterplot = figure;
            tempName=EEGch;
            tempName{end+1}='RT';
            multiCorr2GroupSubplot(6,6,[PSDtemp RTtemp],size(PSDtemp,2)+1,tempName,Param)
            LuFontStandard;

            papersizePX=[0 0 6*6 6*6];
            set(gcf, 'PaperUnits', 'centimeters');
            set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));

            % saveas(gcf,[SaveTemp Param.Corr BName{iBOI} BName2{iBOI} 'PSDAcc_' allGroupNames],'pdf');
            saveas(gcf,[SaveTemp Param.Corr BName{iBOI} BName2{iBOI} 'PSDRT_' allGroupNames],'png');
            saveas(gcf,[SaveTemp Param.Corr BName{iBOI} BName2{iBOI} 'PSDRT_' allGroupNames '.eps'],'epsc');

            %% Plot PSD-Acc and PSD-RT: R values and P Values topoplots
            %%% Calculate correlations & pvalues 
            [PSDAcc_rSpear(:, iBOI, iTrialType), PSDAcc_pSpear(:, iBOI, iTrialType)] = corr(PSDtemp, Acctemp, 'type', 'spearman', 'rows', 'pairwise');
            [PSDRT_rSpear(:, iBOI, iTrialType), PSDRT_pSpear(:, iBOI, iTrialType)] = corr(PSDtemp, RTtemp, 'type', 'spearman', 'rows', 'pairwise');

            %% % Plot PSD-Acc Correlation
            figure(figureAcc);
            %%%% Plot R value correlation topoplot PSD-Acc(fig3b)
            subplotLU(2,size(BOI,2),1,iBOI);
            topoplot(PSDAcc_rSpear(:,iBOI,iTrialType), ChanEEGLab,'colormap',colorMapPN,'maplimits',[-1 1]); % Changed from 'maplimits',[-1 1]
            if iBOI==1
                yt=text(-0.55,0,'EEG-Behavior R','horizontalalignment','center','verticalalignment','bottom','fontsize',10,'rotation',90);
            end
            if iBOI==size(BOI,2)
                b=colorbar('southoutside');set(b,'position',[0.52 0.93 0.2 0.01],'xtick',[-1 1],'xticklabel',{'-1' '1'},'xlim',[-1 1]);
                xlabel(b,'PSD-Acc Correlation','verticalalignment','top')
            end

            %%%% Plot PSD-Acc P-value topoplot (fig3b supplement)
            subplotLU(2,size(BOI,2),2,iBOI);
            topoplot(log10(PSDAcc_pSpear(:,iBOI,iTrialType)), ChanEEGLab,'colormap',colorMapPN,'maplimits',[-4 4]);
            xlabel(BandName{iBOI});
            text(0,-0.55,[BandName{iBOI} ' (' BandHzName{iBOI} ')'],'horizontalalignment','center','verticalalignment','top','fontsize',10)
            if iBOI==1
                a=ylabel(allGroupNames);
                %                 set(a,'position',[0.01 0.5 0.03 0.4],'verticalalignment','middle')
                set(a,'verticalalignment','middle')
                yt=text(-0.55,0,'P EEG-Behavior','horizontalalignment','center','verticalalignment','bottom','fontsize',10,'rotation',90);
            end

            if iBOI==size(BOI,2)
                %                   c=colorbar('southoutside');set(c,'position',[0.52 0.46 0.2 0.03],'xtick',[-4 0],'xticklabel',{'10e-4' '10e-0'},'xlim',[-4 0]);
                %                   xlabel(c,'P values','verticalalignment','top')
                c=colorbar('southoutside');
                set(gca,'xtick',[],'ytick',[])
                set(c,'position',[0.52 0.51 0.2 0.01],'Limits',[PLim(1) 0],'ticks',[PLim(1) 0],'ticklabels',PLab);
                xlabel(c,'P values','verticalalignment','top')
            end

            %% % Plot PSD-RT Correlation
            figure(figureRT);
            %%%% Plot PSD-RT R value (correlation) heatmap (fig3c)
            subplotLU(2,size(BOI,2),1,iBOI);
            topoplot(PSDRT_rSpear(:,iBOI,iTrialType), ChanEEGLab,'colormap',colorMapPN,'maplimits',[-1 1]); % Changed from 'maplimits',[-1 1]
            if iBOI==1
                yt=text(-0.55,0,'EEG-Behavior R','horizontalalignment','center','verticalalignment','bottom','fontsize',10,'rotation',90);
            end
            if iBOI==size(BOI,2)
                b=colorbar('southoutside');set(b,'position',[0.52 0.93 0.2 0.01],'xtick',[-1 1],'xticklabel',{'-1' '1'},'xlim',[-1 1]);
                xlabel(b,'PSD-RT Correlation','verticalalignment','top')
            end

            %%%% Plot PSD-RT P-value fig3c supp
            subplotLU(2,size(BOI,2),2,iBOI);
            topoplot(log10(PSDRT_pSpear(:,iBOI,iTrialType)), ChanEEGLab,'colormap',colorMapPN,'maplimits',[-4 4]);
            xlabel(BandName{iBOI});
            text(0,-0.55,[BandName{iBOI} ' (' BandHzName{iBOI} ')'],'horizontalalignment','center','verticalalignment','top','fontsize',10)
            if iBOI==1
                a=ylabel(allGroupNames);
                %                 set(a,'position',[0.01 0.5 0.03 0.4],'verticalalignment','middle')
                set(a,'verticalalignment','middle')
                yt=text(-0.55,0,'P EEG-Behavior','horizontalalignment','center','verticalalignment','bottom','fontsize',10,'rotation',90);
            end

            %%%%% add colorbar
            if iBOI==size(BOI,2)
                %                   c=colorbar('southoutside');set(c,'position',[0.52 0.46 0.2 0.03],'xtick',[-4 0],'xticklabel',{'10e-4' '10e-0'},'xlim',[-4 0]);
                %                   xlabel(c,'P values','verticalalignment','top')
                c=colorbar('southoutside');
                set(gca,'xtick',[],'ytick',[])
                set(c,'position',[0.52 0.51 0.2 0.01],'Limits',[PLim(1) 0],'ticks',[PLim(1) 0],'ticklabels',PLab);
                xlabel(c,'P values','verticalalignment','top')
            end
        
        end

        %%% Apply formatting and save figures for heatmaps
        LuFontStandard;
        papersizePX = [0 0 6*size(BOI,2) 6*2+3];

        set(figureAcc, 'PaperUnits', 'centimeters');
        set(figureAcc, 'PaperPosition', papersizePX, 'PaperSize', papersizePX(3:4));
        saveas(figureAcc, [SaveTemp 'SpearmanPSDAcc_' allGroupNames], 'png');

        set(figureRT, 'PaperUnits', 'centimeters');
        set(figureRT, 'PaperPosition', papersizePX, 'PaperSize', papersizePX(3:4));
        saveas(figureRT, [SaveTemp 'SpearmanPSDRT_' allGroupNames], 'png');

        close all;
   
        %% Perform FDR correction for PSD_Acc and PSD_RT   *** stats for fig3b and fig3c 
        % Define the subset of channels and bands
        figstatsChSubset = [1, 32, 16]; % Indices for Fp1, Cz, Oz
        numBandsToInclude = 5; % First 5 bands

        % Initialize container for p-values
        PSDAcc_allPValues = [];
        PSDRT_allPValues = [];

        %%% Extract p-values and r-values for specified channels and bands
        for iBOI = 1:numBandsToInclude
            for iCh = figstatsChSubset
                % Get the p-value and r-value for the current channel and band
                PSDAcc_pValue = PSDAcc_pSpear(iCh, iBOI, iTrialType);
                PSDAcc_rValue = PSDAcc_rSpear(iCh, iBOI, iTrialType);

                PSDRT_pValue = PSDRT_pSpear(iCh, iBOI, iTrialType);
                PSDRT_rValue = PSDRT_rSpear(iCh, iBOI, iTrialType);

                % Append to the list of all p-values with associated metadata
                PSDAcc_allPValues = [PSDAcc_allPValues; PSDAcc_rValue, PSDAcc_pValue, iBOI, iCh];
                PSDRT_allPValues = [PSDRT_allPValues; PSDRT_rValue, PSDRT_pValue, iBOI, iCh];
            end
        end

        %%% Perform FDR correction
        % Extract raw p-values for FDR correction
        PSDAcc_pValuesForFDR = PSDAcc_allPValues(:, 2);
        PSDRT_pValuesForFDR = PSDRT_allPValues(:, 2);

        % Perform FDR correction
        [PSDAcc_fdrAdjustedPValues,Q_Acc,aPrioriProb_Acc,R_squared_Acc] = mafdr(PSDAcc_pValuesForFDR, 'BHFDR', true);
        [PSDRT_fdrAdjustedPValues,Q_RT,aPrioriProb_RT,R_squared_RT] = mafdr(PSDRT_pValuesForFDR, 'BHFDR', true);

        % Append FDR-adjusted p-values
        PSDAcc_allPValues = [PSDAcc_allPValues, PSDAcc_fdrAdjustedPValues];
        PSDRT_allPValues = [PSDRT_allPValues, PSDRT_fdrAdjustedPValues];

        %%% Display and Save Results for PSD_Acc
        outputFile_Acc = fullfile(SaveTemp, 'PSD_Acc_FDR_corrected_stats.txt');
        fid_Acc = fopen(outputFile_Acc, 'w');

        fprintf(fid_Acc, 'Channel\tBand\tR-value\tRaw_P-value\tFDR Corrected P-value\n');
        fprintf('FDR-corrected p-values for PSD-Acc:\n');
        % Display total number of comparisons
        numComparisons = size(PSDAcc_allPValues, 1);
        fprintf('Total number of comparisons used for FDR correction: %d\n', numComparisons);
        fprintf('Channel\tBand\tR-value\tRaw_P-value\tFDR Corrected P-value\n');

        for i = 1:size(PSDAcc_allPValues, 1)
            chName = ChanEEGLab(PSDAcc_allPValues(i, 4)).labels;
            band = BandName{PSDAcc_allPValues(i, 3)};
            rValue = PSDAcc_allPValues(i, 1);
            rawP = PSDAcc_allPValues(i, 2);
            fdrP = PSDAcc_allPValues(i, 5);

            fprintf('%s\t%s\t%.4f\t%.4f\t%.4f\n', chName, band, rValue, rawP, fdrP);
            fprintf(fid_Acc, '%s\t%s\t%.4f\t%.4f\t%.4f\n', chName, band, rValue, rawP, fdrP);
        end

        fclose(fid_Acc);
        fprintf('PSD-Acc results saved to %s\n', outputFile_Acc);

        %%% Display and Save Results for PSD_RT
        outputFile_RT = fullfile(SaveTemp, 'PSD_RT_FDR_corrected_stats.txt');
        fid_RT = fopen(outputFile_RT, 'w');

        fprintf(fid_RT, 'Ch\tBand\tR-value\tRaw_P-value\tFDR Corrected P-value\n');
        fprintf('FDR-corrected p-values for PSD-RT:\n');
        % Display total number of comparisons
        numComparisons = size(PSDRT_allPValues, 1);
        fprintf('Total number of comparisons used for FDR correction: %d\n', numComparisons);
        fprintf('Ch\tBand\tR-value\tRaw_P-value\tFDR Corrected P-value\n');

        for i = 1:size(PSDRT_allPValues, 1)
            chName = ChanEEGLab(PSDRT_allPValues(i, 4)).labels;
            band = BandName{PSDRT_allPValues(i, 3)};
            rValue = PSDRT_allPValues(i, 1);
            rawP = PSDRT_allPValues(i, 2);
            fdrP = PSDRT_allPValues(i, 5);

            fprintf('%s\t%s\t%.4f\t%.4f\t%.4f\n', chName, band, rValue, rawP, fdrP);
            fprintf(fid_RT, '%s\t%s\t%.4f\t%.4f\t%.4f\n', chName, band, rValue, rawP, fdrP);
        end

        fclose(fid_RT);
        fprintf('PSD-RT results saved to %s\n', outputFile_RT);
    end
end

%% PSD-Acc - Loop of correlation scatter plot? (future fig3a+b?)
FlickerColor=[31 125 184; 150 27 27; 219 129 50]/255; %40, Random, Light % blue, red, gold
for iFreqFunc=3%1:length(FunGroupName)
    clear MapGroup diffTmap diffMap Tdata;
    clear diffTPmap diffRPmap diffTmap

    for iTrialType=3%1:length(TrialType)
        % iCom=1;
        SaveTemp=[SubSavePSD TrialTypeName{iTrialType} '\'];
        SaveTemp=[SaveTemp FunGroupName{iFreqFunc} '\' ];
        mkdir(SaveTemp)

        BOI=[5;15];
        BOI=[1 4 8 13 30 39.5 43;4 8 13 30 37 41.5 100];
        BName={'Delta','Theta','Alpha','Beta','Gamma-1','Gamma-E','Gamma-2'};
        BName2={'1-4 Hz','4-8 Hz','8-13Hz','13-30 Hz','30-37 Hz','39-41 Hz','43-100 Hz'};

        PSDBehaviorGroup1 = 1;
        PSDBehaviorGroup2 = 3; 
        PSDBehaviorGroup3 = 6; % 6 = LightRT
        allGroupNames = [GroupName{PSDBehaviorGroup1} GroupName{PSDBehaviorGroup2} GroupName{PSDBehaviorGroup3}];
        dataName = [TrialTypeName{iTrialType} FreqFuncNames{iFreqFunc} allGroupNames];

        IncludedSubj=[SubjG{PSDBehaviorGroup1}(:);SubjG{PSDBehaviorGroup2}(:);SubjG{PSDBehaviorGroup3}(:)];  % 1 3 6 = 40, Random, Light
        FlickerID=[zeros(size(SubjG{PSDBehaviorGroup1}(:)))+1;zeros(size(SubjG{PSDBehaviorGroup2}(:)))+2;zeros(size(SubjG{PSDBehaviorGroup3}(:)))+3];
        Param.Corr='Spearman';  %%%Type of correlation, see Matlab function corr for more details
        Param.Pth=0.05; %%%threshold of Pvalue
        Param.EdgeColor=colorMapPN; %%%Color map for correlation link
        Param.NodeColor=repmat([0.8 0.8 0.8],6,1); %%%Color of Nodes for correlation link plot
        Param.Clim=[-1 1];    %%%Color Limit for Correlation
        Param.Title='Pool All Sample'; %%Any title for label the figure
        Param.MarkerSize=8; %%%MarkerSize of scatter
        Param.SubjIDColor=FlickerColor;
        Param.SubjID=FlickerID;


        %% Scatterplot of correlation - indiv behavior and PSD
        for iBOI=1:size(BOI,2)
            NeedI=find(Fplot>=BOI(1,iBOI)&Fplot<=BOI(2,iBOI));
            figure;

            %          for iCh=1:length(ChanEEGLab)
            if iFreqFunc==3
                PSDtemp=squeeze(FunGroup{iFreqFunc}(LogPSD(IncludedSubj,NeedI,EEGchInd,iTrialType),[],2));
            else
                PSDtemp=squeeze(FunGroup{iFreqFunc}(LogPSD(IncludedSubj,NeedI,EEGchInd,iTrialType),2));
            end
            Acctemp=Acc(IncludedSubj);

            tempName=EEGch;
            tempName{end+1}='Acc';
            multiCorr2GroupSubplot(6,6,[PSDtemp Acctemp],size(PSDtemp,2)+1,tempName,Param)
            LuFontStandard;

            papersizePX=[0 0 6*6 6*6];
            set(gcf, 'PaperUnits', 'centimeters');
            set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));

            % saveas(gcf,[SaveTemp Param.Corr BName{iBOI} BName2{iBOI} 'PSDAcc_' allGroupNames],'pdf');
            saveas(gcf,[SaveTemp Param.Corr BName{iBOI} BName2{iBOI} 'PSDAcc_' allGroupNames],'png');
            saveas(gcf,[SaveTemp Param.Corr BName{iBOI} BName2{iBOI} 'PSDAcc_' allGroupNames '.eps'],'epsc');
            close all
        end
    end
    %
    %     close all
    %
end

%% Preallocate ChannelPairName Cell Array - *** Start of fig4e and fig4f scatterplots and histograms
signifChPairNameCellArray = cell(3,3,3,3,2);
signifChPairNameCellArray{3,3,3,3,2} = [];
WPLIBResults = struct();
%% WPLI-Behavior {FC-B} scatter plot for each channel pair. Alphas Only 7/18/24 MKA *** fig4e fig4f
for iFreqFunc=3%[1:2]%1:length(FunGroupName) % mean = 1, median =2, peak = 3
    %% Set up
    clear MapGroup diffTmap diffMap Tdata;
    clear diffTPmap diffRPmap diffTmap
    globalPvalAlpha = 0.1;
    Param.Corr='Spearman';  %%%Type of correlation, see Matlab function corr for more details, 'Pearson' or 'Spearman'
    todayDate = datestr(now, 'yymmdd');
    %% Iterate by trial type
    for iTrialType=3 %:length(TrialType) % 1=Hit, 2=Miss, 3=HitANDMiss
        SaveTemp=[SubSaveWPLI TrialTypeName{iTrialType} '\' todayDate '_' FunGroupName{iFreqFunc} '_' Param.Corr 'P' num2str(globalPvalAlpha) '\'];

        %% Define Bands of Interest
        BOI=[8 8 10; 13 10 13];
        BName={'Alpha','LowAlpha','HighAlpha'};
        % Gamma-E is  40
        BName2={'8-13Hz','8-10Hz','10-13Hz',};

        %% Define Group Set Names
        GroupSetsName{1} ='40andL'; % 40 and Light
        GroupSetsName{2}='40andR'; % 40 and Random
        GroupSetsName{3}='All3Groups'; % All three groups

        %% Loop thru all groups
        for iGroupSet = 1:2%:length(GroupSetsName)
            if iGroupSet == 1
                %% 40&LightRT
                % included subs
                LightRTSubjG{1}=[];
                [~,LightRTSubjG{1},~]=intersect(SubjID,FlickerSubj{6}); % Light
                IncludedSubj=[SubjG{1}(:);LightRTSubjG{1}(:)];
                FlickerID=[zeros(size(SubjG{1}(:)))+1;zeros(size(LightRTSubjG{1}(:)))+2];

                Param.SubjIDColor=FlickerColor;

            elseif iGroupSet == 2
                %% 40andR
                IncludedSubj=[SubjG{1}(:);SubjG{3}(:)];
                FlickerID=[zeros(size(SubjG{1}(:)))+1;zeros(size(SubjG{3}(:)))+2];
                GFlickerRandColor=[31 125 184; 150 27 27]/255;

                Param.SubjIDColor=GFlickerRandColor;

            elseif iGroupSet == 3
                %% All3Groups
                LightRTSubjG{1}=[];
                [~,LightRTSubjG{1},~]=intersect(SubjID,FlickerSubj{6}); % Light
                IncludedSubj2=[SubjG{1}(:);LightRTSubjG{1}(:);SubjG{3}(:)];
                FlickerID=[zeros(size(SubjG{1}(:)))+1;zeros(size(LightRTSubjG{1}(:)))+2;zeros(size(SubjG{3}(:)))+3];

                Param.SubjIDColor=FlickerColor;

            end
            %% Set up folder FC-B (General)
            SaveTempGroup = [SaveTemp GroupSetsName{iGroupSet} '\'];
            SaveTempAcc=[SaveTempGroup FunGroupName{iFreqFunc} '\WPLI-Acc ChPair Scatterp\' ];
            SaveTempRT = [SaveTempGroup FunGroupName{iFreqFunc} '\WPLI-RT ChPair Scatterp\' ];
            mkdir(SaveTempAcc)
            mkdir(SaveTempRT)
            %% Parameters - FC-B Scatter (GENERAL)

            Param.Pth=0.1; %%%threshold of Pvalue
            Param.EdgeColor=colorMapPN; %%%Color map for correlation link
            Param.NodeColor=repmat([0.8 0.8 0.8],6,1); %%%Color of Nodes for correlation link plot
            Param.Clim=[-1 1];    %%%Color Limit for Correlation
            Param.Title='Pool All Sample'; %%Any title for label the figure
            Param.MarkerSize=8; %%%MarkerSize of scatter
            Param.SubjID=FlickerID;
            %% Scatterplot of correlation - FC-B  General - *** fig4e fig4f scatterplots
            for iBOI=2:size(BOI,2) %iBOI1,2,3 = ALpha,LowAlpha,HighAlpha
                %% Initialize temp variables: Acctemp, RTtemp WPLItemp
                NeedI=find(Fplot>=BOI(1,iBOI)&Fplot<=BOI(2,iBOI));
                figure;

                if iFreqFunc==3
                    WPLItemp=squeeze(FreqFunc{iFreqFunc}(WPLIall(IncludedSubj,NeedI,EEGchInd,EEGchInd,iTrialType),[],2));
                else
                    WPLItemp=squeeze(FreqFunc{iFreqFunc}(WPLIall(IncludedSubj,NeedI,EEGchInd,EEGchInd,iTrialType),2));
                end

                Acctemp=Acc(IncludedSubj); % all accuracies of included subjects
                RTtemp=SubjsAvgRT(IncludedSubj);
                Param.SubjID=FlickerID;

                %% Include only subjects with high accuracy (>80%) - Accuracy cuts
                highAccThreshold = 0.8;
                highAccSubjectsIndex = Acctemp>highAccThreshold;
                % highAccTemp = Acctemp(Acctemp>highAccThreshold);
                % highAccRTtemp = RTtemp(Acctemp>highAccThreshold);
                % highAccWPLItemp = WPLItemp(Acctemp>highAccThreshold);
                Acctemp = Acctemp(highAccSubjectsIndex);
                RTtemp = RTtemp(highAccSubjectsIndex);
                WPLItemp = WPLItemp(highAccSubjectsIndex,:,:);
                Param.SubjID = Param.SubjID(highAccSubjectsIndex);

                %% Initialize channel pair variables and set alpha-pval threhold to save plots
                chPairRAcc = NaN([32 32]);
                chPairPvalAcc = NaN([32 32]);
                chPairRRT = NaN([32 32]);
                chPairPvalRT = NaN([32 32]);
                alphaSavePlot = globalPvalAlpha;  % If p<alphaSavePlot, then that plot gets saved.
                %% Scatterplot iteration fig4
                for iCh = 1:31
                    for jCh = iCh+1:32
                        %% WPLI-Acc Scatter Plot
                        WPLIsingleChPair = squeeze(WPLItemp(:,iCh,jCh));
                        channelPairNAme = [EEGch{iCh} EEGch{jCh}];
                        tempName={channelPairNAme};
                        tempName{end+1}='Acc';

                        [chPairRAcc(iCh,jCh), chPairPvalAcc(iCh,jCh)] = multiCorr2GroupSubplot(1,1,[WPLIsingleChPair Acctemp],size(WPLIsingleChPair,2)+1,tempName,Param);
                        %% Plot significant scatterplots only (fig4)
                        % if chPairPvalAcc(iCh,jCh)<alphaSavePlot
                        %     ylabel("Accuracy")
                        %     xlabel("WPLI Amplitude")
                        %     set(gca,'ylim',[0.8 1], 'ytick',[0.8:0.1:1]);
                        %     papersizePX=[0 0 2.5 2.5];  % Changed from [0 0 6*6 6*6] on 7/18/24
                        %     set(gcf, 'PaperUnits', 'centimeters');
                        %     set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
                        %     LuFontStandard;
                        %     plotname =  ['WPLIAccScatterplot' GroupSetsName{iGroupSet}  '_' channelPairNAme];
                        %     % saveas(gcf,[SaveTemp Param.Corr BName{iBOI} BName2{iBOI} plotname],'pdf');
                        %     saveas(gcf,[SaveTempAcc Param.Corr BName{iBOI} BName2{iBOI} plotname],'png');
                        %     saveas(gcf,[SaveTempAcc Param.Corr BName{iBOI} plotname '.svg'],'svg');
                        %     % saveas(gcf,[SaveTemp Param.Corr BName{iBOI} BName2{iBOI} plotname '.eps'],'epsc');
                        % end
                        close all

                        %% WPLI-RT Scatter plot
                        tempName{end}='RT';

                        [chPairRRT(iCh,jCh), chPairPvalRT(iCh,jCh)] = multiCorr2GroupSubplot(1,1,[WPLIsingleChPair RTtemp],size(WPLIsingleChPair,2)+1,tempName,Param);
                        %% plot significant scatterplots only (fig4 MS)
                        % if chPairPvalRT(iCh,jCh)<alphaSavePlot
                        %     LuFontStandard;
                        %     ylabel("Reaction Time")
                        %     xlabel("WPLI Amplitude")
                        %     % papersizeWPLIBScat = [6 6];
                        %     % paperpositionWPLIBScat = [0 0 papersizeWPLIBScat];
                        %     % papersizePX=[1 1 5 5];  % Changed from [0 0 6*6 6*6] on 7/18/24
                        %     set(gca,'ytick',[0.3:0.1:0.6]);
                        %     set(gcf, 'PaperUnits', 'centimeters');
                        %     % set(gcf,'PaperSize',[5 1]);
                        %     ti = get(gca, 'TightInset');  % Get margins needed for labels, titles, etc.
                        %     % set(gcf, 'PaperPosition', [ti(1) ti(2) 15-ti(1)-ti(3) 10-ti(2)-ti(4)]);
                        %     set(gcf,'PaperPosition',[0 0 2 2]);
                        % 
                        %     plotname =  ['WPLIRTScatter' GroupSetsName{iGroupSet}  '_' channelPairNAme];
                        %     % saveas(gcf,[SaveTemp Param.Corr BName{iBOI} BName2{iBOI} plotname],'pdf');
                        %     saveas(gcf,[SaveTempRT Param.Corr BName{iBOI} BName2{iBOI} plotname],'png');
                        %     saveas(gcf,[SaveTempRT Param.Corr BName{iBOI} BName2{iBOI} plotname],'svg');
                        %     % saveas(gcf,[SaveTemp Param.Corr BName{iBOI} BName2{iBOI} plotname '.eps'],'epsc');
                        % end
                        close all
                    end
                end

                %% Save results in multdimensional structure!
                WPLIBResults(iFreqFunc, iTrialType, iBOI, iGroupSet).chPairRAcc = chPairRAcc;
                WPLIBResults(iFreqFunc, iTrialType, iBOI, iGroupSet).chPairPvalAcc = chPairPvalAcc;
                WPLIBResults(iFreqFunc, iTrialType, iBOI, iGroupSet).chPairRRT = chPairRRT;
                WPLIBResults(iFreqFunc, iTrialType, iBOI, iGroupSet).chPairPvalRT = chPairPvalRT;
                WPLIBResults(iFreqFunc, iTrialType, iBOI, iGroupSet).Param = Param;
                WPLIBResults(iFreqFunc, iTrialType, iBOI, iGroupSet).globalPvalAlpha = globalPvalAlpha;

                %% Save r values and associated p-values to a file
                rpfilename = [SaveTempGroup BName{iBOI} BName2{iBOI} '_' TrialTypeName{iTrialType} '_' FunGroupName{iFreqFunc} '_ACCRT_RPvalues'];
                save(rpfilename, 'chPairRAcc', 'chPairPvalAcc', 'chPairRRT', 'chPairPvalRT','iBOI', 'Param', 'globalPvalAlpha')

                %% Plot r and r^squared histograms for WPLI-ACC and WPLI-RT GENERAL (not used?)
                % subplot(4,1,1);
                % % histogram(chPairRAcc(~isnan(chPairRAcc)))
                % histogram(chPairRAcc)
                % title(['ChPair WPLI-Acc ' Param.Corr ' r'])
                % ylabel('Ch Pair Count')
                % xlabel(['r - nChPairs=' num2str(sum(~isnan(chPairRAcc(:))))])
                % 
                % subplot(4,1,2);
                % chPairRSquaredAcc = chPairRAcc.^2;
                % histogram(chPairRSquaredAcc)
                % title(['WPLI-Acc ' Param.Corr ' r^2'])
                % ylabel('Ch Pair Count')
                % xlabel(['r^2 - nChPairs=' num2str(sum(~isnan(chPairRSquaredAcc(:))))])
                % 
                % subplot(4,1,3);
                % % histogram(chPairRRT(~isnan(chPairRRT)))
                % histogram(chPairRRT)
                % title(['WPLI-RT ' Param.Corr ' r'])
                % ylabel('Ch Pair Count')
                % xlabel(['r - nChPairs=' num2str(sum(~isnan(chPairRRT(:))))])
                % 
                % subplot(4,1,4);
                % chPairRSquaredRT = chPairRRT.^2;
                % histogram(chPairRSquaredRT)
                % title(['WPLI-RT ' Param.Corr ' r^2'])
                % ylabel('Ch Pair Count')
                % xlabel(['r^2 - nChPairs=' num2str(sum(~isnan(chPairRSquaredRT(:))))])
                % 
                % sgtitle([GroupSetsName{iGroupSet} '-' TrialTypeName{iTrialType} FunGroupName{iFreqFunc} BName{iBOI} BName2{iBOI} Param.Corr])
                % saveas(gcf, [SaveTempGroup TrialTypeName{iTrialType} FunGroupName{iFreqFunc} BName{iBOI} BName2{iBOI} '.png'])

                %% Get significant R values
                % Get indices of channel pairs with pvals less than alpha (alpha=0.05)
                alpha = globalPvalAlpha;

                % Get significant Acc
                signifAccPval = chPairPvalAcc > 0 & chPairPvalAcc < alpha; % create logical mask, exclude 0's which correspond to no channel pair or a redundent pair
                signifchPairRAcc = chPairRAcc(signifAccPval);
                signifChPairNameCellArray{iFreqFunc,iTrialType,iGroupSet,iBOI,1} = getListofSignifChPairs(signifAccPval,EEGch);

                % Get significant RT
                signifRtPval = chPairPvalRT > 0 & chPairPvalRT < alpha; % create logical mask
                signifchPairRRT = chPairRRT(signifRtPval);
                signifChPairNameCellArray{iFreqFunc,iTrialType,iGroupSet,iBOI,2} = getListofSignifChPairs(signifRtPval,EEGch);

                if isequal(BName{iBOI},'HighAlpha')
                    HighAlphaRAcc = chPairRAcc;
                    HighAlphaPAcc = chPairPvalAcc;
                    HighAlphaRRT = chPairRRT;
                    HighAlphaPRT = chPairPvalRT;
                    HighAlphaRSigAcc = signifchPairRAcc;
                    HighAlphaRSigRT = signifchPairRRT;
                    AlphaData.HighAlpha.RAcc = chPairRAcc;
                    AlphaData.HighAlpha.PAcc = chPairPvalAcc;
                    AlphaData.HighAlpha.RRT = chPairRRT;
                    AlphaData.HighAlpha.PRT = chPairPvalRT;
                    AlphaData.HighAlpha.RSigAcc = signifchPairRAcc;
                    AlphaData.HighAlpha.RSigRT = signifchPairRRT;
                elseif isequal(BName{iBOI},'LowAlpha')
                    LowAlphaRAcc = chPairRAcc;
                    LowAlphaPAcc = chPairPvalAcc;
                    LowAlphaRRT = chPairRRT;
                    LowAlphaPRT = chPairPvalRT;
                    LowAlphaRSigAcc = signifchPairRAcc;
                    LowAlphaRSigRT = signifchPairRRT;
                    AlphaData.LowAlpha.RAcc = chPairRAcc;
                    AlphaData.LowAlpha.PAcc = chPairPvalAcc;
                    AlphaData.LowAlpha.RRT = chPairRRT;
                    AlphaData.LowAlpha.PRT = chPairPvalRT;
                    AlphaData.LowAlpha.RSigAcc = signifchPairRAcc;
                    AlphaData.LowAlpha.RSigRT = signifchPairRRT;
                elseif isequal(BName{iBOI},'Alpha')
                    AllAlphaRAcc = chPairRAcc;
                    AllAlphaPAcc = chPairPvalAcc;
                    AllAlphaRRT = chPairRRT;
                    AllAlphaPRT = chPairPvalRT;
                    AllAlphaRSigAcc = signifchPairRAcc;
                    AllAlphaRSigRT = signifchPairRRT;
                    AlphaData.AllAlpha.RAcc = chPairRAcc;
                    AlphaData.AllAlpha.PAcc = chPairPvalAcc;
                    AlphaData.AllAlpha.RRT = chPairRRT;
                    AlphaData.AllAlpha.PRT = chPairPvalRT;
                    AlphaData.AllAlpha.RSigAcc = signifchPairRAcc;
                    AlphaData.AllAlpha.RSigRT = signifchPairRRT;
                end

                %% r & r&2 distribution plots (Significant Only) not used for MS?
                % r Distribution for Significant WPLI-Acc Channel pairs
                % subplot(4,1,1);
                % histogram(signifchPairRAcc)        % histogram(chPairRAcc(~isnan(chPairRAcc)))
                % title(['ChPair WPLI-Acc ' Param.Corr ' r (p<' num2str(alpha) ')'])
                % ylabel('Ch Pair Count')
                % xlabel(['r - nChPairs=' num2str(numel(signifchPairRAcc))])
                % 
                % % r^2 distribution for Significant WPLI-Acc Channel pairs
                % subplot(4,1,2);
                % chPairRSquaredAcc = chPairRAcc.^2;
                % signifchPairRSquaredAcc=chPairRSquaredAcc(signifAccPval);
                % histogram(signifchPairRSquaredAcc)
                % title(['WPLI-Acc r^2 (p<' num2str(alpha) ')'])
                % ylabel('Ch Pair Count')
                % xlabel(['r^2 - nChPairs=' num2str(numel(signifchPairRSquaredAcc))])
                % 
                % % r Distribution for Significant WPLI-RT Channel pairs
                % subplot(4,1,3);
                % histogram(signifchPairRRT)      % histogram(chPairRRT(~isnan(chPairRRT)))
                % title(['WPLI-RT  r (p<' num2str(alpha) ')'])
                % ylabel('Ch Pair Count')
                % xlabel(['r - nChPairs=' num2str(numel(signifchPairRRT))])
                % 
                % % r^2 Distribution for Significant WPLI-RT Channel pairs
                % subplot(4,1,4);
                % chPairRSquaredRT = chPairRRT.^2;
                % signifchPairRSquaredRT = chPairRSquaredRT(signifRtPval);
                % histogram(signifchPairRSquaredRT)
                % title(['WPLI-RT r^2 (p<' num2str(alpha) ')'])
                % ylabel('Ch Pair Count')
                % xlabel(['r^2 - nChPairs=' num2str(numel(signifchPairRSquaredRT))])
                % 
                % sgtitle([GroupSetsName{iGroupSet}  '-'  TrialTypeName{iTrialType} ' ' FunGroupName{iFreqFunc} BName{iBOI} BName2{iBOI} ' ' Param.Corr ' p<' num2str(alpha)])
                % saveas(gcf, [SaveTempGroup TrialTypeName{iTrialType} FunGroupName{iFreqFunc} BName{iBOI} BName2{iBOI} '_' Param.Corr 'p' num2str(alpha) '.png'])

                %% Get the top quartile of R values for significant WPLI-Acc channel pairs
                % % Flatten the 2D arrays into a 1D array
                % chPairRAcc_flat = signifchPairRAcc(:);
                % chPairRRT_flat = signifchPairRRT(:);
                %
                % % Sort the 1D array in descending order
                % chPairRAcc_sorted = sort(chPairRAcc_flat, 'descend');
                % chPairRRT_sorted = sort(chPairRRT_flat, 'descend');
                %
                % % Determine the index threshold for the top quartile
                % n = numel(chPairRAcc_sorted);
                % AccQuartile_threshold = ceil(n / 4);
                %
                % n = numel(chPairRRT_sorted);
                % RTQuartile_threshold = ceil(n / 4);
                %
                % % Find the value threshold for the top quartile
                % Acc_threshold_value = chPairRAcc_sorted(AccQuartile_threshold);
                % RT_threshold_value = chPairRRT_sorted(RTQuartile_threshold);
                %
                % % Create the logical mask
                % Acc_topQuartile_mask = signifchPairRAcc >= Acc_threshold_value;
                % RT_topQuartile_mask = signifchPairRRT >= RT_threshold_value;
                %
                % % Display the result
                % % disp('Logical mask for the top quartile values:');
                % % disp(topQuartile_mask);

                %% Plot the top quartile
                % % r Distribution for Significant WPLI-Acc Channel pairs
                % subplot(4,1,1);
                % topQuartile_signifchPairRAcc = signifchPairRAcc(Acc_topQuartile_mask);
                % histogram(topQuartile_signifchPairRAcc)        % histogram(chPairRAcc(~isnan(chPairRAcc)))
                % title(['ChPair WPLI-Acc r (p<' num2str(alpha) ')'])
                % ylabel('Ch Pair Count')
                % xlabel(['r - nChPairs=' num2str(numel(topQuartile_signifchPairRAcc))])
                %
                %
                % % r^2 distribution for Significant WPLI-Acc Channel pairs
                % subplot(4,1,2);
                % chPairRSquaredAcc = chPairRAcc.^2;
                % signifchPairRSquaredAcc=chPairRSquaredAcc(signifAccPval);
                % histogram(signifchPairRSquaredAcc(Acc_topQuartile_mask))
                % title(['WPLI-Acc r^2 (p<' num2str(alpha) ')'])
                % ylabel('Ch Pair Count')
                % xlabel(['r^2 - nChPairs=' num2str(numel(signifchPairRSquaredAcc(Acc_topQuartile_mask)))])
                %
                % % r Distribution for Significant WPLI-RT Channel pairs
                % subplot(4,1,3);
                % histogram(signifchPairRRT(RT_topQuartile_mask))      % histogram(chPairRRT(~isnan(chPairRRT)))
                % title(['WPLI-RT  r (p<' num2str(alpha) ')'])
                % ylabel('Ch Pair Count')
                % xlabel(['r - nChPairs=' num2str(numel(signifchPairRRT(RT_topQuartile_mask)))])
                %
                % % r^2 Distribution for Significant WPLI-RT Channel pairs
                % subplot(4,1,4);
                % chPairRSquaredRT = chPairRRT.^2;
                % signifchPairRSquaredRT = chPairRSquaredRT(signifRtPval);
                % histogram(signifchPairRSquaredRT(RT_topQuartile_mask))
                % title(['WPLI-RT r^2 (p<' num2str(alpha) ')'])
                % ylabel('Ch Pair Count')
                % xlabel(['r^2 - nChPairs=' num2str(numel(signifchPairRSquaredRT(RT_topQuartile_mask)))])
                %
                % sgtitle([GroupSetsName{iGroupSet}  '-' TrialTypeName{iTrialType} FunGroupName{iFreqFun} BName{iBOI} BName2{iBOI} ' p<' num2str(alpha) ' TopQuartile'])
                % saveas(gcf, [SaveTempGroup TrialTypeName{iTrialType} FunGroupName{iFreqFun} BName{iBOI} BName2{iBOI} '_p' num2str(alpha) '_TopQuartile.png'])
            end

            %% Low Alpha vs High Alpha histogrma and ttest comparison
            % Open a text file for writing
            fileName = [SaveTempGroup TrialTypeName{iTrialType} FunGroupName{iFreqFunc} '_ttest_results.txt'];
            fileID = fopen(fileName, 'w');
  

            % Remove NaN values and convert to vectors
            LowAlphaRAcc_non_nan = LowAlphaRAcc(~isnan(LowAlphaRAcc));
            HighAlphaRAcc_non_nan = HighAlphaRAcc(~isnan(HighAlphaRAcc));
            LowAlphaRRT_non_nan = LowAlphaRRT(~isnan(LowAlphaRRT));
            HighAlphaRRT_non_nan = HighAlphaRRT(~isnan(HighAlphaRRT));

            % Perform single-sample t-tests (t-Test for Mean Equal to Zero)
            [h1, p1, ci1, stats1] = ttest(LowAlphaRAcc_non_nan);
            [h2, p2, ci2, stats2] = ttest(HighAlphaRAcc_non_nan);
            [h3, p3, ci3, stats3] = ttest(LowAlphaRRT_non_nan);
            [h4, p4, ci4, stats4] = ttest(HighAlphaRRT_non_nan);

            % Perform paired t-tests
            [h5, p5, ci5, stats5] = ttest(LowAlphaRAcc_non_nan, HighAlphaRAcc_non_nan);
            [h6, p6, ci6, stats6] = ttest(LowAlphaRRT_non_nan, HighAlphaRRT_non_nan);

            % Perform unpaired t-tests (Welch's t-test)
            [h7, p7, ci7, stats7] = ttest2(LowAlphaRSigAcc, HighAlphaRSigAcc, 'Vartype', 'unequal'); % This will perform Welch's t-test, which does not assume equal variances between the two groups.
            [h8, p8, ci8, stats8] = ttest2(LowAlphaRSigRT, HighAlphaRSigRT, 'Vartype', 'unequal');


            % Write to file
            fprintf(fileID, 'Single-sample t-tests (t-Test for Mean Equal to Zero):\n');
            fprintf(fileID, 'LowAlphaRAcc: p-value = %.16f, t-statistic = %.4f, df = %.2f, CI = [%.4f, %.4f]\n', p1, stats1.tstat, stats1.df, ci1(1), ci1(2));
            fprintf(fileID, 'HighAlphaRAcc: p-value = %.16f, t-statistic = %.4f, df = %.2f, CI = [%.4f, %.4f]\n', p2, stats2.tstat, stats2.df, ci2(1), ci2(2));
            fprintf(fileID, 'LowAlphaRRT: p-value = %.16f, t-statistic = %.4f, df = %.2f, CI = [%.4f, %.4f]\n', p3, stats3.tstat, stats3.df, ci3(1), ci3(2));
            fprintf(fileID, 'HighAlphaRRT: p-value = %.16f, t-statistic = %.4f, df = %.2f, CI = [%.4f, %.4f]\n', p4, stats4.tstat, stats4.df, ci4(1), ci4(2));

            fprintf(fileID, '\nPaired t-tests (All channel pairs):\n');
            fprintf(fileID, 'LowAlphaRAcc vs HighAlphaRAcc: p-value = %.16f, t-statistic = %.4f, df = %.2f, CI = [%.4f, %.4f]\n', p5, stats5.tstat, stats5.df, ci5(1), ci5(2));
            fprintf(fileID, 'LowAlphaRRT vs HighAlphaRRT: p-value = %.16f, t-statistic = %.4f, df = %.2f, CI = [%.4f, %.4f]\n', p6, stats6.tstat, stats6.df, ci6(1), ci6(2));

            fprintf(fileID, '\nUnpaired t-tests (Welch''s t-test on Significant Channel Pairs p<%.4f):\n', alpha);
            fprintf(fileID, 'LowAlphaRSigAcc vs HighAlphaRSigAcc: p-value = %.16f, t-statistic = %.4f, df = %.2f, CI = [%.4f, %.4f]\n', p7, stats7.tstat, stats7.df, ci7(1), ci7(2));
            fprintf(fileID, 'LowAlphaRSigRT vs HighAlphaRSigRT: p-value = %.16f, t-statistic = %.4f, df = %.2f, CI = [%.4f, %.4f]\n', p8, stats8.tstat, stats8.df, ci8(1), ci8(2));

            % Close the file
            fclose(fileID);

            % Print to command window
            fprintf('Single-sample t-tests (t-Test for Mean Equal to Zero):\n');
            fprintf('LowAlphaRAcc: p-value = %.16f, t-statistic = %.4f, df = %.2f, CI = [%.4f, %.4f]\n', p1, stats1.tstat, stats1.df, ci1(1), ci1(2));
            fprintf('HighAlphaRAcc: p-value = %.16f, t-statistic = %.4f, df = %.2f, CI = [%.4f, %.4f]\n', p2, stats2.tstat, stats2.df, ci2(1), ci2(2));
            fprintf('LowAlphaRRT: p-value = %.16f, t-statistic = %.4f, df = %.2f, CI = [%.4f, %.4f]\n', p3, stats3.tstat, stats3.df, ci3(1), ci3(2));
            fprintf('HighAlphaRRT: p-value = %.16f, t-statistic = %.4f, df = %.2f, CI = [%.4f, %.4f]\n', p4, stats4.tstat, stats4.df, ci4(1), ci4(2));

            fprintf('\nPaired t-tests (All channel pairs):\n');
            fprintf('LowAlphaRAcc vs HighAlphaRAcc: p-value = %.16f, t-statistic = %.4f, df = %.2f, CI = [%.4f, %.4f]\n', p5, stats5.tstat, stats5.df, ci5(1), ci5(2));
            fprintf('LowAlphaRRT vs HighAlphaRRT: p-value = %.16f, t-statistic = %.4f, df = %.2f, CI = [%.4f, %.4f]\n', p6, stats6.tstat, stats6.df, ci6(1), ci6(2));

            fprintf('\nUnpaired t-tests (Welch''s t-test on Significant Channel Pairs p<%.4f):\n', alpha);
            fprintf('LowAlphaRSigAcc vs HighAlphaRSigAcc: p-value = %.16f, t-statistic = %.4f, df = %.2f, CI = [%.4f, %.4f]\n', p7, stats7.tstat, stats7.df, ci7(1), ci7(2));
            fprintf('LowAlphaRSigRT vs HighAlphaRSigRT: p-value = %.16f, t-statistic = %.4f, df = %.2f, CI = [%.4f, %.4f]\n', p8, stats8.tstat, stats8.df, ci8(1), ci8(2));

            %% Create overlapping histograms for LowAlphaRAcc and HighAlphaRAcc
            figure;

            % Font size for text annotations
            fontSize = 6;  % You can adjust this value to make the font smaller or larger
            legendFontSize = 5;
            yOffset = -0.25;
            binWidth = 0.05;  % Bin size for histograms

            darkgreen = "#77AC30";
            grey = '#727272';
            lowAlphaColor = darkgreen;
            highAlphaColor = grey;

            % 1st subplot: LowAlphaRAcc vs HighAlphaRAcc
            subplot(2, 2, 1);
            histogram(LowAlphaRAcc_non_nan, 'Normalization', 'count', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'BinWidth', binWidth, 'FaceColor', darkgreen);
            hold on;
            histogram(HighAlphaRAcc_non_nan , 'Normalization', 'count', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'BinWidth', binWidth, 'FaceColor', grey);
            hold off;
            title(['WPLI-Acc ' Param.Corr ' r']);
            ylabel('Ch Pair Count');
            xlabel(['r - Paired: t(' num2str(stats5.df) ')=' num2str(stats5.tstat, '%.2f') ', p=' num2str(p5, '%.3f')] , 'FontSize', fontSize);
            % legend('LowAlphaRAcc', 'HighAlphaRAcc');
            lgd1 = legend({['LowAlphaRAcc: t(' num2str(stats1.df) ')=' num2str(stats1.tstat, '%.2f') ', p=' num2str(p1, '%.3f')], ...
                ['HighAlphaRAcc: t(' num2str(stats2.df) ')=' num2str(stats2.tstat, '%.2f') ', p=' num2str(p2, '%.3f')], ...
                ['Paired: t(' num2str(stats5.df) ')=' num2str(stats5.tstat, '%.2f') ', p=' num2str(p5, '%.3f')]}, 'FontSize', legendFontSize);
            set(lgd1, 'Color', 'none');


            % 2nd subplot: LowAlphaRRT vs HighAlphaRRT
            subplot(2, 2, 2);
            histogram(LowAlphaRRT_non_nan , 'Normalization', 'count', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'BinWidth', binWidth, 'FaceColor', darkgreen);
            hold on;
            histogram(HighAlphaRRT_non_nan , 'Normalization', 'count', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'BinWidth', binWidth, 'FaceColor', grey);
            hold off;
            title(['WPLI-RT ' Param.Corr ' r']);
            ylabel('Ch Pair Count');
            xlabel(['r  - Paired: t(' num2str(stats6.df) ')=' num2str(stats6.tstat, '%.2f') ', p=' num2str(p6, '%.3f')], 'FontSize', fontSize);
            % legend('LowAlphaRRT', 'HighAlphaRRT');
            lgd2= legend({['LowAlphaRRT: t(' num2str(stats3.df) ')=' num2str(stats3.tstat, '%.2f') ', p=' num2str(p3, '%.3f')], ...
                ['HighAlphaRRT: t(' num2str(stats4.df) ')=' num2str(stats4.tstat, '%.2f') ', p=' num2str(p4, '%.3f')], ...
                ['Paired: t(' num2str(stats6.df) ')=' num2str(stats6.tstat, '%.2f') ', p=' num2str(p6, '%.3f')]}, 'FontSize', legendFontSize);
            set(lgd2, 'Color', 'none');


            % 3rd subplot: LowAlphaRSigAcc vs HighAlphaRSigAcc
            subplot(2, 2, 3);
            histogram(LowAlphaRSigAcc, 'Normalization', 'count', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'BinWidth', binWidth, 'FaceColor', darkgreen);
            hold on;
            histogram(HighAlphaRSigAcc, 'Normalization', 'count', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'BinWidth', binWidth, 'FaceColor', grey);
            hold off;
            title(['WPLI-Acc ' Param.Corr ' r (p<' num2str(alpha) ')']);
            ylabel('Ch Pair Count');
            xlabel(['r - Welch: t(' num2str(stats7.df) ')=' num2str(stats7.tstat, '%.2f') ', p=' num2str(p7, '%.3f')] , 'FontSize', fontSize);
            lgd3 = legend(['LowAlphaRSigAcc nChPairs=' num2str(numel(LowAlphaRSigAcc))], ['HighAlphaRSigAcc nChPairs=' num2str(numel(HighAlphaRSigAcc))], 'FontSize', legendFontSize);
            set(lgd3, 'Color', 'none');

            % 4th subplot: LowAlphaRSigRT vs HighAlphaRSigRT
            subplot(2, 2, 4);
            histogram(LowAlphaRSigRT, 'Normalization', 'count', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'BinWidth', binWidth, 'FaceColor', darkgreen);
            hold on;
            histogram(HighAlphaRSigRT, 'Normalization', 'count', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'BinWidth', binWidth, 'FaceColor', grey);
            hold off;
            title(['WPLI-RT ' Param.Corr ' r (p<' num2str(alpha) ')']);
            ylabel('Ch Pair Count');
            xlabel(['r - Welch: t(' num2str(stats8.df) ')=' num2str(stats8.tstat, '%.2f') ', p=' num2str(p8, '%.3f')], 'FontSize', fontSize);
            lgd4 = legend(['LowALphaRSigRT nChPairs=' num2str(numel(LowAlphaRSigRT))], ['HighALphaRSigRT nChPairs=' num2str(numel(HighAlphaRSigRT))], 'FontSize', legendFontSize);
            set(lgd4, 'Color', 'none');


            % Add a title over all subplots
            sgtitle(['Low vs High Alpha ' TrialTypeName{iTrialType} ' ' FunGroupName{iFreqFunc}]);

            % Save the figure as a single file
            saveas(gcf, [SaveTempGroup TrialTypeName{iTrialType} FunGroupName{iFreqFunc} 'LA-HA_Histograms.png']);
            saveas(gcf, [SaveTempGroup TrialTypeName{iTrialType} FunGroupName{iFreqFunc} 'LA-HA_Histograms.svg'], 'svg');


            %% All Alpha vs Low Alpha histogrma and ttest comparison
            % Open a text file for writing
            fileID = fopen([SaveTempGroup TrialTypeName{iTrialType} FunGroupName{iFreqFunc} '_AA-LA-ttest_results.txt'], 'w');

            % Remove NaN values and convert to vectors
            LowAlphaRAcc_non_nan = LowAlphaRAcc(~isnan(LowAlphaRAcc));
            AllAlphaRAcc_non_nan = AllAlphaRAcc(~isnan(AllAlphaRAcc));
            LowAlphaRRT_non_nan = LowAlphaRRT(~isnan(LowAlphaRRT));
            AllAlphaRRT_non_nan = AllAlphaRRT(~isnan(AllAlphaRRT));

            % Perform single-sample t-tests (t-Test for Mean Equal to Zero)
            [h1, p1, ci1, stats1] = ttest(LowAlphaRAcc_non_nan);
            [h2, p2, ci2, stats2] = ttest(AllAlphaRAcc_non_nan);
            [h3, p3, ci3, stats3] = ttest(LowAlphaRRT_non_nan);
            [h4, p4, ci4, stats4] = ttest(AllAlphaRRT_non_nan);

            % Perform paired t-tests
            [h5, p5, ci5, stats5] = ttest(LowAlphaRAcc_non_nan, AllAlphaRAcc_non_nan);
            [h6, p6, ci6, stats6] = ttest(LowAlphaRRT_non_nan, AllAlphaRRT_non_nan);

            % Perform unpaired t-tests (Welch's t-test)
            [h7, p7, ci7, stats7] = ttest2(LowAlphaRSigAcc, AllAlphaRSigAcc, 'Vartype', 'unequal'); % This will perform Welch's t-test, which does not assume equal variances between the two groups.
            [h8, p8, ci8, stats8] = ttest2(LowAlphaRSigRT, AllAlphaRSigRT, 'Vartype', 'unequal');

            % Write results to file
            fprintf(fileID, [TrialTypeName{iTrialType} ' ' FunGroupName{iFreqFunc}]);
            fprintf(fileID, 'Single-sample t-tests (t-Test for Mean Equal to Zero):\n');
            fprintf(fileID, 'LowAlphaRAcc: p-value = %.4f, t-statistic = %.4f, df = %.2f, CI = [%.4f, %.4f]\n', p1, stats1.tstat, stats1.df, ci1(1), ci1(2));
            fprintf(fileID, 'AllAlphaRAcc: p-value = %.4f, t-statistic = %.4f, df = %.2f, CI = [%.4f, %.4f]\n', p2, stats2.tstat, stats2.df, ci2(1), ci2(2));
            fprintf(fileID, 'LowAlphaRRT: p-value = %.4f, t-statistic = %.4f, df = %.2f, CI = [%.4f, %.4f]\n', p3, stats3.tstat, stats3.df, ci3(1), ci3(2));
            fprintf(fileID, 'AllAlphaRRT: p-value = %.4f, t-statistic = %.4f, df = %.2f, CI = [%.4f, %.4f]\n', p4, stats4.tstat, stats4.df, ci4(1), ci4(2));

            fprintf(fileID, '\nPaired t-tests (All channel pairs):\n');
            fprintf(fileID, 'LowAlphaRAcc vs AllAlphaRAcc: p-value = %.4f, t-statistic = %.4f, df = %.2f, CI = [%.4f, %.4f]\n', p5, stats5.tstat, stats5.df, ci5(1), ci5(2));
            fprintf(fileID, 'LowAlphaRRT vs AllAlphaRRT: p-value = %.4f, t-statistic = %.4f, df = %.2f, CI = [%.4f, %.4f]\n', p6, stats6.tstat, stats6.df, ci6(1), ci6(2));

            fprintf(fileID, '\nUnpaired t-tests (Welch''s t-test on Significant Channel Pairs p<%.4f):\n',alpha);
            fprintf(fileID, 'LowAlphaRSigAcc vs AllAlphaRSigAcc: p-value = %.4f, t-statistic = %.4f, df = %.2f, CI = [%.4f, %.4f]\n', p7, stats7.tstat, stats7.df, ci7(1), ci7(2));
            fprintf(fileID, 'LowALphaRSigRT vs AllAlphaRSigRT: p-value = %.4f, t-statistic = %.4f, df = %.2f, CI = [%.4f, %.4f]\n', p8, stats8.tstat, stats8.df, ci8(1), ci8(2));

            % Close the file
            fclose(fileID);

            %% Create overlapping histograms for LowAlphaRAcc and AllAlphaRAcc
            figure;

            % Font size for text annotations
            fontSize = 6;  % You can adjust this value to make the font smaller or larger
            legendFontSize = 5;
            yOffset = -0.25;
            binWidth = 0.05;  % Bin size for histograms

            % 1st subplot: LowAlphaRAcc vs AllAlphaRAcc
            subplot(2, 2, 1);
            histogram(LowAlphaRAcc_non_nan, 'Normalization', 'count', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'BinWidth', binWidth);
            hold on;
            histogram(AllAlphaRAcc_non_nan , 'Normalization', 'count', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'BinWidth', binWidth, 'FaceColor', 'green');
            hold off;
            title(['WPLI-Acc ' Param.Corr ' r']);
            ylabel('Ch Pair Count');
            xlabel(['r - Paired: t(' num2str(stats5.df) ')=' num2str(stats5.tstat, '%.2f') ', p=' num2str(p5, '%.3f')] , 'FontSize', fontSize);
            % legend('LowAlphaRAcc', 'AllAlphaRAcc');
            lgd1 = legend({['LowAlphaRAcc: t(' num2str(stats1.df) ')=' num2str(stats1.tstat, '%.2f') ', p=' num2str(p1, '%.3f')], ...
                ['AllAlphaRAcc: t(' num2str(stats2.df) ')=' num2str(stats2.tstat, '%.2f') ', p=' num2str(p2, '%.3f')], ...
                ['Paired: t(' num2str(stats5.df) ')=' num2str(stats5.tstat, '%.2f') ', p=' num2str(p5, '%.3f')]}, 'FontSize', legendFontSize);
            set(lgd1, 'Color', 'none');

            % 2nd subplot: LowAlphaRRT vs AllAlphaRRT
            subplot(2, 2, 2);
            histogram(LowAlphaRRT_non_nan , 'Normalization', 'count', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'BinWidth', binWidth);
            hold on;
            histogram(AllAlphaRRT_non_nan , 'Normalization', 'count', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'BinWidth', binWidth, 'FaceColor', 'green');
            hold off;
            title(['WPLI-RT ' Param.Corr ' r']);
            ylabel('Ch Pair Count');
            xlabel(['r  - Paired: t(' num2str(stats6.df) ')=' num2str(stats6.tstat, '%.2f') ', p=' num2str(p6, '%.3f')], 'FontSize', fontSize);
            % legend('LowAlphaRRT', 'AllAlphaRRT');
            lgd2= legend({['LowAlphaRRT: t(' num2str(stats3.df) ')=' num2str(stats3.tstat, '%.2f') ', p=' num2str(p3, '%.3f')], ...
                ['AllAlphaRRT: t(' num2str(stats4.df) ')=' num2str(stats4.tstat, '%.2f') ', p=' num2str(p4, '%.3f')], ...
                ['Paired: t(' num2str(stats6.df) ')=' num2str(stats6.tstat, '%.2f') ', p=' num2str(p6, '%.3f')]}, 'FontSize', legendFontSize);
            set(lgd2, 'Color', 'none');


            % 3rd subplot: LowAlphaRSigAcc vs AllAlphaRSigAcc
            subplot(2, 2, 3);
            histogram(LowAlphaRSigAcc, 'Normalization', 'count', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'BinWidth', binWidth);
            hold on;
            histogram(AllAlphaRSigAcc, 'Normalization', 'count', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'BinWidth', binWidth, 'FaceColor', 'green');
            hold off;
            title(['WPLI-Acc ' Param.Corr ' r (p<' num2str(alpha) ')']);
            ylabel('Ch Pair Count');
            xlabel(['r - Welch: t(' num2str(stats7.df) ')=' num2str(stats7.tstat, '%.2f') ', p=' num2str(p7, '%.3f')] , 'FontSize', fontSize);
            lgd3 = legend(['LowAlphaRSigAcc nChPairs=' num2str(numel(LowAlphaRSigAcc))], ['AllAlphaRSigAcc nChPairs=' num2str(numel(AllAlphaRSigAcc))], 'FontSize', legendFontSize);
            set(lgd3, 'Color', 'none');

            % 4th subplot: LowAlphaRSigRT vs AllAlphaRSigRT
            subplot(2, 2, 4);
            histogram(LowAlphaRSigRT, 'Normalization', 'count', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'BinWidth', binWidth);
            hold on;
            histogram(AllAlphaRSigRT, 'Normalization', 'count', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'BinWidth', binWidth, 'FaceColor', 'green');
            hold off;
            title(['WPLI-RT ' Param.Corr ' r (p<' num2str(alpha) ')']);
            ylabel('Ch Pair Count');
            xlabel(['r - Welch: t(' num2str(stats8.df) ')=' num2str(stats8.tstat, '%.2f') ', p=' num2str(p8, '%.3f')], 'FontSize', fontSize);
            lgd4 = legend(['LowALphaRSigRT nChPairs=' num2str(numel(LowAlphaRSigRT))], ['AllAlphaRSigRT nChPairs=' num2str(numel(AllAlphaRSigRT))], 'FontSize', legendFontSize);
            set(lgd4, 'Color', 'none');

            % Add a title over all subplots
            sgtitle(['Low vs All Alpha ' TrialTypeName{iTrialType} ' ' FunGroupName{iFreqFunc}]);

            % Save the figure as a single file
            saveas(gcf, [SaveTempGroup TrialTypeName{iTrialType} FunGroupName{iFreqFunc} 'AA-LA_Histograms.png']);
            saveas(gcf, [SaveTempGroup TrialTypeName{iTrialType} FunGroupName{iFreqFunc} 'AA-LA_Histograms.svg']);

            %% All Alpha vs High Alpha histogrma and ttest comparison
            % Open a text file for writing
            fileID = fopen([SaveTempGroup TrialTypeName{iTrialType} FunGroupName{iFreqFunc} '_HA-AA-ttest_results.txt'], 'w');

            % Remove NaN values and convert to vectors
            HighAlphaRAcc_non_nan = HighAlphaRAcc(~isnan(HighAlphaRAcc));
            AllAlphaRAcc_non_nan = AllAlphaRAcc(~isnan(AllAlphaRAcc));
            HighAlphaRRT_non_nan = HighAlphaRRT(~isnan(HighAlphaRRT));
            AllAlphaRRT_non_nan = AllAlphaRRT(~isnan(AllAlphaRRT));

            % Perform single-sample t-tests (t-Test for Mean Equal to Zero)
            [h1, p1, ci1, stats1] = ttest(HighAlphaRAcc_non_nan);
            [h2, p2, ci2, stats2] = ttest(AllAlphaRAcc_non_nan);
            [h3, p3, ci3, stats3] = ttest(HighAlphaRRT_non_nan);
            [h4, p4, ci4, stats4] = ttest(AllAlphaRRT_non_nan);

            % Perform paired t-tests
            [h5, p5, ci5, stats5] = ttest(HighAlphaRAcc_non_nan, AllAlphaRAcc_non_nan);
            [h6, p6, ci6, stats6] = ttest(HighAlphaRRT_non_nan, AllAlphaRRT_non_nan);

            % Perform unpaired t-tests (Welch's t-test)
            [h7, p7, ci7, stats7] = ttest2(HighAlphaRSigAcc, AllAlphaRSigAcc, 'Vartype', 'unequal'); % This will perform Welch's t-test, which does not assume equal variances between the two groups.
            [h8, p8, ci8, stats8] = ttest2(HighAlphaRSigRT, AllAlphaRSigRT, 'Vartype', 'unequal');

            % Write results to file
            fprintf(fileID, [TrialTypeName{iTrialType} ' ' FunGroupName{iFreqFunc}]);
            fprintf(fileID, 'Single-sample t-tests (t-Test for Mean Equal to Zero):\n');
            fprintf(fileID, 'HighAlphaRAcc: p-value = %.4f, t-statistic = %.4f, df = %.2f, CI = [%.4f, %.4f]\n', p1, stats1.tstat, stats1.df, ci1(1), ci1(2));
            fprintf(fileID, 'AllAlphaRAcc: p-value = %.4f, t-statistic = %.4f, df = %.2f, CI = [%.4f, %.4f]\n', p2, stats2.tstat, stats2.df, ci2(1), ci2(2));
            fprintf(fileID, 'HighAlphaRRT: p-value = %.4f, t-statistic = %.4f, df = %.2f, CI = [%.4f, %.4f]\n', p3, stats3.tstat, stats3.df, ci3(1), ci3(2));
            fprintf(fileID, 'AllAlphaRRT: p-value = %.4f, t-statistic = %.4f, df = %.2f, CI = [%.4f, %.4f]\n', p4, stats4.tstat, stats4.df, ci4(1), ci4(2));

            fprintf(fileID, '\nPaired t-tests (All channel pairs):\n');
            fprintf(fileID, 'HighAlphaRAcc vs AllAlphaRAcc: p-value = %.4f, t-statistic = %.4f, df = %.2f, CI = [%.4f, %.4f]\n', p5, stats5.tstat, stats5.df, ci5(1), ci5(2));
            fprintf(fileID, 'HighAlphaRRT vs AllAlphaRRT: p-value = %.4f, t-statistic = %.4f, df = %.2f, CI = [%.4f, %.4f]\n', p6, stats6.tstat, stats6.df, ci6(1), ci6(2));

            fprintf(fileID, '\nUnpaired t-tests (Welch''s t-test on Significant Channel Pairs p<%.4f):\n',alpha);
            fprintf(fileID, 'HighAlphaRSigAcc vs AllAlphaRSigAcc: p-value = %.4f, t-statistic = %.4f, df = %.2f, CI = [%.4f, %.4f]\n', p7, stats7.tstat, stats7.df, ci7(1), ci7(2));
            fprintf(fileID, 'HighAlphaRSigRT vs AllAlphaRSigRT: p-value = %.4f, t-statistic = %.4f, df = %.2f, CI = [%.4f, %.4f]\n', p8, stats8.tstat, stats8.df, ci8(1), ci8(2));

            % Close the file
            fclose(fileID);

            %% Create overlapping histograms for HighAlphaRAcc and AllAlphaRAcc
            figure;

            % Font size for text annotations
            fontSize = 6;  % You can adjust this value to make the font smaller or larger
            legendFontSize = 5;
            yOffset = -0.25;
            binWidth = 0.05;  % Bin size for histograms

            % 1st subplot: HighAlphaRAcc vs AllAlphaRAcc
            subplot(2, 2, 1);
            histogram(HighAlphaRAcc_non_nan, 'Normalization', 'count', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'BinWidth', binWidth, 'FaceColor', '#D95319');
            hold on;
            histogram(AllAlphaRAcc_non_nan , 'Normalization', 'count', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'BinWidth', binWidth, 'FaceColor', 'green');
            hold off;
            title(['WPLI-Acc ' Param.Corr ' r']);
            ylabel('Ch Pair Count');
            xlabel(['r - Paired: t(' num2str(stats5.df) ')=' num2str(stats5.tstat, '%.2f') ', p=' num2str(p5, '%.3f')] , 'FontSize', fontSize);
            % legend('HighAlphaRAcc', 'AllAlphaRAcc');
            lgd1 = legend({['HighAlphaRAcc: t(' num2str(stats1.df) ')=' num2str(stats1.tstat, '%.2f') ', p=' num2str(p1, '%.3f')], ...
                ['AllAlphaRAcc: t(' num2str(stats2.df) ')=' num2str(stats2.tstat, '%.2f') ', p=' num2str(p2, '%.3f')], ...
                ['Paired: t(' num2str(stats5.df) ')=' num2str(stats5.tstat, '%.2f') ', p=' num2str(p5, '%.3f')]}, 'FontSize', legendFontSize);
            set(lgd1, 'Color', 'none');

            % 2nd subplot: HighAlphaRRT vs AllAlphaRRT
            subplot(2, 2, 2);
            histogram(HighAlphaRRT_non_nan , 'Normalization', 'count', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'BinWidth', binWidth, 'FaceColor', '#D95319');
            hold on;
            histogram(AllAlphaRRT_non_nan , 'Normalization', 'count', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'BinWidth', binWidth, 'FaceColor', 'green');
            hold off;
            title(['WPLI-RT ' Param.Corr ' r']);
            ylabel('Ch Pair Count');
            xlabel(['r  - Paired: t(' num2str(stats6.df) ')=' num2str(stats6.tstat, '%.2f') ', p=' num2str(p6, '%.3f')], 'FontSize', fontSize);
            % legend('HighAlphaRRT', 'AllAlphaRRT');
            lgd2= legend({['HighAlphaRRT: t(' num2str(stats3.df) ')=' num2str(stats3.tstat, '%.2f') ', p=' num2str(p3, '%.3f')], ...
                ['AllAlphaRRT: t(' num2str(stats4.df) ')=' num2str(stats4.tstat, '%.2f') ', p=' num2str(p4, '%.3f')], ...
                ['Paired: t(' num2str(stats6.df) ')=' num2str(stats6.tstat, '%.2f') ', p=' num2str(p6, '%.3f')]}, 'FontSize', legendFontSize);
            set(lgd2, 'Color', 'none');


            % 3rd subplot: HighAlphaRSigAcc vs AllAlphaRSigAcc
            subplot(2, 2, 3);
            histogram(HighAlphaRSigAcc, 'Normalization', 'count', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'BinWidth', binWidth, 'FaceColor', '#D95319');
            hold on;
            histogram(AllAlphaRSigAcc, 'Normalization', 'count', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'BinWidth', binWidth, 'FaceColor', 'green');
            hold off;
            title(['WPLI-Acc ' Param.Corr ' r (p<' num2str(alpha) ')']);
            ylabel('Ch Pair Count');
            xlabel(['r - Welch: t(' num2str(stats7.df) ')=' num2str(stats7.tstat, '%.2f') ', p=' num2str(p7, '%.3f')] , 'FontSize', fontSize);
            lgd3 = legend(['HighAlphaRSigAcc nChPairs=' num2str(numel(HighAlphaRSigAcc))], ['AllAlphaRSigAcc nChPairs=' num2str(numel(AllAlphaRSigAcc))], 'FontSize', legendFontSize);
            set(lgd3, 'Color', 'none');

            % 4th subplot: HighAlphaRSigRT vs AllAlphaRSigRT
            subplot(2, 2, 4);
            histogram(HighAlphaRSigRT, 'Normalization', 'count', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'BinWidth', binWidth, 'FaceColor', '#D95319');
            hold on;
            histogram(AllAlphaRSigRT, 'Normalization', 'count', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'BinWidth', binWidth, 'FaceColor', 'green');
            hold off;
            title(['WPLI-RT ' Param.Corr ' r (p<' num2str(alpha) ')']);
            ylabel('Ch Pair Count');
            xlabel(['r - Welch: t(' num2str(stats8.df) ')=' num2str(stats8.tstat, '%.2f') ', p=' num2str(p8, '%.3f')], 'FontSize', fontSize);
            lgd4 = legend(['HighAlphaRSigRT nChPairs=' num2str(numel(HighAlphaRSigRT))], ['AllAlphaRSigRT nChPairs=' num2str(numel(AllAlphaRSigRT))], 'FontSize', legendFontSize);
            set(lgd4, 'Color', 'none');

            % Add a title over all subplots
            sgtitle(['High vs All Alpha - ' TrialTypeName{iTrialType} ' ' FunGroupName{iFreqFunc}]);

            % Save the figure as a single file
            saveas(gcf, [SaveTempGroup TrialTypeName{iTrialType} FunGroupName{iFreqFunc} 'AA-HA_Histograms.png']);


        end
    end
end
%% Fig4E and fig4F stats - chi-squared OLD
% Total channel pairs
total_pairs = 496;

% Updated Accuracy table data
accuracy_counts = [10, 4; 14, 21]; % lower alpha, upper alpha for 40vLight and 40vRandom
reaction_time_counts = [6, 4; 14, 17]; % lower alpha, upper alpha for 40vLight and 40vRandom

% Extract groups
acc_40vLight = accuracy_counts(:, 1); % 40vLight
acc_40vRandom = accuracy_counts(:, 2); % 40vRandom

rt_40vLight = reaction_time_counts(:, 1); % 40vLight
rt_40vRandom = reaction_time_counts(:, 2); % 40vRandom

% Perform chi-squared tests
[chi_stat_acc_light, p_val_acc_light] = chi_squared_test(acc_40vLight);
[chi_stat_acc_random, p_val_acc_random] = chi_squared_test(acc_40vRandom);

[chi_stat_rt_light, p_val_rt_light] = chi_squared_test(rt_40vLight);
[chi_stat_rt_random, p_val_rt_random] = chi_squared_test(rt_40vRandom);

% Display the results
disp('Chi-squared test for Accuracy (40vLight):');
disp(['Chi-stat: ', num2str(chi_stat_acc_light)]);
disp(['P-value: ', num2str(p_val_acc_light)]);

disp('Chi-squared test for Accuracy (40vRandom):');
disp(['Chi-stat: ', num2str(chi_stat_acc_random)]);
disp(['P-value: ', num2str(p_val_acc_random)]);

disp('Chi-squared test for Reaction Time (40vLight):');
disp(['Chi-stat: ', num2str(chi_stat_rt_light)]);
disp(['P-value: ', num2str(p_val_rt_light)]);

disp('Chi-squared test for Reaction Time (40vRandom):');
disp(['Chi-stat: ', num2str(chi_stat_rt_random)]);
disp(['P-value: ', num2str(p_val_rt_random)]);

%% Fig4F and Fi4G stats - chisquare MS
% Define observed frequency tables

% 40 Hz & Light - Accuracy (ACC)
obs_Acc_Light = [10, 4; 2, 10];

% 40 Hz & Random - Accuracy (ACC)
obs_Acc_Random = [14, 20; 1, 1];

% 40 Hz & Light - Reaction Time (RT)
obs_RT_Light = [6, 5; 4, 10]; 

% 40 Hz & Random - Reaction Time (RT)
obs_RT_Random = [14, 1; 17, 0];

% Perform chi-square tests and calculate effect sizes
[h_Acc_Light, p_Acc_Light, stats_Acc_Light, v_Acc_Light] = chi2test_effect(obs_Acc_Light);
[h_Acc_Random, p_Acc_Random, stats_Acc_Random, v_Acc_Random] = chi2test_effect(obs_Acc_Random);
[h_RT_Light, p_RT_Light, stats_RT_Light, v_RT_Light] = chi2test_effect(obs_RT_Light);
[h_RT_Random, p_RT_Random, stats_RT_Random, v_RT_Random] = chi2test_effect(obs_RT_Random);

% Display results
disp('Chi-square test results and effect sizes (Cramér''s V):');
disp(['p-value (Accuracy, 40 & Light): ', num2str(p_Acc_Light), ', Cramér''s V: ', num2str(v_Acc_Light)]);
disp(['p-value (Accuracy, 40 & Random): ', num2str(p_Acc_Random), ', Cramér''s V: ', num2str(v_Acc_Random)]);
disp(['p-value (RT, 40 & Light): ', num2str(p_RT_Light), ', Cramér''s V: ', num2str(v_RT_Light)]);
disp(['p-value (RT, 40 & Random): ', num2str(p_RT_Random), ', Cramér''s V: ', num2str(v_RT_Random)]);

%%  Save entire workspace at end of script - This creates a very large file and takes a long time to save.  Please comment out if not needed.
tic
saveDate = datestr(datetime, 'yy-mm-dd_HHMMSSFFF');
S4EndSaveFileName =['S4EndData_' saveDate '_' scriptName '.mat'];
S4End_Save_Path = [SavePath S4EndSaveFileName];
save(S4End_Save_Path,'-v7.3')  % Check this
toc
%% Create output text
endDateAndTime = char(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss'));
endDuration = toc(S4Start);
SummaryTextOutput = sprintf('%s\nStep4 completed on: %s\nRun duration: %g seconds.', SummaryTextOutput, endDateAndTime, endDuration);

%Display the output message in the command window
disp(SummaryTextOutput);

%Create a .txt file with the output message
filePath = fullfile(SaveFolder, [scriptName, endDuration '_CompletionSummary.txt']);

fileID = fopen(filePath, 'w');
if fileID == -1
    error('Failed to open or create the file: %s', filePath);
else
    fprintf(fileID, '%s', SummaryTextOutput);
    fclose(fileID);
    fprintf('File created and written: %s\n', filePath);
end

%% Functions - must be at the end of the script
%% Chi Squared test OLD
function [chi2, p, V, dof, grandTotal] = chi_squared_test(observed)
% chi_squared_test performs a chi-squared test on a 2x2 contingency table.
%
% Inputs:
%   observed - A 2x2 matrix of observed frequencies.
%
% Outputs:
%   chi2       - Chi-squared statistic.
%   p          - p-value of the test.
%   V          - Cramér's V (effect size).
%   dof        - Degrees of freedom.
%   grandTotal - Total number of observations in the table.
%
% Example:
%   observed = [93, 403; 32, 464];
%   [chi2, p, V, dof, grandTotal] = chi_squared_test(observed);

    % Ensure the input is a 2x2 matrix
    if ~isequal(size(observed), [2, 2])
        error('Input observed must be a 2x2 matrix.');
    end

    % Calculate row totals, column totals, and grand total
    rowTotals = sum(observed, 2);
    colTotals = sum(observed, 1);
    grandTotal = sum(observed(:));

    % Calculate expected frequencies
    expected = (rowTotals * colTotals) / grandTotal;

    % Compute chi-squared statistic
    chi2 = sum((observed - expected).^2 ./ expected, 'all');

    % Degrees of freedom
    dof = (size(observed, 1) - 1) * (size(observed, 2) - 1);

    % p-value from chi-squared distribution
    p = 1 - chi2cdf(chi2, dof);

    % Calculate Cramér's V (effect size)
    k = min(size(observed)); % Smaller dimension of the table
    V = sqrt(chi2 / (grandTotal * (k - 1)));

    % Display results
    fprintf('Chi-squared Test Results:\n');
    fprintf('Chi-squared (X^2): %.2f\n', chi2);
    fprintf('Degrees of Freedom (dof): %d\n', dof);
    fprintf('Total Observations (grandTotal): %d\n', grandTotal);
    fprintf('p-value: %.4f\n', p);
    fprintf('Cramér''s V (Effect size): %.4f\n', V);
end
%% Chi Squared test MS
function [hypothesis_test_results, p, stats, V] = chi2test_effect(obs)
    % Compute expected values
    expected = sum(obs,2) * sum(obs,1) / sum(obs(:));
    
    % Compute Chi-square statistic
    chi2_stat = sum((obs - expected).^2 ./ expected, 'all');
    
    % Degrees of freedom
    df = (size(obs,1)-1) * (size(obs,2)-1);
    
    % Compute p-value
    p = 1 - chi2cdf(chi2_stat, df);
    
    % Hypothesis test result
    hypothesis_test_results = p < 0.05;  % Reject H0 if p < 0.05
    
    % Compute Cramér's V
    N = sum(obs(:)); % Total sample size
    k = min(size(obs)); % Smallest dimension (rows or columns)
    V = sqrt(chi2_stat / (N * (k - 1))); 
    
    % Store stats
    stats.chi2 = chi2_stat;
    stats.df = df;
end
%% Function to apply FDR correction
function [adjusted_p_values] = fdr_correction(p_values)
    [sorted_p, sort_idx] = sort(p_values); % Sort p-values
    n = length(p_values); % Number of comparisons
    fdr_thresholds = (1:n) / n * 0.05; % FDR thresholds at alpha = 0.05
    adjusted_p = min(1, cummin((sorted_p ./ fdr_thresholds) .* n)); % Adjusted p-values
    adjusted_p_values = zeros(size(p_values));
    adjusted_p_values(sort_idx) = adjusted_p; % Rearrange to original order
end