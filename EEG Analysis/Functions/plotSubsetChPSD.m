%% plotSubsetChPSDAllParticipants.m Plot PSD of Frontal Channels
%   Plots PSDs for a subset of channels (e.g frontal channels only) for
%   each subject. Matty Attokaren 4/17/24

frontalHalfChs = [1:6,25:31]; % 13 Chns, All F chns, in front half of scalp
frontal3rowsChs = [1:4, 27:31]; % 9 chns, no FC channels Fp1AF3F7F3F4F8AF4Fp2Fz
occipital_Channels = [14:18]; %#ok<NBRAK2> % 5 occiptal channels, PO3, O1, Oz, O2, PO4
extractedChannels = frontal3rowsChs;  % The subset of channels that are extracted.

filenameSuffix = '_Frontal_PSD_PostICA_BandVLines';  % CHANGE THIS TO BE UNIQUE

% epochBinRange = [-4000.0 1200.0];
%% Dataset order
nDatasets = size(ALLEEG,2);  % number of datasets loaded into EEGLAB
if nDatasets == 0
    error('No datasets loaded.  Please load datasets into EEGLAB and try again.')
end
lastDataset = nDatasets; % set as nDatasets to end at the last participant loaded in EEGLAB
firstDataset = 1; % set as 1 to start from the first participant loaded in EEGLAB
if firstDataset>nDatasets  % incase firstdata set is hardcoded
    firstDataset = 1;
end
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',firstDataset,'study',0);

%% Iterate thru all datasets
for iDataset = firstDataset:lastDataset
    %% Create new EEG structure with just the channels of interest
    EEG_subset = pop_select( EEG, 'channel',extractedChannels); % Extract chns from EEG
    % EEG_frontal = pop_epochbin( EEG_frontal , epochBinRange,  'pre');

    %% Plot
    figure; pop_spectopo(EEG_subset, 1, [0 max(EEG_subset.times)], 'EEG' , 'freqrange',[1 61],'electrodes','off');
    chLabels = [EEG_subset.chanlocs.labels];  % Gets channel labels e.g [PO3O1OzO2PO4]
    title([EEG_subset.setname ' - ' chLabels])
    % [spectopo_outputs] = pop_spectopo( EEG, dataflag, timerange, ...
    % process, 'key', 'val',...); % returns spectopo() outputs

    BOIticks = [1 4 8 13 30 40 60];
    set(gca,'xtick',BOIticks);
    ax=gca;
    ax.XGrid = 'on';
    ax.YGrid = 'off';

    %% Save Figure
    saveas(gcf, ['IndivPSDs\' EEG_subset.setname filenameSuffix], 'jpg'); %
    % gcf -> get current Figure!

    %% Check if last file, ifnot load next set. Display completion message.
    if iDataset ~= lastDataset
        % Set next dataset as current dataset, if not at the end
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, iDataset,'retrieve',iDataset+1,'study',0);
    end

    disp(['Plot completed for: ' EEG_subset.filename]);
    close
end