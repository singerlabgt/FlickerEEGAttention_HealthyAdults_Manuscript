function [outputArg1] = calculateWPLI(DataA,DataB)
    %calculateWPLI Calculations WPLI from DataA and DataB
    %   DataA corresponds to SubA and Chi; DataB corresponds to SubB and Chj

    %%  Calculate for pair of channels
    % if sum(sum(isnan(tempSig1)))>1||sum(sum(isnan(tempSig2)))>1
    %     continue
    % end
    % tic
    nTrials = size(DataA,2); % get number of Trials.  # of Trials should be the same for DataA and DataB
    for iTrial=1:nTrials  % loop by trial
        if nTrials>1
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
    save([SaveTemp 'Ch' num2str(iCh) 'Ch' num2str(jCh) '.mat'],'TrialSpec','ValidIndex','psdParameter');


    for iTrialType=3%1:length(TrialType)  % group: hit, miss, etc.
        TrialI=[];
        parfor j=1:length(TrialType{iTrialType})
            TrialI=union(TrialI,find(TrialTypeTemp==TrialType{iTrialType}(j)));
        end
        TrialI=intersect(TrialI,ValidIndex);
        if ~isempty(TrialI)
            %                    CohGroup{iGroup,iFile}{iCh,jCh}=Coh_TrialIndex(TrialSpec1,TrialI);
            CohGroup{iTrialType,iFile}{iCh,jCh}=crossspec_TrialIndex(TrialSpec,TrialI);
        end
    end