subj = 48;
isHit = (FileStruct(subj).dotsynch(:,3)==1);  % get logical index for all trials that are Hits (misses and premature hits will = 0)
colorchangeTimes = FileStruct(subj).dotsynch(:,1); % get time of color change for Hit trials only        
subjRTtimes = FileStruct(subj).dotsynch(:,2); % get RT time for Hit trials (this should ignore premature hits)
subjRTduration =  subjRTtimes - colorchangeTimes; % subtract RT time from color change time to get duration of RT (AKA the reaction time)
subjRTdurInSecs = subjRTduration/512; % sample rate is generally 512 samples per second
indivSubjsAvgRT = mean(subjRTdurInSecs);
AvgSubjsAvgRT = mean(SubjsAvgRT)