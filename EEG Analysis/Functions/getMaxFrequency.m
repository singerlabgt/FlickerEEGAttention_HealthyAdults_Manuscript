BOI=[1 4 8 13 30 39.5 43;4 8 13 30 37 41.5 100];
        BName={'Delta','Theta','Alpha','Beta','Gamma-1','Gamma-E','Gamma-2'};
        BName2={'1-4 Hz','4-8 Hz','8-13Hz','13-30 Hz','30-37 Hz','39-41 Hz','43-100 Hz'};
iBOI=6; %currently set to Gamma-E when iBOI=6
NeedI=find(Fplot>=BOI(1,iBOI)&Fplot<=BOI(2,iBOI));
%%when running not the peak remove the second to last part with []
[varia, indices]=FunGroup{iFun}(LogPSD(SubjG{iG},NeedI,EEGchInd,iCom),[],2);
maxFrequencies=Fplot(NeedI(indices)); % should give you the frequency for each channel where the max occurs for each participant 

