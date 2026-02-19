function data_table=RateHist_LME1(DataGroup,SavePathFile)

%%%%This is designed for ratehistogram plot, repeat or non-repeat measures;
%%%%output{1} a is for figure legend
%%%%output{2} Outstats is the statistics

Rfolder='Y:\singer\LuZhang\Project6-EEG\Analysis\RFolder\'; % R code - change path, copy code to new folder

Session=1:size(DataGroup{1},2);
% DataGroup=DataPlot;

DataAll=zeros(size(DataGroup{1},1),size(DataGroup{1},2),length(DataGroup));
GroupAll=DataAll;
SessionAll=DataAll;

for iG=1:length(DataGroup)
    DataAll(:,:,iG)=DataGroup{iG};
    GroupAll(:,:,iG)=iG;
end

SubjAll=repmat([1:size(DataAll,1)]',1,length(Session),length(DataGroup));
SessionAll=repmat(Session(:)',size(DataAll,1),1,length(DataGroup));


tbl = array2table([DataAll(:),SubjAll(:),SessionAll(:),GroupAll(:)],'VariableNames',{'LFP','Subj','Fre','RTGroup'});
writetable(tbl,[Rfolder 'tempData.csv'])  % 
RunRcode([Rfolder 'LME1RTgroupImpact.R'],'C:\Program Files\R\R-4.0.2\bin'); % Install R 4.0.2 - Must put into this folder
% above runs R code in R folder

%l1=csvread([Rfolder 'tempResult.csv']);
lme=readcell([Rfolder 'tempResult.csv']);

variable_names = lme(1,:);
data = lme(2:end,:);
data_table = cell2table(data, 'VariableNames', variable_names);

% delete([Rfolder 'tempResult.csv']);
delete([Rfolder 'tempData.csv']);
movefile([Rfolder 'tempResult.csv'],SavePathFile);