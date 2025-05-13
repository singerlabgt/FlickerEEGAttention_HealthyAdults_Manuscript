import os
import numpy as np
import pandas as pd
from plots_for_AttentionTask import *
from MaxStreak import *
from StatsTest_Attention import *
# from Window_Size_Violin_by_group_for_onefolder import *
from Acc_RT_Overlay_for_onefolder import *
import datetime

### Script to plot all figures
x = datetime.datetime.now()
runTime = x.strftime("%y%m%d-%H%M%S")

outputFolder = os.path.abspath(os.path.join(os.path.dirname(__file__), 'Output and Figures' "\\"))

subjects_dataframe = pd.read_pickle("240112-095859_AttnTaskDF.pkl")  # df = data_frame
Acc_Cut_dataframe = cutLowAccSubs(subjects_dataframe, 0.8)

folder_name = f'{runTime}_AttentionTaskFigures'

os.makedirs(folder_name, exist_ok=True)
plot_filename = 'AttentionTaskData.xlsx'
plot_path = os.path.join(folder_name, plot_filename)
subjects_dataframe.to_excel(plot_path)

AccVPlot(Acc_Cut_dataframe, folder_name)
RtVPlot(Acc_Cut_dataframe, folder_name)
RtvVPlot(Acc_Cut_dataframe, folder_name)
RtvAccScatter(Acc_Cut_dataframe, folder_name)
RtAccScatter(Acc_Cut_dataframe, folder_name)
RTvsAcc_heatmap(Acc_Cut_dataframe, folder_name)

AllMissedStreakNotGrouped(Acc_Cut_dataframe, folder_name)
HighestMissedStreakNotGrouped(Acc_Cut_dataframe, folder_name)
HMSGroupAndScattervAcc(Acc_Cut_dataframe, folder_name)
AllMissedStreaksGrouped(Acc_Cut_dataframe, folder_name)
AllStatsTests(Acc_Cut_dataframe, folder_name)

# create_three_violins(window_size=10, folname=folder_name)
# create_three_violins(window_size=25, folname=folder_name)
# create_three_violins(window_size=40, folname=folder_name)
# create_and_overlay_accuracy_rt_plots(window_size=10, folname=folder_name)
# create_and_overlay_accuracy_rt_plots(window_size=25, folname=folder_name)
# create_and_overlay_accuracy_rt_plots(window_size=40, folname=folder_name)
