import os
import numpy as np
import pandas as pd
from MaxStreak import *
from StatsTest_Attention import *
import datetime
from plots_for_AttentionTask import *

from plots_for_Fig1_illustrator import *

### Script to plot all graphs for Figure 1 - publication format

x = datetime.datetime.now()
runTime = x.strftime("%y%m%d-%H%M%S")

# define folder path to save outputs
outputFolder = os.path.abspath(os.path.join(os.path.dirname(__file__), 'Output and Figures' "\\"))
folder_name = f'{runTime}_Fig1'
folder_name = os.path.abspath(os.path.join(os.path.dirname(__file__), 'Outputs' "\\" f'{runTime}_Fig1'))
os.makedirs(folder_name, exist_ok=True)

# create dataframe
df = create_dataframe()
dataframe_pkl = create_Attention_Task_pickle(df, folder_name, "rawAttnTaskDf.pkl")
# dataframe_pkl = "231201-094026AttnTaskDF.pkl" # Hard code a previously made pickle!
df = pd.read_pickle(dataframe_pkl)  # df = data_frame - example dataframe_pkl: "231201-094026AttnTaskDF.pkl"

accuracyCutoff = 0.8  # the minimum attention task accuracy percentage needed to be included in dataframe and
# subsequent plots
cutDf = cutLowAccSubs(df, accuracyCutoff)
cutDfMeans = add_means_to_df(cutDf)

# Save copy of dataframe as an Excel spreadsheet for easy viewing
df_filename = runTime + '_AttentionTaskData.xlsx'
plot_path = os.path.join(folder_name, df_filename)
cutDfMeans.to_excel(plot_path)
create_Attention_Task_pickle(df, folder_name, "cutAttnTaskDfMeans.pkl")  # saves cut df with means to a pickle
stats_test_Fig1(cutDfMeans, folder_name)

# %% Plots to plot!
# Panel E - Accuracy, Reaction Time and RTV Violin Plots
acc_violin_plot_Fig1(cutDfMeans, folder_name)
RTV_violin_plot_Fig1(cutDfMeans, folder_name)

# Panel F - Accuracy vs Reaction Time Scatterplot
RTvsAcc_scatter_Fig1(cutDfMeans, folder_name)
# RTVvACC(cutDfMeans, folder_name)

# Panel F (right) - Heatmaps
RTvsAcc_heatmap_Fig1(cutDfMeans, folder_name)

### Complete version of plots (Non-fig1 versions, with p-values)
# acc_violin_plot(cutDf, folder_name)
# RTV_violinplot(cutDf, folder_name)
RTvsAcc_scatter(cutDf, folder_name)
RTvsAcc_heatmap(cutDf, folder_name)

print("The 'plot_fig1_graphs.py' script has completed running. " + runTime)
