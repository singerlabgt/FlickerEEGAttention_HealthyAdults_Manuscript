import os
import pandas as pd
import numpy as np
import datetime

"""
Creates a python data frame from .csv data files produced by PsychoPy behavior attention tasks.
Resulting data frame contains subjects' accuracy and reaction times. 
"""

x = datetime.datetime.now()
runTime = x.strftime("%y%m%d-%H%M%S")


def create_dataframe(master_path='S10-87_AttnTaskCsvPaths.csv', NewCohortDate="2023"):
    starting_dataframe = pd.read_csv(master_path)  # reads .csv file and saves to dataframe variable
    output_dataframe = starting_dataframe
    csvList = output_dataframe["File"]
    sessions = output_dataframe["Session"]

    # Initialize arrays
    rt_arrays = []
    accuracy_arrays = []
    crt_arrays = []
    stimulus_counts = []  # New array to store stimulus counts
    avgAccArr = []
    dates_array = []
    OldOrNewEEG_array = []

    folderName = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
    dataFolderName = os.path.join(folderName, 'Data Subset')

    for rowInCsvList in range(len(csvList)):
        folder = sessions[rowInCsvList]
        if "," in csvList[rowInCsvList]:  # if there are multiple data.csv files for this subject
            print(csvList[rowInCsvList])
            data_frames = []
            split_list = csvList[rowInCsvList].split(",")
            for csvFile in range(len(split_list)):
                # fileSessionName = os.path.abspath(
                fileSessionName = os.path.join(dataFolderName, split_list[csvFile].strip())
                # os.path.
                raw_data = pd.read_csv(fileSessionName, index_col=0)
                data_frames.append(raw_data)
                if csvFile == 0:
                    date = data_frames[csvFile].date.iloc[0]
                    dates_array.append(date)
            raw_data = pd.concat(data_frames)  # Concatenate DataFrames

            # Rest of your processing for this case

        else:
            # fileSessionName = os.path.abspath(
            #     os.path.join(os.path.dirname(__file__), 'Data Subset', csvList[rowInCsvList]))
            fileSessionName = os.path.join(dataFolderName, csvList[rowInCsvList])
            raw_data = pd.read_csv(fileSessionName, index_col=0)
            date = raw_data.date.iloc[0]
            dates_array.append(date)

            # rest of your processing for this case
        try:
            rt_arrays.append(raw_data["reactionKey.rt"])
            accuracy_arrays.append(raw_data["reactionKey.corr"])
            crt_arrays.append(raw_data["reactionKey.started"])
            stimulus_counts.append(len(raw_data))  # Add the number of stimuli to the array
        except:
            continue

        if date[0:4] >= NewCohortDate:
            OldOrNewEEG_array.append("New")
        else:
            OldOrNewEEG_array.append("Old")
    for entry2 in range(len(accuracy_arrays)):
        avgAccRaw = accuracy_arrays[entry2]
        avgAcc = avgAccRaw.mean()
        avgAccArr.append(avgAcc)

    output_dataframe.insert(len(output_dataframe.columns), "Date", dates_array)
    output_dataframe.insert(len(output_dataframe.columns), "Cohort", OldOrNewEEG_array)
    output_dataframe.insert(len(output_dataframe.columns), "Reaction_Times", rt_arrays)
    output_dataframe.insert(len(output_dataframe.columns), "Accuracies", accuracy_arrays)
    output_dataframe.insert(len(output_dataframe.columns), "AccuraciesAVG", avgAccArr)
    output_dataframe.insert(len(output_dataframe.columns), "Cumulative_Reaction_Times", crt_arrays)
    output_dataframe.insert(len(output_dataframe.columns), "Stimulus_Count", stimulus_counts)
    # Add the new array to DataFrame

    return output_dataframe


def cutLowAccSubs(datfra, accThreshold=0.8):
    # accThreshold = 0.5  # threshold of accuracy that participants must be above to be included in the plot
    lowAccuracyIndexes = datfra.index[
        datfra["AccuraciesAVG"] < accThreshold].tolist()  # finds indices of rows with low acc
    aboveThresdf = datfra.drop(labels=lowAccuracyIndexes, axis=0)  # removes low acc rows
    return aboveThresdf


def add_means_to_df(data_frame): # took out the folder name here since not sure why we would do pkl here when we have a seperate function for it
    rtMeans = []
    accMeans = []
    rt_stds = []
    rt_CVs = []

    # for subject in range(0, len(data_frame)):
    #     rtMeans += [data_frame["Reaction_Times"][subject].mean()]
    #     accMeans += [data_frame["Accuracies"][subject].mean()]
    for iSubject, Subject in data_frame.iterrows():  # how to find mean of array of series without for loop?
        rtMeans += [Subject["Reaction_Times"].mean()]
        subxMean = Subject["Reaction_Times"].mean()
        accMeans += [Subject["Accuracies"].mean()]
        subxStd = Subject["Reaction_Times"].std()
        rt_stds += [subxStd]
        subx_rt_CV = (subxStd / subxMean)
        rt_CVs.append(subx_rt_CV)

    data_frame["Avg_RT"] = rtMeans
    data_frame["AccuraciesAVG"] = accMeans
    data_frame["RTV"] = rt_CVs

    # grandMean = np.mean(rtMeans)
    # print("The mean is", grandMean)
    # plot_filename = runTime + "AttnTaskDF.pkl"
    # plot_path = os.path.join(folder_name, plot_filename)
    # data_frame.to_pickle(plot_path)
    return data_frame


def find_mean_rt(data_frame):  # don't think we use this
    rtMeans = []
    for subject in range(0, len(data_frame["Reaction_Times"])):
        rtMeans += [data_frame["Reaction_Times"][subject].mean()]

    grandMean = np.mean(rtMeans)
    print("The mean is", grandMean)
    return grandMean


# "Pickles" and saves dataframe as a .pkl file which can be unpickled later
# Pickling is necessary since dataframe is a heterogeneous table that does not export to .csv or excel easily
def create_Attention_Task_pickle(dataframe: object, folder_name, pickle_name="AttnTaskDF.pkl") -> object:  # exm
    pkl_filename = "./" + runTime + "_" + pickle_name
    pickle_path = os.path.join(folder_name, pkl_filename)
    dataframe.to_pickle(pickle_path)
    return pickle_path

# only comment out if you need to generate a new pkl without using the original code - not neccesary for running the Attention Task Script
# folder_name = f'{runTime}_AttentionTaskFigures'

# os.makedirs(folder_name, exist_ok=True)
# create_Attention_Task_pickle(add_means_to_df(create_dataframe()), folder_name)

