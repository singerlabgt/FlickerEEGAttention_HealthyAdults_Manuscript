"""
This module contains functions that take in an Attention Task behavior dataframe and plots figures such as scatter, violin and heatmap plots for metrics such as accuracy (Acc), reaction time (RT) and reaction time variability (RTV).
"""

import matplotlib.pyplot as plt

plt.rcParams['svg.fonttype'] = 'none'  # allows the font to be readable when saving plot figures as .svg for illustrator
from matplotlib import pyplot as plt
import seaborn as sns
# from dataframe_creation import create_dataframe
# from dataframe_creation import get_date
from textwrap3 import wrap
import scipy.stats as stats
import datetime
from dataframe_creation import *

x = datetime.datetime.now()
runTime = x.strftime("%y%m%d-%H%M%S")

outputFolder = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'Output and Figures' "\\"))

font = 'Arial'
axis_font_size = 8
legend_font_size = 6
violin_fig_size = (1.5, 1.5)
scatter_fig_size = (2, 2)
heatmap_fig_size = (2, 6)


def acc_violin_plot_Fig1(data_frame, folder_name):  # plots violins for Grouped Accuracy and RT and not grouped Accuracy
    # session = input("For which session type would you like to make figures? Enter EEG or MRI: ")
    session = "EEG"
    violin_df = data_frame
    violin_df["Stimulation"] = violin_df["Stimulation"].str.strip()
    violin_df = violin_df.groupby("Session")
    if session == "EEG":
        violin_df = violin_df.get_group("EEG")
    elif session == "MRI":
        violin_df = violin_df.get_group("MRI")
    else:
        print("Invalid session type. Try again by rerunning the script.")

    grouped = violin_df.groupby("Stimulation")
    gamma_group = grouped.get_group("40Hz")
    random_group = grouped.get_group("Random")
    # newgamma_group = grouped.get_group("New 40Hz")
    light_group = grouped.get_group("Light")
    # combinedgamma_group = pd.concat([newgamma_group, gamma_group])
    rt_dict = {"Stimulation": [], "Reaction Time": []}
    acc_dict = {"Stimulation": [], "Accuracy": []}

    stimulation_values = violin_df["Stimulation"].unique()
    stimulation_indices = np.arange(len(stimulation_values))
    stimulation_counts = violin_df["Stimulation"].value_counts()

    for subject in random_group["Reaction_Times"]:
        rt_dict["Stimulation"] += ["Random"]
        rt_dict["Reaction Time"] += [subject.mean()]
    for subject in gamma_group["Reaction_Times"].dropna():
        rt_dict["Stimulation"] += ["40Hz"]
        rt_dict["Reaction Time"] += [subject.mean()]
    for subject in light_group["Reaction_Times"].dropna():
        rt_dict["Stimulation"] += ["Light"]
        rt_dict["Reaction Time"] += [subject.mean()]
    # for subject in newgamma_group["Reaction_Times"]:
    #    rt_dict["Stimulation"] += ["New 40Hz"]
    #    rt_dict["Reaction Time"] += [subject.mean()]
    # for subject in combinedgamma_group["Reaction_Times"]:
    #   rt_dict["Stimulation"] += ["Random"]
    #   rt_dict["Reaction Time"] += [subject.mean()]
    for subject in random_group["Accuracies"].dropna():
        acc_dict["Stimulation"] += ["Random"]
        acc_dict["Accuracy"] += [subject.mean() * 100]
    for subject in gamma_group["Accuracies"].dropna():
        acc_dict["Stimulation"] += ["40Hz"]
        acc_dict["Accuracy"] += [subject.mean() * 100]
    for subject in light_group["Accuracies"].dropna():
        acc_dict["Stimulation"] += ["Light"]
        acc_dict["Accuracy"] += [subject.mean() * 100]
    # for subject in newgamma_group["Accuracies"]:
    #    acc_dict["Stimulation"] += ["New 40Hz"]
    #    acc_dict["Accuracy"] += [subject.mean()]
    # for subject in combinedgamma_group["Accuracies"]:
    #  acc_dict["Stimulation"] += ["New 40Hz"]
    #  acc_dict["Accuracy"] += [subject.mean()]
    # gamma_acc_percents = np.array(gamma_acc_means) * 100
    # random_acc_percents = np.array(random_acc_means) * 100
    # # newgamma_acc_percents = np.array(newgamma_acc_means) * 100
    # light_acc_percents = np.array(light_acc_means) * 100
    rt_frame = pd.DataFrame(rt_dict)
    acc_frame = pd.DataFrame(acc_dict)
    # rt_frame=rt_frame.dropna()
    # acc_frame=acc_frame.dropna()
    lowest_acc = min(acc_frame.Accuracy)
    y_minacc = lowest_acc - 0.025
    gamma_rt_means = []
    random_rt_means = []
    light_rt_means = []
    gamma_acc_means = []
    random_acc_means = []
    light_acc_means = []
    gamma_rt_stds = []
    random_rt_stds = []
    light_rt_stds = []
    gamma_rt_CVs = []
    random_rt_CVs = []
    light_rt_CVs = []
    AllAccuracies_array = []
    Subjects_array = []
    for subject in random_group["Reaction_Times"]:
        random_rt_means += [subject.mean()]
    for subject in gamma_group["Reaction_Times"]:
        gamma_rt_means += [subject.mean()]
    for subject in light_group["Reaction_Times"]:
        light_rt_means += [subject.mean()]
    for subject in random_group["Accuracies"]:
        random_acc_means += [subject.mean()]
    for subject in light_group["Accuracies"]:
        light_acc_means += [subject.mean()]
    for subject in gamma_group["Accuracies"]:
        gamma_acc_means += [subject.mean()]

    for subject in gamma_group["Reaction_Times"]:
        subxMean = subject.mean()
        subxStd = subject.std()
        gamma_rt_stds += [subxStd]
        subx_rt_CV = [subxStd / subxMean]
        gamma_rt_CVs += subx_rt_CV

    for subject in random_group["Reaction_Times"]:
        subxMean = subject.mean()
        subxStd = subject.std()
        random_rt_stds += [subxStd]
        subx_rt_CV = [subxStd / subxMean]
        random_rt_CVs += subx_rt_CV

    for subject in light_group["Reaction_Times"]:
        subxMean = subject.mean()
        subxStd = subject.std()
        light_rt_stds += [subxStd]
        subx_rt_CV = [subxStd / subxMean]
        light_rt_CVs += subx_rt_CV

    for subject in violin_df["Accuracies"]:
        AllAccuracies_array += [subject.mean() * 100]
    for subject in violin_df["Subject"]:
        Subjects_array.append(subject)

    # fig, (ax, ax0) = plt.subplots(nrows=1, ncols=2, figsize=(6, 4))
    # sns.violinplot(y="Reaction Time", x="Stimulation", ax=ax, data=rt_frame, saturation=0.7, linewidth=0.6)
    # sns.stripplot(y="Reaction Time", x="Stimulation", data=rt_frame, color="gray", edgecolor="black", ax=ax)
    # sns.violinplot(y="Accuracy", x="Stimulation", ax=ax0, data=acc_frame, saturation=0.7, linewidth=0.6)
    # sns.stripplot(y="Accuracy", x="Stimulation", data=acc_frame, color="gray", edgecolor="black", ax=ax0)
    # title = plt.suptitle("\n".join(wrap("Reaction Time/Accuracy Subject Distribution Based on Stimulus Type for " + session + " Dot Task", 60)))
    # ax0.set_ylim(0.3, 1)
    # plt.subplots_adjust(left=0.09, right=0.91, wspace=0.25, top=0.85)
    #
    # plt.savefig(plot_path + "\\_" + runTime + "-" + session + " Grouped Dot Task Subject Performance Distribution")
    # plt.show()

    NotGroupedAccuracies = pd.DataFrame(zip(Subjects_array, AllAccuracies_array),
                                        columns=["Subject", "Accuracies"])

    # stats tests
    GvRRT_RS = round(stats.ranksums(gamma_rt_means, random_rt_means).pvalue,
                     3)  # Rank Sums for All Gamma versus Random Accuracy
    GvLRT_RS = round(stats.ranksums(gamma_rt_means, light_rt_means).pvalue,
                     3)  # Rank Sums for All Gamma versus Light Accuracy
    RvLRT_RS = round(stats.ranksums(random_rt_means, light_rt_means).pvalue,
                     3)  # Rank Sums for Random versus Light Accuracy
    GvRAcc_RS = round(stats.ranksums(gamma_acc_means, random_acc_means).pvalue,
                      6)  # Rank Sums for All Gamma versus Random Accuracy
    GvLAcc_RS = round(stats.ranksums(gamma_acc_means, light_acc_means).pvalue,
                      3)  # Rank Sums for All Gamma versus Light Accuracy
    RvLAcc_RS = round(stats.ranksums(random_acc_means, light_acc_means).pvalue,
                      3)  # Rank Sums for Random versus Light Accuracy

    fig, (ax1) = plt.subplots(nrows=1, ncols=1, figsize=violin_fig_size)
    # Set the default font to Arial
    plt.rcParams['font.family'] = font
    # plt.tight_layout()
    colors = ["goldenrod", "#ff2c2c", "#0583D2"]
    sns.violinplot(y="Reaction Time", x="Stimulation", ax=ax1, data=rt_frame, saturation=0.7, linewidth=0.6,
                   palette=colors, order=["Light", "Random", "40Hz"], inner=None)
    sns.swarmplot(y="Reaction Time", x="Stimulation", data=rt_frame, color="black", edgecolor="black", ax=ax1, size=2)
    sns.boxplot(data=rt_frame, x="Stimulation", y="Reaction Time", orient="x", palette=colors, fliersize=3,
                order=["Light", "Random", "40Hz"], width=0.25)
    # title1 = plt.suptitle(
    #     "\n".join(wrap("RT for " + session + " Attention Task", 60)),
    #     y=0.97)
    ax1.set_ylim(0.2, 0.7)
    # Add number of occurrences to x-axis labels
    order = ["Light", "Random", "40Hz"]
    group_counts = violin_df["Stimulation"].value_counts()
    ax1.set_ylabel("Reaction Times (s)", fontsize=axis_font_size)

    ax1.set_xticklabels(
        [f"{stimulation.get_text()}" for stimulation in
         ax1.get_xticklabels()],
        fontsize=axis_font_size)  # Add stim counts:  \n (n={group_counts[stimulation.get_text()]})
    box = ax1.get_position()
    ax1.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
    caption = (
        f'Rank Sums p-value: \n Gamma v. Light:{GvLRT_RS}, Gamma v. Random:{GvRRT_RS}')
    # ax1.annotate(caption, xy=(0.5, -0.15), xycoords='axes fraction', ha='center', fontsize=6)
    ax1.spines[['right', 'top']].set_visible(False)

    # Set fonts to something that may be editable in illustrator
    # Set the default font to Arial
    plt.rcParams['font.family'] = font

    # fig.set_size_inches(3, 2)  # set size of figure
    os.makedirs(folder_name, exist_ok=True)
    plot_filename = runTime + "-AvgRT_for_" + session + "_Attention_Task_Fig1"
    plot_path = os.path.join(folder_name, plot_filename)
    plt.savefig(plot_path)
    plt.savefig(plot_path + '.svg')
    plt.close()
    # plt.show()

    fig2, (ax2) = plt.subplots(nrows=1, ncols=1, figsize=violin_fig_size)
    ax2.set_ylim(ymin=float(y_minacc), ymax=100)
    sns.violinplot(y="Accuracy", x="Stimulation", ax=ax2, data=acc_frame, saturation=0.7, linewidth=0.6, palette=colors,
                   order=["Light", "Random", "40Hz"], inner=None)
    sns.swarmplot(y="Accuracy", x="Stimulation", data=acc_frame, color="black", edgecolor="black", ax=ax2, size=2)
    sns.boxplot(data=acc_frame, x="Stimulation", y="Accuracy", orient="x", palette=colors, fliersize=3,
                order=["Light", "Random", "40Hz"], width=0.25)
    # title2 = plt.suptitle(
    #     "\n".join(wrap("Accuracies for " + session + " Attention Task", 60)),
    #     y=0.97)

    group_counts = violin_df["Stimulation"].value_counts()
    ax2.set_ylabel("Accuracy (%)")
    ax2.set_xticklabels(
        [f"{stimulation.get_text()}" for stimulation in
         ax2.get_xticklabels()],
        fontsize=axis_font_size)  # Add stim counts:  \n (n={group_counts[stimulation.get_text()]})
    box = ax2.get_position()
    ax2.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
    caption = (f'Rank Sums p-value: \n Gamma v. Light:{GvLAcc_RS}, Gamma v. Random:{GvRAcc_RS}')
    # ax2.annotate(caption, xy=(0.5, -0.15), xycoords='axes fraction', ha='center', fontsize=6)
    ax2.spines[['right', 'top']].set_visible(False)

    # Set fonts to something that may be editable in illustrator
    # Set the default font to Arial
    plt.rcParams['font.family'] = font

    # fig2.set_size_inches(3, 2)  # set size of figure
    os.makedirs(folder_name, exist_ok=True)
    plot_filename = runTime + "-" + session + "_AvgAcc_AttnTask_Fig1"
    plot_path = os.path.join(folder_name, plot_filename)
    plt.savefig(plot_path)
    plt.savefig(plot_path + '.svg')
    plt.close()
    # plt.show()

    print("RT and Accuracy violin plots plotted")


def RTV_violin_plot_Fig1(data_frame, folder_name):  # Violin Plot showing RTV by Stimulation group
    data_frame["Stimulation"] = data_frame["Stimulation"].str.strip()
    rt_means = []
    rt_stds = []
    rt_CVs = []
    # acc_means = []

    for subject in data_frame["Reaction_Times"]:
        subxMean = np.nanmean(subject)  # mean for subject
        rt_means.append(subxMean)
        subxStd = np.nanstd(subject)
        rt_stds.append(subxStd)
        subx_rt_CV = subxStd / subxMean if subxMean != 0 else 0  # CV = std/mean (avoid division by zero)
        rt_CVs.append(subx_rt_CV)
    data_frame["RT_std"] = rt_stds
    data_frame["RT_CV"] = rt_CVs

    stimulation_values = data_frame["Stimulation"].unique()
    stimulation_indices = np.arange(len(stimulation_values))

    grouped = data_frame.groupby("Stimulation")
    gamma_group = grouped.get_group("40Hz")
    random_group = grouped.get_group("Random")
    # newgamma_group = grouped.get_group("New 40Hz")
    light_group = grouped.get_group("Light")
    # combinedgamma_group = pd.concat([grouped.get_group("New 40Hz"), grouped.get_group("40Hz")])

    gamma_rt_means = []
    random_rt_means = []
    newgamma_rt_means = []
    light_rt_means = []
    gamma_acc_means = []
    random_acc_means = []
    newgamma_acc_means = []
    light_acc_means = []
    allGammaRtMeans = []
    allGammaAccMeans = []
    gamma_rt_stds = []
    random_rt_stds = []
    gamma_rt_CVs = []
    random_rt_CVs = []
    newgamma_rt_stds = []
    light_rt_stds = []
    newgamma_rt_CVs = []
    light_rt_CVs = []
    allGamma_rt_stds = []
    allGamma_rt_CVs = []

    for subject in random_group["Reaction_Times"]:
        random_rt_means += [subject.mean()]
    for subject in gamma_group["Reaction_Times"]:
        gamma_rt_means += [subject.mean()]
        allGammaRtMeans += [subject.mean()]
    for subject in light_group["Reaction_Times"]:
        light_rt_means += [subject.mean()]
    for subject in random_group["Accuracies"]:
        random_acc_means += [subject.mean()]
    for subject in gamma_group["Accuracies"]:
        gamma_acc_means += [subject.mean()]
        allGammaAccMeans += [subject.mean()]
    for subject in light_group["Accuracies"]:
        light_acc_means += [subject.mean()]

    for subject in gamma_group["Reaction_Times"]:
        subxMean = subject.mean()
        subxStd = subject.std()
        gamma_rt_stds += [subxStd]
        subx_rt_CV = [subxStd / subxMean]
        gamma_rt_CVs += subx_rt_CV

    for subject in random_group["Reaction_Times"]:
        subxMean = subject.mean()
        subxStd = subject.std()
        random_rt_stds += [subxStd]
        subx_rt_CV = [subxStd / subxMean]
        random_rt_CVs += subx_rt_CV

    for subject in light_group["Reaction_Times"]:
        subxMean = subject.mean()
        subxStd = subject.std()
        light_rt_stds += [subxStd]
        subx_rt_CV = [subxStd / subxMean]
        light_rt_CVs += subx_rt_CV
    # GvRRTV_RS = round(stats.ranksums(gamma_rt_CVs, random_rt_CVs).pvalue,
    #                   2)  # Rank Sums for All Gamma versus Random Accuracy
    # GvLRTV_RS = round(stats.ranksums(gamma_rt_CVs, light_rt_CVs).pvalue,
    #                   2)  # Rank Sums for All Gamma versus Light Accuracy
    # RvLRTV_RS = round(stats.ranksums(random_rt_CVs, light_rt_CVs).pvalue,
    #                   2)  # Rank Sums for Random versus Light Accuracy

    subframe = data_frame[["Subject", "Stimulation", "RT_CV"]]
    # data_frame["Reaction_Times"] = pd.to_numeric(data_frame["Reaction_Times"], errors="coerce")

    # Making plots with two violins per plot, one for each stimulus type
    fig, (ax1) = plt.subplots(nrows=1, ncols=1, figsize=violin_fig_size)
    # Set the default font to Arial
    plt.rcParams['font.family'] = font
    colors = ["goldenrod", "#ff2c2c", "#0583D2"]

    sns.violinplot(y="RT_CV", x="Stimulation", ax=ax1, data=subframe, saturation=0.7, linewidth=0.6, palette=colors,
                   order=["Light", "Random", "40Hz"], inner=None)
    sns.swarmplot(y="RT_CV", x="Stimulation", data=subframe, color="black", edgecolor="black", ax=ax1, size=2)
    sns.boxplot(data=subframe, x="Stimulation", y="RT_CV", orient="x", palette=colors, fliersize=3,
                order=["Light", "Random", "40Hz"], width=0.25)
    means = data_frame.groupby("Stimulation")["RT_CV"].mean()
    x_vals = np.arange(len(means))
    group_counts = data_frame["Stimulation"].value_counts()
    # for ipoint in range(len(subframe)):
    #     if ipoint==26 or ipoint==34 or ipoint==36:
    #         continue
    #     else:
    #         ax1.text(x=subframe["Stimulation"][ipoint], y=subframe["RT_CV"][ipoint], s=subframe["Subject"][ipoint], horizontalalignment='center', size='small', color='black')
    ax1.set_ylabel("RT CV", fontsize=axis_font_size)
    ax1.set_xlabel(" ", fontsize=axis_font_size)
    ax1.set_xticklabels(
        [f"{stimulation.get_text()}" for stimulation in
         ax1.get_xticklabels()],
        fontsize=axis_font_size)  # Add stim counts:  \n (n={group_counts[stimulation.get_text()]})
    box = ax1.get_position()
    ax1.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
    ax1.spines[['right', 'top']].set_visible(False)
    # caption = (
    #     f'Rank Sums p-value: \n Gamma v. Light:{GvLRTV_RS}, Gamma v. Random:{GvRRTV_RS}')
    # ax1.annotate(caption, xy=(0.5, -0.45), xycoords='axes fraction', ha='center', fontsize=2)
    # title1 = plt.suptitle("\n".join(
    #     wrap("RTV", 60)), y=0.97)
    ax1.set_ylim(0.0, 0.4)

    # Set fonts to something that may be editable in illustrator
    # Set the default font to Arial
    plt.rcParams['font.family'] = font

    # plt.subplots_adjust(left=0.09, right=0.91, wspace=2)
    os.makedirs(folder_name, exist_ok=True)
    # fig.set_size_inches(3, 2)  # set size of figure
    plot_filename = runTime + "_RTV_Attention_Task_Fig1"
    plot_path = os.path.join(folder_name, plot_filename)
    plt.savefig(plot_path)
    plt.savefig(plot_path + '.svg')
    # plt.tight_layout()
    plt.close()


def RTvsAcc_scatter_Fig1(data_frame, folder_name):  # Plots RT versus Acc in a Scatterplot
    data_frame["Stimulation"] = data_frame["Stimulation"].str.strip()
    grouped = data_frame.groupby("Stimulation")
    gamma_group = grouped.get_group("40Hz")
    random_group = grouped.get_group("Random")
    # newgamma_group = grouped.get_group("New 40Hz")
    light_group = grouped.get_group("Light")
    # combinedgamma_group = pd.concat([newgamma_group, gamma_group])

    # counts the occurrences of each stimulation for each group
    stimulation_values = data_frame["Stimulation"].unique()
    stimulation_indices = np.arange(len(stimulation_values))
    stimulation_counts = data_frame["Stimulation"].value_counts()

    # intiate Variables
    gamma_rt_means = []
    random_rt_means = []
    gamma_rt_stds = []
    random_rt_stds = []
    gamma_rt_CVs = []
    random_rt_CVs = []
    gamma_acc_means = []
    random_acc_means = []
    newgamma_rt_means = []
    light_rt_means = []
    newgamma_rt_stds = []
    light_rt_stds = []
    newgamma_rt_CVs = []
    light_rt_CVs = []
    newgamma_acc_means = []
    light_acc_means = []
    combinedgamma_rt_means = []
    combinedgamma_rt_stds = []
    combinedgamma_rt_CVs = []
    combinedgamma_acc_means = []

    # Corrected mean, std, and CV calculations for gamma group
    for participant_RTs in gamma_group["Reaction_Times"]:
        subxMean = participant_RTs.mean()
        gamma_rt_means += [subxMean]
        subxStd = participant_RTs.std()
        gamma_rt_stds += [subxStd]
        subx_rt_CV = [subxStd / subxMean]
        gamma_rt_CVs += subx_rt_CV

    for participant_hits in gamma_group["Accuracies"]:
        gamma_acc_means += [participant_hits.mean()]

    # Corrected mean, std, and CV calculations for random group
    for random_subject_RTs in random_group["Reaction_Times"]:
        subxMean = random_subject_RTs.mean()
        random_rt_means += [subxMean]
        subxStd = random_subject_RTs.std()
        random_rt_stds += [subxStd]
        subx_rt_CV = [subxStd / subxMean]
        random_rt_CVs += subx_rt_CV

    for random_subject_hits in random_group["Accuracies"]:
        random_acc_means += [random_subject_hits.mean()]

    # Corrected mean, std, and CV calculations for light group
    for light_subject_RTs in light_group["Reaction_Times"]:
        subxMean = light_subject_RTs.mean()
        light_rt_means += [subxMean]
        subxStd = light_subject_RTs.std()
        light_rt_stds += [subxStd]
        subx_rt_CV = [subxStd / subxMean]
        light_rt_CVs += subx_rt_CV

    for x in light_group["Accuracies"]:
        light_acc_means += [x.mean()]

        # Corrected mean, std, and CV calculations for newgamma group
        # for x in newgamma_group["Reaction_Times"]:
        subxMean = x.mean()
        newgamma_rt_means += [subxMean]
        subxStd = x.std()
        newgamma_rt_stds += [subxStd]
        subx_rt_CV = [subxStd / subxMean]
        newgamma_rt_CVs += subx_rt_CV

        # for x in newgamma_group["Accuracies"]:
        newgamma_acc_means += [x.mean()]
        # for x in combinedgamma_group["Reaction_Times"]:
        subxMean = x.mean()
        combinedgamma_rt_means += [subxMean]
        subxStd = x.std()
        combinedgamma_rt_stds += [subxStd]
        subx_rt_CV = [subxStd / subxMean]
        combinedgamma_rt_CVs += subx_rt_CV

        # for x in combinedgamma_group["Accuracies"]:
        combinedgamma_acc_means += [x.mean()]

    gamma_acc_percents = np.array(gamma_acc_means) * 100
    random_acc_percents = np.array(random_acc_means) * 100
    # newgamma_acc_percents = np.array(newgamma_acc_means) * 100
    light_acc_percents = np.array(light_acc_means) * 100
    # combinedgamma_acc_percents = np.array(combinedgamma_acc_means) * 100

    # Corrected mean, std, and CV calculations for newgamma group
    gamma_acc_percents = np.array(gamma_acc_means) * 100
    random_acc_percents = np.array(random_acc_means) * 100
    # newgamma_acc_percents = np.array(newgamma_acc_means) * 100
    light_acc_percents = np.array(light_acc_means) * 100
    # combinedgamma_acc_percents = np.array(combinedgamma_acc_means)*100

    colors = ["goldenrod", "#ff2c2c", "#0583D2"]

    # Set fonts to something that may be editable in illustrator
    # Set the default font to Arial
    plt.rcParams['font.family'] = font

    fig, (ax) = plt.subplots(nrows=1, ncols=1, figsize=(1.5, 1.5))

    plt.scatter(light_acc_percents, light_rt_means, color=colors[0], s=3)
    plt.scatter(random_acc_percents, random_rt_means, color=colors[1], s=3)
    plt.scatter(gamma_acc_percents, gamma_rt_means, color=colors[2], s=3)
    # plt.scatter(newgamma_acc_percents, newgamma_rt_means)
    # plt.scatter(combinedgamma_acc_percents, combinedgamma_rt_means)

    # Creating custom legend labels with counts

    random_legend_label = f"Random (n={stimulation_counts.get('Random', 0)})"
    gamma_legend_label = f"40Hz (n={stimulation_counts.get('40Hz', 0)})"
    light_legend_label = f"Light (n={stimulation_counts.get('Light', 0)})"
    # newgamma_legend_label = f"New 40Hz (n={stimulation_counts.get('New 40Hz', 0)})"
    # combinedgamma_legend_label = f"40Hz (n={stimulation_counts.get( '40Hz', 0)+stimulation_counts.get('New 40Hz',0)})"

    # Adding legends for the scatter points with custom labels
    # plt.legend([random_legend_label, gamma_legend_label, light_legend_label, newgamma_legend_label],loc='lower left')

    plt.legend([light_legend_label, random_legend_label, gamma_legend_label], loc='lower left',
               fontsize=legend_font_size)

    ar, br = np.polyfit(random_acc_percents, random_rt_means, 1)
    plt.plot(random_acc_percents, ar * random_acc_percents + br, color=colors[1])
    ag, bg = np.polyfit(gamma_acc_percents, gamma_rt_means, 1)
    plt.plot(gamma_acc_percents, ag * gamma_acc_percents + bg, color=colors[2])
    ar, br = np.polyfit(light_acc_percents, light_rt_means, 1)
    plt.plot(light_acc_percents, ar * light_acc_percents + br, color=colors[0])

    # Spearman_AccvRtLight = round(stats.spearmanr(light_rt_means, light_acc_means).statistic, 3)
    # SpearmanSquared_AccvRtLight = round(
    #     stats.spearmanr(light_rt_means, light_acc_means).statistic * stats.spearmanr(light_rt_means,
    #                                                                                  light_acc_means).statistic, 3)
    # Spearman_AccvRtRandom = round(stats.spearmanr(random_rt_means, random_acc_means).statistic, 3)
    # SpearmanSquared_AccvRtRandom = round(
    #     stats.spearmanr(random_rt_means, random_acc_means).statistic * stats.spearmanr(random_rt_means,
    #                                                                                    random_acc_means).statistic, 3)
    # Spearman_AccvRtGamma = round(stats.spearmanr(gamma_rt_means, gamma_acc_means).statistic, 3)
    # SpearmanSquared_AccvRtGamma = round(
    #     stats.spearmanr(gamma_rt_means, gamma_acc_means).statistic * stats.spearmanr(gamma_rt_means,
    #                                                                                  gamma_acc_means).statistic, 3)
    plt.xlabel("Accuracy (%)", fontsize=axis_font_size)
    plt.ylabel("Reaction Time (s)", fontsize=axis_font_size)
    # plt.title('Accuracy vs RT for EEG Attention Task')
    # caption = (f'Light: r= {Spearman_AccvRtLight} r\u00b2={SpearmanSquared_AccvRtLight} '
    #            f' Random: r= {Spearman_AccvRtRandom} r\u00b2={SpearmanSquared_AccvRtRandom}'
    #            f' Gamma: r= {Spearman_AccvRtGamma} r\u00b2={SpearmanSquared_AccvRtGamma}')
    # plt.annotate(caption, xy=(0.5, -0.11), xycoords='axes fraction', ha='center', fontsize=10)
    # fig, ax = plt.gcf()

    ax.spines[['right', 'top']].set_visible(False)
    fig.set_size_inches(2.5, 2)
    # plt.show()
    os.makedirs(folder_name, exist_ok=True)
    plot_filename = runTime + "- ACCvsRTscatter_Fig1"
    plot_path = os.path.join(folder_name, plot_filename)
    plt.savefig(plot_path)
    plt.savefig(plot_path + '.svg')
    plt.close()


def RTvsAcc_heatmap_Fig1(data_frame, folder_name):  # Plots RT versus Acc in a Scatterplot
    data_frame["Stimulation"] = data_frame["Stimulation"].str.strip()
    grouped = data_frame.groupby("Stimulation")
    gamma_group = grouped.get_group("40Hz")
    random_group = grouped.get_group("Random")
    light_group = grouped.get_group("Light")

    # counts the occurrences of each stimulation for each group
    stimulation_values = data_frame["Stimulation"].unique()
    stimulation_indices = np.arange(len(stimulation_values))
    stimulation_counts = data_frame["Stimulation"].value_counts()

    # initiate Variables
    gamma_rt_means = []
    random_rt_means = []
    gamma_rt_stds = []
    random_rt_stds = []
    gamma_rt_CVs = []
    random_rt_CVs = []
    gamma_acc_means = []
    random_acc_means = []
    light_rt_means = []
    light_rt_stds = []
    light_rt_CVs = []
    light_acc_means = []

    # Corrected mean, std, and CV calculations for gamma group
    for subject in gamma_group["Reaction_Times"]:
        subxMean = subject.mean()
        gamma_rt_means += [subxMean]
        subxStd = subject.std()
        gamma_rt_stds += [subxStd]
        subx_rt_CV = [subxStd / subxMean]
        gamma_rt_CVs += subx_rt_CV

    for subject in gamma_group["Accuracies"]:
        gamma_acc_means += [subject.mean()]

    # Corrected mean, std, and CV calculations for random group
    for subject in random_group["Reaction_Times"]:
        subxMean = subject.mean()
        random_rt_means += [subxMean]
        subxStd = subject.std()
        random_rt_stds += [subxStd]
        subx_rt_CV = [subxStd / subxMean]
        random_rt_CVs += subx_rt_CV

    for subject in random_group["Accuracies"]:
        random_acc_means += [subject.mean()]

    # Corrected mean, std, and CV calculations for light group
    for subject in light_group["Reaction_Times"]:
        subxMean = subject.mean()
        light_rt_means += [subxMean]
        subxStd = subject.std()
        light_rt_stds += [subxStd]
        subx_rt_CV = [subxStd / subxMean]
        light_rt_CVs += subx_rt_CV

    for subject in light_group["Accuracies"]:
        light_acc_means += [subject.mean()]

    gamma_acc_percents = np.array(gamma_acc_means) * 100
    random_acc_percents = np.array(random_acc_means) * 100
    light_acc_percents = np.array(light_acc_means) * 100

    # AccuracyBins = 10
    # RTBins = 0
    # bins = 10
    accEdges = [0.80, 0.825, 0.85, 0.875, 0.90, 0.925, 0.95, 0.975, 1.00]  # x-axis
    rtEdges = [0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6]  # y-axis
    bins = (accEdges, rtEdges)
    colorMap = 'viridis'  # blueToYellow: cmap='viridis' or 'gray'
    lightHeatmap, xedges, yedges = np.histogram2d(light_acc_means, light_rt_means, bins)
    randomHeatmap, _, _ = np.histogram2d(random_acc_means, random_rt_means, bins)
    gammaHeatmap, _, _ = np.histogram2d(gamma_acc_means, gamma_rt_means, bins)
    extent = [accEdges[0], xedges[-1], yedges[0], yedges[-1]]

    plt.clf()
    plt.close()
    fig, axs = plt.subplots(layout='constrained')
    axs.spines[:].set_visible(False)  # hides all spines
    axs.set_xlabel('Accuracy', labelpad=20)
    # plt.xlabel("Accuracy", labelpad=-10)
    # plt.ylabel("Reaction Time (s)")
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)

    heatmap_xticks = (0.8, 0.9, 1.0)

    # Plot the first heatmap
    plt.subplot(131)
    plt.imshow(lightHeatmap.T, cmap=colorMap, vmin=0, vmax=6, origin='lower', extent=extent)
    plt.title("Light")
    plt.ylabel("Reaction Time (s)")
    plt.xticks(heatmap_xticks)

    # Plot the second heatmap
    plt.subplot(132)
    plt.imshow(randomHeatmap.T, cmap=colorMap, vmin=0, vmax=6, origin='lower', extent=extent)
    plt.title("Random")
    plt.xticks(heatmap_xticks)


    # Plot the third heatmap
    plt.subplot(133)
    im = plt.imshow(gammaHeatmap.T, cmap=colorMap, vmin=0, vmax=6, origin='lower', extent=extent)
    plt.xticks(heatmap_xticks)
    plt.title("40Hz")

    # Adjust layout to prevent overlapping
    # plt.tight_layout()
    fig.colorbar(im, ax=axs, orientation="horizontal", pad=0, shrink=0.6, location='top')
    fig.set_size_inches(3.6, 2)  # set size of figure

    # Set fonts to something that may be editable in illustrator
    # Set the default font to Arial
    plt.rcParams['font.family'] = font

    os.makedirs(folder_name, exist_ok=True)
    # Save as a png (read to open and view on PC)
    plot_filename = runTime + "- ACCvsRTHeatmap_Fig1"
    plot_path = os.path.join(folder_name, plot_filename)
    plt.savefig(plot_path)
    plt.savefig(plot_path + '.svg')  # Save as svg - vector graphic easily to edit in Adobe Illustrator

    plt.close()
    # Plot
