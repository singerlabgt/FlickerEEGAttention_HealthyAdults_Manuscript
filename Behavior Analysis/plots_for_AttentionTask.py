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
from StatsTest_Attention import *

x = datetime.datetime.now()
runTime = x.strftime("%y%m%d-%H%M%S")

outputFolder = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'Output and Figures' "\\"))


def AccVPlot(subjects_dataframe, folder_name):  # plots violins for Grouped Accuracy and RT and not grouped Accuracy
    acc_frame=pd.DataFrame(zip(subjects_dataframe["Stimulation"], (subjects_dataframe["AccuraciesAVG"]*100)),
                                          columns=["Stimulation","Accuracy"])
    lowest_acc = min(acc_frame.Accuracy)
    y_minacc = lowest_acc - 0.025
    colors = ["goldenrod", "#ff2c2c", "#0583D2"]

    ranksums_fdr_dataframe=RankSums_FDR_dataframe(subjects_dataframe)
    GammavLightAcc_CorrectedRankSum=round(ranksums_fdr_dataframe["40HzvsLight"]["Accuracy"],3)
    GammavRandomAcc_CorrectedRankSum=round(ranksums_fdr_dataframe["40HzvsRandom"]["Accuracy"], 3)
    fig1, (ax1) = plt.subplots(nrows=1, ncols=1, figsize=(8, 7))
    ax1.set_ylim(ymin=float(y_minacc), ymax=100)
    sns.violinplot(y="Accuracy", x="Stimulation", ax=ax1, data=acc_frame, saturation=0.7, linewidth=0.6,
                   palette=colors,
                   order=["Light", "Random", "40Hz"], inner=None)
    sns.swarmplot(y="Accuracy", x="Stimulation", data=acc_frame, color="black", edgecolor="black", ax=ax1, size=3)
    sns.boxplot(data=acc_frame, x="Stimulation", y="Accuracy", orient="x",  palette=colors, fliersize=3,
                order=["Light", "Random", "40Hz"], width=0.25)
    title1 = plt.suptitle(
        "\n".join(wrap("Accuracies for EEG Attention Task", 60)),
        y=0.97)
    # plt.subplots_adjust(left=0.09, right=0.91, wspace=0.25)
    group_counts = subjects_dataframe["Stimulation"].value_counts()
    ax1.set_ylabel("Accuracy (%)")
    ax1.set_xticklabels(
        [f"{stimulation.get_text()} (n={group_counts[stimulation.get_text()]})" for stimulation in
         ax1.get_xticklabels()])
    box = ax1.get_position()
    ax1.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
    caption = (f'FDR Corrected Rank Sums p-value: \n Gamma v. Light:{GammavLightAcc_CorrectedRankSum}, Gamma v. Random:{GammavRandomAcc_CorrectedRankSum}')
    ax1.annotate(caption, xy=(0.5, -0.15), xycoords='axes fraction', ha='center', fontsize=10)
    ax1.legend()
    # handles, labels = plt.get_legend_handles_labels()
    # hadnles.append(mpathces.Patch(color='none', label=caption))
    # ax1.legend(handles=handles)
    # ax1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
    #           fancybox=True, shadow=True, handles=0, labels=caption, ncol=5)
    # Save the figure
    os.makedirs(folder_name, exist_ok=True)
    plot_filename = runTime + "-AccVPlot"
    plot_path = os.path.join(folder_name, plot_filename)
    plt.savefig(plot_path)
    plt.savefig(plot_path + '.svg')
    plt.close()
    # plt.show()
    print("RT and Accuracy violin plots plotted")

def RtVPlot(subjects_dataframe, folder_name):  # Violin Plot showing RT by Stimulation group
    rt_frame = pd.DataFrame(zip(subjects_dataframe["Stimulation"], subjects_dataframe["Avg_RT"]),
                            columns=["Stimulation", "Reaction Time"])
    ranksums_fdr_dataframe = RankSums_FDR_dataframe(subjects_dataframe)
    GammavLightRT_CorrectedRankSum = round(ranksums_fdr_dataframe["40HzvsLight"]["RT"], 3)
    GammavRandomRT_CorrectedRankSum = round(ranksums_fdr_dataframe["40HzvsRandom"]["RT"], 3)
    fig, (ax1) = plt.subplots(nrows=1, ncols=1, figsize=(8, 7))
    colors = ["goldenrod", "#ff2c2c", "#0583D2"]
    sns.violinplot(y="Reaction Time", x="Stimulation", ax=ax1, data=rt_frame, saturation=0.7, linewidth=0.6,
                   palette=colors, order=["Light", "Random", "40Hz"], inner=None)
    sns.swarmplot(y="Reaction Time", x="Stimulation", data=rt_frame, color="black", edgecolor="black", ax=ax1, size=4)
    sns.boxplot(data=rt_frame, x="Stimulation", y="Reaction Time", orient="x", palette=colors, fliersize=3,
                order=["Light", "Random", "40Hz"], width=0.25)
    title1 = plt.suptitle(
        "\n".join(wrap("RT for EEG Attention Task", 60)),
        y=0.97)
    ax1.set_ylim(0.2, 0.7)
    # Add number of occurrences to x-axis labels
    order = ["Light", "Random", "40Hz"]
    group_counts = subjects_dataframe["Stimulation"].value_counts()
    ax1.set_ylabel("Reaction Times (s)")

    ax1.set_xticklabels(
        [f"{stimulation.get_text()} \n (n={group_counts[stimulation.get_text()]})" for stimulation in
         ax1.get_xticklabels()])
    box = ax1.get_position()
    ax1.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
    caption = (
        f'Rank Sums p-value: \n Gamma v. Light:{GammavLightRT_CorrectedRankSum}, Gamma v. Random:{GammavRandomRT_CorrectedRankSum}')
    ax1.annotate(caption, xy=(0.5, -0.18), xycoords='axes fraction', ha='center', fontsize=10)
    os.makedirs(folder_name, exist_ok=True)
    plot_filename = runTime + "-RtVPlot"
    plot_path = os.path.join(folder_name, plot_filename)
    plt.savefig(plot_path)
    plt.savefig(plot_path + '.svg')
    plt.close()
    # plt.show()


def RtvVPlot(subjects_dataframe, folder_name):  # Violin Plot showing RTV by Stimulation group
    rtv_frame = pd.DataFrame(zip(subjects_dataframe["Stimulation"], subjects_dataframe["RTV"]),columns=["Stimulation", "RT_CV"])
    ranksums_fdr_dataframe = RankSums_FDR_dataframe(subjects_dataframe)
    GammavLightRTV_CorrectedRankSum = round(ranksums_fdr_dataframe["40HzvsLight"]["RT"], 3)
    GammavRandomRTV_CorrectedRankSum = round(ranksums_fdr_dataframe["40HzvsRandom"]["RT"], 3)
    fig, (ax1) = plt.subplots(nrows=1, ncols=1, figsize=(9, 7))
    colors = ["goldenrod", "#ff2c2c", "#0583D2"]
    sns.violinplot(y="RT_CV", x="Stimulation", ax=ax1, data=rtv_frame, saturation=0.7, linewidth=0.6,
                   palette=colors,
                   order=["Light", "Random", "40Hz"], inner=None)
    sns.swarmplot(y="RT_CV", x="Stimulation", data=rtv_frame, color="black", edgecolor="black", ax=ax1, size=4)
    sns.boxplot(data=rtv_frame, x="Stimulation", y="RT_CV", orient="x",  palette=colors, fliersize=3,
                order=["Light", "Random", "40Hz"], width=0.25)

    group_counts = rtv_frame["Stimulation"].value_counts()

    ax1.set_ylabel("Coefficient of Reaction Time Variability")
    ax1.set_xticklabels(
        [f"{stimulation.get_text()} (n={group_counts[stimulation.get_text()]})" for stimulation in
         ax1.get_xticklabels()])
    box = ax1.get_position()
    ax1.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
    caption = (f'Rank Sums p-value: \n Gamma v. Light:{GammavLightRTV_CorrectedRankSum}, Gamma v. Random:{GammavRandomRTV_CorrectedRankSum}')
    ax1.annotate(caption, xy=(0.5, -0.15), xycoords='axes fraction', ha='center', fontsize=10)
    title1 = plt.suptitle("\n".join(
        wrap("RT Variability for EEG Attention Task", 60)), y=0.97)
    ax1.set_ylim(0.0, 0.4)
    # plt.subplots_adjust(left=0.09, right=0.91, wspace=2)
    os.makedirs(folder_name, exist_ok=True)
    plot_filename = runTime + "RtvVPlot"
    plot_path = os.path.join(folder_name, plot_filename)
    plt.savefig(plot_path)
    plt.savefig(plot_path + '.svg')
    # plt.tight_layout()
    plt.close()
    # plt.show()

def RtAccScatter(subjects_dataframe, folder_name):  # Plots RT versus Acc in a Scatterplot
    subjects_dataframe["Stimulation"] = subjects_dataframe["Stimulation"].str.strip()
    grouped = subjects_dataframe.groupby("Stimulation")
    gamma_group = grouped.get_group("40Hz")
    random_group = grouped.get_group("Random")
    light_group = grouped.get_group("Light")

    light_acc_percents= light_group["AccuraciesAVG"]*100
    light_rt_means= light_group["Avg_RT"]
    random_acc_percents=random_group["AccuraciesAVG"]*100
    random_rt_means= random_group["Avg_RT"]
    gamma_acc_percents=gamma_group["AccuraciesAVG"]*100
    gamma_rt_means= gamma_group["Avg_RT"]
    # counts the occurrences of each stimulation for each group
    stimulation_values = subjects_dataframe["Stimulation"].unique()
    stimulation_indices = np.arange(len(stimulation_values))
    stimulation_counts = subjects_dataframe["Stimulation"].value_counts()

    colors = ["goldenrod", "#ff2c2c", "#0583D2"]

    plt.scatter(light_acc_percents, light_rt_means, color=colors[0])
    plt.scatter(random_acc_percents, random_rt_means, color=colors[1])
    plt.scatter(gamma_acc_percents, gamma_rt_means, color=colors[2])

    # Creating custom legend labels with counts

    random_legend_label = f"Random (n={stimulation_counts.get('Random', 0)})"
    gamma_legend_label = f"40Hz (n={stimulation_counts.get('40Hz', 0)})"
    light_legend_label = f"Light (n={stimulation_counts.get('Light', 0)})"

    # Adding legends for the scatter points with custom labels
    plt.legend([light_legend_label, random_legend_label, gamma_legend_label], loc='lower left')

    ar, br = np.polyfit(random_acc_percents, random_rt_means, 1)
    plt.plot(random_acc_percents, ar * random_acc_percents + br, color=colors[1])
    ag, bg = np.polyfit(gamma_acc_percents, gamma_rt_means, 1)
    plt.plot(gamma_acc_percents, ag * gamma_acc_percents + bg, color=colors[2])
    ar, br = np.polyfit(light_acc_percents, light_rt_means, 1)
    plt.plot(light_acc_percents, ar * light_acc_percents + br, color=colors[0])

    spearman_df = Spearman_dataframe(subjects_dataframe)
    Spearman_Light_AccvRT = round(spearman_df["Accuracy v RT"]["Light"], 3)
    Spearman_Gamma_AccvRT = round(spearman_df["Accuracy v RT"]["Gamma"], 3)
    Spearman_Random_AccvRT = round(spearman_df["Accuracy v RT"]["Random"], 3)
    plt.xlabel("Accuracy (%)")
    plt.ylabel("Reaction Time (s)")
    plt.title('Accuracy vs RT for EEG Attention Task')
    caption = (f'Light: r= {Spearman_Light_AccvRT} r\u00b2={round(Spearman_Light_AccvRT*Spearman_Light_AccvRT,3)} '
               f' Random: r= {Spearman_Random_AccvRT} r\u00b2={Spearman_Random_AccvRT*Spearman_Random_AccvRT}'
               f' Gamma: r= {Spearman_Gamma_AccvRT} r\u00b2={round(Spearman_Gamma_AccvRT*Spearman_Gamma_AccvRT,3)}')
    plt.annotate(caption, xy=(0.5, -0.11), xycoords='axes fraction', ha='center', fontsize=10)
    fig = plt.gcf()
    fig.set_size_inches(8, 7)
    # plt.show()
    os.makedirs(folder_name, exist_ok=True)
    plot_filename = runTime + "-RTAccScatter"
    plot_path = os.path.join(folder_name, plot_filename)
    plt.savefig(plot_path)
    plt.savefig(plot_path + '.svg')
    plt.close()


def RtvAccScatter(subjects_dataframe, folder_name):  # scatterplot of Reaction Time Variability versus Accuracy
    subjects_dataframe["Stimulation"] = subjects_dataframe["Stimulation"].str.strip()  # makes sure that spaces are removed from the end of stimulation time to ensure leading or lagging spaces are removed and not excluded from the data
    grouped = subjects_dataframe.groupby("Stimulation")  # grouping participants by their stimulation type
    gamma_group = grouped.get_group("40Hz")
    random_group = grouped.get_group("Random")
    light_group = grouped.get_group("Light")

    # counts the occurences of each stimuatlion for each group
    stimulation_values = subjects_dataframe["Stimulation"].unique()
    stimulation_indices = np.arange(len(stimulation_values))
    stimulation_counts = subjects_dataframe["Stimulation"].value_counts()

    light_acc_percents = light_group["AccuraciesAVG"] * 100
    light_rt_CVs = light_group["RTV"]
    random_acc_percents = random_group["AccuraciesAVG"] * 100
    random_rt_CVs = random_group["RTV"]
    gamma_acc_percents = gamma_group["AccuraciesAVG"] * 100
    gamma_rt_CVs = gamma_group["RTV"]

    # colors that will be used for the scatterplot, also sets the colors by the order of the legend, so make sure you order in the order of you legend or order you legend in the order of this
    colors = ["goldenrod", "#ff2c2c", "#0583D2"]

    plt.scatter(light_acc_percents, light_rt_CVs, color=colors[0])  # order here doesn't matter much either, just make sure color=colors[] correlates to the color you want
    plt.scatter(random_acc_percents, random_rt_CVs, color=colors[1])
    plt.scatter(gamma_acc_percents, gamma_rt_CVs, color=colors[2])

    # Creating custom legend labels with counts, order here doesn't matter too much
    random_legend_label = f"Random (n={stimulation_counts.get('Random', 0)})"
    gamma_legend_label = f"40Hz (n={stimulation_counts.get('40Hz', 0)})"
    light_legend_label = f"Light (n={stimulation_counts.get('Light', 0)})"
    # newgamma_legend_label = f"New 40Hz (n={stimulation_counts.get('New 40Hz', 0)})"
    # combinedgamma_legend_label = f"40Hz (n={stimulation_counts.get( '40Hz', 0)+stimulation_counts.get('New 40Hz',0)})"

    # Adding legends for the scatter points with custom labels, Order Matters Here: if label is not correlating to color graphed, change the order here to adapt to the order of colors
    plt.legend([light_legend_label, random_legend_label, gamma_legend_label], loc='lower left')

    ar, br = np.polyfit(random_acc_percents, random_rt_CVs, 1)
    plt.plot(random_acc_percents, ar * random_acc_percents + br,
             color=colors[1])  # make sure colors coordinate here too
    ag, bg = np.polyfit(gamma_acc_percents, gamma_rt_CVs, 1)
    plt.plot(gamma_acc_percents, ag * gamma_acc_percents + bg, color=colors[2])
    ar, br = np.polyfit(light_acc_percents, light_rt_CVs, 1)
    plt.plot(light_acc_percents, ar * light_acc_percents + br, color=colors[0])

    spearman_df= Spearman_dataframe(subjects_dataframe)
    Spearman_Light_AccvRTV = round(spearman_df["Accuracy v RTV"]["Light"], 3)
    Spearman_Gamma_AccvRTV = round(spearman_df["Accuracy v RTV"]["Gamma"], 3)
    Spearman_Random_AccvRTV= round(spearman_df["Accuracy v RTV"]["Random"], 3)

    plt.xlabel("Accuracy (%)")
    plt.ylabel("Coefficient of Variability")
    plt.title('Accuracy vs RT Variability for EEG Attention Task')
    caption = (f'Light: r= {Spearman_Light_AccvRTV} r\u00b2={Spearman_Light_AccvRTV*Spearman_Light_AccvRTV} '
               f' Random: r= {Spearman_Random_AccvRTV} r\u00b2={round(Spearman_Random_AccvRTV*Spearman_Random_AccvRTV, 3)}'
               f' Gamma: r= {Spearman_Gamma_AccvRTV} r\u00b2={Spearman_Gamma_AccvRTV*Spearman_Gamma_AccvRTV}')
    plt.annotate(caption, xy=(0.5, -0.11), xycoords='axes fraction', ha='center', fontsize=10)
    fig = plt.gcf()
    fig.set_size_inches(8, 7)
    # plt.show()
    os.makedirs(folder_name, exist_ok=True)
    plot_filename = runTime + "-RTVAccScatter"
    plot_path = os.path.join(folder_name, plot_filename)
    plt.savefig(plot_path)
    plt.savefig(plot_path + '.svg')
    plt.close()


def RTvsAcc_heatmap(subjects_dataframe, folder_name):  # Plots RT versus Acc in a Scatterplot
    subjects_dataframe["Stimulation"] = subjects_dataframe["Stimulation"].str.strip()
    grouped = subjects_dataframe.groupby("Stimulation")
    gamma_group = grouped.get_group("40Hz")
    random_group = grouped.get_group("Random")
    light_group = grouped.get_group("Light")

    # counts the occurrences of each stimulation for each group
    stimulation_values = subjects_dataframe["Stimulation"].unique()
    stimulation_indices = np.arange(len(stimulation_values))
    stimulation_counts = subjects_dataframe["Stimulation"].value_counts()

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

    # Plot the first heatmap
    plt.subplot(131)
    plt.imshow(lightHeatmap.T, cmap=colorMap, vmin=0, vmax=6, origin='lower', extent=extent)
    plt.title("Light")
    plt.ylabel("Reaction Time (s)")

    # Plot the second heatmap
    plt.subplot(132)
    plt.imshow(randomHeatmap.T, cmap=colorMap, vmin=0, vmax=6, origin='lower', extent=extent)
    plt.title("Random")

    # Plot the third heatmap
    plt.subplot(133)
    im = plt.imshow(gammaHeatmap.T, cmap=colorMap, vmin=0, vmax=6, origin='lower', extent=extent)
    plt.title("40Hz")

    # Adjust layout to prevent overlapping
    # plt.tight_layout()
    fig.colorbar(im, ax=axs, orientation="horizontal", pad=0, shrink=0.6, location='top')
    fig.set_size_inches(7.2, 4)  # set size of figure

    # Set fonts to something that may be editable in illustrator
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42

    os.makedirs(folder_name, exist_ok=True)
    # Save as a png (read to open and view on PC)
    plot_filename = runTime + "-ACCvsRTHeatmap"
    plot_path = os.path.join(folder_name, plot_filename)
    plt.savefig(plot_path)
    plt.savefig(plot_path + '.svg')  # Save as svg - vector graphic easily to edit in Adobe Illustrator

    # Save as svg - vector graphic easily to edit in Adobe Illustrator
    plot_filename_svg = runTime + "-ACCvsRTHeatmap.svg"
    plot_path_svg = os.path.join(folder_name, plot_filename_svg)
    plt.savefig(plot_path_svg)
    plt.close()


def plotHeatmapPlotly(rxnTimes, Accuracies):
    fig = px.density_heatmap(df, x="total_bill", y="tip", nbinsx=20, nbinsy=20, color_continuous_scale="Viridis")
    fig.show()


def plot2DHeatmapMeshGrid(rxnTimes, Accuracies):
    # generate 2 2d grids for the x & y bounds
    y, x = np.meshgrid(np.linspace(-3, 3, 100), np.linspace(-3, 3, 100))

    z = (1 - x / 2. + x ** 5 + y ** 3) * np.exp(-x ** 2 - y ** 2)
    # x and y are bounds, so z should be the value *inside* those bounds.
    # Therefore, remove the last value from the z array.
    z = z[:-1, :-1]
    z_min, z_max = -np.abs(z).max(), np.abs(z).max()

    fig, ax = plt.subplots()

    c = ax.pcolormesh(x, y, z, cmap='RdBu', vmin=z_min, vmax=z_max)
    ax.set_title('pcolormesh')
    # set the limits of the plot to the limits of the data
    ax.axis([x.min(), x.max(), y.min(), y.max()])
    fig.colorbar(c, ax=ax)

    plt.show()


def plot2DHeatmapImshow(AvgRTs, AvgAccs):
    heatmap, xedges, yedges = np.histogram2d(AvgRTs, AvgAccs, bins=2)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

    plt.clf()
    plt.imshow(heatmap.T, extent=extent, origin='lower')
    plt.show()
