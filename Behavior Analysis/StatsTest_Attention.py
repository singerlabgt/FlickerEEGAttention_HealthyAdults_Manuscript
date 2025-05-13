import os
import pandas as pd
import scipy.stats as stats
import datetime
import numpy as np

# from AttentionTaskFiguresSM import *

RunTime = datetime.datetime.now()
RunTime = RunTime.strftime("%y%m%d-%H%M%S")


def RankSums_NotCorrected_dataframe(subjects_dataframe):
    subjects_dataframe["Stimulation"] = subjects_dataframe["Stimulation"].str.strip()
    grouped = subjects_dataframe.groupby("Stimulation")
    gamma_group = grouped.get_group("40Hz")
    random_group = grouped.get_group("Random")
    light_group = grouped.get_group("Light")
    GammaAccMeans=gamma_group["AccuraciesAVG"]
    RandomAccMeans=random_group["AccuraciesAVG"]
    LightAccMeans=light_group["AccuraciesAVG"]
    GammaRTMeans=gamma_group["Avg_RT"]
    RandomRTMeans=random_group["Avg_RT"]
    LightRTMeans=light_group["Avg_RT"]
    GammaRTV=gamma_group["RTV"]
    RandomRTV=random_group["RTV"]
    LightRTV=light_group["RTV"]

    RankSum_GammavRandomAcc = stats.ranksums(GammaAccMeans, RandomAccMeans)
    RankSum_GammavRandomRT = stats.ranksums(GammaRTMeans, RandomRTMeans)
    RankSum_GammavRandomRTV = stats.ranksums(GammaRTV, RandomRTV)
    RankSum_GammavLightAcc=stats.ranksums(GammaAccMeans, LightAccMeans)
    RankSum_GammavLightRT=stats.ranksums(GammaRTMeans, LightRTMeans)
    RankSum_GammavLightRTV=stats.ranksums(GammaRTV, LightRTV)

    acc_pvals = [RankSum_GammavRandomAcc[1], RankSum_GammavLightAcc[1]]
    rt_pvals = [RankSum_GammavRandomRT[1], RankSum_GammavLightRT[1]]
    rtv_pvals = [RankSum_GammavRandomRTV[1], RankSum_GammavLightRTV[1]]

    colHeaders = ['40HzvsRandom', '40HzvsLight']
    rowHeaders = ['Accuracy', 'RT', 'RTV']
    no_fdrcorrection_frame = pd.DataFrame(np.array([acc_pvals, rt_pvals, rtv_pvals]), columns=colHeaders)
    no_fdrcorrection_frame.index = rowHeaders

    return no_fdrcorrection_frame


def RankSums_FDR_dataframe(subjects_dataframe):
    subjects_dataframe["Stimulation"] = subjects_dataframe["Stimulation"].str.strip()
    grouped = subjects_dataframe.groupby("Stimulation")
    gamma_group = grouped.get_group("40Hz")
    random_group = grouped.get_group("Random")
    light_group = grouped.get_group("Light")
    GammaAccMeans=gamma_group["AccuraciesAVG"]
    RandomAccMeans=random_group["AccuraciesAVG"]
    LightAccMeans=light_group["AccuraciesAVG"]
    GammaRTMeans=gamma_group["Avg_RT"]
    RandomRTMeans=random_group["Avg_RT"]
    LightRTMeans=light_group["Avg_RT"]
    GammaRTV=gamma_group["RTV"]
    RandomRTV=random_group["RTV"]
    LightRTV=light_group["RTV"]

    RankSum_GammavRandomAcc = stats.ranksums(GammaAccMeans, RandomAccMeans)
    RankSum_GammavRandomRT = stats.ranksums(GammaRTMeans, RandomRTMeans)
    RankSum_GammavRandomRTV = stats.ranksums(GammaRTV, RandomRTV)
    RankSum_GammavLightAcc=stats.ranksums(GammaAccMeans, LightAccMeans)
    RankSum_GammavLightRT=stats.ranksums(GammaRTMeans, LightRTMeans)
    RankSum_GammavLightRTV=stats.ranksums(GammaRTV, LightRTV)

    #FDR Correction:
    acc_pvals = [RankSum_GammavRandomAcc[1], RankSum_GammavLightAcc[1]]
    acc_fdr_pvals = stats.false_discovery_control(acc_pvals)
    rt_pvals = [RankSum_GammavRandomRT[1], RankSum_GammavLightRT[1]]
    rt_fdr_pvals = stats.false_discovery_control(rt_pvals)
    rtv_pvals = [RankSum_GammavRandomRTV[1], RankSum_GammavLightRTV[1]]
    rtv_fdr_pvals = stats.false_discovery_control(rtv_pvals)

    colHeaders = ['40HzvsRandom', '40HzvsLight']
    rowHeaders = ['Accuracy', 'RT', 'RTV']
    no_fdrcorrection_frame = pd.DataFrame(np.array([acc_pvals, rt_pvals, rtv_pvals]), columns=colHeaders)
    no_fdrcorrection_frame.index = rowHeaders
    fdr_frame = pd.DataFrame(np.array([acc_fdr_pvals, rt_fdr_pvals, rtv_fdr_pvals]), columns=colHeaders)
    fdr_frame.index = rowHeaders

    return fdr_frame


def Normality_dataframe(subjects_dataframe):
    subjects_dataframe["Stimulation"] = subjects_dataframe["Stimulation"].str.strip()
    grouped = subjects_dataframe.groupby("Stimulation")
    gamma_group = grouped.get_group("40Hz")
    random_group = grouped.get_group("Random")
    light_group = grouped.get_group("Light")
    GammaAccMeans = gamma_group["AccuraciesAVG"]
    RandomAccMeans = random_group["AccuraciesAVG"]
    LightAccMeans = light_group["AccuraciesAVG"]
    GammaRTMeans = gamma_group["Avg_RT"]
    RandomRTMeans = random_group["Avg_RT"]
    LightRTMeans = light_group["Avg_RT"]
    GammaRTV = gamma_group["RTV"]
    RandomRTV = random_group["RTV"]
    LightRTV = light_group["RTV"]

    gammaNormality = []
    randomNormality=[]
    lightNormality=[]

    GammaAcc_Normality=stats.normaltest(GammaAccMeans)
    RandomAcc_Normality=stats.normaltest(RandomAccMeans)
    LightAcc_Normality=stats.normaltest(LightAccMeans)
    GammaRT_Normality=stats.normaltest(GammaRTMeans)
    RandomRT_Normality=stats.normaltest(RandomRTMeans)
    LightRT_Normality=stats.normaltest(LightRTMeans)
    GammaRTV_Normality=stats.normaltest(GammaRTV)
    RandomRTV_Normality=stats.normaltest(RandomRTV)
    LightRTV_Normality=stats.normaltest(LightRTV)

    gammaNormality.append(GammaAcc_Normality.pvalue)
    gammaNormality.append(GammaRT_Normality.pvalue)
    gammaNormality.append(GammaRTV_Normality.pvalue)

    randomNormality.append(RandomAcc_Normality.pvalue)
    randomNormality.append(RandomRT_Normality.pvalue)
    randomNormality.append(RandomRTV_Normality.pvalue)

    lightNormality.append(LightAcc_Normality.pvalue)
    lightNormality.append(LightRT_Normality.pvalue)
    lightNormality.append(LightRTV_Normality.pvalue)

    colHeaders=['Accuracy', 'RT', 'RTV']
    rowHeaders=['Gamma', 'Random', 'Light']
    normality_frame=pd.DataFrame(np.array([gammaNormality, randomNormality, lightNormality]), columns=colHeaders)
    normality_frame.index = rowHeaders

    return normality_frame


def Spearman_dataframe(subjects_dataframe):
    subjects_dataframe["Stimulation"] = subjects_dataframe["Stimulation"].str.strip()
    grouped = subjects_dataframe.groupby("Stimulation")
    gamma_group = grouped.get_group("40Hz")
    random_group = grouped.get_group("Random")
    light_group = grouped.get_group("Light")

    GammaAccMeans = gamma_group["AccuraciesAVG"]
    RandomAccMeans = random_group["AccuraciesAVG"]
    LightAccMeans = light_group["AccuraciesAVG"]
    GammaRTMeans = gamma_group["Avg_RT"]
    RandomRTMeans = random_group["Avg_RT"]
    LightRTMeans = light_group["Avg_RT"]
    GammaRTV = gamma_group["RTV"]
    RandomRTV = random_group["RTV"]
    LightRTV = light_group["RTV"]

    Spearman_AccvRtLight = stats.spearmanr(LightRTMeans, LightAccMeans)
    Spearman_AccvRtRandom = stats.spearmanr(RandomRTMeans, RandomAccMeans)
    Spearman_AccvRtGamma = stats.spearmanr(GammaRTMeans, GammaAccMeans)

    Spearman_AccvRtvLight = stats.spearmanr(LightRTV, LightAccMeans)
    Spearman_AccvRtvRandom = stats.spearmanr(RandomRTV, RandomAccMeans)
    Spearman_AccvRtvGamma = stats.spearmanr(GammaRTV, GammaAccMeans)

    gammaSpearman = []
    randomSpearman = []
    lightSpearman = []

    gammaSpearman.append(Spearman_AccvRtGamma.pvalue)
    gammaSpearman.append(Spearman_AccvRtvGamma.pvalue)
    randomSpearman.append(Spearman_AccvRtRandom.pvalue)
    randomSpearman.append(Spearman_AccvRtvRandom.pvalue)
    lightSpearman.append(Spearman_AccvRtLight.pvalue)
    lightSpearman.append(Spearman_AccvRtvLight.pvalue)

    colHeaders = ['Accuracy v RT', 'Accuracy v RTV']
    rowHeaders = ['Gamma', 'Random', 'Light']
    spearman_frame = pd.DataFrame(np.array([gammaSpearman, randomSpearman, lightSpearman]), columns=colHeaders)
    spearman_frame.index = rowHeaders

    return spearman_frame


def AllStatsTests(subjects_dataframe, folder_name):
    os.makedirs(folder_name, exist_ok=True)
    filename = f'StatsTestOutput {RunTime}.txt'
    plot_path = os.path.join(folder_name, filename)
    txtfile = open(plot_path, 'x')

    spearman_df=Spearman_dataframe(subjects_dataframe)
    no_correction_df=RankSums_NotCorrected_dataframe(subjects_dataframe)
    fdr_df=RankSums_FDR_dataframe(subjects_dataframe)
    normality_df=Normality_dataframe(subjects_dataframe)
    print("Rank Sums", file=txtfile)
    print("Uncorrected p-values: ", file=txtfile)
    print(no_correction_df, file=txtfile)
    print("\nFDR Corrected p-values: ", file=txtfile)
    print(fdr_df, file=txtfile)
    print("\nSpearman", file=txtfile)
    print(spearman_df, file=txtfile)
    print("\nNormality", file=txtfile)
    print(normality_df, file=txtfile)
    txtfile.close()