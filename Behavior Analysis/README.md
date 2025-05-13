# SensoryFx_Attention_Task_Analysis
Attention Task Analysis for Young Adult EEG Flicker Study (Matty) 

## Project Background: 

The goal of this project is to understand the influence of gamma stimulation (40Hz) on healthy young adults. There are three experimental groups: random, gamma, and light groups that help in testing this. These groups are adminstered audiovisual stimulation correlating to each group in order to understand any impacts gamma audiovisual stimulation induces. 

## Overview: 
For this specific code, the end goal is to create figures that showcase key information from the Attention Task. This is the period of 1 hour stimulation in which the participant is completing the dot task. These figures share key information about the accuracy, reaction time, and reaction time variability in order to compare the three stimulation groups. 

### Project Files Description: 
**File 1 - dataframe_creation: **

This specific piece of code generates the dataframe that hosues all of the key information used in future code. A key factor of the code that can be changed is the dates for old and new cohorts in order to specify changes in the gamma group.  

**File 2 - plot_creation_for_AttentionTask:**

This file houses a list of functions that are called later by File 7. In this file are: 

acc_violin_plot: 

![231204-131332-EEG Attention Task RT](https://github.com/singerlabgt/SensoryFx_Attention_Task_Analysis/assets/133933250/81d63718-a17a-486e-9b0f-cf2816b36acb)
![231204-131332-EEG Attention Task Acc](https://github.com/singerlabgt/SensoryFx_Attention_Task_Analysis/assets/133933250/638b907d-5177-42a2-b213-d1326a1b33e7)
![231204-131332-EEG Not Grouped Attention Task Acc](https://github.com/singerlabgt/SensoryFx_Attention_Task_Analysis/assets/133933250/8ac88bcc-749e-49cd-b53b-a693aeba3661)

RTVvACC: 

![231204-131332- ACCvsRTVscatter](https://github.com/singerlabgt/SensoryFx_Attention_Task_Analysis/assets/133933250/5730b34d-2d72-46e4-a9c7-6334ce66f63d)

RTV_violinplot:  

![ Grouped Dot Task Subject Performance Distribution RT](https://github.com/singerlabgt/SensoryFx_Attention_Task_Analysis/assets/133933250/2f7f6d59-edb3-4576-a1dd-16d45e426ed9)

RTvsAcc_scatter: 

![231204-131332- ACCvsRTscatter](https://github.com/singerlabgt/SensoryFx_Attention_Task_Analysis/assets/133933250/4b4b4d9a-218d-4a07-a24f-ecc865ce7240)



**File 3 - Acc_RT_Overlay_for_onefolder**


**File 4 - Window_Size_Violin_by_group_for_onefolder**


**File 5 - MaxStreak**

This file houses a list of functions that are called later by File 7. In this file are: 

AllMissedStreakNotGrouped:

![231204-131332- AllMissedStreaksNotGrouped](https://github.com/singerlabgt/SensoryFx_Attention_Task_Analysis/assets/133933250/530c79aa-6212-47f7-ad81-e7b25bc1d94a)

HighestMissedStreakNotGrouped: 

![231204-131332- HighestMissedStreaksNotGrouped](https://github.com/singerlabgt/SensoryFx_Attention_Task_Analysis/assets/133933250/db7f0965-b04d-44b3-a5e1-598832fc89c9)

HMSGroupAndScattervAcc:

![231204-131332- HighestMissedStreakGrouped](https://github.com/singerlabgt/SensoryFx_Attention_Task_Analysis/assets/133933250/cf788f2d-9970-40f4-80dd-eeb6cf6eb6dc)
![231204-131332- ACCvsMissedStreakScatter](https://github.com/singerlabgt/SensoryFx_Attention_Task_Analysis/assets/133933250/873c9d2c-eb3e-4ad6-b61c-f7da40bb9cf2)

also the dataframe for Highed Missed Streak gets made into a csv file in the folder 

AllMissedStreaksGrouped: 

![231204-131332- AllMissedStreaksGrouped](https://github.com/singerlabgt/SensoryFx_Attention_Task_Analysis/assets/133933250/ffc71610-22ac-47c8-b213-7f78dcd5c2a4)

also the dataframe for All Missed Hits gets made into a csv file in the folder 

**File 6 - StatsTest_Attention**

This file runs a complete stat analysis and generates the product in a text file.

**File 7 - AllFigs:** 

This is the file that calls all of the various functions listed above and generates all of the figures in one folder. This is the file you run when actually generating the figures and where you can comment out certain graphs you don't want made. 

## Getting Started: 
### Dependencies: 
packages: matplotlib, numpy, pandas, seaborn, scipy, datetime, textwrap3, itertools, colorbrewer, os
python version: run using 3.10 - can work for other versions 
### Installation:
install all of the files into one folder in Project_FlickerHealthyYoungAdults\EEG Data Analysis\Attention Task and then make clone the repository in a sepearte folder, this is because of the location of the DataSubset, you could also copy over the datasubset from this file and then create this in a new folder where two folders up the data subset is there

When running the code, some key changes that may need to be made include the cohort date is this changes fro dividing old from new, the pkl file used to ensure you have the proper data, the accuracy cut if desired

### Executing Program: 
to run this code you just need to run AllFigs and this file will call the other scripts to create all of the figures seen above 

### Example Figures 

### Authors
Matty. Sindhu, Shangze