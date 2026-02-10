# FlickerEEGAttention_HealthyAdults_Manuscript
## Project Background
- EEG and behavior analysis code for the flicker study with healthy young adults.
## Overview
- Summary Visual (if applicable)
- Brief description of the goals of the project (what this code does)
## Project Files Description
### Analysis Files
#### Behavior Analysis
- __.py
#### EEG Preprocessing
- S0_ConvertBDFtoSet_and_Preprocess_for_Syncing.m
	- Converts .bdf EEG file into a .set file
- S1_Sync_EEGSet_to_AttentionCsv.m
- S2_LoadEEG_CheckSync.m
- S3_CompletePreprocessing.m
#### EEG Analysis & Figures
- S4_ALLSubjects_NoCut_NoNotch_2_100hz.m

### Data
#### raw data
- location: insert OpenNeuro link
	- format: BIDS
#### processed data
- location: 
- format: stored in the Singer Lab data format as .mat
- analysis code and results:
	-  location: this folder, /ad.gatech.edu/bme/labs/singer/YourName/ProjectName
	- sub folders: -doc/ - documentation files results/ - output figures and intermediary data structures requirements.txt - software requirements for the project project-name/ - analysis code related to the project

## Getting Started
### Dependencies
- Windows 11
- Python 3.11
- MATLAB 2023a
- EEGLAB v2025.0
	- EEGLAB extensions
		- AMICA1.7: v1.7
		- Biosig3.8.4
		- EEG-BIDS10.2
		- ERPLAB9.10
		- Fieldtrip-lite20240111
		- ICLabel
		- LIMO4.1.2
		- bva-io1.73
		- clean_rawdata
		- dipfit
		- firfilt
		- zapline-plus1.2.1`
### Installation/Set-up
- How/where to download your program
- Any modifications needed to be made to files/folders
	- e.g. How to set-up folder structure.  Where should sample data be located?
### Executing program
- place channel locations file, "32BioSemiChannelCoordinates.ced", in a folder calls "Functions"
#### Preprocessing steps
1. Open MATLAB and Set current folder to be the folder with S0_ConvertBDFtoSet_and_Preprocess_for_Syncing.m
2. open S0_ConvertBDFtoSet_and_Preprocess_for_Syncing.m
3. Open eeglab
4. Load all raw set files into eeglab
	- ![Load raw set files](images/readme-draft-1.png)
5. Select the first dataset
	- ![Select first dataset](images/readme-draft-2.png)
6. If preprocessing BDF files,
- code blocks for commands
#### (Include example images of plots/outputs)
## Authors
Contributors names and contact info
## related papers
(name et al., Title, Journal, Year)
Attokaren et al., 40 Hz Audiovisual Stimulation Improves Sustained Attention and Related Brain Oscillations
