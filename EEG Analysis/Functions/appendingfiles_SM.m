%edit this line for how many files you are merging ([1 2 ...])
EEG = pop_mergeset( ALLEEG, [1  2 ], 0);
%edit the name you want the file to save as: this will save over other files if you do not change the name, set to save only a .set file right now but if you wanted a .fdt and .set remove the 'save mode' 'one file' part
EEG = pop_saveset( EEG, 'filename','S050_F23.bdf','filepath','C:\Users\smettupalli8\Box\Project_FlickerHealthyYoungAdults\EEG Data Analysis\EEGDataAnalysis2024\Sindhu_EEGSynchingRepository\SensoryFx_EEG_Analysis\1_AllDatasets','savemode', 'onefile');