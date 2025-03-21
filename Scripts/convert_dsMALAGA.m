clear; clc; close all;
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;  % Initialize EEGLAB

% Define dataset directory
data_path = 'D:\master\Iniciacion_Investigacion\Datasets\dsMalaga\Raw_Data\Raw_Data\';

% List of subjects
subjects = {'0009'};

% Loop through each subject
for i = 1:length(subjects)
    subjID = subjects{i};
    vhdr_file = ['S_', subjID, '_E_0007_0001.vhdr'];
    
    % Load EEG data
    EEG = pop_loadbv(data_path, vhdr_file);
    
    % Add dataset to EEGLAB
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0, 'setname', ['Subject_' subjID], 'gui', 'off'); 
    
    % Save dataset
    EEG = pop_saveset(EEG, 'filename', ['S_00', subjID, '.set'], 'filepath', data_path);
end

disp('All subjects processed successfully!');
