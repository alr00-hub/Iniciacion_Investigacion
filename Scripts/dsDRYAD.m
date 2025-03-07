%% Script for dsDRYAD - Separate .set Files per Song

%% Basic Setup
clear; clc;
n_trials = 1;
n_subs = 1;
path_to_ds = 'D:\master\Iniciacion_Investigacion\Datasets\dsDRYAD\diliBach_4dryad_CND\diliBach_4dryad_CND';
file_base_name = 'dataSub';

%% Process each subject
for sub_idx = 1:n_subs
    
    % Converting sub number to string
    subject_id = sprintf('%d', sub_idx);

    % Load subject data
    mat_file = [file_base_name subject_id '.mat'];
    subject_data = load(fullfile(path_to_ds, mat_file));

    % Create separate folder for each subject
    subject_folder = fullfile(path_to_ds, ['sub-' subject_id]);
    if ~exist(subject_folder, 'dir')
        mkdir(subject_folder);
    end

    % Extract variables
    EEG_Trials = subject_data.eeg.data; % Trials EEG struct
    Fs = subject_data.eeg.fs;               % Sampling frequency
    chanlocs = subject_data.eeg.chanlocs; % Channels lcoation structure
    n_channels = size(chanlocs, 2); % Number of channels

    %% Save each trial separately as a .set file
    for trial_idx = 1:n_trials
        
        % Create separate folder for each song
        trial_folder = fullfile(subject_folder, ['trial-' num2str(trial_idx, '%02d')]);
        if ~exist("trial_folder", 'dir')
            mkdir(trial_folder);
        end
        
        % Extract EEG data for the song
        EEG_Trial = EEG_Trials{trial_idx}';
        
        % Create EEG structure
        EEG = pop_importdata('dataformat', 'matlab', 'nbchan', n_channels, 'data', EEG_Trial, ...
            'chanlocs', chanlocs, 'setname', ['sub-' subject_id 'ses-' num2str(trial_idx)], ...
            'srate', Fs, 'pnts', size(EEG_Trial, 2), 'xmin', 0);
        EEG = eeg_checkset(EEG);

        % Save .set file
        output_filename = ['sub-' subject_id '-ses' num2str(trial_idx) '.set'];
        pop_saveset(EEG, 'filename', output_filename, 'filepath', trial_folder);

    end

end