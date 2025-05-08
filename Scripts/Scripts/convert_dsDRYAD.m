%% Script for dsDRYAD - Aggregate 3 Trials per Song as Epochs with Events

%% Basic Setup
clear; clc;
n_trials = 3; % Each song appears 3 times
n_songs = 10; % 10 unique songs
n_subs = 20;   % Number of subjects
path_to_ds = 'D:\master\Iniciacion_Investigacion\Datasets\dsDRYAD\diliBach_4dryad_CND\diliBach_4dryad_CND';
file_base_name = 'dataSub';

%% Process each subject
for sub_idx = 1:n_subs
    
    display(['Processing subject ' num2str(sub_idx)]);
    % Convert subject number to string
    subject_id = sprintf('%d', sub_idx);
    subject_str = sprintf('%03d', sub_idx);

    % Load subject data
    mat_file = [file_base_name subject_id '.mat'];
    subject_data = load(fullfile(path_to_ds, mat_file));

    % Create separate folder for each subject
    subject_folder = fullfile(path_to_ds, ['sub-' subject_str]);
    if ~exist(subject_folder, 'dir')
        mkdir(subject_folder);
    end

    % Extract variables
    EEG_Trials = subject_data.eeg.data; % 1x30 cell array
    Fs = subject_data.eeg.fs;           % Sampling frequency
    chanlocs = subject_data.eeg.chanlocs; % Channel locations
    padding = subject_data.eeg.paddingStartSample; % Mentioned to be 512 in the description (1s)
    n_channels = size(chanlocs, 2); % Number of channels
    
    %% Save 3 trials per song as epochs in a single .set file
    for song_idx = 1:n_songs
        
        song_str = sprintf('%02d', song_idx);
        % Create folder for each song
        song_folder = fullfile(subject_folder, ['song-' song_str]);
        if ~exist(song_folder, 'dir')
            mkdir(song_folder);
        end
        
        % Initialize EEG data matrix (Channels x Samples x Trials) and
        % event related structures
        EEG_Data = []; % Structure for the 3D EEG data 
        event_latencies = []; % Latencies of each event
        event_types = {}; % Types of each event
        event_positions = [];  % These are the latencies of each event if it were continuous data
        
        sample_offset = 0; % Counter

        % Extract trials for this song (song appears once in each 10-trial block)
        for trial_rep = 1:n_trials
            trial_index = song_idx + (trial_rep - 1) * 10;
            trial_data = EEG_Trials{trial_index}'; % Transposed because originally it was (samples x channels)
            
            % Store as 3D EEG data (n_channels x n_samples x epoch)
            EEG_Data(:, :, trial_rep) = trial_data;
            
            % Define event latency and position
            event_latencies = [event_latencies, sample_offset + (padding + 1)];  % Set the event latency to 513 (paddingStartSample + 1)
            event_types = [event_types, {sprintf('Song %d/Trial %d', song_idx, trial_rep)}]; % Unique name for each Song/Trial
            
            sample_offset = sample_offset + size(trial_data, 2);  % Update the sample offset for the next trial
        end

        % Create EEGLAB EEG structure in .set format
        EEG = pop_importdata('dataformat', 'array', 'nbchan', n_channels, 'data', EEG_Data, ...
            'chanlocs', chanlocs, 'setname', ['sub-' subject_id '_song-' num2str(song_idx)], ...
            'srate', Fs, 'pnts', size(trial_data, 2), 'xmin', 0);
        
        % Create events and epochs in EEG structure
        EEG.event = [];
        EEG.epoch = [];
    
        % Create events and fill EEG.epoch for each trial
        for i = 1:n_trials
            % Fill event struct
            EEG.event(i).type = event_types{i};
            EEG.event(i).latency = event_latencies(i);
            EEG.event(i).urevent = i;
            EEG.event(i).epoch = i;
            EEG.epoch(i).event = i;

            % Fill epoch struct
            EEG.epoch(i).eventlatency = padding+1; % The latency is the same for each epoch
            EEG.epoch(i).eventposition = event_latencies(i); % The latency is different in continuous data
            EEG.epoch(i).eventtype = event_types{i};
            EEG.epoch(i).eventurevent = i;
        end

        
        %% Check EEG structure and save
        
        EEG = eeg_checkset(EEG);  % Ensure EEG is valid
        output_filename = ['sub-' subject_id '_song-' num2str(song_idx) '.set'];
        pop_saveset(EEG, 'filename', output_filename, 'filepath', song_folder); % Save .set file

    end

end
