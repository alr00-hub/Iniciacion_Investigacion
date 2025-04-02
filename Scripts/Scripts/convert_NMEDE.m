%% Basic Setup
clear; clc;

n_songs = 2; % 1 song but 2 variants
n_subs = 24; % Number of subjects
n_channels = 129;
path_to_ds = 'D:\master\GlobusPC\NMED-E'; % Directory of the dataset
file_base_name = 'RawEEG_S'; % Base name for the subject files
spf_file = 'GSN129.sfp'; % Chanlocs File

%% Process each subject
for sub_idx = 1:n_subs
    
    display(['Processing subject ' num2str(sub_idx)]);
    % Convert subject number to string
    subject_id = sprintf('%02d', sub_idx);

    % Load subject data
    mat_file = [file_base_name subject_id '.mat'];
    subject_data = load(fullfile(path_to_ds, mat_file));

    % Create separate folder for each subject
    subject_folder = fullfile(path_to_ds, ['sub-' subject_id]);
    if ~exist(subject_folder, 'dir')
        mkdir(subject_folder);
    end

    % Extract variables
    EEG_Data = subject_data.X; % Channels x Samples
    Fs = subject_data.fs;   % Sampling frequency (1000 Hz)
    DIN_1 = subject_data.DIN_1; % Event information
    
    % Parse the DIN triggers and onsets
    [allTriggers, allOnsets] = parseDIN(DIN_1);
    event_latencies = [];
    event_types = {};

    % Loop through songs
    for song_idx = 1:n_songs
        
        % Get start and end of each song (excluding the ratings)
        song_start_idx = find(allTriggers == 21);
        song_start_idx = song_start_idx(song_idx);
        song_end_idx = find(allTriggers == (11 + song_idx));
        

        precise_onset_arr = find(allTriggers == 128 & allOnsets > allOnsets(song_start_idx));
        precise_onset_baseline = allOnsets(precise_onset_arr(1)) - 1000;
        precise_onset_song = allOnsets(precise_onset_arr(2)) - 1000;

        event_latencies = [event_latencies, precise_onset_baseline];
        event_latencies = [event_latencies, precise_onset_song];
        event_latencies = [event_latencies, allOnsets(song_end_idx)];
        event_types = [event_types, {sprintf("Baseline Trial %d", song_idx)}];
        event_types = [event_types, {sprintf("Trial %d", song_idx)}];
        event_types = [event_types, {sprintf("Ratings Trial %d", song_idx)}];

    end

    chanlocs = readlocs(fullfile(path_to_ds, spf_file));
    EEG = pop_importdata('dataformat', 'array', 'nbchan', n_channels, 'data', EEG_Data, ...
            'chanlocs', chanlocs, 'setname', ['sub-' subject_id], ...
            'srate', Fs, 'pnts', size(EEG_Data, 2), 'xmin', 0);

    EEG.event = [];
    
    for i = 1:length(event_latencies)
        EEG.event(i).type = event_types{i};
        EEG.event(i).latency = event_latencies(i);
        EEG.event(i).epoch = 1;
        EEG.epoch(i).event = i;
    end
  
    EEG = eeg_checkset(EEG);  % Ensure EEG is valid
    output_filename = ['sub-' subject_id '.set'];
    pop_saveset(EEG, 'filename', output_filename, 'filepath', subject_folder); % Save .set file
end
