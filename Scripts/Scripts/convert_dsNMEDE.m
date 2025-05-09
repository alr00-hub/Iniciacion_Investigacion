%% Script for dsNMEDE conversion
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
    EEG_Data = subject_data.X; % Channels x Samples
    Fs = subject_data.fs;   % Sampling frequency (1000 Hz)
    DIN_1 = subject_data.DIN_1; % Event information
    
    % Parse the DIN triggers and onsets
    [allTriggers, allOnsets] = parseDIN(DIN_1);
    event_latencies = [];
    event_types = {};

    fst_song = 0;
    for i = 1:length(allTriggers)
        if allTriggers(i) == 22, fst_song = 1; break; end
        if allTriggers(i) == 23, fst_song = 2; break; end
    end

    % Loop through songs
    for song_idx = 1:n_songs
        
        song_str = sprintf('%02d', song_idx);

        song_folder = fullfile(subject_folder, ['song-' song_str]);
        if ~exist(song_folder, 'dir')
            mkdir(song_folder);
        end
        
        baseline_start_idx = find(allTriggers == 21);
        if fst_song == 1, baseline_start_idx = baseline_start_idx(song_idx); end
        if fst_song == 2, baseline_start_idx = baseline_start_idx(1 + n_songs - song_idx); end

        % Get start and end of each song (excluding the ratings)
        song_start_idx = find(allTriggers == 21 + song_idx);
        %song_start_idx = song_start_idx(song_idx);
        song_end_idx = find(allTriggers == (11 + song_idx));
        
        % Find the onset times for baseline and ratings trial
        precise_onset_arr = find(allTriggers == 128 & allOnsets < allOnsets(song_end_idx) & allOnsets > allOnsets(baseline_start_idx));
        precise_onset_baseline = allOnsets(precise_onset_arr(1)) - 1000; % Baseline trial onset (with offset)
        precise_onset_song = allOnsets(precise_onset_arr(2)) - 1000;     % Song trial onset (with offset)
        
        % Set up the latencies for the events
        event_latencies = [precise_onset_baseline, precise_onset_song];
        event_types = {sprintf("Baseline Trial %d", song_idx), sprintf("Trial %d", song_idx)};

        % Find the onset time for Ratings Trial start (event 11 + song_idx)
        ratings_start_idx = find(allTriggers == (11 + song_idx));
        ratings_start_onset = allOnsets(ratings_start_idx);
        event_latencies = [event_latencies, ratings_start_onset];
        event_types = [event_types, {sprintf("Ratings Trial %d", song_idx)}];

        % Create EEG struct for each song (with only relevant data)
        chanlocs = readlocs(fullfile(path_to_ds, spf_file));
        EEG = pop_importdata('dataformat', 'array', 'nbchan', n_channels, 'data', EEG_Data, ...
                'chanlocs', chanlocs, 'setname', ['sub-' subject_id '_song_' num2str(song_idx)], ...
                'srate', Fs, 'pnts', size(EEG_Data, 2), 'xmin', 0);

        % Filter EEG data for the selected time range
        % Select data from baseline trial start to ratings trial start
        data_start_idx = find(EEG.times >= precise_onset_baseline, 1, 'first');
        data_end_idx = find(EEG.times >= ratings_start_onset, 1, 'first');
        
        EEG.data = EEG.data(:, data_start_idx:data_end_idx); % Crop EEG data to the desired segment
        EEG.times = EEG.times(data_start_idx:data_end_idx);  % Adjust the times accordingly
        
        % Create the event structure for this segment of data
        EEG.event = [];
        for i = 1:length(event_latencies)
            event_idx = find(EEG.times >= event_latencies(i), 1, 'first');
            EEG.event(i).type = event_types{i};
            EEG.event(i).latency = event_idx;
            EEG.event(i).epoch = 1;
            EEG.epoch(i).event = i;
        end
        
        EEG = eeg_checkset(EEG);  % Ensure EEG is valid
        
        % Save the EEG data for the specific song and segment
        output_filename = ['sub-' subject_str '_song-' song_str '.set'];
        pop_saveset(EEG, 'filename', output_filename, 'filepath', song_folder); % Save .set file
        
        % Clear event latencies and event types for the next song
        event_latencies = [];
        event_types = {};
    end
end
