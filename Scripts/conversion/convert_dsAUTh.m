%% Script for dsAUTh conversion

%% Basic Setup
clear; clc;
n_songs = 30;
n_subs = 20;
path_to_ds = 'D:\master\Iniciacion_Investigacion\Datasets\dsARISTOTLE\';
file_base_name = 'music_listening_experiment_s';

%% Process each subject
for sub_idx = 1:n_subs
    
    % Converting sub number to 2-digit string
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
    EEG_Rest = subject_data.EEG_Rest;   % Resting EEG
    EEG_Songs = subject_data.EEG_Songs; % Songs EEG
    Fs = subject_data.Fs;               % Sampling frequency
    chan_info = subject_data.sensor_info; % Contains labels and 2D locs
    n_channels = size(chan_info.loc, 1); % Number of channels
    song_ratings = subject_data.song_ratings; % Song ratings
    
    % Prepare channel locations struct
    chanlocs = struct(...
        'labels', chan_info.labels, ...
        'X', num2cell(chan_info.loc(:, 1)), ...
        'Y', num2cell(chan_info.loc(:, 2)), ...
        'Z', num2cell(zeros(size(chan_info.loc, 1), 1)) ... % Assuming Z=0 for simplicity
    );
    
    %% Save each song separately as a .set file
    for song_idx = 1:n_songs
        
        song_str = sprintf('%02d', song_idx);
        % Create separate folder for each song
        song_folder = fullfile(subject_folder, ['song-' song_str]);
        if ~exist(song_folder, 'dir')
            mkdir(song_folder);
        end
        
        % Extract EEG data for the song
        EEG_song = squeeze(EEG_Songs(song_idx, :, :));
        
        % Create EEG structure
        EEG = pop_importdata('dataformat', 'matlab', 'nbchan', n_channels, 'data', EEG_song, ...
            'chanlocs', chanlocs, 'setname', ['sub-' subject_str '_song-' song_str], ...
            'srate', Fs, 'pnts', size(EEG_song, 2), 'xmin', 0);
        
        % Add event marker with song rating
        EEG.event = struct();
        EEG.event(1).type = ['song:' song_str '/rating:' num2str(song_ratings(song_idx))];
        EEG.event(1).latency = 1; % Mark event at the beginning
        EEG = eeg_checkset(EEG);

           
        % Save .set file
        output_filename = ['sub-' subject_str '_song-' song_str '.set'];
        pop_saveset(EEG, 'filename', output_filename, 'filepath', song_folder);
    end
    
    fprintf('Conversion completed for %s\n', subject_id);
end
