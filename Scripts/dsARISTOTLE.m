%%%%%%%%%%%%%%%     Script for dsARISTOTLE      %%%%%%%%%%%%%%%

%% Basic Setup
clear; clc; 
n_songs = 30;
path_to_ds = 'D:\master\Iniciacion_Investigacion\Datasets\dsARISTOTLE\';

%% Load .mat file
sub1 = load([path_to_ds 'music_listening_experiment_s01.mat']);

%% Extract variables
EEG_Rest = sub1.EEG_Rest;   % (n_channels × n_time_points)
EEG_Songs = sub1.EEG_Songs; % (n_songs × n_channels × n_time_points)
Fs = sub1.Fs;               % Sampling frequency
chan_info = sub1.sensor_info;
n_channels = size(chan_info.loc, 1); % Number of channels
song_ratings = sub1.song_ratings; % Ratings of the songs

%% Concatenate EEG data
EEG_concat = EEG_Rest; % Start with resting EEG
event_latencies = [1]; % Store event latencies (first event is at time 0)

for song_idx = 1:n_songs
    event_latencies = [event_latencies, size(EEG_concat, 2) + 1]; % Save event time
    EEG_concat = [EEG_concat, squeeze(EEG_Songs(song_idx, :, :))]; % Update concat 
end

%% Prepare the chanlocs struct
% chan_info.loc only has 2D (X, Y), so we add Z = 0 to convert to 3D


chanlocs = struct('labels', [], 'X', [], 'Y', [], 'Z', []); % Initialize struct

% Load the data into the chanlocs struct
for i = 1:n_channels
    chanlocs(i).labels = chan_info.labels{i}    % Labels
    chanlocs(i).X = chan_info.loc(i, 1);        % X coord
    chanlocs(i).Y = chan_info.loc(i, 2);        % Y coord
    chanlocs(i).Z = 0;                          % Z coord, assuming 0 for simplicity (2D plane)
end



%% Creating the EEG .set file
EEG = pop_importdata('dataformat', 'matlab', 'nbchan', n_channels, 'data', EEG_concat, ...
    'chanlocs', chanlocs, 'setname', 'kk', 'srate', Fs, 'pnts', 0, 'xmin', 0);

%% Adding events
EEG.event = [];

EEG.event(1).type = 'rest_start'; % First event is resting
EEG.event(1).latency = event_latencies(1); % Load first latency (0) 

% Loop across al songs
for song_idx = 1:n_songs
    EEG.event(end+1).type = ['song:' num2str(song_idx) '/rating:' num2str(song_ratings(song_idx))]; % Event name: song/rating
    EEG.event(end).latency = event_latencies(song_idx + 1); % Add latency
end

EEG = eeg_checkset(EEG);  % Check consistency

% Save the dataset
pop_saveset(EEG, 'filename', 'Subject1_AllEEG.set', 'filepath', path_to_ds);
