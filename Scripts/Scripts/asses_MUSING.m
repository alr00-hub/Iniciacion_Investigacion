%% Script to assess quality of MUSING dataset

clear; clc;
n_songs = 12; % Number of unique songs
n_subs = 20;   % Number of subjects
path_to_ds = 'D:\master\Iniciacion_Investigacion\Datasets\ds003774\sourcedata';
sub_base_name = 'sub-';
song_base_name = 'song-';
addpath(genpath('D:\eeglab\eeglab_current\eeglab2024.2\plugins\clean_rawdata2.10'));
usable_songs_per_subject = zeros(n_subs, 1);
all_bad_channels = []; % Store all bad channels for summary

%% Loop through subjects
for sub_idx = 1:n_subs
    display(['Processing subject ' num2str(sub_idx)]);
    subject_str = [sub_base_name num2str(sub_idx, '%03d')];
    usable_songs = 0;
    figure;
    %% Loop through songs
    for song_idx = 1:n_songs
        display(['Processing song ' num2str(song_idx)]);
        song_str = [song_base_name num2str(song_idx, '%02d')];
        
        %% Build path and file names
        filename = [subject_str '_' song_str '.set'];
        filepath = [path_to_ds '\' subject_str '\eeg\' song_str];

        %% Load the data from EEGLAB file (.set)
        EEG = pop_loadset('filename', filename, 'filepath', filepath);

        %% Downsample if sfreq > 256
        if EEG.srate > 256
            EEG = pop_resample(EEG, 256);
        end

        %% Re-reference data (average)
        EEG = pop_reref(EEG, []);
        
        %% High pass filter at 0.4 Hz
        EEG = pop_eegfiltnew(EEG, 'locutoff', 0.4);

        %% Criteria for bad channels detection
        n_channels = size(EEG.data, 1);

        % Define chanlocs for 3D
        chanlocs = [EEG.chanlocs.X; EEG.chanlocs.Y; EEG.chanlocs.Z]';

        % Compute Euclidean distance between each 2-pair of channels
        distances = nan(n_channels, n_channels);
        for i = 1:n_channels
            for j = 1:n_channels

                % We set the distance between channel with itself to nan
                if i == j
                    dist = nan;
                
                % Compute Euclidean distance
                else
                    dist = sqrt((chanlocs(i, 1) - chanlocs(j, 1))^2 + ...
                                (chanlocs(i, 2) - chanlocs(j, 2))^2 + ...
                                (chanlocs(i, 3) - chanlocs(j, 3))^2);
                end
                distances(i, j) = dist;
            end
        end

        % Set the neighbour threshold to 1/4 of the max distance
        neighbour_threshold = median(max(distances, [], 2, "omitnan"), "omitnan")/4;
        neighbor_matrix = distances <= neighbour_threshold; % Neighbour mask (0 or 1)
        
        % We only have 1 trial per song
        data = EEG.data;

        %% Criteria 1: Detect flat channels
        c1_bad_channels = [];
        % 1st Option -> Using clean_flatlines (default parameters)
        EEG_clean_flat = clean_flatlines(EEG);
        c1_bad_labels = setdiff({EEG.chanlocs.labels}, {EEG_clean_flat.chanlocs.labels});
        c1_bad_channels = find(ismember({EEG.chanlocs.labels}, c1_bad_labels));
        
        % 2nd Option -> Looking for abnormally low std (some fixed value)
        %channel_std = squeeze(std(data, 0, 2));
        %std_threshold = 75; % µV
        %c1_bad_channels = find(channel_std < std_threshold);

        %% Criteria 2: Low correlation with neighbour channels
        c2_bad_channels = [];
        % 1st Option -> Using clean_channels (default parameters) 
        EEG_clean = clean_channels(EEG);
        c2_bad_labels = setdiff({EEG.chanlocs.labels}, {EEG_clean.chanlocs.labels});
        c2_bad_channels = find(ismember({EEG.chanlocs.labels}, c2_bad_labels));
        
        % 2nd Option -> Computing correlation for each channel with the
        % MEDIAN of its neighbours

        % corrs_all = zeros(n_channels,1);
        % for channel = 1:n_channels
        %     neighbors = find(neighbor_matrix(channel, :) == 1);
        %     neighbor_avg = mean(data(neighbors, :), 1);
        %     corrs_all(channel) = corr(data(channel, :)', neighbor_avg');
        % end
        % 
        % % Z-score across all correlations
        % z_corrs = (corrs_all - mean(corrs_all)) / std(corrs_all);
        % c2_bad_channels = find(z_corrs < -1.25);  
        % Using some fixed correlation thershold
        % c2_bad_channels = find(corrs_all < 0.8);

        %% Criteria 3 & 4: Spectral Power Outliers
        [pxx, f] = pwelch(data', [], [], [], EEG.srate);
        low_band = mean(pxx(f >= 1 & f <= 10, :), 1);
        high_band = mean(pxx(f >= 65 & f <= 90, :), 1);
        w = 3;

        detect_outliers = @(x, w) ...
            find(x < quantile(x, 0.25) - w * iqr(x) | x > quantile(x, 0.75) + w * iqr(x));
        
        c3_bad_channels = detect_outliers(low_band, w);
        c4_bad_channels = detect_outliers(high_band, w);

        %% Combine all bad channels

        bad_channels = unique([c1_bad_channels(:); c2_bad_channels(:); c3_bad_channels(:); c4_bad_channels(:)]);

        removed_percentage = (length(bad_channels) / n_channels) * 100;
        if removed_percentage <= 25
            usable_songs = usable_songs + 1;
        end

        
        %% Visualize the bad channels using topoplot

        % Create a binary vector for bad channels (1 for removed, 0 for kept)
        removed_mask = zeros(1, n_channels);
        removed_mask(bad_channels) = 1;  % Mark the bad channels

        % Subplot to create a single summary figure per subject
        subplot(2, 6, song_idx);
        topoplot(removed_mask, EEG.chanlocs, 'style', 'blank', 'electrodes', 'on');
        title(['Song ' num2str(song_idx)], 'Position', [0, 0.5, 1]);

        % Calculate removed percentage per criteria
        pct_c1 = 100 * numel(c1_bad_channels) / n_channels;
        pct_c2 = 100 * numel(c2_bad_channels) / n_channels;
        pct_c3 = 100 * numel(c3_bad_channels) / n_channels;
        pct_c4 = 100 * numel(c4_bad_channels) / n_channels;

        % Show all results
        text(0, -0.5, sprintf('Tot: %.1f%%', removed_percentage), ...
            'Color', 'red', 'FontSize', 9, 'HorizontalAlignment', 'center');
        text(0, -0.6, sprintf('C1: %.1f%%', pct_c1), ...
            'Color', 'k', 'FontSize', 8, 'HorizontalAlignment', 'center');
        text(0, -0.7, sprintf('C2: %.1f%%', pct_c2), ...
            'Color', 'k', 'FontSize', 8, 'HorizontalAlignment', 'center');
        text(0, -0.8, sprintf('C3: %.1f%%', pct_c3), ...
            'Color', 'k', 'FontSize', 8, 'HorizontalAlignment', 'center');
        text(0, -0.9, sprintf('C4: %.1f%%', pct_c4), ...
            'Color', 'k', 'FontSize', 8, 'HorizontalAlignment', 'center');
    
    end

    usable_songs_per_subject(sub_idx) = usable_songs;
end

usable_songs_percentage = (sum(usable_songs_per_subject) / (n_songs * n_subs)) * 100;
usable_subjects = sum(usable_songs_per_subject >= 0.75 * n_songs);
usable_subjects_percentage = (usable_subjects / n_subs) * 100;

%% Display summary results
disp('--------------------------------');
disp('Dataset Quality Summary:');
disp(['Percentage of usable songs: ', num2str(usable_songs_percentage, '%.2f'), '%']);
disp(['Percentage of usable subjects: ', num2str(usable_subjects_percentage, '%.2f'), '%']);
disp('--------------------------------');


