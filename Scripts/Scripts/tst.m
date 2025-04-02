%% Script to assess quality of dsDRYAD dataset

%clear; clc;
n_songs = 12; % Number of unique songs
n_subs = 4;   % Number of subjects
path_to_ds = 'D:\master\Iniciacion_Investigacion\Datasets\ds003774';
sub_base_name = 'sub-';
song_base_name = 'ses-';

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

        filename = [subject_str '_' song_str '_task-MusicListening_run-' num2str(song_idx) '_eeg.set'];
        filepath = [path_to_ds '\' subject_str '\' song_str '\eeg'];

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

                % We set distance between channel with itself to nan
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
        neighbour_threshold = median(max(distances, [], 2, "omitnan"), "omitnan")/2;
        neighbor_matrix = distances <= neighbour_threshold; % Neighbour mask (0 or 1)

        % Compute the average over the 3 trials
        data_avg = EEG.data(:, :);
        %data_avg = EEG.data(:, :, 1);

        % Ciretia 1: Detect flat channels (low standard deviation)
        channel_std = squeeze(std(data_avg, 0, 2));
        threshold_low = median(channel_std) - 1.5 * iqr(channel_std);
        c1_bad_channels = find(channel_std < threshold_low);

        % Criteria 2: Compute correlation with nearest channels
        corr_matrix = corrcoef(data_avg');
        minimum_corr = 1/2; % Minimum corr between neighbour channels
        c2_bad_channels = [];
        
        for channel = 1:n_channels
            neighbors = find(neighbor_matrix(channel, :) == 1);
            avg_corr = length(find((corr_matrix(channel, neighbors) > minimum_corr)));
            
            % If a channel has less than ONE neighbour with high correlation
            if avg_corr < 2
                c2_bad_channels = [c2_bad_channels; channel];
            end

        end

        % Criteria 3 & 4: Detect outliers in power spectrum
        [pxx, f] = pwelch(data_avg', [], [], [], EEG.srate); % Compute FFT
        low_power = mean(pxx(f >= 1 & f <= 10, :), 1);
        high_power = mean(pxx(f >= 65 & f <= 90, :), 1);

        % Outliers based on IQR
        c3_bad_channels = find(abs(low_power - median(low_power)) > 1.5 * iqr(low_power));
        c4_bad_channels = find(abs(high_power - median(high_power)) > 1.5 * iqr(high_power));

        %% Combine all bad channels

        bad_channels = unique([c1_bad_channels(:); c2_bad_channels(:); c3_bad_channels(:); c4_bad_channels(:)]);
        removed_percentage = (length(bad_channels) / n_channels) * 100;
        if removed_percentage <= 25
            usable_songs = usable_songs + 1;
        end

        %display(['Bad channels detected for ' subject_str '/' song_str ': ', num2str(bad_channels')]);
        
        %% Visualize the bad channels using topoplot

        % Create a binary vector for bad channels (1 for removed, 0 for kept)
        removed_mask = zeros(1, n_channels);
        removed_mask(bad_channels) = 1;  % Mark the bad channels

        % Subplot to create a single summary figure per subject
        subplot(2, 6, song_idx);
        topoplot(removed_mask, EEG.chanlocs, 'style', 'blank', 'electrodes', 'on');
        title(['Song ' num2str(song_idx)], 'Position', [0, 0.5, 1]);

        % Display the percentage of channels removed in the plot
        text(0, -0.6, sprintf('%.2f%% of channels flagged as bad', removed_percentage), ...
            'Color', 'red', 'FontSize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    
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
