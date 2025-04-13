%% Script to assess quality of dsDRYAD dataset

%clear; clc;
n_songs = 10; % Number of unique songs
n_subs = 3;   % Number of subjects
path_to_ds = 'D:\master\Iniciacion_Investigacion\Datasets\dsDRYAD\diliBach_4dryad_CND\diliBach_4dryad_CND';
sub_base_name = 'sub-';
song_base_name = 'song-';

usable_songs_per_subject = zeros(n_subs, 1);
all_bad_channels = []; % Store all bad channels for summary

%% Loop through subjects
for sub_idx = 3:n_subs
    display(['Processing subject ' num2str(sub_idx)]);
    subject_str = [sub_base_name num2str(sub_idx)];
    usable_songs = 0;
    figure;
    %% Loop through songs
    for song_idx = 1:n_songs
        display(['Processing song ' num2str(song_idx)]);
        song_str = [song_base_name num2str(song_idx)];
        
        %% Build path and file names
        filename = [subject_str '_' song_str '.set'];
        filepath = [path_to_ds '\' subject_str '\' song_str];

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
        %data_avg = squeeze(mean(EEG.data, 3));
        data_avg = EEG.data(:, :, 1);

        % Criterio 1: Detección de canales planos (desviación estándar baja)
        channel_std = squeeze(std(data_avg, 0, 2));
        
        % Establecer umbral fijo
        fixed_threshold = 75;  % Ajusta según el dataset
        c1_bad_channels = find(channel_std < fixed_threshold);

        % Criteria 2 (Modified): Correlate each channel with the AVERAGE of its neighbors
        c2_bad_channels = [];
        
        corrs_all = zeros(n_channels,1);
        for channel = 1:n_channels
            neighbors = find(neighbor_matrix(channel, :) == 1);
            neighbor_avg = median(data_avg(neighbors, :), 1);
            corrs_all(channel) = corr(data_avg(channel, :)', neighbor_avg');
        end
        
        % Z-score across all correlations
        %z_corrs = (corrs_all - mean(corrs_all)) / std(corrs_all);
        %c2_bad_channels = find(z_corrs < -0.5);  % flag channels with z-score below
        c2_bad_channels = find(corrs_all < 0.7);  % flag channels with z-score below

        %% === Criteria 3 & 4: Spectral Power Outliers ===
        [pxx, f] = pwelch(data_avg', [], [], [], EEG.srate);
        low_band = mean(pxx(f >= 1 & f <= 10, :), 1);
        high_band = mean(pxx(f >= 65 & f <= 90, :), 1);

        detect_outliers = @(x, w) ...
            find(x < quantile(x, 0.25) - w * iqr(x) | x > quantile(x, 0.75) + w * iqr(x));
        w = 3;

        c3_bad_channels = detect_outliers(low_band, w);
        c4_bad_channels = detect_outliers(high_band, w);

        %% Combine all bad channels

        bad_channels = unique([c1_bad_channels(:); c2_bad_channels(:); c3_bad_channels(:); c4_bad_channels(:)]);
        %bad_channels = unique([c2_bad_channels(:)]);
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
        subplot(2, 5, song_idx);
        topoplot(removed_mask, EEG.chanlocs, 'style', 'blank', 'electrodes', 'on');
        title(['Song ' num2str(song_idx)], 'Position', [0, 0.5, 1]);

        pct_c1 = 100 * numel(c1_bad_channels) / n_channels;
        pct_c2 = 100 * numel(c2_bad_channels) / n_channels;
        pct_c3 = 100 * numel(c3_bad_channels) / n_channels;
        pct_c4 = 100 * numel(c4_bad_channels) / n_channels;
 
        % Mostrar los porcentajes debajo del topoplot
        text(0, -0.55, sprintf('Tot: %.1f%%', removed_percentage), ...
            'Color', 'red', 'FontSize', 9, 'HorizontalAlignment', 'center');
        text(0, -0.65, sprintf('C1: %.1f%%', pct_c1), ...
            'Color', 'k', 'FontSize', 8, 'HorizontalAlignment', 'center');
        text(0, -0.72, sprintf('C2: %.1f%%', pct_c2), ...
            'Color', 'k', 'FontSize', 8, 'HorizontalAlignment', 'center');
        text(0, -0.79, sprintf('C3: %.1f%%', pct_c3), ...
            'Color', 'k', 'FontSize', 8, 'HorizontalAlignment', 'center');
        text(0, -0.86, sprintf('C4: %.1f%%', pct_c4), ...
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


