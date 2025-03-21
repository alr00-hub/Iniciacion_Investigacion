%% Script to assess quality of dsDRYAD dataset

clear; clc;
n_trials = 3; % Each song appears 3 times
n_songs = 1; % Number of unique songs
n_subs = 1;   % Number of subjects
path_to_ds = 'D:\master\Iniciacion_Investigacion\Datasets\dsDRYAD\diliBach_4dryad_CND\diliBach_4dryad_CND';
sub_base_name = 'sub-';
song_base_name = 'song-';

%% Loop through subjects
for sub_idx = 1:n_subs
    display(['Processing subject ' num2str(sub_idx)]);
    subject_str = [sub_base_name num2str(sub_idx)];

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


        % Degine chanlocs for 3D
        chanlocs = [EEG.chanlocs.X; EEG.chanlocs.Y; EEG.chanlocs.Z]';

        % Compute euclidean distance between each 2-pair of channels
        distances = nan(n_channels, n_channels);
        for i = 1:n_channels
            for j = i+1:n_channels
                % Calcular distancia euclidiana entre el canal i y el canal j
                dist = sqrt((chanlocs(i, 1) - chanlocs(j, 1))^2 + ...
                            (chanlocs(i, 2) - chanlocs(j, 2))^2 + ...
                            (chanlocs(i, 3) - chanlocs(j, 3))^2);
                distances(i, j) = dist;
            end
        end

        neighbour_threshold = mean(max(distances, [], 2, "omitnan"), "omitnan")/4;
        neighbor_matrix = distances <= neighbour_threshold;


        % Ciretia 1: Detect flat channels (low standard deviation)
        channel_std = squeeze(std(EEG.data, 0, 2));
        threshold_low = median(channel_std) - 1.5 * iqr(channel_std);
        c1_bad_channels = find(channel_std < threshold_low);

        % Criteria 2: Compute correlation with nearest channels
        data_avg = squeeze(mean(EEG.data, 3));
        corr_matrix = corrcoef(data_avg');
        minimum_corr = 1/3; % Minimum corr between neighbour channels
        c2_bad_channels = [];
        
        for channel = 1:n_channels
            neighbors = find(neighbor_matrix(channel, :) == 1);
            avg_corr = mean(corr_matrix(channel, neighbors));
            
            if avg_corr < minimum_corr
                c2_bad_channels = [c2_bad_channels; channel];
            end

        end

        % Step 3: Detect outliers in power spectrum
        [pxx, f] = pwelch(data_avg', [], [], [], EEG.srate);
        low_freq_idx = f >= 0.5 & f <= 10;
        high_freq_idx = f >= 65 & f <= 90;

        low_power = mean(pxx(low_freq_idx, :), 1);
        high_power = mean(pxx(high_freq_idx, :), 1);

        % Outliers based on IQR
        threshold_low_power = median(low_power) + 1.5 * iqr(low_power);
        threshold_high_power = median(high_power) + 1.5 * iqr(high_power);
        
        bad_channels_low_power = find(low_power > threshold_low_power);
        bad_channels_high_power = find(high_power > threshold_high_power);

        %% Combine all bad channels
        c1_bad_channels = c1_bad_channels(:); % Convertir a columna
        c2_bad_channels = c2_bad_channels(:);
        bad_channels_low_power = bad_channels_low_power(:);
        bad_channels_high_power = bad_channels_high_power(:);
        
        bad_channels = unique([c1_bad_channels; c2_bad_channels; bad_channels_low_power; bad_channels_high_power]);
        %bad_channels = unique([bad_channels_flat; bad_channels_corr; bad_channels_low_power; bad_channels_high_power]);
        
        display(['Bad channels detected: ', num2str(bad_channels')]);
    end
end
