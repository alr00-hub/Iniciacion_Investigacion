% EEG Quality Assessment Script for dsDRYAD dataset
clear; clc;
addpath(genpath('D:\eeglab\eeglab_current\eeglab2024.2\plugins\clean_rawdata2.10'));  % Update path as needed

% General Parameters
n_songs = 10;
n_subjects = 3;
path_to_ds = 'D:\master\Iniciacion_Investigacion\Datasets\dsDRYAD\diliBach_4dryad_CND\diliBach_4dryad_CND';
sub_base = 'sub-';
song_base = 'song-';
usable_songs_per_subject = zeros(n_subjects, 1);

%% Loop through subjects
for sub_idx = 3:n_subjects
    disp(['Processing subject ' num2str(sub_idx)]);
    subject_str = [sub_base num2str(sub_idx)];
    usable_songs = 0;
    figure;

    %% Loop through songs
    for song_idx = 1:n_songs
        disp(['Processing song ' num2str(song_idx)]);
        song_str = [song_base num2str(song_idx)];
        filename = [subject_str '_' song_str '.set'];
        filepath = fullfile(path_to_ds, subject_str, song_str);

        % Load EEG data
        EEG = pop_loadset('filename', filename, 'filepath', filepath);

        % Preprocessing
        if EEG.srate > 256
            EEG = pop_resample(EEG, 256);
        end
        EEG = pop_reref(EEG, []);
        EEG = pop_eegfiltnew(EEG, 'locutoff', 0.4); % High-pass filter

        % Compute average data
        data_avg = squeeze(mean(EEG.data, 3));
        n_channels = size(data_avg, 1);

        %% === Criterion 1: Flat channels (low STD) ===
        channel_std = squeeze(std(data_avg, 0, 2));
        std_threshold = 75; % µV
        c1_bad_channels = find(channel_std < std_threshold);

        %% === Criterion 2: Bad correlation via clean_rawdata (RANSAC) ===
        EEG_event = pop_selectevent(EEG, 'event', 1);
        % Aplicar clean_rawdata al EEG promedio
        %[cleaned_sig, removed_chan] = clean_channels(EEG_event);
        EEG_clean = clean_artifacts(EEG_event);
        
        % Identify bad channels based on RANSAC method
        c2_bad_labels = setdiff({EEG.chanlocs.labels}, {EEG_clean.chanlocs.labels});
        c2_bad_channels = find(ismember({EEG.chanlocs.labels}, c2_bad_labels));

        % Detectar canales malos

        %c2_bad_channels = find(removed_chan);

        %% === Criteria 3 & 4: Spectral Power Outliers ===
        [pxx, f] = pwelch(data_avg', [], [], [], EEG.srate);
        low_band = mean(pxx(f >= 1 & f <= 10, :), 1);
        high_band = mean(pxx(f >= 65 & f <= 90, :), 1);

        detect_outliers = @(x, w) ...
            find(x < quantile(x, 0.25) - w * iqr(x) | x > quantile(x, 0.75) + w * iqr(x));
        w = 3;

        c3_bad_channels = detect_outliers(low_band, w);
        c4_bad_channels = detect_outliers(high_band, w);

        %% === Combine all bad channels ===
        bad_channels = unique([c1_bad_channels(:); c2_bad_channels(:); c3_bad_channels(:); c4_bad_channels(:)]);
        %bad_channels = unique([c2_bad_channels(:)]);
        removed_percentage = (length(bad_channels) / n_channels) * 100;

        if removed_percentage <= 25
            usable_songs = usable_songs + 1;
        end

        %% === Topographical visualization of bad channels ===
        subplot(2, 5, song_idx);
        removed_mask = zeros(1, n_channels);
        removed_mask(bad_channels) = 1;
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

%% === Dataset Quality Summary ===
usable_songs_pct = (sum(usable_songs_per_subject) / (n_songs * n_subjects)) * 100;
usable_subjects = sum(usable_songs_per_subject >= 0.75 * n_songs);
usable_subjects_pct = (usable_subjects / n_subjects) * 100;

disp('--------------------------------');
disp('Dataset Quality Summary:');
disp(['% of usable songs: ', num2str(usable_songs_pct, '%.2f'), '%']);
disp(['% of usable subjects: ', num2str(usable_subjects_pct, '%.2f'), '%']);
disp('--------------------------------');
