%% Script to assess quality of MUSING dataset

clear; clc;
n_songs = 12; % Number of unique songs
n_subs = 1;   % Number of subjects
path_to_ds = 'D:\master\Iniciacion_Investigacion\Datasets\ds003774\sourcedata';
sub_base_name = 'sub-';
song_base_name = 'song-';
addpath(genpath('D:\eeglab\eeglab_current\eeglab2024.2\plugins\clean_rawdata2.10'));
usable_songs_per_subject = zeros(n_subs, 1);
all_bad_channels = []; % Store all bad channels for summary

for sub_idx = 1:n_subs
    display(['Processing subject ' num2str(sub_idx)]);
    subject_str = [sub_base_name num2str(sub_idx, '%03d')];
    usable_songs = 0;

    fig_topos = figure('Visible','off','Units','normalized','Position',[0 0 1 1]);
    fig_c1c2 = figure('Visible','off','Units','normalized','Position',[0 0 1 1]);
    fig_c3 = figure('Visible','off','Units','normalized','Position',[0 0 1 1]);
    fig_c4 = figure('Visible','off','Units','normalized','Position',[0 0 1 1]);

    for song_idx = 1:n_songs
        display(['Processing song ' num2str(song_idx)]);
        song_str = [song_base_name num2str(song_idx, '%02d')];
        filename = [subject_str '_' song_str '.set'];
        filepath = [path_to_ds '\' subject_str '\eeg\' song_str];
        EEG = pop_loadset('filename', filename, 'filepath', filepath);

        if EEG.srate > 256
            EEG = pop_resample(EEG, 256);
        end

        EEG = pop_reref(EEG, []);
        EEG = pop_eegfiltnew(EEG, 'locutoff', 0.4);

        n_channels = size(EEG.data, 1);
        chanlocs = [EEG.chanlocs.X; EEG.chanlocs.Y; EEG.chanlocs.Z]';
        distances = nan(n_channels, n_channels);
        for i = 1:n_channels
            for j = 1:n_channels
                if i == j
                    dist = nan;
                else
                    dist = sqrt((chanlocs(i, 1) - chanlocs(j, 1))^2 + ...
                                (chanlocs(i, 2) - chanlocs(j, 2))^2 + ...
                                (chanlocs(i, 3) - chanlocs(j, 3))^2);
                end
                distances(i, j) = dist;
            end
        end

        neighbour_threshold = median(max(distances, [], 2, "omitnan"), "omitnan")/4;
        neighbor_matrix = distances <= neighbour_threshold;
        data = EEG.data;

        % Criteria 1: Detect flat channels
        EEG_clean_flat = clean_flatlines(EEG);
        c1_bad_labels = setdiff({EEG.chanlocs.labels}, {EEG_clean_flat.chanlocs.labels});
        c1_bad_channels = find(ismember({EEG.chanlocs.labels}, c1_bad_labels));

        % Criteria 2: Low correlation with neighbour channels
        % EEG_clean = clean_channels(EEG);
        % c2_bad_labels = setdiff({EEG.chanlocs.labels}, {EEG_clean.chanlocs.labels});
        % c2_bad_channels = find(ismember({EEG.chanlocs.labels}, c2_bad_labels));

        corrs_all = zeros(n_channels,1);
        for channel = 1:n_channels
            neighbors = find(neighbor_matrix(channel, :) == 1);
            neighbor_avg = mean(data(neighbors, :), 1);
            corrs_all(channel) = corr(data(channel, :)', neighbor_avg');
        end
        c2_bad_channels = find(corrs_all < 0.8);
        % Criteria 3 & 4: Spectral Power Outliers
        [pxx, f] = pwelch(data', [], [], [], EEG.srate);
        low_band = mean(pxx(f >= 1 & f <= 10, :), 1);
        high_band = mean(pxx(f >= 65 & f <= 90, :), 1);
        w = 3;
        detect_outliers = @(x, w) ...
            find(x < quantile(x, 0.25) - w * iqr(x) | x > quantile(x, 0.75) + w * iqr(x));
        c3_bad_channels = detect_outliers(low_band, w);
        c4_bad_channels = detect_outliers(high_band, w);

        bad_channels = unique([c1_bad_channels(:); c2_bad_channels(:); c3_bad_channels(:); c4_bad_channels(:)]);
        removed_percentage = (length(bad_channels) / n_channels) * 100;
        if removed_percentage <= 25
            usable_songs = usable_songs + 1;
        end

        % Visualization
        % Topoplot agrupado
        set(0, 'CurrentFigure', fig_topos);
        set(fig_topos, 'Visible', 'off');
        subplot(4, 3, song_idx); hold on;
        removed_mask = zeros(1, n_channels);
        removed_mask(bad_channels) = 1;
        topoplot(removed_mask, EEG.chanlocs, 'style', 'blank', 'electrodes', 'on');
        title(['Song ' num2str(song_idx)]);

        % Calculate removed percentage per criteria
        pct_c1 = 100 * numel(c1_bad_channels) / n_channels;
        pct_c2 = 100 * numel(c2_bad_channels) / n_channels;
        pct_c3 = 100 * numel(c3_bad_channels) / n_channels;
        pct_c4 = 100 * numel(c4_bad_channels) / n_channels;

        % Show all results
        text(0, -0.5, sprintf('Tot: %.1f%%', removed_percentage), ...
            'Color', 'red', 'FontSize', 9, 'HorizontalAlignment', 'center');
        text(0, -0.65, sprintf('C1: %.1f%%', pct_c1), ...
            'Color', 'k', 'FontSize', 8, 'HorizontalAlignment', 'center');
        text(0, -0.8, sprintf('C2: %.1f%%', pct_c2), ...
            'Color', 'k', 'FontSize', 8, 'HorizontalAlignment', 'center');
        text(0, -0.95, sprintf('C3: %.1f%%', pct_c3), ...
            'Color', 'k', 'FontSize', 8, 'HorizontalAlignment', 'center');
        text(0, -1.1, sprintf('C4: %.1f%%', pct_c4), ...
            'Color', 'k', 'FontSize', 8, 'HorizontalAlignment', 'center');
    
        % Series temporales C1 y C2
        set(0, 'CurrentFigure', fig_c1c2);
        set(fig_c1c2, 'Visible', 'off');
        subplot(6, 2, song_idx); hold on;
        for ch = c1_bad_channels
            plot(EEG.times, data(ch,:) + 300*ch, 'r');
        end
        for ch = c2_bad_channels
            plot(EEG.times, data(ch,:) + 300*ch, 'b');
        end
        title(['Song ' num2str(song_idx)]);
        
        
    
        % Potencia espectral C3
        set(0, 'CurrentFigure', fig_c3);
        set(fig_c3, 'Visible', 'off');
        subplot(4, 3, song_idx); hold on;
        bar(low_band, 'FaceColor', 'red'); hold on;
        scatter(c3_bad_channels, low_band(c3_bad_channels), 10, 'black', 'filled');
        title(['Song ' num2str(song_idx)]);
        
    
        % Potencia espectral C4
        set(0, 'CurrentFigure', fig_c4);
        set(fig_c4, 'Visible', 'off');
        subplot(4, 3, song_idx); hold on;
        bar(high_band, 'FaceColor', 'red'); hold on;
        scatter(c4_bad_channels, high_band(c4_bad_channels), 10, 'black', 'filled');
        title(['Song ' num2str(song_idx)]);
       

    end
    
    sgtitle(fig_topos, 'Canales eliminados por criterio');
    sgtitle(fig_c1c2, 'Criterio 1: Flat (rojo), Criterio 2: Corr (azul)');
    sgtitle(fig_c3, 'Criterio 3: Potencia 1–10 Hz');
    sgtitle(fig_c4, 'Criterio 4: Potencia 65–90 Hz');
    usable_songs_per_subject(sub_idx) = usable_songs;

    % Exportar todas las figuras en orden
    filename_pdf = ['Subject_' num2str(sub_idx, '%03d') '_Quality.pdf'];
    exportgraphics(fig_topos, filename_pdf, 'ContentType','vector');
    exportgraphics(fig_c1c2, filename_pdf, 'Append', true);
    exportgraphics(fig_c3, filename_pdf, 'Append', true);
    exportgraphics(fig_c4, filename_pdf, 'Append', true);
    
    % Cierre de figuras
    close(fig_topos); close(fig_c1c2); close(fig_c3); close(fig_c4);
end

usable_songs_percentage = (sum(usable_songs_per_subject) / (n_songs * n_subs)) * 100;
usable_subjects = sum(usable_songs_per_subject >= 0.75 * n_songs);
usable_subjects_percentage = (usable_subjects / n_subs) * 100;

disp('--------------------------------');
disp('Dataset Quality Summary:');
disp(['Percentage of usable songs: ', num2str(usable_songs_percentage, '%.2f'), '%']);
disp(['Percentage of usable subjects: ', num2str(usable_subjects_percentage, '%.2f'), '%']);
disp('--------------------------------');
