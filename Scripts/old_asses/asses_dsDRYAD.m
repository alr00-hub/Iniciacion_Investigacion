%% Script to assess quality of DRYAD dataset

clear; clc;
n_songs = 10; % Number of unique songs
n_subs = 5;   % Number of subjects
path_to_ds = 'D:\master\Iniciacion_Investigacion\Datasets\dsDRYAD\diliBach_4dryad_CND\diliBach_4dryad_CND';
sub_base_name = 'sub-';
song_base_name = 'song-';
addpath(genpath('D:\eeglab\eeglab_current\eeglab2024.2\plugins\clean_rawdata2.10'));
usable_songs_per_subject = zeros(n_subs, 1);
all_bad_channels = []; % Store all bad channels for summary
export = true;

for sub_idx = 5:n_subs
    display(['Processing subject ' num2str(sub_idx)]);
    subject_str = [sub_base_name num2str(sub_idx, '%d')];
    usable_songs = 0;

    if export
        export_title_page(sub_idx);
        fig_topos = figure('Visible','off','Units','normalized','Position',[1 1 1 1]);
        fig_c1c2 = figure('Visible','off','Units','normalized','Position',[1 1 1 1]);
        fig_c3 = figure('Visible','off','Units','normalized','Position',[1 1 1 1]);
        fig_c4 = figure('Visible','off','Units','normalized','Position',[1 1 1 1]);
    end

    for song_idx = 1:n_songs
        display(['Processing song ' num2str(song_idx)]);
        song_str = [song_base_name num2str(song_idx, '%d')];
        filename = [subject_str '_' song_str '.set'];
        filepath = [path_to_ds '\' subject_str '\' song_str];
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
        data = squeeze(mean(EEG.data, 3));

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
            neighbor_median = median(data(neighbors, :), 1);
            corrs_all(channel) = corr(data(channel, :)', neighbor_median');
        end
        %z_corrs = (corrs_all - mean(corrs_all)) / std(corrs_all);
        %c2_bad_channels = find(z_corrs < -1.25);  
        c2_bad_channels = find(corrs_all < 0.75);
        % Criteria 3 & 4: Spectral Power Outliers
        [pxx, f] = pwelch(data', [], [], [], EEG.srate);
        low_freq_idx = f >= 1 & f <= 10;
        high_freq_idx = f >= 65 & f <= 90;
        low_band_data = pxx(low_freq_idx, :);
        high_band_data = pxx(high_freq_idx, :);
        low_band = mean(low_band_data, 1);
        high_band = mean(high_band_data, 1);
        w = 3;
        detect_outliers = @(x, w) ...
            find(x < quantile(x, 0.25) - w * iqr(x) | x > quantile(x, 0.75) + w * iqr(x));
        c3_bad_channels = detect_outliers(low_band, w);
        c4_bad_channels = detect_outliers(high_band, w);

        bad_channels = unique([c1_bad_channels(:); c2_bad_channels(:); c3_bad_channels(:); c4_bad_channels(:)]);
        removed_percentage = (length(bad_channels) / n_channels) * 100;
        if removed_percentage <= 20
            usable_songs = usable_songs + 1;
        end

        % Visualization
        if export
            create_figures(sub_idx, song_idx, EEG, data, n_channels, removed_percentage, ...
            neighbor_matrix, c1_bad_channels, c2_bad_channels, c3_bad_channels, c4_bad_channels, ...
            bad_channels, low_band, high_band, fig_topos, fig_c1c2, fig_c3, fig_c4);
        end

    end
    
    usable_songs_per_subject(sub_idx) = usable_songs;
    if export
        export_report_to_pdf(sub_idx, fig_topos, fig_c1c2, fig_c3, fig_c4);
    end
end

usable_songs_percentage = (sum(usable_songs_per_subject) / (n_songs * n_subs)) * 100;
usable_subjects = sum(usable_songs_per_subject >= 0.75 * n_songs);
usable_subjects_percentage = (usable_subjects / n_subs) * 100;

disp('--------------------------------');
disp('Dataset Quality Summary:');
disp(['Percentage of usable songs: ', num2str(usable_songs_percentage, '%.2f'), '%']);
disp(['Percentage of usable subjects: ', num2str(usable_subjects_percentage, '%.2f'), '%']);
disp('--------------------------------');


function export_title_page(sub_idx)
    filename_pdf = ['Subject_' num2str(sub_idx, '%03d') '_Quality.pdf'];
    fig_title = figure('Visible','off','Units','normalized','Position',[0 0 1 1]);
    axis off;
    text(0.5, 0.6, sprintf('dsDRYAD'), ...
        'FontSize', 34, 'HorizontalAlignment', 'center');
    text(0.5, 0.4, sprintf('Quality Report\nSubject %03d', sub_idx), ...
        'FontSize', 24, 'HorizontalAlignment', 'center');
    text(0.5, 0.2, datestr(now, 'mmmm dd, yyyy'), ...
        'FontSize', 14, 'HorizontalAlignment', 'center');
    exportgraphics(fig_title, filename_pdf, 'Append', false, 'ContentType', 'image', 'Resolution', 300);
    close(fig_title);
end

function create_figures(sub_idx, song_idx, EEG, data, n_channels, removed_percentage, ...
    neighbor_matrix, c1_bad_channels, c2_bad_channels, c3_bad_channels, c4_bad_channels, ...
    bad_channels, low_band, high_band, fig_topos, fig_c1c2, fig_c3, fig_c4)

    filename_pdf = ['Subject_' num2str(sub_idx, '%03d') '_Quality.pdf'];

    % Grouped Topoplot
    set(0, 'CurrentFigure', fig_topos);
    set(fig_topos, 'Visible', 'off');
    subplot(2, 5, song_idx); hold on;
    removed_mask = zeros(1, n_channels);
    removed_mask(bad_channels) = 1;
    topoplot(removed_mask, EEG.chanlocs, 'style', 'blank', 'electrodes', 'on');

    % Calculate removed percentage per criteria
    pct_c1 = 100 * numel(c1_bad_channels) / n_channels;
    pct_c2 = 100 * numel(c2_bad_channels) / n_channels;
    pct_c3 = 100 * numel(c3_bad_channels) / n_channels;
    pct_c4 = 100 * numel(c4_bad_channels) / n_channels;

    % Show all results
    text(0, -0.55, sprintf('Song: %d - Total Removed: %.1f%%', song_idx, removed_percentage), ...
        'Color', 'red', 'FontSize', 9, 'HorizontalAlignment', 'center');
    text(0, -0.75, sprintf('C1: %.1f%% || C2: %.1f%% || C3: %.1f%% || C4: %.1f%%', ...
        pct_c1, pct_c2, pct_c3, pct_c4), 'Color', 'k', 'FontSize', 8, ...
        'HorizontalAlignment', 'center');

    % Temporal series C1 y C2
    set(0, 'CurrentFigure', fig_c1c2);
    set(fig_c1c2, 'Visible', 'off');
    subplot(10, 1, song_idx); hold on;
    if ~isempty(c1_bad_channels)
        for ch = c1_bad_channels
            plot(EEG.times(:), squeeze(data(ch,:)), 'b');
        end
    end
    
    if ~isempty(c2_bad_channels)
        for ch = c2_bad_channels
            plot(EEG.times(:), squeeze(data(ch,:)), 'r');
        end
    end
    title(['Song ' num2str(song_idx)]);

    % C2 vs neighbors
    n_neighbors = sum(neighbor_matrix, 2);
    peripheral_threshold = prctile(n_neighbors, 10);
    peripheral_channels = find(n_neighbors <= peripheral_threshold);
    peripheral_c2_bad = intersect(c2_bad_channels, peripheral_channels);
    peripheral_c2_good = setdiff(peripheral_channels, c2_bad_channels);
    n_bad = min(length(peripheral_c2_bad), 6);
    n_good = min(length(peripheral_c2_good), 6);
    peripheral_c2_bad = peripheral_c2_bad(1:n_bad);
    peripheral_c2_good = peripheral_c2_good(1:n_good);
    if ~isempty(peripheral_c2_bad)
        fig_c2_neighbors_bad = figure('Visible','off','Units','normalized','Position',[1 1 1 1]);
        set(fig_c2_neighbors_bad, 'Visible', 'off');

        for i = 1:n_bad
            ch = peripheral_c2_bad(i);
            neighbors = find(neighbor_matrix(ch, :) == 1);
            subplot(n_bad, 1, i); hold on;

            offset = 0;
            yticks_vals = [];
            yticklabels_str = {};

            % Neighbors (blue)
            for nb = neighbors
                plot(EEG.times, data(nb,:) + offset, 'b');
                yticks_vals(end+1) = offset;
                yticklabels_str{end+1} = EEG.chanlocs(nb).labels;
                offset = offset + 1000;
            end
            % Bad channel (red)
            plot(EEG.times, data(ch,:) + offset, 'r');
            yticks_vals(end+1) = offset;
            yticklabels_str{end+1} = EEG.chanlocs(ch).labels;

            % Labels
            set(gca, 'YTick', yticks_vals, 'YTickLabel', yticklabels_str, 'YTickLabelRotation', 0);
            title(['Channel ' EEG.chanlocs(ch).labels ' (number of neighbours: ' num2str(numel(neighbors)) ')']);
            xlabel('Time (ms)');
            ylabel('Channels');
        end

        sgtitle(fig_c2_neighbors_bad, ['C2 check:Time series of a bad channel (red) and its neighbours (blue) - Song ' num2str(song_idx)]);
        exportgraphics(fig_c2_neighbors_bad, filename_pdf, 'Append', true, 'ContentType', 'image', 'Resolution', 300);
    end

    if ~isempty(peripheral_c2_good)
        fig_c2_neighbors_good = figure('Visible','off','Units','normalized','Position',[1 1 1 1]);
        set(fig_c2_neighbors_good, 'Visible', 'off');

        for i = 1:n_good
            ch = peripheral_c2_good(i);
            neighbors = find(neighbor_matrix(ch, :) == 1);
            subplot(n_good, 1, i); hold on;

            offset = 0;
            yticks_vals = [];
            yticklabels_str = {};

            % Neighbors (blue)
            for nb = neighbors
                plot(EEG.times, data(nb,:) + offset, 'b');
                yticks_vals(end+1) = offset;
                yticklabels_str{end+1} = EEG.chanlocs(nb).labels;
                offset = offset + 1000;
            end
            % Good channel (red)
            plot(EEG.times, data(ch,:) + offset, 'r');
            yticks_vals(end+1) = offset;
            yticklabels_str{end+1} = EEG.chanlocs(ch).labels;

            % Labels
            set(gca, 'YTick', yticks_vals, 'YTickLabel', yticklabels_str, 'YTickLabelRotation', 0);
            title(['Channel ' EEG.chanlocs(ch).labels ' (number of neighbours: ' num2str(numel(neighbors)) ')']);
            xlabel('Time (ms)');
            ylabel('Channels');
        end

        sgtitle(fig_c2_neighbors_good, ['C2 check:Time series of a good channel (red) and its neighbours (blue) - Song ' num2str(song_idx)]);
        exportgraphics(fig_c2_neighbors_good, filename_pdf, 'Append', true, 'ContentType', 'image', 'Resolution', 300);
    end

    % Spectral Power C3
    if ~isempty(c3_bad_channels)
        set(0, 'CurrentFigure', fig_c3);
        set(fig_c3, 'Visible', 'off');
        subplot(2, 5, song_idx); hold on;
        bar(low_band, 'FaceColor', 'white'); hold on;
        scatter(c3_bad_channels, low_band(c3_bad_channels), 10, 'red', 'filled');
        title(['Song ' num2str(song_idx)]);
    end

   if ~isempty(c4_bad_channels) 
        % Spectral Power C4
        set(0, 'CurrentFigure', fig_c4);
        set(fig_c4, 'Visible', 'off');
        subplot(2, 5, song_idx); hold on;
        bar(high_band, 'FaceColor', 'white'); hold on;
        scatter(c4_bad_channels, high_band(c4_bad_channels), 10, 'red', 'filled');
        title(['Song ' num2str(song_idx)]);
   end
    

   
end

function export_report_to_pdf(sub_idx, fig_topos, fig_c1c2, fig_c3, fig_c4)
        
    filename_pdf = ['Subject_' num2str(sub_idx, '%03d') '_Quality.pdf']; 
    sgtitle(fig_topos, 'Removed Channels by each Criteria');
    sgtitle(fig_c1c2, 'C1 check: Flat (blue) and C2 check: Low Corr (red)');
    sgtitle(fig_c3, 'C3 check: Low Power Outlier (red dot)');
    sgtitle(fig_c4, 'C4 check: High Power Outlier (red dot)');

    % Export figures
    exportgraphics(fig_topos, filename_pdf, 'Append', true, 'ContentType', 'image', 'Resolution', 300);
    exportgraphics(fig_c1c2, filename_pdf, 'Append', true, 'ContentType', 'image', 'Resolution', 300);
    exportgraphics(fig_c3, filename_pdf, 'Append', true, 'ContentType', 'image', 'Resolution', 300);
    exportgraphics(fig_c4, filename_pdf, 'Append', true, 'ContentType', 'image', 'Resolution', 300);
    
    % Close all
    close(fig_topos); close(fig_c1c2); close(fig_c3); close(fig_c4);

end
