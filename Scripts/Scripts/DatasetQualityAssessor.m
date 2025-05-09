classdef DatasetQualityAssessor
    methods(Static)
        function run(ds_name, n_songs, n_subs, path_to_ds, sub_base_name, song_base_name, export)
            % Adding path for clean_rawdata module
            addpath(genpath('D:\eeglab\eeglab_current\eeglab2024.2\plugins\clean_rawdata2.10'));
            usable_songs_per_subject = zeros(n_subs, 1);
            
            % Loop over subjects
            for sub_idx = 1:n_subs
                disp(['Processing subject ' num2str(sub_idx)]);
                usable_songs = 0;

                % Export title to pdf if flag
                if export
                    DatasetQualityAssessor.export_title_page(ds_name, sub_idx);
                    figs = DatasetQualityAssessor.init_figures();
                end

                % Loop over songs
                for song_idx = 1:n_songs
                    disp(['Processing song ' num2str(song_idx)]);
                    EEG = DatasetQualityAssessor.load_and_preprocess(path_to_ds, sub_base_name, song_base_name, sub_idx, song_idx);
                    [neighbor_matrix, data] = DatasetQualityAssessor.compute_neighbors(EEG);
                    [c1, c2, corrs] = DatasetQualityAssessor.criteria_flat_and_corr(EEG, data, neighbor_matrix);
                    [c3, c4, low_band, high_band] = DatasetQualityAssessor.criteria_power(EEG, data);
                    bad_channels = unique([c1(:); c2(:); c3(:); c4(:)]);
                    removed_percentage = 100 * numel(bad_channels) / size(EEG.data,1);

                    if removed_percentage <= 20
                        usable_songs = usable_songs + 1;
                    end
                    
                    % Create figures if flag
                    if export
                        DatasetQualityAssessor.create_figures(ds_name, sub_idx, song_idx, EEG, data, removed_percentage, ...
                            corrs, neighbor_matrix, c1, c2, c3, c4, bad_channels, low_band, high_band, figs);
                    end
                end

                usable_songs_per_subject(sub_idx) = usable_songs;

                % Export previous figures to pdf if flag
                if export
                    DatasetQualityAssessor.export_report_to_pdf(ds_name, sub_idx, figs);
                end
            end

            DatasetQualityAssessor.report_quality(usable_songs_per_subject, n_songs, n_subs);
        end

        function EEG = load_and_preprocess(path_to_ds, sub_base, song_base, sub_idx, song_idx)
            % Set variables
            subject_str = [sub_base num2str(sub_idx, '%03d')];
            song_str = [song_base num2str(song_idx, '%02d')];
            filename = [subject_str '_' song_str '.set'];
            filepath = [path_to_ds '\' subject_str '\' song_str];
            EEG = pop_loadset('filename', filename, 'filepath', filepath);

            % Downsample to 256 
            if EEG.srate > 256
                EEG = pop_resample(EEG, 256);
            end
            EEG = pop_reref(EEG, []); % Average reference
            EEG = pop_eegfiltnew(EEG, 'locutoff', 0.4); % High-pass filter (0.4 Hz)
        end

        function [neighbor_matrix, data] = compute_neighbors(EEG)
            % Set variables
            n_channels = size(EEG.data, 1);
            chanlocs = [EEG.chanlocs.X; EEG.chanlocs.Y; EEG.chanlocs.Z]';
            distances = nan(n_channels, n_channels);

            % Create distance matrix
            for i = 1:n_channels
                for j = 1:n_channels
                    if i == j
                        dist = nan; % A channel is not its own neighbor
                    else
                        dist = sqrt((chanlocs(i, 1) - chanlocs(j, 1))^2 + ...
                                    (chanlocs(i, 2) - chanlocs(j, 2))^2 + ...
                                    (chanlocs(i, 3) - chanlocs(j, 3))^2);
                    end
                    distances(i, j) = dist;
                end
            end
            neighbor_threshold = median(max(distances, [], 2, 'omitnan'), 'omitnan') / 4; % Divide space into 4 regions
            neighbor_matrix = distances <= neighbor_threshold; 
            data = squeeze(mean(EEG.data, 3)); % If the ds has only 1 epoch it returns the original data
        end

        function [c1, c2, corrs] = criteria_flat_and_corr(EEG, data, neighbor_matrix)
            EEG_clean_flat = clean_flatlines(EEG); % Using default clean_flatlines
            c1 = find(~ismember({EEG.chanlocs.labels}, {EEG_clean_flat.chanlocs.labels}));

            % Set variables
            n_channels = size(data,1);
            corrs = zeros(n_channels,1);
            for ch = 1:n_channels
                neighbors = find(neighbor_matrix(ch, :));
                neighbor_median = median(data(neighbors, :), 1); % Median of all neighbors
                corrs(ch) = corr(data(ch, :)', neighbor_median'); % Correlation with that median
            end
            c2 = find(corrs < 0.75); % If corr < 0.75 then labeled as bad
        end

        function [c3, c4, low_band, high_band] = criteria_power(EEG, data)
            % Compute FFT using Welch's method
            [pxx, f] = pwelch(data', [], [], [], EEG.srate);
            low_idx = f >= 1 & f <= 10; % low_freq band threshold
            high_idx = f >= 65 & f <= 90; % high_freq band threshold
            low_band = mean(pxx(low_idx, :), 1);
            high_band = mean(pxx(high_idx, :), 1);
            detect_outliers = @(x, w) find(x < quantile(x, 0.25) - w * iqr(x) | x > quantile(x, 0.75) + w * iqr(x));
            w = 3;
            c3 = detect_outliers(low_band, w); % low freq outlier
            c4 = detect_outliers(high_band, w);% high freq outlier
        end

        function export_title_page(ds_name, sub_idx)
            % Export figures as png with high resolution to pdf
            filename_pdf = [ds_name '_Subject_' num2str(sub_idx, '%03d') '_QualityReport.pdf'];
            fig = figure('Visible','off','Units','normalized','Position',[0 0 1 1]);
            axis off;
            text(0.5, 0.6, ds_name, 'FontSize', 34, 'HorizontalAlignment', 'center');
            text(0.5, 0.4, sprintf('Quality Report\nSubject %03d', sub_idx), 'FontSize', 24, 'HorizontalAlignment', 'center');
            text(0.5, 0.2, datestr(now, 'mmmm dd, yyyy'), 'FontSize', 14, 'HorizontalAlignment', 'center');
            exportgraphics(fig, filename_pdf, 'Append', false, 'ContentType', 'image', 'Resolution', 300);
            close(fig);
        end

        function figs = init_figures()
            % Set figures
            f = @(tag) figure('Visible','off','Units','normalized','Position',[1 1 1 1], 'Tag', tag);
            figs.topos = f('topos');
            figs.c1c2 = f('c1c2');
            figs.c3 = f('c3');
            figs.c4 = f('c4');
        end

        function create_figures(ds_name, sub_idx, song_idx, EEG, data, removed_pct, corrs, neighbor_matrix, ...
                c1, c2, c3, c4, bad, low_band, high_band, figs)
            % Create figures with subplots for each criteria
            filename_pdf = [ds_name '_Subject_' num2str(sub_idx, '%03d') '_QualityReport.pdf'];
            n_channels = size(data,1);
            set(0, 'CurrentFigure', figs.topos);
            subplot(4, 3, song_idx);
            mask = zeros(1, n_channels); mask(bad) = 1;
            topoplot(mask, EEG.chanlocs, 'style', 'blank', 'electrodes', 'on');
            text(0, -0.55, sprintf('Song: %d - Total Removed: %.1f%%', song_idx, removed_pct), ...
                'Color', 'red', 'FontSize', 9, 'HorizontalAlignment', 'center');
            text(0, -0.75, sprintf('C1: %.1f%% || C2: %.1f%% || C3: %.1f%% || C4: %.1f%%', ...
                100*numel(c1)/n_channels, 100*numel(c2)/n_channels, 100*numel(c3)/n_channels, 100*numel(c4)/n_channels), ...
                'Color', 'k', 'FontSize', 8, 'HorizontalAlignment', 'center');

            set(0, 'CurrentFigure', figs.c1c2);
            subplot(12,1,song_idx); hold on;
            for ch = c1, plot(EEG.times, data(ch,:), 'b'); end
            for ch = c2, plot(EEG.times, data(ch,:), 'r'); end
            title(['Song ' num2str(song_idx)]);

            if ~isempty(c3)
                set(0, 'CurrentFigure', figs.c3);
                subplot(4,3,song_idx); bar(low_band, 'FaceColor', 'white'); hold on;
                scatter(c3, low_band(c3), 10, 'red', 'filled');
                title(['Song ' num2str(song_idx)]);
            end

            if ~isempty(c4)
                set(0, 'CurrentFigure', figs.c4);
                subplot(4,3,song_idx); bar(high_band, 'FaceColor', 'white'); hold on;
                scatter(c4, high_band(c4), 10, 'red', 'filled');
                title(['Song ' num2str(song_idx)]);
            end

            DatasetQualityAssessor.plot_c2_neighbors(song_idx, EEG, data, corrs, neighbor_matrix, c2, filename_pdf);
        end

        function plot_c2_neighbors(song_idx, EEG, data, corrs, neighbor_matrix, c2_bad_channels, filename_pdf)
            % Set variables
            n_neighbors = sum(neighbor_matrix, 2); % Number of neighbors for each channel
            threshold = prctile(n_neighbors, 10);
            peripheral = find(n_neighbors <= threshold); % Take peripheral (i.e. low number of neighbors)
            bad = intersect(c2_bad_channels, peripheral); % Those peripheral who were labeled as bad by c2
            good = setdiff(peripheral, c2_bad_channels); % Those peripheral who were labeled as good by c2

            % Plot some of them (up to 6 timeseries)
            DatasetQualityAssessor.plot_neighbor_figures('bad', song_idx, EEG, data, corrs, neighbor_matrix, bad, filename_pdf);
            DatasetQualityAssessor.plot_neighbor_figures('good', song_idx, EEG, data, corrs, neighbor_matrix, good, filename_pdf);
        end

        function plot_neighbor_figures(type, song_idx, EEG, data, corrs, neighbor_matrix, channels, filename_pdf)
            n = min(6, length(channels));
            if n == 0, return; end
            fig = figure('Visible','off','Units','normalized','Position',[1 1 1 1]);
            for i = 1:n
                ch = channels(i);
                neighbors = find(neighbor_matrix(ch, :));
                subplot(n,1,i); hold on;
                offset = 0; ytick_vals = []; ytick_labels = {};
                for nb = neighbors
                    plot(EEG.times, data(nb,:) + offset, 'b'); % Paint it blue
                    ytick_vals(end+1) = offset;
                    ytick_labels{end+1} = EEG.chanlocs(nb).labels;
                    offset = offset + 1000;
                end
                switch upper(type) 
                    case 'GOOD'
                        plot(EEG.times, data(ch,:) + offset, 'g'); % Paint it green
                    case 'BAD'
                        plot(EEG.times, data(ch,:) + offset, 'r'); % Paint it red
                end
                ytick_vals(end+1) = offset;
                ytick_labels{end+1} = EEG.chanlocs(ch).labels;
                set(gca, 'YTick', ytick_vals, 'YTickLabel', ytick_labels);
                title(['Channel ' EEG.chanlocs(ch).labels ' - Correlation with neighbors: ' num2str(corrs(ch)) ' (number of neighbours: ' num2str(numel(neighbors)) ')']);
                xlabel('Time (ms)');
                ylabel('Channels');
            end
            sgtitle(fig, ['C2 check: Some ' upper(type) ' peipheral channels vs Neighbors - Song ' num2str(song_idx)]);
            exportgraphics(fig, filename_pdf, 'Append', true, 'ContentType', 'image', 'Resolution', 300);
            close(fig);
        end

        function export_report_to_pdf(ds_name, sub_idx, figs)
            filename_pdf = [ds_name '_Subject_' num2str(sub_idx, '%03d') '_QualityReport.pdf'];
            sgtitle(figs.topos, 'Removed Channels by each Criteria');
            sgtitle(figs.c1c2, 'C1 check: Flat (blue) and C2 check: Low Corr (red)');
            sgtitle(figs.c3, 'C3 check: Low Power Outlier (1 Hz - 10 Hz)');
            sgtitle(figs.c4, 'C4 check: High Power Outlier (65 Hz - 90 Hz)');
            exportgraphics(figs.topos, filename_pdf, 'Append', true, 'ContentType', 'image', 'Resolution', 300);
            exportgraphics(figs.c1c2, filename_pdf, 'Append', true, 'ContentType', 'image', 'Resolution', 300);
            exportgraphics(figs.c3, filename_pdf, 'Append', true, 'ContentType', 'image', 'Resolution', 300);
            exportgraphics(figs.c4, filename_pdf, 'Append', true, 'ContentType', 'image', 'Resolution', 300);
            close(figs.topos); close(figs.c1c2); close(figs.c3); close(figs.c4);
        end

        function report_quality(usable_songs_per_subject, n_songs, n_subs)
            usable_songs_percentage = 100 * sum(usable_songs_per_subject) / (n_songs * n_subs);
            usable_subjects = sum(usable_songs_per_subject >= 0.75 * n_songs);
            usable_subjects_percentage = 100 * usable_subjects / n_subs;
            disp('--------------------------------');
            disp(['Percentage of usable songs: ', num2str(usable_songs_percentage, '%.2f'), '%']);
            disp(['Percentage of usable subjects: ', num2str(usable_subjects_percentage, '%.2f'), '%']);
            disp('--------------------------------');
        end
    end
end
