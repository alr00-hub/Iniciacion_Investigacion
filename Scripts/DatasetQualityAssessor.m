classdef DatasetQualityAssessor
    methods(Static)
        function run(ds_name, n_songs, n_subs, path_to_ds, sub_base_name, song_base_name, subplot_sizes, export)
            % Quality struct
            quality_data = struct();
            
            % Adding path for clean_rawdata module
            addpath(genpath('D:\eeglab\eeglab_current\eeglab2024.2\plugins\clean_rawdata2.10'));
            usable_songs_per_subject = zeros(n_subs, 1);
            low_band_bad_windows = zeros(n_subs, 1);
            high_band_bad_windows = zeros(n_subs, 1);
            total_windows = zeros(n_subs, 1);

            % Loop over subjects
            for sub_idx = 1:n_subs
                disp(['Processing subject ' num2str(sub_idx)]);
                % Init struct for subject
                quality_data(sub_idx).subject_id = sub_idx;
                
                usable_songs = 0;
                lb_bad_windows = 0;
                hb_bad_windows = 0;
                n_total_windows = 0;

                % Export title to pdf if flag
                if export
                    DatasetQualityAssessor.export_title_page(ds_name, sub_idx);
                    figs = DatasetQualityAssessor.init_figures();
                end

                % Loop over songs
                for song_idx = 1:n_songs
                    disp(['Processing song ' num2str(song_idx)]);
                    EEG = DatasetQualityAssessor.load_and_preprocess(path_to_ds, sub_base_name, song_base_name, sub_idx, song_idx);
                    %neighbor_matrix = DatasetQualityAssessor.compute_neighbors(EEG);
                    neighbor_matrix = NaN;
                    [c1_bad_channels, c2_bad_channels, corrs] = DatasetQualityAssessor.criteria_flat_and_corr(EEG, neighbor_matrix);
                    [c3_bad_channels, c4_bad_channels, c3_bad_windows, c4_bad_windows, low_band, high_band,] = DatasetQualityAssessor.criteria_power(EEG);
                    bad_channels = unique([c1_bad_channels(:); c2_bad_channels(:); c3_bad_channels(:); c4_bad_channels(:)]);
                    removed_percentage = 100 * numel(bad_channels) / size(EEG.data,1);

                    if removed_percentage <= 20
                        usable_songs = usable_songs + 1;
                    end
                    lb_bad_windows = lb_bad_windows + numel(c3_bad_windows);
                    hb_bad_windows = hb_bad_windows + numel(c4_bad_windows);
                    n_total_windows = n_total_windows + floor(size(EEG.data, 2)/EEG.srate);

                    % Save quality data
                    quality_data(sub_idx).songs(song_idx).song_id = song_idx;
                    quality_data(sub_idx).songs(song_idx).total_channels = size(EEG.data, 1);
                    quality_data(sub_idx).songs(song_idx).c1_bad_channels = c1_bad_channels;
                    quality_data(sub_idx).songs(song_idx).c2_bad_channels = c2_bad_channels;
                    quality_data(sub_idx).songs(song_idx).c3_bad_channels = c3_bad_channels;
                    quality_data(sub_idx).songs(song_idx).c4_bad_channels = c4_bad_channels;
                    quality_data(sub_idx).songs(song_idx).total_windows = floor(size(EEG.data, 2)/EEG.srate);
                    quality_data(sub_idx).songs(song_idx).c3_bad_windows = c3_bad_windows;
                    quality_data(sub_idx).songs(song_idx).c4_bad_windows = c4_bad_windows;

                    % Create figures if flag
                    if export
                        DatasetQualityAssessor.create_figures(ds_name, sub_idx, song_idx, EEG, removed_percentage, ...
                            corrs, neighbor_matrix, c1_bad_channels, c2_bad_channels, c3_bad_channels, c4_bad_channels, ...
                            bad_channels, low_band, high_band, c3_bad_windows, c4_bad_windows, subplot_sizes, figs);
                    end
                end

                usable_songs_per_subject(sub_idx) = usable_songs;
                low_band_bad_windows(sub_idx) = lb_bad_windows;
                high_band_bad_windows(sub_idx) = hb_bad_windows;
                total_windows(sub_idx) = n_total_windows;

                % Export previous figures to pdf if flag
                if export
                    DatasetQualityAssessor.export_report_to_pdf(ds_name, sub_idx, figs);
                end
            end

            DatasetQualityAssessor.report_quality(quality_data, ds_name ,usable_songs_per_subject, low_band_bad_windows, high_band_bad_windows, total_windows, n_songs, n_subs);
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

        function [neighbor_matrix] = compute_neighbors(EEG)
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
        end

       function [c1_bad_channels, c2_bad_channels, corrs] = criteria_flat_and_corr(EEG, ~)
            EEG_clean_flat = clean_flatlines(EEG); % Using default clean_flatlines
            c1_bad_channels = find(~ismember({EEG.chanlocs.labels}, {EEG_clean_flat.chanlocs.labels}));
            EEG.data = squeeze(mean(EEG.data, 3)); % If the ds has only 1 epoch it returns the original data
            %data = EEG.data;
            EEG_clean = clean_channels_nolocs(EEG);
            c2_bad_labels = setdiff({EEG.chanlocs.labels}, {EEG_clean.chanlocs.labels});
            c2_bad_channels = find(ismember({EEG.chanlocs.labels}, c2_bad_labels));
            corrs = NaN; % For now, we don't use this variable
            % ----------------------- C2 variant -----------------------
            % % Set variables
            % n_channels = size(data,1);
            % n_samples = size(data,2);
            % 
            % % Sliding window: size 20% of samples, step 50%
            % win_len = floor(0.2 * n_samples);
            % step = floor(0.5 * win_len);
            % n_windows = floor((n_samples - win_len) / step) + 1;
            % 
            % corrs = zeros(n_channels, 1);
            % for ch = 1:n_channels
            %     neighbors = find(neighbor_matrix(ch, :));
            %     corr_vals = zeros(n_windows, 1);
            %     for w = 1:n_windows
            %         idx = (w-1)*step + (1:win_len);
            %         seg_ch = data(ch, idx)';
            %         seg_neighbors = median(data(neighbors, idx), 1)';
            %         corr_vals(w) = corr(seg_ch, seg_neighbors);
            %     end
            % 
            %     corrs(ch) = median(corr_vals);
            % end
            % 
            % c2_bad_channels = find(corrs < 0.75); 
            % -- ------------------------------------------------------
        end

        function [c3_bad_channels, c4_bad_channels, c3_bad_windows, c4_bad_windows, low_band, high_band] = criteria_power(EEG)
        
            % Define number of windows and their length
            n_channels = size(EEG.data, 1);
            n_samples = size(EEG.data, 2);
            win_len = EEG.srate; % 1 second
            n_windows = floor(n_samples/win_len);

            low_band = zeros(n_windows, n_channels);
            high_band = zeros(n_windows, n_channels);
            for win = 1:n_windows
                idx = win_len * (win-1) + (1:win_len);
                % Sanity check. This should not happen for our ds
                if idx(end) > n_samples, break; end
                seg = EEG.data(:, idx);
                [pxx, f] = pwelch(seg', [], [], [], EEG.srate);
                pxx_db = 10 * log10(pxx); % to dB
                low_idx = f >= 1 & f <= 10;
                high_idx = f >= 65 & f <= 90;
                low_band(win, :) = mean(pxx_db(low_idx, :), 1);
                high_band(win, :) = mean(pxx_db(high_idx, :), 1);
            end
            
            % We compute median over time windows
            median_low_power = median(low_band, 1);
            median_high_power = median(high_band, 1);

            % Detect outliers based on IQR
            detect_outliers = @(x, w) find(x < quantile(x, 0.25) - w * iqr(x) | ...
                                            x > quantile(x, 0.75) + w * iqr(x));

            w = 3; % Whisker length
            chan_threshold = 0.1;
            win_threshold = 0.05;

            % We compute the outlier windows & channels, for individual channels only 
            [c3_win_idx, c3_chan_idx] = detect_outliers(low_band, w);
            [c4_win_idx, c4_chan_idx] = detect_outliers(high_band, w);

            % We compute the outliers for median power over windows, over all channels
            c3_bad_chan_median = detect_outliers(median_low_power, w);
            c4_bad_chan_median = detect_outliers(median_high_power, w);

            % Computing bad channels considering just the number of outlier
            % windows
            c3_win_outlier_per_chan = accumarray(c3_chan_idx, 1, [n_channels, 1]);
            c4_win_outlier_per_chan = accumarray(c4_chan_idx, 1, [n_channels, 1]);
            c3_bad_channels = find((c3_win_outlier_per_chan/n_windows) > chan_threshold);
            c4_bad_channels = find((c4_win_outlier_per_chan/n_windows) > chan_threshold);

            % Merge bad channels: those with many outlier windows (relative to themselves)
            % and those that are outliers in terms of median power across all channels
            c3_bad_channels = unique(union(c3_bad_channels, c3_bad_chan_median));
            c4_bad_channels = unique(union(c4_bad_channels, c4_bad_chan_median));

            % Identify bad windows: windows with a high proportion of outlier channels
            % (Note: we do not apply a median-based criterion across windows as we do for channels)
            c3_chan_outlier_per_win = accumarray(c3_win_idx, 1, [n_windows, 1]);
            c4_chan_outlier_per_win = accumarray(c4_win_idx, 1, [n_windows, 1]);
            c3_bad_windows = find((c3_chan_outlier_per_win/n_channels) > win_threshold);
            c4_bad_windows = find((c4_chan_outlier_per_win/n_channels) > win_threshold);

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
            f = @(tag) figure('Visible','on','Units','normalized','Position',[1 1 1 1], 'Tag', tag);
            figs.topos = f('topos');
            figs.c1c2 = f('c1c2');
            figs.c3_bad_channels = f('c3_bad_channels');
            figs.c4_bad_channels = f('c4_bad_channels');
        end

        function create_figures(~, ~, song_idx, EEG, removed_pct, ~, ~, ...
            c1_bad_channels, c2_bad_channels, c3_bad_channels, c4_bad_channels, ...
            bad, low_band, high_band, c3_bad_windows, c4_bad_windows, ...
            subplot_sizes, figs)

            % Create figures with subplots for each criteria
            %filename_pdf = [ds_name '_Subject_' num2str(sub_idx, '%03d') '_QualityReport.pdf'];

            % --- General Parameters ---
            n_channels = size(EEG.data,1);
            mid = floor(length(EEG.times)/2);
            srate = EEG.srate;
            win = round(srate); % 1 second
            idx = mid - floor(win/2) : mid + floor(win/2);
            normalize = @(x) (x - mean(x)) ./ std(x);

            % --- TOPO Plot (Bad Channels Overview) ---
            set(0, 'CurrentFigure', figs.topos);
            subplot(subplot_sizes(1), subplot_sizes(2), song_idx);

            mask = zeros(1, n_channels); 
            mask(bad) = 1;

            topoplot(mask, EEG.chanlocs, 'style', 'blank', 'electrodes', 'on');
            text(0, -0.55, sprintf('Song: %d - Total Removed: %.1f%%', song_idx, removed_pct), ...
                'Color', 'red', 'FontSize', 9, 'HorizontalAlignment', 'center');
            text(0, -0.75, sprintf('C1: %.1f%% || C2: %.1f%% || C3: %.1f%% || C4: %.1f%%', ...
                100*numel(c1_bad_channels)/n_channels, ...
                100*numel(c2_bad_channels)/n_channels, ...
                100*numel(c3_bad_channels)/n_channels, ...
                100*numel(c4_bad_channels)/n_channels), ...
                'Color', 'k', 'FontSize', 8, 'HorizontalAlignment', 'center');

            % --- C1 & C2 Plot (Normalized Time-Series) ---
            set(0, 'CurrentFigure', figs.c1c2);
            subplot(subplot_sizes(1), subplot_sizes(2), song_idx); 
            hold on;

            offset = 0; 
            ytick_vals = []; 
            ytick_labels = {};

            for ch = c1_bad_channels
                plot(EEG.times(idx), normalize(EEG.data(ch, idx)) + offset, 'b');
                ytick_vals(end+1) = offset;
                ytick_labels{end+1} = EEG.chanlocs(ch).labels;
                offset = offset + 3;
            end

            for ch = c2_bad_channels
                plot(EEG.times(idx), normalize(EEG.data(ch, idx)) + offset, 'r');
                ytick_vals(end+1) = offset;
                ytick_labels{end+1} = EEG.chanlocs(ch).labels;
                offset = offset + 3;
            end

            set(gca, 'YTick', ytick_vals, 'YTickLabel', ytick_labels);
            xlabel('Time (ms)');
            ylabel('Channels (Normalized)');
            title(['Song ' num2str(song_idx) ' - 1s centered']);

            % --- C3 Plot (Low Band PSD Boxplot) ---
            set(0, 'CurrentFigure', figs.c3_bad_channels);
            subplot(subplot_sizes(1), subplot_sizes(2), song_idx);

            
            boxplot(low_band, 'Widths', 0.9, 'OutlierSize', 2, 'Whisker', 3);
            ax = gca;
            ax.XTickLabelRotation = 90;
            ax.FontSize = ceil(4 * 64/size(EEG.data, 1));

            title(['Song ' num2str(song_idx) ' - Bad 1s windows (counting all channels): ' num2str(numel(c3_bad_windows))], 'FontSize', 9);
            xlabel('Channel', 'FontSize', 9);
            ylabel('PSD (dB)', 'FontSize', 9);

            h = findobj(gca, 'Tag', 'Box');
            numBoxes = length(h);
            cmap = turbo(numBoxes);

            for k = 1:numBoxes
                patch(get(h(k), 'XData'), get(h(k), 'YData'), ...
                    cmap(k, :), 'FaceAlpha', 0.6, 'EdgeColor', 'black');
            end

            if ~isempty(c3_bad_channels)
                new_labels = string(ax.XTickLabel); 
                new_labels(c3_bad_channels) = "(bad) " + new_labels(c3_bad_channels);
                ax.XTickLabel = cellstr(new_labels);
            end

            % --- C4 Plot (High Band PSD Boxplot) ---
            set(0, 'CurrentFigure', figs.c4_bad_channels);
            subplot(subplot_sizes(1), subplot_sizes(2), song_idx);

            boxplot(high_band, 'Widths', 0.9, 'OutlierSize', 2, 'Whisker', 3);
            ax = gca;
            ax.XTickLabelRotation = 90;
            ax.FontSize = ceil(4 * 64/size(EEG.data, 1));
            title(['Song ' num2str(song_idx) ' - Bad 1s windows (counting all channels): ' num2str(numel(c4_bad_windows))], 'FontSize', 9);
            xlabel('Channel', 'FontSize', 9);
            ylabel('PSD (dB)', 'FontSize', 9);

            h = findobj(gca, 'Tag', 'Box');
            numBoxes = length(h);
            cmap = turbo(numBoxes);

            for k = 1:numBoxes
                patch(get(h(k), 'XData'), get(h(k), 'YData'), ...
                    cmap(k, :), 'FaceAlpha', 0.6, 'EdgeColor', 'black');
            end

            if ~isempty(c4_bad_channels)
                new_labels = string(ax.XTickLabel); 
                new_labels(c4_bad_channels) = "(bad) " + new_labels(c4_bad_channels);
                ax.XTickLabel = cellstr(new_labels);
            end

            %DatasetQualityAssessor.plot_c2_neighbors(song_idx, EEG, data, corrs, neighbor_matrix, c2_bad_channels, filename_pdf);
        end


        function plot_c2_neighbors(song_idx, EEG, corrs, neighbor_matrix, c2_bad_channels, filename_pdf)
            % Set variables
            n_neighbors = sum(neighbor_matrix, 2); % Number of neighbors for each channel
            threshold = prctile(n_neighbors, 10);
            peripheral = find(n_neighbors <= threshold); % Take peripheral (i.e. low number of neighbors)
            bad = intersect(c2_bad_channels, peripheral); % Those peripheral who were labeled as bad by c2_bad_channels
            good = setdiff(peripheral, c2_bad_channels); % Those peripheral who were labeled as good by c2_bad_channels

            % Plot some of them (up to 6 timeseries)
            DatasetQualityAssessor.plot_neighbor_figures('bad', song_idx, EEG, corrs, neighbor_matrix, bad, filename_pdf);
            DatasetQualityAssessor.plot_neighbor_figures('good', song_idx, EEG, corrs, neighbor_matrix, good, filename_pdf);
        end

        function plot_neighbor_figures(type, song_idx, EEG, corrs, neighbor_matrix, channels, filename_pdf)
            n = min(6, length(channels));
            if n == 0, return; end
            fig = figure('Visible','off','Units','normalized','Position',[1 1 1 1]);
            mid = floor(length(EEG.times)/2);
            srate = EEG.srate;
            win = round(srate); % 1 second
            idx = mid - floor(win/2) : mid + floor(win/2);
            
            % Normalize window for better signal contrast in plot
            normalize = @(x) (x - mean(x)) ./ std(x); 

            for i = 1:n
                ch = channels(i);
                neighbors = find(neighbor_matrix(ch, :));
                subplot(n,1,i); hold on;
                offset = 0; ytick_vals = []; ytick_labels = {};
                for nb = neighbors
                    plot(EEG.times(idx), normalize(EEG.data(nb,idx)) + offset, 'b'); % Paint it blue
                    ytick_vals(end+1) = offset;
                    ytick_labels{end+1} = EEG.chanlocs(nb).labels;
                    offset = offset + 3; % Small offset because data is normalized
                end
                switch upper(type) 
                    case 'GOOD'
                        plot(EEG.times(idx), normalize(EEG.data(ch,idx)) + offset, 'g'); % Paint it green
                    case 'BAD'
                        plot(EEG.times(idx), normalize(EEG.data(ch,idx)) + offset, 'r'); % Paint it red
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
            sgtitle(figs.c3_bad_channels, 'C3 check: Low Power Outlier (1 Hz - 10 Hz)');
            sgtitle(figs.c4_bad_channels, 'C4 check: High Power Outlier (65 Hz - 90 Hz)');
            exportgraphics(figs.topos, filename_pdf, 'Append', true, 'ContentType', 'image', 'Resolution', 300);
            exportgraphics(figs.c1c2, filename_pdf, 'Append', true, 'ContentType', 'image', 'Resolution', 300);
            exportgraphics(figs.c3_bad_channels, filename_pdf, 'Append', true, 'ContentType', 'image', 'Resolution', 600);
            exportgraphics(figs.c4_bad_channels, filename_pdf, 'Append', true, 'ContentType', 'image', 'Resolution', 600);
            close(figs.topos); close(figs.c1c2); close(figs.c3_bad_channels); close(figs.c4_bad_channels);
        end

        function report_quality(quality_data, ds_name, usable_songs_per_subject, low_band_bad_windows, high_band_bad_windows, total_windows, n_songs, n_subs)
            % Export quality data to a .mat file
            save([ds_name '_quality_data2.mat'], 'quality_data');
            usable_songs_percentage = 100 * sum(usable_songs_per_subject) / (n_songs * n_subs);
            usable_subjects = sum(usable_songs_per_subject >= 0.75 * n_songs);
            usable_subjects_percentage = 100 * usable_subjects / n_subs;
            disp('--------------------------------');
            disp(['Percentage of usable songs: ', num2str(usable_songs_percentage, '%.2f'), '%']);
            disp(['Percentage of usable subjects: ', num2str(usable_subjects_percentage, '%.2f'), '%']);
            disp(['Total bad section (1 Hz - 10 Hz): ', num2str(sum(low_band_bad_windows), '%.2f'), ' seconds ~ ' num2str(100*sum(low_band_bad_windows)/sum(total_windows), '%.2f'), '%']);
            disp(['Total bad section (65 Hz - 90 Hz): ', num2str(sum(high_band_bad_windows), '%.2f'), ' seconds ~ ' num2str(100*sum(high_band_bad_windows)/sum(total_windows), '%.2f'), '%']);
            disp('--------------------------------');

        end

    end
end
