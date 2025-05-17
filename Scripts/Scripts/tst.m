clear;
n_chan = 64; n_win = 10;
matrix = zeros(n_win, n_chan);

for ch = 1:n_chan 
    win_sig = linspace(10, 100, 10) + ch;
    matrix(:, ch) = win_sig;
end
w = 3;
matrix(10, 3) = 1e3;
matrix(9, 6) = 2e3;
matrix(5, 63) = 3e3;
matrix(4, 63) = 4e3;
iqr_pinga = iqr(matrix);
q75 = quantile(matrix, 0.75);
superior_outlier_matrix = q75 + w * iqr_pinga;


detect_outliers = @(x, w) find(x < quantile(x, 0.25) - w * iqr(x) | x > quantile(x, 0.75) + w * iqr(x));

[win_idx, chan_idx] = detect_outliers(matrix, w);

outlier_win_per_chan = accumarray(chan_idx, 1, [n_chan, 1])';
total_bad_windows = sum(outlier_win_per_chan);
c3 = find((outlier_win_per_chan/n_win) >= 0.1);

[idx_win, idx_chan] = find(matrix < 34);
n_found_chan = accumarray(idx_win, 1, [n_win, 1]);
n_found_chan