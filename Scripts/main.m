%% Script to run the DatasetQualityAssessor
% Change the parameters below to run the DatasetQualityAssessor for each dataset
ds_name = 'dsNMEDE';
n_songs = 2;
n_subs = 24;
path_to_ds = 'D:\master\GlobusPC\NMED-E';
sub_base_name = 'sub-';
song_base_name = 'song-';
subplot_sizes = [2 ,2];
export = true;
DatasetQualityAssessor.run(ds_name, n_songs, n_subs, path_to_ds, sub_base_name, song_base_name, subplot_sizes, export)