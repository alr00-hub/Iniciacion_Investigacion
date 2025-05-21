%% Script to assess quality of DRYAD dataset
ds_name = 'dsRIUMA';
n_songs = 16;
n_subs = 6;
path_to_ds = 'D:\master\Iniciacion_Investigacion\Datasets\dsMalaga\ConvertedData2';
sub_base_name = 'sub-';
song_base_name = 'song-';
subplot_sizes = [5 ,2];
export = false;
DatasetQualityAssessor.run(ds_name, n_songs, n_subs, path_to_ds, sub_base_name, song_base_name, subplot_sizes, export)