%% Script to assess quality of DRYAD dataset
clear; clc;
ds_name = 'dsAUTh';
n_songs = 30;
n_subs = 5;
path_to_ds = 'D:\master\Iniciacion_Investigacion\Datasets\dsARISTOTLE';
sub_base_name = 'sub-';
song_base_name = 'ses-';
export = true;
DatasetQualityAssessor.run(ds_name, n_songs, n_subs, path_to_ds, sub_base_name, song_base_name, export)