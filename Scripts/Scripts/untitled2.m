%% Script to assess quality of DRYAD dataset
clear; clc;
ds_name = 'testing';
n_songs = 2;
n_subs = 7;
path_to_ds = 'D:\master\Iniciacion_Investigacion\Datasets\ds003774\sourcedata';
sub_base_name = 'sub-';
song_base_name = 'song-';
export = true;
DatasetQualityAssessor.run(ds_name, n_songs, n_subs, path_to_ds, sub_base_name, song_base_name, export)