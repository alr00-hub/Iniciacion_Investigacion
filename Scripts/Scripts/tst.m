clear; clc;
n_songs = 12; % Número de canciones
n_subs = 20;   % Número de sujetos
path_to_ds = 'D:\master\Iniciacion_Investigacion\Datasets\dsDRYAD\diliBach_4dryad_CND\diliBach_4dryad_CND';
% Root folder (change if needed)
rootDir = path_to_ds;

% Get list of all subfolders matching sub-XXX format
subFolders = dir(fullfile(rootDir, 'sub-*'));

for i = 1:length(subFolders)
    if subFolders(i).isdir
        eegPath = fullfile(rootDir, subFolders(i).name, 'eeg');
        if exist(eegPath, 'dir')
            % Get all song-XX directories inside the eeg folder
            songFolders = dir(fullfile(eegPath, 'song-*'));
            for j = 1:length(songFolders)
                songFolderPath = fullfile(eegPath, songFolders(j).name);
                if songFolders(j).isdir
                    fprintf('Removing folder: %s\n', songFolderPath);
                    rmdir(songFolderPath, 's');  % 's' removes folder contents
                end
            end
        end
    end
end
