%% Script for dsRIUMA conversion

clear; clc;
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;  % Initialize EEGLAB

% Define dataset directory
data_path = 'D:\master\Iniciacion_Investigacion\Datasets\dsMalaga\Raw_Data\Raw_Data\';
output_path = 'D:\master\Iniciacion_Investigacion\Datasets\dsMalaga\ConvertedData';
mkdir(output_path);
% List of subjects
subjects = [9, 10, 13, 14, 15, 17];

n_songs = 16;

% Loop through each subject
for i = 1:length(subjects)
    subjID = sprintf('%04d', subjects(i));
    subject_str = sprintf('%03d', subjects(i));
    
    subject_folder = fullfile(output_path, ['sub-' subject_str]);
    if ~exist(subject_folder, 'dir')
        mkdir(subject_folder);
    end
    vhdr_file = ['S_' subjID '_E_0007_0001.vhdr'];
    
    % Load EEG data
    EEG = pop_loadbv(data_path, vhdr_file);
    % Add dataset to EEGLAB
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0, 'setname', ['sub-' subject_str], 'gui', 'off'); 
    
    for song_idx = 1:n_songs
        song_str = sprintf('%02d',song_idx);
        song_folder = fullfile(subject_folder, ['song-' song_str]);
        if ~exist(song_folder, 'dir')
            mkdir(song_folder);
        end
        event_name = sprintf('S %2d', song_idx);
        song_eeg = pop_rmdat( EEG, {event_name}, [0 30], 0 );
        pop_saveset(song_eeg, 'filename', ['sub-' subject_str '_song-' song_str '.set'], 'filepath', song_folder);
    end
    
end

disp('All subjects processed successfully!');