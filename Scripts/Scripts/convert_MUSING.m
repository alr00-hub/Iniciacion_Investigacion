%% Script for dsDRYAD - Aggregate 3 Trials per Song as Epochs with Events

%% Basic Setup
clear; clc;
n_songs = 12; % Número de canciones
n_subs = 20;   % Número de sujetos
path_to_ds = 'D:\master\Iniciacion_Investigacion\Datasets\ds003774\sourcedata';
sub_base = 'sub-';

%% Process each subject
for sub_idx = 1:n_subs
    
    display(['Processing subject ' num2str(sub_idx)]);
    % Convertir el número del sujeto a una cadena
    subject_str = [sub_base sprintf('%03d', sub_idx)];
    
    % Cargar datos del sujeto
    filename = [subject_str '_task-ListeningandResponse_eeg.set'];
    filepath = fullfile(path_to_ds, subject_str, 'eeg');
    EEG = pop_loadset('filename', filename, 'filepath', filepath);

    %% Indices de los eventos
    stm_indices = find(strcmp({EEG.event.type}, 'stm+'));  % Indices de 'stm+'
    fxnd_indices = find(strcmp({EEG.event.type}, 'fxnd')); % Indices de 'fxnd'

    %% Procesar cada canción
    for song_idx = 1:n_songs
        % Verificar que hay suficientes eventos
        if song_idx > length(stm_indices) || song_idx > length(fxnd_indices)
            warning('No suficientes eventos para la canción %d en el sujeto %s. Omitiendo.', song_idx, subject_str);
            continue;  % Si no hay suficientes eventos, saltar esta canción
        end
        
        % Obtener el índice de inicio (stm+) y fin (fxnd) de la canción
        song_start_idx = stm_indices(song_idx);
        song_end_idx = fxnd_indices(song_idx);
        
        % Obtener las latencias de los eventos
        start_latency = EEG.event(song_start_idx).latency;
        end_latency = EEG.event(song_end_idx).latency;

        % Definir el intervalo para la extracción de datos
        start_time = (start_latency - EEG.event(1).latency) / EEG.srate;  % Convertir a segundos
        end_time = (end_latency - EEG.event(1).latency) / EEG.srate;      % Convertir a segundos

        % Extraer los datos entre los eventos stm+ y fxnd para la canción
        song_EEG = pop_select(EEG, 'time', [start_time, end_time]);  % Seleccionar el segmento de tiempo

        % Asegurarse de que la estructura EEG sea válida
        song_EEG = eeg_checkset(song_EEG);
        
        converted_filepath = fullfile(path_to_ds, subject_str);
        % Crear la carpeta para guardar el archivo
        song_folder = fullfile(converted_filepath, ['song-' num2str(song_idx, '%02d')]);
        if ~exist(song_folder, 'dir')
            mkdir(song_folder);
        end
        
        % Guardar el archivo .set correspondiente para esta canción
        output_filename = [subject_str '_song-' num2str(song_idx, '%02d') '.set'];
        pop_saveset(song_EEG, 'filename', output_filename, 'filepath', song_folder);
        
        display(['Saved: ' fullfile(song_folder, output_filename)]);
    end

end
