# Datasets M/EEG

# 1. (EEG) EEG Database for musical genres detection

Al igual que el anterior, su fecha de publicación es muy reciente. Consta de 6 sujetos, que escuchan 16 canciones de
varios generos musicales. Las porciones escuchadas de las canciones se escogieron de manera arbitraria, teniendo en
cuenta la parte más "pegadiza" de cada una. También cuenta con valoraciones subjetivas como el grado de familiardad y
disfrute por canción.

## Info. adicional

Fecha: 24/01/2025;  
Peso total: 4.5GB;  
Los datos son guardados en formato .eeg;  
Pocos sujetos;  
Enlace: https://riuma.uma.es/xmlui/handle/10630/36947

# 2. (EEG) Music Listening- Genre EEG dataset

Este es el que se uso en el TFG. Cuenta con 20 sujetos que han estado escuchando 12 canciones de distinto género y han dado 
valoraciones de familiaridad y disfrute con notas del 1 al 5. Lo usaremos porque también sirve como una métrica de comparación adicional
con los otros.

## Info. adicional

Fecha: 23/08/2021;  
Peso total: 10.12GB;
Los datos están en formato .set;  
La calidad no es muy buena;  
Enlace: https://openneuro.org/datasets/ds003774/versions/1.0.2;  

# 3. (EEG) Naturalistic Music EEG Dataset - Elgar (NMED-E)

El dataset contiene grabaciones EEG de 24 sujetos durante sesiones de escucha de música junto con una serie de valoraciones
subjetivas evaluadas del 1 al 9. Entre las valoraciones estan familiardad, "arousal" y nivel de interés. El principal problema
que he visto es que la velocidad de descarga es baja y no he encontrado una manera facil de descargar todos los datos juntos. Hay varios
datasets parecidos a estos de la misma fuente (Universidad Stanford), algunos no tan recientes: NMED-M (2021), NMED-RP (2018), 
NMED-T (2017), NMED-H (2014 - 2016).

## Info. adicional

Fecha: 18/04/2021 (Creado) - 13/05/2024 (Modificado);  
Peso total: 20GB aprox.;  
Los datos son guardados en formato .mat;  
Cuenta con documentación del formato de datos más detallada;  
Enlace: https://purl.stanford.edu/pp371jh5722  

# 4. (EEG) Naturalistic Music EEG Dataset - Minimalism (NMED-M)

Se menciona anteriormente, viene también de la Universidad de Stanford y tiene datos similares al anterior. En este caso son 30 sujetos,
escuchando otras 5 piezas distintas. En una escucha posterior, se pidió a los participantes que evaluen las piezas en base a diferentes
factores como difrute, interés, etc. En este caso sí que hay un fichero comprimido .zip de todos los datos raw medidos.

## Info. adicional

Fecha: 11/05/2021 (Creado) - 05/12/2022 (Modificado);  
Peso total: 34GB;  
Los datos son guardados en formato .mat;  
Cuenta con documentación del formtado de datos más detallada;  
Enlace: https://purl.stanford.edu/kt396gb0630  

# 5. (EEG) On the estimate of music appraisal from surface EEG: a dynamic-network approach based on cross-sensor PAC measurements

En este caso se publican los datos en un repositorio de GitHub. Contienen grabaciones EEG de 20 sujetos durante la escucha de música de 
30 canciones distintas, 80 segundos cada una. En este caso, cada participante tuvo una playlist personalizada, siguiendo una
distribución uniforme en cuanto a las puntuaciones. Los 80 segundos de cada canción se eligieron específicamente, centrándose sobre todo
en el estribillo y la parte más llamativa.  

## Info. adicional

Fecha: 02/06/2021;  
Peso total: 700MB;  
Los datos son guardados en formato .mat;  
Al parecer la tasa de muestreo es bastante baja, 125 (Hz)  
Enlace: https://github.com/AdamosDA/Music-EEG?tab=readme-ov-file

# 6. (M/EEG) Tracking predictions in naturalistic music listening using MEG and computational models of music

Este es de los pocos dataset reciente relacionado con la escucha de música que contiene datos MEG. Cuenta con 35 sujetos y un total de 19 
canciones con duración media de 2 minutos. Los sujetos estaban sentados y podían elegir cuando empezará la canción a través de un botón. 
También se les presento un punto fijo en la pantalla para evitar movimientos de la cabeza u oculares. 

## Info. adicional

Fecha: 08/06/2022 (Creado) - 10/01/2023 (Modificado);  
Peso total: 333.3GB;  
Los datos son guardados en formato .ds (CTF);  
Solo para los datos MEG, cada sujeto ocupa casi 5GB;  
Enlace: https://data.ru.nl/collections/di/dccn/DSC_3018045.02_116  

# 7. (EEG) Cortical encoding of melodic expectations in human temporal cortex

Un total de 20 sujetos, mitad son músicos, escuchando piezas de piano de Bach. 

## Info. adicional

Fecha: 22/01/2020;  
Peso total: 5GB;  
Los datos son guardados en formato .mat;  
Enlace: https://datadryad.org/dataset/doi:10.5061/dryad.g1jwstqmh#citations;  


# Otros enlaces:
https://openneuro.org/datasets/ds005403/versions/1.0.1    
https://openneuro.org/datasets/ds002893/versions/2.0.0  
https://openneuro.org/datasets/ds002721/versions/1.0.3  
