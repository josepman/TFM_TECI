clear all
close all
seed = 123;

% Crea una estructura de parámetros.
sP = simtb_create_sP;

sP.out_path = '/Users/hose/Desktop/TFM_TECI/simulated_data/';
sP.prefix = 'task'; 
sP.saveNII_FLAG = 0;
sP.verbose_display = 0;
sP.M = 10;   % 10 sujetos
sP.nC = 30;  % 30 componentes cada uno
sP.nT = 300;   % Duracion del experimento
sP.nV = 200;   % número de voxels

%% BLOCK DESIGN - TAREAS
% Ajusta los parámetros de los bloques de tareas
sP.TC_block_n = 2;       % Experimento simple con 2 condiciones (ON-OFF y derecha-izquierda)
sP.TC_block_length = 20; % ON durante 20 TRs
sP.TC_block_ISI = 30;    % OFF durante 30 TRs
sP.TC_block_amp = zeros(sP.nC, sP.TC_block_n);    
sP.TC_block_same_FLAG = 1; % 1 = block structure same for all subjects, 0 = order randomized
sP.TC_event_n = 0;      % No events

%     Right-hand
% Equivalen a las componentes 19, 23 y 27
sP.TC_block_amp(19,1) = 1;      % Fuertemente correlacionadas con la tarea 1
sP.C_block_amp(19,2) = -1;      % Pero todo lo contrario con la tarea 2
sP.TC_block_amp(23,1) = 2;      % Distintos niveles de relacion
sP.TC_block_amp(23,2) = -2;      
sP.TC_block_amp(27,1) = 0.5;
sP.TC_block_amp(27,2) = -0.5;
%     Left-hand
% Equivalen a las componentes 20, 22 y 28
sP.TC_block_amp(20,1) = -1;    % Correlacionadas negativamente con la tarea 
sP.TC_block_amp(20,2) = 1;    % Correlacionadas negativamente con la tarea 
sP.TC_block_amp(22,1) = -2;
sP.TC_block_amp(22,2) = 2;
sP.TC_block_amp(28,1) = -0.75;     
sP.TC_block_amp(28,2) = 0.75;  

% El resto de componentes las dejamos sin relación aparente


%% UNIQUE EVENTS (NO)
sP.TC_unique_FLAG = 0; % 1 = include unique events
%TC_unique_prob = 0.2*ones(1,nC); % [1 x nC] prob of unique event at each TR
%TC_unique_amp = ones(M,nC); % [M x nC] matrix of amplitude of unique events
% [...]
% SPATIAL COMPONENTS AND SPATIAL VARIABILITY
sP.SM_source_ID = [1:sP.nC];

%% LABEL COMPONENTS
% Here, we label components or component groups that may be used later
% Frontal: 1 second temporal delay from bilateral frontal
comp_TC = find(sP.SM_source_ID == 1);
comp_T1 = find(sP.SM_source_ID == 2);
comp_T2 = find(sP.SM_source_ID == 3);

% Frontal: 1 second temporal delay from bilateral frontal
comp_F1 = find(sP.SM_source_ID == 4);
comp_F2 = find(sP.SM_source_ID == 5);

% Medial Frontal: has lower baseline intensity (signal dropout)
comp_MF = find(sP.SM_source_ID == 6);

% Precuneus: activation only to targets
comp_P = find(sP.SM_source_ID == 7);

% DMN: negative activation to task events
comp_DMN = find(sP.SM_source_ID == 8);

comp_9 = find(sP.SM_source_ID == 9);
comp_10 = find(sP.SM_source_ID == 10);
comp_11 = find(sP.SM_source_ID == 11);
comp_12 = find(sP.SM_source_ID == 12);
comp_13 = find(sP.SM_source_ID == 13);

% CSF and white matter: unaffected by task, but has signal amplitude differences
comp_CSF1 = find(sP.SM_source_ID == 14);
comp_CSF2 = find(sP.SM_source_ID == 15);
comp_WM1 = find(sP.SM_source_ID == 16);
comp_WM2 = find(sP.SM_source_ID == 17);

% Dorsal Attention Network: activation to novels more than targets
comp_DAN = find(sP.SM_source_ID == 18);

% (Sensori)Motor: activation to targets and novels 
comp_M1_L = find(sP.SM_source_ID == 19);
comp_M1_R = find(sP.SM_source_ID == 20);
comp_M2_L = find(sP.SM_source_ID == 22);
comp_M2_R = find(sP.SM_source_ID == 23);

% Bilateral frontal: positive activation to for targets and novels
comp_BF = find(sP.SM_source_ID == 24);

comp_21 = find(sP.SM_source_ID == 21);
comp_25 = find(sP.SM_source_ID == 25);
comp_26 = find(sP.SM_source_ID == 26);

% Auditory:
comp_AUD1 = find(sP.SM_source_ID == 27);
comp_AUD2 = find(sP.SM_source_ID == 28);

% Hippocampus: activation only to novels
comp_H1 = find(sP.SM_source_ID == 29);
comp_H2 = find(sP.SM_source_ID == 30);

% compile list of all defined components of interest
sP.complist = [comp_AUD1, comp_AUD2, comp_DMN, comp_BF, comp_F1, comp_F2 ...
    comp_P, comp_DAN, comp_H1, comp_H2, comp_TC, comp_T1, comp_T2...
    comp_CSF1, comp_CSF2, comp_WM1, comp_WM2, comp_MF ...
    comp_M1_L, comp_M1_R, comp_M2_L, comp_M2_R];

sP.SM_present(:,:) = ones(sP.M, sP.nC);

%% VARIABILIDAD
% Y ahora introduce la variabilidad para cada componente de cada sujeto.
% Por simplicidad, la variabilidad estará espacialmente lineada para cada sujeto, aunque 
% una variabilidad con distribución normal representarían unos datos más reales.

% Source 19 and 20: offset in y position, linearly spaced between -5 to +5 voxels
sP.compIND = find(sP.SM_source_ID == 19); % get the component index for the source
sP.SM_translate_y(:, sP.compIND) = linspace(-5, 5, sP.M)';
sP.compIND = find(sP.SM_source_ID == 20); % get the component index for the source
sP.SM_translate_y(:, sP.compIND) = linspace(-5, 5, sP.M)';

% Source 22: offset in x position, linearly spaced between -6 to +6 voxels
sP.compIND = find(sP.SM_source_ID == 22);
sP.SM_translate_x(:, sP.compIND) = linspace(-6, 6, sP.M)';
% Source 23: offset in x position, linearly spaced between -6 to +6 voxels
sP.compIND = find(sP.SM_source_ID == 23);
sP.SM_translate_x(:, sP.compIND) = linspace(-6, 6, sP.M)';

% En este último variamos la dispersión en lugar de la localizacion
sP.compIND = find(sP.SM_source_ID == 27);
sP.SM_spread(:, sP.compIND) = linspace(0.5, 2.5, sP.M)';
% Source 28: rotation, linearly spaced between -30 to +30 degrees
compIND = find(sP.SM_source_ID == 28);
sP.SM_spread(:, sP.compIND) = linspace(0.5, 2.5, sP.M)';

% Check que los parámetros se han actualizado correctamente
%simtb_figure_params(sP,{'SM_translate_y','SM_translate_x','SM_theta','SM_spread'});


%% EVENTS
sP.TC_event_n = 4;  % spikes in CSF no related to the task
sP.TC_event_same_FLAG = 0;  %different in each subject
% event probabilities (0.7 standards, 0.075 targets and novels, 0.05 CSF spikes)
sP.TC_event_prob = [0.7, 0.075, 0.075, 0.05]; 
 % initialize [nC x TC_event_n] matrix of task-modulation amplitudes
sP.TC_event_amp = zeros(sP.nC, sP.TC_event_n);

 % event type 1: standard tone
sP.TC_event_amp([comp_M1_L comp_M1_R comp_M2_L comp_M2_R], 1) = 1.5; % very strong
sP.TC_event_amp([comp_AUD1 comp_F1 comp_AUD2], 1) = 0.9; % moderate task-modulation
sP.TC_event_amp([comp_BF comp_F2 comp_DAN], 1) = 0.6; % mild
sP.TC_event_amp([comp_DMN comp_10], 1) =-0.5; % negative weak

 % event type 2: target tone
sP.TC_event_amp([comp_AUD1 comp_AUD2], 2) = 1.2; % strong
sP.TC_event_amp([comp_BF comp_F1 comp_F2 comp_21], 2) = 1.0; % moderate
sP.TC_event_amp([comp_DAN], 2) = 0.8; % mild
sP.TC_event_amp([comp_P comp_MF comp_9 comp_T2 comp_11], 2) = 0.5; % weak
sP.TC_event_amp([comp_M1_L comp_M1_R comp_M2_L comp_M2_R], 2) = 1.6; % very_strong
sP.TC_event_amp([comp_DMN], 2) =-0.3; % negative weak

% event type 3: novel tone
sP.TC_event_amp([comp_AUD1 comp_AUD2 ], 3) = 1.0; % moderate
sP.TC_event_amp([comp_BF comp_F1 comp_13 comp_T1], 3) = 0.5; % weak
sP.TC_event_amp([comp_DAN comp_25 comp_26], 3) = 1.2; % strong
sP.TC_event_amp([comp_H1 comp_F2], 3) = 0.7; % mild
sP.TC_event_amp([comp_M1_L comp_M1_R comp_M2_L comp_M2_R], 3) = 1.5; % very strong
sP.TC_event_amp([comp_DMN comp_12 comp_H2 comp_TC], 3) = -0.3; % negative weak

% event type 4: 'spikes' in CSF (not related to task)
sP.TC_event_amp([comp_CSF1 comp_CSF2 comp_WM1 comp_WM2], 4) = 1.0; % moderate


%% MODELO HEMODINÁMICO PARA GENERAR LOS TC
% Ponemos el mismo para todos (convolucion con HRF), exceptos para los tejidos CSF y WM
sP.TC_source_type = ones(1,sP.nC); % convolution with HRF for most components
sP.TC_source_type([comp_CSF1, comp_CSF2]) = 3; % spike model for CSF

sP.TC_source_params = cell(sP.M, sP.nC); % initialize the cell structure
sP.P(1) = 5; % delay of response (relative to onset)
sP.P(2) = 13; % delay of undershoot (relative to onset)
sP.P(3) = 1; % dispersion of response
sP.P(4) = 1; % dispersion of undershoot
sP.P(5) = 6; % ratio of response to undershoot
sP.P(6) = 0; % onset (seconds)
sP.P(7) = 32; % length of kernel (seconds)
[sP.TC_source_params{:}] = deal(sP.P);

%% DATA BASELINE Y TIPOS DE TEJIDO
sP.D_baseline = 800*ones(1, sP.M); % [1 x M] vector of baseline signal intensity
% FLAG to include different tissue types (distinct baselines in the data)
sP.D_TT_FLAG = 1; % if 0, baseline intensity is constant
sP.D_TT_level = [1.15, 0.8, 1, 1.2]; % TT fractional intensities

%% RUIDO
sP.D_noise_FLAG = 1;    % Para añadir rician noise a los datos
minCNR = 1;
maxCNR = 2;
sP.D_CNR = rand(1, sP.M)*(maxCNR-minCNR) + minCNR;  %CNR distribuido uniformemente entre [1-2] entre los sujetos

%% (SIN) MOVIMIENTO 
sP.D_motion_FLAG = 0; % 1=motion, 0=no motion
%D_motion_TRANSmax = 0.02; % max translation, proportion of entire image
%D_motion_ROTmax = 5; % max rotation, in degrees
%D_motion_deviates = ones(M,3); % proportion of max each subject moves
%D_motion_deviates(1,:) = 0.5; % Subject 1 moves half as much

%% % CAMBIO SEÑAL PICO A PICO
% con el objetivo de que se parezcan más a datos reales
sP.D_pSC = 3 + 0.25*randn(sP.M, sP.nC); % [M x nC] matrix of percent signal changes
sP.D_pSC(:,comp_CSF1) = 1.2*sP.D_pSC(:,comp_CSF1);
sP.D_pSC(:,comp_CSF2) = 1.2*sP.D_pSC(:,comp_CSF2);
sP.D_pSC(:,comp_WM1) = 0.5*sP.D_pSC(:,comp_WM1);
sP.D_pSC(:,comp_WM2) = 0.5*sP.D_pSC(:,comp_WM2);

simtb_main(sP);