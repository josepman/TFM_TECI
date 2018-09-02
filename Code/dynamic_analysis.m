%% DYNAMIC ANALYSIS
% http://commdetect.weebly.com/ - manual for scripts

clear all
close all

ROI_timeseries = csvread('/Users/hose/Desktop/TFM_TECI/Code/ts_con_CSF.csv',1,1);

nNodes = size(ROI_timeseries,2);
w = 20;                       % window size
overlapping = 0.8;
l = size(ROI_timeseries,1);   % total length of timeseries

% A partir de las series temporales de cada nodo, genera las ventanas
% temporales de tamaño w y solapamiento temporal = overlapping
sliding_windows = gen_sliding_window(ROI_timeseries, w, overlapping);

% Create the connectivity matrix
metric_rule = 'corr';  %'corr', 'dist_corr', 'partial_corr', 'cross_corr', 'MI'
connectivity_matrix = gen_connectivity_matrix(sliding_windows, metric_rule);

% NULL-MODELS
for i=1:size(ROI_timeseries,2);
    tseries_null_model(:,i) = tseries_rand(ROI_timeseries(:,i));                          % Randomize the timeserie of each ROI (before to create matrix conn)
end
tseries_null_model_sw = gen_sliding_window(tseries_null_model, w, overlapping);
tseries_null_model_conn = gen_connectivity_matrix(tseries_null_model_sw, metric_rule);

for i=1:size(connectivity_matrix,3)
     nodal_null_model(:,:,i) = randmio_und_connected(connectivity_matrix(:,:,i), 50);
end

temporal_null_model_conn = temporal_rand(connectivity_matrix);                   % Randomize layers of the multilayer network

% gamma = 2;
% omega = 1;
% conn_mat = cell(1,size(connectivity_matrix,3));
% for i=1:size(connectivity_matrix,3)
%     conn_mat{1,i} = connectivity_matrix(:,:,i);
% end
% [B,twom] = multicat(conn_mat,gamma,omega);                                         % NG null-model for categorical multi-layer networks
% 
models = {connectivity_matrix,tseries_null_model_conn,temporal_null_model_conn };


%[G,Gb] = graphs_from_sliding_window(connectivity_matrix, thr); %G=weighted, Gb=binary
for j=1:length(models);
    connectivity_matrix = models{3};
    % Get the communities
    n = 100;
    final_partition = zeros(nNodes,size(connectivity_matrix,3));
    final_Q = zeros(1,size(connectivity_matrix,3));
    thr_nodal = cell(1,size(connectivity_matrix,3));
    quality = zeros(1,size(connectivity_matrix,3));
    consensus = zeros(nNodes,size(connectivity_matrix,3));
    consensus_simm = zeros(1,size(connectivity_matrix,3));
    agreement_nmi = zeros(1,size(connectivity_matrix,3));
    agreement_z = zeros(1,size(connectivity_matrix,3));
    SWP = zeros(1,size(connectivity_matrix,3));
    delta_c = zeros(1,size(connectivity_matrix,3));
    delta_L = zeros(1,size(connectivity_matrix,3));
    persistence_LMC = cell(1,size(connectivity_matrix,3));

    for i=1:size(connectivity_matrix,3)
        % final partition is created from consensus clustering. See
        % community_detection script
        [final_partition(:,i), Q, thr_nodal{1,i}, quality(i), consensus(:,i), consensus_simm(i)] = community_detection(connectivity_matrix(:,:,i), n, nNodes);

        %Similarity between both methods
        agreement_nmi(i) = normalized_mutual_information(consensus(:,i), final_partition(:,i), 'adjusted');
        agreement_z(i) = zrand(consensus(:,i), final_partition(:,i));

        [SWP(i), delta_C(i), delta_L(i)] = small_world_propensity(connectivity_matrix(:,:,i)); %'O' = onela, 'B' = barrat, 'bin' = binary

        persistence_LMC{1,i} = sig_lmc(final_partition(:,1), connectivity_matrix(:,:,1));
    end

    %%%%%
    % Improvement of multilayer partitions of categorical networks
    final_partition = postprocess_categorical_multilayer(final_partition,size(final_partition,2));

    % reorders nodes and layers to emphasize persistent structure in an unordered multilayer partition
    [final_partition_sorted,s1,s2]=sort_categorical(final_partition); 
    pers=categorical_persistence(final_partition_sorted)
    %%%%


    % agreement between windows:
    agreement_w_nmi = zeros(size(connectivity_matrix,3),size(connectivity_matrix,3));
    agreement_w_z = agreement_w_nmi;
    for i=1:size(connectivity_matrix,3)
        for j=1:size(connectivity_matrix,3)
            agreement_w_nmi(i,j) = normalized_mutual_information(final_partition(:,i), final_partition(:,j), 'adjusted');
            agreement_w_z(i,j) = zrand(final_partition(:,i), final_partition(:,j));
        end
    end

    % Coeff of changes of comm. per node
    flexib = flexibility(final_partition', 'cat');

    % Cohesion (if a node change, other also change):
    [Cij,node_cohesion,node_disjoint,node_flexibility,strength_cohesion,commChanges,commCohesion,commDisjoint,commIndex] = calc_node_cohesion(final_partition);
    % Promiscuity (fraction of all communities in which a node participates at
    % least once)
    P = promiscuity(final_partition');

    % Persistance
    persist = persistence(final_partition');

    % El resto de métricas se pueden extraer con la función bct_measures:
    thr = 0.5;    % para option='binary'
    option = 'weighted';
    
    G{i} = bct_measures(connectivity_matrix, option, thr);
end

