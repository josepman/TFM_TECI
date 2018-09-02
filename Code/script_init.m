
for i=1:30;
    if(i < 10);
        subject = 'CorrMatrix_sample_0%d.mat';
    else
        subject = 'CorrMatrix_sample_%d.mat';
    end;
    roi_subject = sprintf(subject,i);
    subject_xx = load(roi_subject);
    formatSpec = 'Pval_Matix_subject_%d.csv';
    str = sprintf(formatSpec,i);
    csvwrite(str,subject_xx.PValMatrix);
end



prueba = arrayToContactSeq(dcor_mat,0);
prueba2 = arrayToContactSeq(dcor_matB(:,:,1),0);
y = betweennessCentrality(prueba,0);  % one for matrix
networksFromContacts(prueba,0)
correlationnet = networksFromContacts(prueba,0);
sm = temporalSmallWorldness(prueba,0,nNodes);
reachability = reachabilityGraph(prueba,0, 1:(size(ROISignals,1)/size_w)-1);


for i=1:size(slid_wind,3);
   figure(10),
   subplot(3,2,i), heatmap(correlationnet(:,:,i));
   %figure(2),
   %subplot(5,1,i), clustergram(correlationnet(:,:,i));
end


%% GRAPHS
Results = struct;
Results.G = {};  %G is a cell containing all the graphs. rows = subjects, cols = windows
Results.Gb = {};   % The same for binary connectivity graphs
Results.G_degree = {};
Results.Gb_degree = {};

Results.params = {};
%Results.params{ROIs} = 
n_subjects = 30;
w = 20;        % Window size
overlapping = 0.5;
metric_rule = 'dist_corr';
thr = 0.5;

for n=1:n_subjects;
    if(n < 10)
        path_file = 'C:\\Users\\gangsta\\Desktop\\Signals\\ROISignals_sample_0%d.mat';
        roi_subject = sprintf(path_file,n)
    else
        path_file = 'C:\\Users\\gangsta\\Desktop\\Signals\\ROISignals_sample_%d.mat';
        roi_subject = sprintf(path_file,n)
    end;
    ROI_timeseries = load(roi_subject);
    sliding_windows = gen_sliding_window(ROI_timeseries.ROISignals, w, overlapping);
    conn_mat = gen_connectivity_matrix(sliding_windows, metric_rule);
    
    Results.ROI_ts{n} = ROI_timeseries;
    Results.sliding_window{n} = sliding_windows;
    Results.connectivity_mat{n} = conn_mat;
    %[Results.G,Results.Gb] = graphs_from_sliding_window(conn_mat);
    %Empieza a recoger metricas
    %Exporta plots y tablas con los resultados
end


% Once created the adjacency matrix, we can create the graph
G = graph(BMatrix);
plot(G,'XData',BrainRegions.x,'YData',BrainRegions.y,'ZData',BrainRegions.z)

% or a weigthed graph
%Plot the graph, labeling the edges with their weights, and making the width of the edges proportional to their weights. Use a rescaled version of the edge weights to determine the width of each edge, such that the widest line has a width of 5.
G_w = graph(CorrMatrix);
LWidths = 5*G_w.Edges.Weight/max(G_w.Edges.Weight);
plot(G_w,'XData',BrainRegions.x,'YData',BrainRegions.y,'ZData',BrainRegions.z,'EdgeLabel',G_w.Edges.Weight,'LineWidth',LWidths)

CorrMatrix2 = CorrMatrix;
CorrMatrix2(abs(CorrMatrix2)<0.5) = 0;
G_w = graph(CorrMatrix2);
LWidths = 5*G_w.Edges.Weight/max(G_w.Edges.Weight);
plot(G_w,'XData',BrainRegions.x,'YData',BrainRegions.y,'ZData',BrainRegions.z,'EdgeLabel',G_w.Edges.Weight,'LineWidth',LWidths)

G = G_w;
deg = degree(G);
neigh = neighbors(G);
nears = nearest(G,nodeX,distance);  %if weighted, weight are used as distances

