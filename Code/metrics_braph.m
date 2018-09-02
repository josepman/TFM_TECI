graph_bu = GraphBU(CorrMatrix,'threshold',0.1);   
graph_wu = GraphWU(abs(CorrMatrix));
graph = graph_wu;
% A partir de aqui, todo es igual
N = graph.nodenumber();
between = graph.betweenness();   % one per node
closen = graph.closeness();      % one per node

[clus_graph_avg,clus_per_node] = graph.cluster();
deg = graph.degree();                         % todos conectados con todos
avdg =  graph.measure(Graph.DEGREEAV);        % average value
graph.distance();                             % is a matrix
ecc = graph.eccentricity();    % one per node
rad_min = min(ecc);            % 1 value
diam_max = max(ecc);            % 1 value
rad = graph.radius();            % 1 value
diam = graph.diameter();            % 1 value

global_eff = graph.geff();      % one per node
[local_eff_nodes, local_eff_graph] = graph.leff();   %loc_eff_graph = one per node, loc_eff_nodes = 1 value
particip = graph.participation();   % WU must be positive.  one per node
pth_length = graph.pl();                % path length.      one per node
carac_pth = graph.measure(Graph.CPL);   % characteristic path length = 1 value
[community_structure,modularity_graph] = graph.structure();   % only for positive values. comm=one per node, mod = 1 value
transit = graph.transitivity();     % one value
tr = graph.triangles();             % one per node
test_zscore = graph.zscore();       % one per node

% Strength only for weighted graphs
str = graph.strength();                     % one per node     = sum(CorrMatrix)
avstrng = graph.measure(Graph.STRENGTHAV);  % average = 1 value

