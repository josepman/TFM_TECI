function G = bct_measures(CorrMatrix, option, thr, XYZ)
% From BrainRegions.csv
% BrainRegions.Properties.VariableNames{4} = 'x';
% BrainRegions.Properties.VariableNames{5} = 'y';
% BrainRegions.Properties.VariableNames{6} = 'z';
% coord = [BrainRegions.x,BrainRegions.y, BrainRegions.z];


    if option=='weighted'
            G.graph = CorrMatrix;
    elseif option=='binary'
            thr = 0.5;
            CorrMatrix(abs(CorrMatrix) < thr) = 0;
            CorrMatrix(abs(CorrMatrix) > thr) = 1;
            G.graph = CorrMatrix;
    end
    
    for i=1:size(G.graph,3)
        %% GENERAL MEASURES
        [G.dens(i),G.nodes(i),G.edges(i)] = density_und(G.graph(:,:,i));
        %[G.N,G.E] = rentian_scaling(G.graph, XYZ, 6);    % Topo-measure of the cost-efficiency integration. See value of n (num. of partitions)
                                               % Outputs:
                                            %	N       nx1 vector of the number of nodes in each of the n partitions.
                                            %	E       nx1 vector of the number of edges crossing the boundary of each
                                            %           partition.

        % Small-World
    %     if (C>Cr & abs(L-Lr)<0.3)
    %         G.sw = C/L;
    %     else
    %         G.sw = 'NA';
    %     end

    % Community detection
        %   The optimal community structure is a subdivision of the network into
        %   groups of nodes which have a high number of within-group connections
        %   and a low number of between group connections. 
        gamma = 1;
        %M0 = ; %initial vector community affiliation
        B = 'negative_asym';
        [G.community_afiliation(:,i), G.optimized_B(i)]=community_louvain(abs(G.graph(:,:,i)),gamma);    % negative weights are allowed.   Louvain --> non-overlapping nodes
                                                  % gamma > 1 == smaller modules
                                                  % gamma < 1 == larger modules
                                                  % gamma = 1 == classic modularity
                                                  % (default)
                                                  % B={modularity, potts,
                                                  % negative_sym, negative_asym}.
                                                  % See Rubinov and Sporns (2011) Neuroimage 56:2068-79 for a discussion
        [G.comm_struc(:,i), G.maximized_mod(i)] = modularity_und(G.graph(:,:,i));        %   The modularity is a statistic that quantifies the degree to which the
                                                  %   network may be subdivided into such clearly delineated groups.


        % Degree
        G.local_degree(:,i) = degrees_und(G.graph(:,:,i))';
        G.avg_degree(:,i) = G.local_degree(:,i)/length(G.graph(:,:,i))';
        G.node_str(:,i) = strengths_und(G.graph(:,:,i))';                % = sum(G)
        G.avg_str(:,i) = G.node_str(:,i)/length(G.graph(:,:,i))';

        % Centrality measures:
        G.Z = module_degree_zscore(G.graph(:,:,i), G.community_afiliation(:,i));
        G.large_eig = eigenvector_centrality_und(G.graph(:,:,i));
        %G.r = pagerank_centrality(G.graph, G.D, falff);
        [G.Erange,G.eta,G.Eshort,G.fs] = erange(G.graph(:,:,i));
        G.particip(:,i) = participation_coef(G.graph(:,:,i), G.community_afiliation(:,i));   % 0=undirected, default

        % Similarity
        G.match_ind{i} = matching_ind_und(G.graph(:,:,i));
        % D = distance matrix, char_path_len = landa
                                                                                                %   Outputs:    lambda,         network characteristic path length
                                                                                                %               efficiency,     network global efficiency
                                                                                                %               ecc,            nodal eccentricity
                                                                                                %               radius,         network radius
                                                                                                %               diameter,       network diameter


        %% EXCLUSIVE WEIGTHED MESURES                                        
        if option=='weighted'
            % Degree
            [G.Sp, G.Sn, G.vpos(:,i), G.vneg(:,i)] = strengths_und_sign(G.graph(:,:,i));
            G.Spos = G.Sp';
            G.Sneg = G.Sn';
            % Centrality measures
            G.BC(:,i) = betweenness_wei(G.graph(:,:,i));
            [G.edgeBC_matrix{i}, G.nodeBC_vector(:,i)] = edge_betweenness_wei(G.graph(:,:,i));
            [G.Hpos(:,i), G.Hneg(:,i)] = diversity_coef_sign(G.graph(:,:,i), G.community_afiliation(:,i));          % Requires community afiliation vector Ci. G could contain negative weights
            [G.GWpos(:,i), G.GWneg(:,i)] = gateway_coef_sign(G.graph(:,:,i), G.community_afiliation(:,i), 2);   % centtype=1, node strength; = 2 betwenness centrality
            [G.Ppos(:,i), G.Pneg(:,i)] = participation_coef_sign(G.graph(:,:,i), G.community_afiliation(:,i));

            % Clustering
            G.C(:,i) = clustering_coef_wu(G.graph(:,:,i));
            G.transitv(:,i) = transitivity_wu(G.graph(:,:,i));

            % Core
            G.assort(:,i) = assortativity_wei(G.graph(:,:,i), 0);
            [G.loc_assort_pos(:,i), G.loc_assort_neg(:,i)] = local_assortativity_wu_sign(G.graph(:,:,i));
            G.rich(:,i) = rich_club_wu(G.graph(:,:,i))';       % rich club curve for weighted graph
            %[G.s_core, G.sn] = score_wu(G.graph,s);        % =k-core, with strength of node instead of degree

            % Distance
            G.D{i} = distance_wei(G.graph(:,:,i));                   % lengths of shortest paths between all pairs of nodes.
            G.Eglob(i) = efficiency_wei(G.graph(:,:,i));
            G.Eloc(:,i) = efficiency_wei(G.graph(:,:,i),1);


        %% EXCLUSIVE BINARY MEASURES
        elseif option == 'binary'
            % Centrality measures
            G.BC(:,i) = betweenness_bin(G.graph(:,:,i));
            [G.edgeBC_matrix, G.nodeBC_vector] = edge_betweenness_bin(G.graph(:,:,i));
            [G.coreness,G.kn_central] = kcoreness_centrality_bu(G.graph(:,:,i));
            G.sub_centrality = subgraph_centrality(G.graph(:,:,i));

            % Clustering
            G.C(:,i) = clustering_coef_bu(G.graph(:,:,i));
            [G.comps,G.comp_sizes] = get_components(G.graph(:,:,i));                        % isolated nodes appears as size = 1
            G.transitv = transitivity_bu(G.graph(:,:,i));

            % Core
            G.assort = assortativity_bin(G.graph(:,:,i), 0);
            %[G.graphkcore, G.kn_core, G.peelorder, G.peellevel] = kcore_bu(G.graph, k);      %k = min. degree of node for the core
            [G.rich, G.Nk, G.Ek] = rich_club_bu(G.graph(:,:,i));                    %k = max. level at calculated

            % Distance
            G.D = distance_bin(G.graph(:,:,i));                                             % lengths of shortest paths between all pairs of nodes.
            G.Eglob = efficiency_bin(G.graph(:,:,i));
            G.Eloc = efficiency_bin(G.graph(:,:,i),1);

            % Similarity
            [G.edge_overlap_mat, G.edge_overlap_vec, G.deg_ij] = edge_nei_bu(G.graph(:,:,i));
            G.grapht = gtom(G.graph(:,:,i),100);
        end;

        % Distance
        [G.char_path_len(i), G.global_efficiency(i), G.ecc{i}, G.rad, G.diam] = charpath(G.D{i}); 
        G.radius = G.rad';
        G.diameter = G.diam';

    end
end
    
