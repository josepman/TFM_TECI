function [final_partition, final_Q, thr_nodal, quality, consensus, consensus_simm] = community_detection(G, n, n_regions)
    
    % n = number of iterations of the optmization Louvain algorithm
    % n_regions = number of regions
    comm = zeros(n_regions,n);
    modularity = zeros(1,n);
    for i=1:n
        %Ci = modularity_louvain_und(G.graph);
        [Ci Q] = modularity_louvain_und(G);
        comm(:,i) = Ci;
        modularity(i) = Q;
        n_comm(i) = max(Ci); 
        %[Ci_h Q_h] = modularity_louvain_und(G.graph,1);
    end

    % obtain a consensus partition from the n iterations
    [final_partition, final_Q, thr_nodal, quality] = consensus_iterative(comm');
    final_partition = final_partition(1,:)';
    final_Q = Q(1);
    
    [consensus, consensus_simm, pairwise_simm] = consensus_similarity(comm');
    consensus = consensus'; 
end