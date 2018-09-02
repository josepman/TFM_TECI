function [G,Gb] = graphs_from_sliding_window(connectivity_matrix, thr)

% It returns 2 cells (weigthed and binary) with a graphs object for each sliding-window connectivity matrix.                                                 

if ~exist('thr','var') || isempty(thr);
    thr = 0.5;
end

conn_bin = connectivity_matrix;
conn_bin(abs(conn_bin)<thr) = 0;
conn_bin(abs(conn_bin)>=thr) = 1;

G = {};
Gb = {};
for i=1:size(connectivity_matrix,3);
    G{i} = graph(connectivity_matrix(:,:,i));
    Gb{i} = graph(conn_bin(:,:,i));
end

end