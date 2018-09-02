function temporal_null_model = temporal_rand(connectivity_matrix)
    v = [1:size(connectivity_matrix,3)];
    shuf = v(randperm(length(v)));
    temporal_null_model = connectivity_matrix(:,:,shuf);
end


