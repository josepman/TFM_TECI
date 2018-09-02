function conn_mat = gen_connectivity_matrix(sliding_window, metric_rule);

    %% GET THE CONNECTIVITY MATRIX
     if metric_rule=='corr'
         for k=1:size(sliding_window,3)
             for i=1:size(sliding_window,2)
                 for j=1:size(sliding_window,2)
                     conn_mat(i,j,k) = corr(sliding_window(:,i,k), sliding_window(:,j,k));
                 end
             end
         end
 
         %binary Connectivity matrix
         conn_matB = conn_mat;
         conn_matB(abs(conn_matB)<0.5)=0;
         conn_matB(abs(conn_matB)>=0.5)=1;
%         
%     elseif metric_rule=='partial_corr'
%         % Generate sparser matrixes than pearson correlation
%         
     elseif metric_rule == 'dist_corr'
        % Get the distance correlation matrix
        for k=1:size(sliding_window,3);
            for i=1:size(sliding_window,2);
                for j=1:size(sliding_window,2);
                    conn_mat(i,j,k) = distcorr(sliding_window(:,i,k), sliding_window(:,j,k));
                end
            end
        end

        %binary Connectivity matrix
        conn_matB = conn_mat;
        conn_matB(abs(conn_matB)<0.5)=0;
        conn_matB(abs(conn_matB)>=0.5)=1;
    
    % elseif metric_rule == 'cross-corr';   %--> to see if there are
    % connectivity in different parts (activity go along a path) --> maybe
    % it is a good idea after community searching?? and check only inside
    % the communities
    
    %elseif metric_rule == 'mutual_information';
        
    % Add your custom functions here if needed
    
    %end

end