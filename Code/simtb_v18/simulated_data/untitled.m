for i=1:size(sliding_windows,3)
    ncols = size(sliding_windows,3);
    nrows = 3;
    figure(1)
    corr_mat = corr(sliding_windows(:,:,i));
    nodes_of_interest = [3,4,5,6,12,16,17,18,21,22];
    
    for j=1:size(sliding_windows,2)
        if any(nodes_of_interest(:) == j)
            subplot(nrows,ncols,i)
            plot(sliding_windows(:,j,i)), hold on
        else
            subplot(nrows,ncols,i+size(sliding_windows,3))
            plot(sliding_windows(:,j,i)), hold on
        end
    end
    
    subplot(nrows,ncols,i+2*size(sliding_windows,3))
    pcolor(corr_mat)
end


close all
figure('Position', [100, 100, 400, 400]);
axes
axis('square')  % [EDITED: or better 'square' !]
for i=1:size(sliding_windows,3)
    ncols = 4;
    nrows = 2*ceil(size(sliding_windows,3)/ncols);
    

    corr_mat = corr(sliding_windows(:,:,i));
    nodes_of_interest = [3,4,5,6,12,16,17,18,21,22];

    k = ceil(i/4);
    for j=1:size(sliding_windows,2)
        if any(nodes_of_interest(:) == j)
            subplot(nrows,ncols, i+ncols*(k-1))
            plot(sliding_windows(:,j,i)), hold on
            plot(task_windows(:,i), 'k--', 'linewidth',2)
        %else
        %    subplot(nrows,ncols,k+4)
        %    plot(sliding_windows(:,j,i)), hold on
        %    plot(task_windows(:,i), 'k--', 'linewidth',2)
        end
    end
    
    subplot(nrows,ncols,i+ncols*k)
    pcolor(corr_mat)
    
end

