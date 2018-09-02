function sliding_windows = gen_sliding_window(ROI_timeseries, w, overlapping)

    % For just a timeserie
    if length(size(ROI_timeseries))==1
        n_windows = ceil((length(ROI_timeseries)-w)/(w*(1-overlapping)));
        sliding_windows = zeros(w, n_windows);
        sliding_windows(:,1) = ROI_timeseries(1:w); %initialize the 3D Matrix
        for i=2:n_windows;
            inc = i*w*(1-overlapping);
            sliding_windows(:,i) = ROI_timeseries(inc+1:inc+w);
        end
    else
    % For matrix of timeseries
    % Create the different sliding-windows for a given subject
        sliding_windows = ROI_timeseries(1:w,:); %initialize the 3D Matrix
        n_windows = ceil((size(ROI_timeseries,1)-w)/(w*(1-overlapping)));
        for i=1:n_windows-1;
            inc = i*w*(1-overlapping);
            sliding_windows(:,:,i+1) = ROI_timeseries(inc+1:inc+w,:);
        end
    end
end

