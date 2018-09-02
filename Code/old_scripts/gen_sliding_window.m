function sliding_windows = gen_sliding_window(ROI_timeseries, w, overlapping)
    % Create the different sliding-windows for a given subject
    sliding_windows = ROI_timeseries(1:w,:); %initialize the 3D Matrix
    for i=1:size(ROI_timeseries,2)-1;
        inc = i*w*overlapping;
        sliding_windows(:,:,i+1) = ROI_timeseries(inc+1:inc+w,:);
    end
end