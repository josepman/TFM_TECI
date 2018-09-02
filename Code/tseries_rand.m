function tseries_null_model = tseries_rand(ROI_timeseries)
    tseries_null_model = zeros(size(ROI_timeseries));
    for i=1:size(ROI_timeseries,2)
        v = [1:size(ROI_timeseries,1)];
        shuf = v(randperm(length(v)));
        tseries_null_model(:,1) = ROI_timeseries(shuf,i);
    end
end

