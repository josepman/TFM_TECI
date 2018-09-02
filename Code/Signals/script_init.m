
for i=1:30;
    if(i < 10);
        subject = 'ROISignals_sample_0%d.mat';
    else
        subject = 'ROISignals_sample_%d.mat';
    end;
    roi_subject = sprintf(subject,i);
    subject_xx = load(roi_subject);
    formatSpec = "ROISignal_subject_%d.csv";
    str = sprintf(formatSpec,i);
    csvwrite(str,subject_xx.ROISignals);
end

% Get the distance correlation matrix
x = subject_xx.ROISignals;
for i=1:size(subject_xx.ROISignals,2);
    for j=1:size(subject_xx.ROISignals,2);
        dcor_mat(i,j) = distcorr(x(:,i),x(:,j));
    end
end
