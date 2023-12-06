function res = eq1_lab(x)
    global lab_filter;
    global lab_wout_filter;
    convtd = (x(1:3,:)*lab_wout_filter')' + repmat(x(4,:), size(lab_wout_filter, 1), 1);
    diffs = lab_filter - convtd;
    euc_dist = sqrt(sum(diffs.^2.0, 2));
    res = sqrt(sum(euc_dist.^2.0));
    return;
end
