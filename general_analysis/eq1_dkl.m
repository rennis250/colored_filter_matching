function res = eq1_dkl(x)
    global dkl_filter;
    global dkl_wout_filter;
    convtd = (x(1:3,:)*dkl_wout_filter')' + repmat(x(4,:), size(dkl_wout_filter, 1), 1);
    diffs = dkl_filter - convtd;
    euc_dist = sqrt(sum(diffs.^2.0, 2));
    res = sqrt(sum(euc_dist.^2.0));
    return;
end
