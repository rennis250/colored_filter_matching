function rmse = conv_transform_to_min(x)
	global the_orig_data;
	global the_filter_data;
	global chroma_or_no;

	if chroma_or_no
	    convtd = (x(1:2,:)*the_orig_data')' + repmat(x(3,:), size(the_orig_data, 1), 1);
	else
	    convtd = (x(1:3,:)*the_orig_data')' + repmat(x(4,:), size(the_orig_data, 1), 1);
	end
    diffs = the_filter_data - convtd;

	% i was previously computing the root-squared-error in a slightly roundabout way,
	% not the root-MEAN-squared-error
    % euc_dist = sqrt(sum(diffs.^2.0, 2));
    % res = sqrt(sum(euc_dist.^2.0));

	% let's actually compute the root-MEAN-squared error (RMSE)
	n = size(the_filter_data, 1);
	err = diffs.^2.0;
	err = sum(err(:));
	rmse = sqrt(err/n);

    return;
end
