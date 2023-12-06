function [rmse_fmin, rmse_dzmura, rmse_rigid, rmse_affinemap, rmse_identity] = compute_forms_of_convergence(orig_data, filter_data, chroma_only)
	n = size(orig_data, 1);

	if chroma_only
		Bfit = affine_convergence(orig_data(:, 2:3), filter_data(:, 2:3), chroma_only);
	else
		Bfit = affine_convergence(orig_data, filter_data, chroma_only);
	end

	if chroma_only
		err = Bfit - filter_data(:, 2:3);
	else
		err = Bfit - filter_data;
	end
	err = err .* err;
	err = sum(err(:));
	rmse_fmin = sqrt(err/n);

	[ret_R, ret_t] = rigid_transform_3D(orig_data, filter_data);
	Bfit = (ret_R*orig_data') + repmat(ret_t, 1, n);
	Bfit = Bfit';

	if chroma_only
		err = Bfit(:, 2:3) - filter_data(:, 2:3);
	else
		err = Bfit - filter_data;
	end
	err = err .* err;
	err = sum(err(:));
	rmse_rigid = sqrt(err/n);

	[A, b] = affinemap(orig_data, filter_data);
	Bfit = A*orig_data' + b;
	Bfit = Bfit';

	global Bfit_affinemap;
	Bfit_affinemap = Bfit;

	if chroma_only
		err = Bfit(:, 2:3) - filter_data(:, 2:3);
	else
		err = Bfit - filter_data;
	end
	err = err .* err;
	err = sum(err(:));
	rmse_affinemap = sqrt(err/n);
    
	% also, do a version where we just compute the error without doing anything.
	% this let's us see how much the various models improve on handling the
	% extra variance.

	if chroma_only
		err = orig_data(:, 2:3) - filter_data(:, 2:3);
	else
		err = orig_data - filter_data;
	end
	err = err .* err;
	err = sum(err(:));
	rmse_identity = sqrt(err/n);
    
    % and lastly, compute a version of the accepted 4 parameter model,
    % with only one scaling factor for all channels and 3 parameters for
    % the translation vector
    
	if chroma_only
		Bfit = convergence_model(orig_data(:, 2:3), filter_data(:, 2:3), chroma_only);
	else
		Bfit = convergence_model(orig_data, filter_data, chroma_only);
	end

	if chroma_only
		err = Bfit - filter_data(:, 2:3);
	else
		err = Bfit - filter_data;
	end
	err = err .* err;
	err = sum(err(:));
	rmse_dzmura = sqrt(err/n);

	return;
end