clear all; close all; clc;

orig_data = csvread('../data/data_from_dzmura/orig_data.csv');
orig_data = [zeros(size(orig_data, 1), 1), orig_data(:, 1), orig_data(:, 2)];

obsvrs = {'CC', 'MD', 'MW', 'OR'};

rmse = {};

for obs = 1:length(obsvrs)
  obs_data = csvread(['../data/data_from_dzmura/' obsvrs{obs} '.csv']);
  obs_data = [zeros(size(obs_data, 1), 1), obs_data(:, 1), obs_data(:, 2)];

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%       DKL           %%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%%%%%%%%%%%
  % make chroma only version of DKL data and compute forms of convergence for it

  [rmse_fmin, rmse_rigid, rmse_affinemap, rmse_identity] = compute_forms_of_convergence(orig_data, obs_data, true);
  rmse = [rmse; {rmse_fmin, rmse_rigid, rmse_affinemap, rmse_identity, 'dkl', 'chroma', obsvrs{obs}}];

  global Bfit_affinemap;
  figure(1);
  clf;
  hold on;
  axis square;
  p1 = plot(obs_data(:, 2), obs_data(:, 3), 'bo');
  p2 = plot(Bfit_affinemap(:, 2), Bfit_affinemap(:, 3), 'ro');
  for x = 1:size(obs_data, 1)
    line([obs_data(x, 2) Bfit_affinemap(x, 2)], [obs_data(x, 3) Bfit_affinemap(x, 3)], 'Color', 'black');
  end
  axis([-1 1 -1 1]);
  legend('Observer settings', 'Best general affine fit');

  fn = ['./output/dzmura_convergence/affine_map_' obsvrs{obs} '.png'];
  print('-dpng', fn);
end

rmse_dzmura = [{'rmse_fmin', 'rmse_rigid', 'rmse_affinemap', 'rmse_identity', 'cspace', 'chroma_or_full', 'obs'}; rmse];
writecell(rmse_dzmura);
