clear all; close all; clc;

expns = {'3'};
exp_names = {'filter'};
filter_exp = 1;
eizo_exps = [filter_exp]; % which exp in the iteration loop is it?
filter_exps = [filter_exp];

avg_bkgd_color = [0.8, 0.8, 0.8];

gRGB_filter_exp = csvread('../calibration/eizo_gamma.csv');

gR = gRGB_filter_exp(1);
gG = gRGB_filter_exp(2);
gB = gRGB_filter_exp(3);

vor_texture = csvread('vor.csv');
lms_absorp_curves = csvread('linss2_10e_1.csv');
lms_absorp_curves = lms_absorp_curves(1:14:end, 2:end);

expc = 1;
expn = expns{expc};

monxyY = csvread('../calibration/eizo_chroma.csv');
DKL2RGB = dkl2rgbFromCalib('../calibration/eizo_chroma.csv');
RGB2LMS = rgb2lmsFromCalib('../calibration/eizo_mon_spectra.csv');

LMS2RGB = inv(RGB2LMS);
monxyz = xyY2XYZ(monxyY);
RGB2DKL_T = inv(DKL2RGB);

fns = dir(['../data/exp' expn '/*.mat']);
imgns = dir(['../images/exp' expn '/*.png']);

% exps x imgs x masks
load('../data/imgstats.mat');

stats_for_curr_exp = squeeze(stats(expc, :, :));

%% first make the average filter match of observers
redtrans_file = 'munsell_red_better.spd';
greentrans_file = 'munsell_green_better.spd';
bluetrans_file = 'munsell_blue_better.spd';
yellowtrans_file = 'munsell_yellow_better.spd';

monitor_spectra = csvread('../calibration/eizo_mon_spectra.csv');
% [vor_filter_map, filter_idxs, redtrans, bluetrans, greentrans, yellowtrans] = ...
%    gen_voronoi_filter_map(monitor_spectra, redtrans_file, bluetrans_file, greentrans_file, yellowtrans_file, false, false);
% more random achrom colors
% [vor_filter_map, filter_idxs, redtrans, bluetrans, greentrans, yellowtrans] = ...
%     gen_voronoi_filter_map(monitor_spectra, redtrans_file, bluetrans_file, greentrans_file, yellowtrans_file, false, true);
% in this case, we make the voronoi map multicolored
% [vor_filter_map, filter_idxs, redtrans, bluetrans, greentrans, yellowtrans] = ...
%   gen_voronoi_filter_map(monitor_spectra, redtrans_file, bluetrans_file, greentrans_file, yellowtrans_file, true, false);
% and here we make a super random colored voronoi map
[vor_filter_map, filter_idxs, redtrans, bluetrans, greentrans, yellowtrans] = ...
  gen_voronoi_filter_map(monitor_spectra, redtrans_file, bluetrans_file, greentrans_file, yellowtrans_file, true, true);

save('vor_rand_col.mat', 'vor_filter_map');

ldrgyvs = [];
goodics = [];
illums = {};
hilos = [];
bodies = [];
ic = 1;
for imgc = 1:length(imgns)
  imgn = imgns(imgc).name;
  fparts = strsplit(imgn, '_');
  if strcmp(fparts{end}, 'lighter.png')
    if strcmp(fparts{6}, 'red') || strcmp(fparts{6}, 'green')
      continue;
    end
    illum = fparts{6};
    body = fparts{9};
    hilo = 1;
  else
    if strcmp(fparts{7}, 'red') || strcmp(fparts{7}, 'green')
      continue;
    end
    illum = fparts{7};
    body = fparts{10};
    hilo = 0;
  end

  oc = 1;
  for obsc = 1:length(fns)
    fn = fns(obsc).name;
    load(['../data/exp' expn '/' fn], 'data', 'subID', 'trialOrder');

    if expc == filter_exp
      if strcmp(subID, '05') || strcmp(subID, 'go') || strcmp(subID, 'rr')

       continue;
      end
    end

    idxs = find(trialOrder == imgc);
    trialc = 1;
    for idx = idxs
      % reproduce obs match
      rg_mix = data(idx, 1);
      by_mix = data(idx, 2);
      ld = data(idx, 3);
      ldrgyvs(ic, oc, trialc, :) = [ld, rg_mix, by_mix];
      trialc = trialc + 1;
    end
    oc = oc + 1;
  end

  goodics = [goodics; imgc];
  illums{ic} = illum;
  hilos(ic) = hilo;
  bodies{ic} = body;
  ic = ic + 1;
end

ldrgyv = squeeze(mean(ldrgyvs, 3)); % obs mean
ldrgyvm = squeeze(mean(ldrgyv, 2)); % pop mean

vor_col_wout_filter = vor_filter_map;
vor_col_wout_filter(isnan(vor_col_wout_filter) | isinf(vor_col_wout_filter)) = 0.0;
vor_wout_filter_rgb = (vor_col_wout_filter*lms_absorp_curves*LMS2RGB').^(1.0/2.2);

vor_wout_filt_rgbs_gc_lin(:, 1) = vor_wout_filter_rgb(:, 1).^gR;
vor_wout_filt_rgbs_gc_lin(:, 2) = vor_wout_filter_rgb(:, 2).^gG;
vor_wout_filt_rgbs_gc_lin(:, 3) = vor_wout_filter_rgb(:, 3).^gB;

vor_wout_filter_dkl = real(rgb2dkl(RGB2DKL_T, vor_wout_filt_rgbs_gc_lin'))';

mask = imread('../images/masks/obj_mask_correct.png');
mask_r = mask(:, :, 1);
mask_g = mask(:, :, 2);
mask_b = mask(:, :, 3);
mask = mask_r;
reject_pxs = mask_r == 127 & mask_g == 127 & mask_b == 127;
reject_pxs = logical(reject_pxs(:));
bkgd_pxs = mask_r == 0 & mask_g == 0 & mask_b == 0;
bkgd_pxs = logical(bkgd_pxs(:));
obj_pxs = ~(reject_pxs | bkgd_pxs);
obj_pxs = obj_pxs(:);
mask(obj_pxs) = 1;
mask = logical(mask(:));

%% then load the corresponding glaven and calculate convergence for it

h = figure(1);
fullfig(h);
% for ic = 1:length(goodics)
for ic = 9
  ldrgyv = ldrgyvm(ic, :);
  ld = ldrgyv(1);
  rg_mix = ldrgyv(2);
  by_mix = ldrgyv(3);

  filter = ld.*(((rg_mix.*redtrans(:,2) + (1 - rg_mix).*greentrans(:,2)) + (by_mix.*bluetrans(:,2) + (1 - by_mix).*yellowtrans(:,2))));
  vor_col2 = vor_filter_map;
  vor_col2(filter_idxs(:), :) = bsxfun(@times, vor_filter_map(filter_idxs(:), :), filter'.^2.2);
  vor_col2(isnan(vor_col2) | isinf(vor_col2)) = 0.0;
  vor_rgb = (vor_col2*lms_absorp_curves*LMS2RGB').^(1.0/2.2);

  vor_filt_rgbs_gc_lin(:, 1) = vor_rgb(:, 1).^gR;
  vor_filt_rgbs_gc_lin(:, 2) = vor_rgb(:, 2).^gG;
  vor_filt_rgbs_gc_lin(:, 3) = vor_rgb(:, 3).^gB;

  vor_filter_dkl = real(rgb2dkl(RGB2DKL_T, vor_filt_rgbs_gc_lin'))';

  [ldgr_d1_vor, ldgr_d3_vor, rggr_d1_vor, rggr_d2_vor, yvgr_d2_vor, yvgr_d3_vor, ld_interp_d1_vor, ld_interp_d3_vor, rg_interp_d1_vor, rg_interp_d2_vor, yv_interp_d2_vor, yv_interp_d3_vor] = divergence_map(vor_wout_filter_dkl(filter_idxs(:), :), vor_filter_dkl(filter_idxs(:), :));
  [fmin_map_filt, aff_filt, rmse_filt, rmse_dz, ~, rmse_aff, rmse_id, dz_map] = compute_forms_of_convergence(vor_wout_filter_dkl(filter_idxs(:), :), vor_filter_dkl(filter_idxs(:), :), false);

  ldrgyv
  [(rmse_id - rmse_filt)/rmse_id, (rmse_id - rmse_aff)/rmse_id, (rmse_id - rmse_dz)/rmse_id]
  fmin_map_filt
  dz_map

  imgn = imgns(goodics(ic)).name;
  imgn
  orig_img_with_filter = im2double(imread(['../images/exp' expn '/' imgn]));
  img_with_filter = gammaCorr(orig_img_with_filter, gRGB_filter_exp);
  img_wout_filter = gammaCorr(im2double(imread(['../images/empty_scene_' illums{ic} '_illum.png'])), gRGB_filter_exp);

  rwof = img_wout_filter(:, :, 1); rwof = rwof(:); r_wout_filter = rwof(obj_pxs(:));
  gwof = img_wout_filter(:, :, 2); gwof = gwof(:); g_wout_filter = gwof(obj_pxs(:));
  bwof = img_wout_filter(:, :, 3); bwof = bwof(:); b_wout_filter = bwof(obj_pxs(:));

  dkl_wout_filter = real(rgb2dkl(RGB2DKL_T, [r_wout_filter(:), g_wout_filter(:), b_wout_filter(:)]'))';

  r = img_with_filter(:, :, 1); r = r(:); r_img = r(obj_pxs(:));
  g = img_with_filter(:, :, 2); g = g(:); g_img = g(obj_pxs(:));
  b = img_with_filter(:, :, 3); b = b(:); b_img = b(obj_pxs(:));

  dkl_img = real(rgb2dkl(RGB2DKL_T, [r_img(:), g_img(:), b_img(:)]'))';

  [ldgr_d1, ldgr_d3, rggr_d1, rggr_d2, yvgr_d2, yvgr_d3, ld_interp_d1, ld_interp_d3, rg_interp_d1, rg_interp_d2, yv_interp_d2, yv_interp_d3] = divergence_map(dkl_wout_filter, dkl_img);
  [fmin_map_glaven, aff_glaven, rmse_glaven, ~, ~, glav_aff, ~, dz_glav] = compute_forms_of_convergence(dkl_wout_filter, dkl_img, false);

  [rmse_glaven, glav_aff]
  fmin_map_glaven

  clf;
  set(gcf, 'color', 'w');

  subplot(2, 2, 1);
  hold on;
  imshow(orig_img_with_filter);

  subplot(2, 2, 2);
  hold on;
  imshow(reshape(vor_rgb, 256, 256, 3));

  subplot(2, 2, 3);
  hold on;
  grid on;
  h1 = streamslice(yvgr_d3, ldgr_d3, yv_interp_d3, ld_interp_d3);
  set(h1, 'Color', [0, 0, 0]);
  % quiver(yvgr_d3_vor, ldgr_d3_vor, yv_interp_d3_vor, ld_interp_d3_vor, 15, 'Color', 'red');
  h2 = streamslice(yvgr_d3_vor, ldgr_d3_vor, yv_interp_d3_vor, ld_interp_d3_vor);
  set(h2, 'Color', [1, 0, 0]);
  % if hilos(ic) == 1
  %   for x = 1:length(the_covms_ld_yv)
  %     if hilos(x) == 0 && strcmp(illums{ic}, illums{x}) && strcmp(bodies{ic}, bodies{x})
  %       break;
  %     end
  %     companion_filt_to_plot = ldrgyvm(x, :);
  %     companion_ld = companion_filt_to_plot(1);
  %     rg_mix = companion_filt_to_plot(2);
  %     by_mix = companion_filt_to_plot(3);
  %   end
  % end
  xlabel('S-(L+M)');
  ylabel('L+M+S');
  axis square;

  subplot(2, 2, 4);
  hold on;
  grid on;
  h1 = streamslice(rggr_d2, yvgr_d2, rg_interp_d2, yv_interp_d2);
  set(h1, 'Color', [0, 0, 0]);
  % quiver(rggr_d2_vor, yvgr_d2_vor, rg_interp_d2_vor, yv_interp_d2_vor, 50, 'Color', 'red');
  h2 = streamslice(rggr_d2_vor, yvgr_d2_vor, rg_interp_d2_vor, yv_interp_d2_vor);
  set(h2, 'Color', [1, 0, 0]);
  axis([-0.3, 0.3, -0.3, 0.3]);
  ylabel('S-(L+M)');
  xlabel('L-M');
  axis square;

  export_fig(['../figures/glaven_filter_conv_comp_' num2str(imgn) '.tif'], '-m2', '-r500');
end
