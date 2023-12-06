clear all; close all; clc;

DKL2RGB = dkl2rgbFromCalib('../calibration/eizo_chroma.csv');
RGB2DKL_T = inv(DKL2RGB);
RGB2LMS = rgb2lmsFromCalib('../calibration/eizo_mon_spectra.csv');
LMS2RGB = inv(RGB2LMS);

gRGB_filter_exp = csvread('../calibration/eizo_gamma.csv');

gR = gRGB_filter_exp(1);
gG = gRGB_filter_exp(2);
gB = gRGB_filter_exp(3);

imgn = '../images/patterned/mitsuba_caustics_rgb_trans_lightest_illum_blue_refl_munsell_blue_lightest.png';

% ldrgyv = [0.3714, 0.3892, 0.8737];
ldrgyv = [0.6923, 0.2980, 0.8878];

% extracted from random color map, but with affinemap function
% fmin_map_filt = [0.8124, -0.0102, -0.0146;
%    -0.0116, 0.8475, -0.0037;
%     0.0039, -0.0006, 0.8490;
%    -0.1732, -0.0111, 0.0066];

% extracted from random color map
% fmin_map_filt = [0.8699, 0.0032, -0.0049;
%     0.0210, 0.8500, 0.0007;
%    -0.1111, -0.0127, 0.8417;
%    -0.1353, 0.0106, -0.0701];

% extracted from original achrom dist
% fmin_map_filt = [0.7440  -28.6767   -4.1774;
%    -0.0677   11.0990   -0.8800;
%     0.3289   -8.3033   -2.2506;
%    -0.1571    0.0614    0.4422];
% fmin_M_filt = fmin_map_filt(1:3, :);
% fmin_t_filt = squeeze(fmin_map_filt(4, :));

% extracted from original achrom with dzmura model
% mt = [0.1704, -0.8324, -0.0412, 0.0421];
% extracted from colored with dzmura model
% mt = [0.1694, -0.8334, -0.0526, 0.0530];
% diff colored filter - light blue
mt = [0.6377, -0.3626, -0.2573, 0.1164];
fmin_M_filt = [mt(1), 0.0, 0.0;
               0.0, mt(1), 0.0;
               0.0, 0.0, mt(1)];
fmin_t_filt = [mt(2), mt(3), mt(4)];

fmin_map_glaven = [0.6325, -0.6668, 0.1296;
   -0.2868, 0.3457, -0.0393;
    0.3818, -0.4217, 0.1448;
   -0.3804, -0.2701, 0.3836];
fmin_M_glaven = fmin_map_glaven(1:3, :);
fmin_t_glaven = squeeze(fmin_map_glaven(4, :));

vor_texture = csvread('vor.csv');
lms_absorp_curves = csvread('linss2_10e_1.csv');
lms_absorp_curves = lms_absorp_curves(1:14:end, 2:end);
monitor_spectra = csvread('../calibration/eizo_mon_spectra.csv');

redtrans_file = '../base_stimuli/spectra/munsell_red_better.spd';
greentrans_file = '../base_stimuli/spectra/munsell_green_better.spd';
bluetrans_file = '../base_stimuli/spectra/munsell_blue_better.spd';
yellowtrans_file = '../base_stimuli/spectra/munsell_yellow_better.spd';

[vor_filter_map, filter_idxs, redtrans, bluetrans, greentrans, yellowtrans] = ...
    gen_voronoi_filter_map(monitor_spectra, redtrans_file, bluetrans_file, greentrans_file, yellowtrans_file, false, false);

% random colored voronoi
load('vor_rand_col.mat', 'vor_filter_map');
vor_filter_map(isnan(vor_filter_map) | isinf(vor_filter_map)) = 0.0;

ld = ldrgyv(1);
rg_mix = ldrgyv(2);
by_mix = ldrgyv(3);
filter = ld.*(((rg_mix.*redtrans(:,2) + (1 - rg_mix).*greentrans(:,2)) + (by_mix.*bluetrans(:,2) + (1 - by_mix).*yellowtrans(:,2))));
vor_col2 = vor_filter_map;
vor_col2(filter_idxs(:), :) = bsxfun(@times, vor_filter_map(filter_idxs(:), :), filter'.^2.2);
vor_col2(isnan(vor_col2) | isinf(vor_col2)) = 0.0;

vor_rgb = (vor_col2*lms_absorp_curves*LMS2RGB').^(1/2.2);
vor_rgb_gc_lin(:, 1) = vor_rgb(:, 1).^gR;
vor_rgb_gc_lin(:, 2) = vor_rgb(:, 2).^gG;
vor_rgb_gc_lin(:, 3) = vor_rgb(:, 3).^gB;
vor_filter_dkl = real(rgb2dkl(RGB2DKL_T, vor_rgb_gc_lin'))';
vor2_filter_rgb = real(dkl2rgb(DKL2RGB, vor_filter_dkl'))';

% make a version without the filter
vor_wout_filt_rgb = (vor_filter_map*lms_absorp_curves*LMS2RGB').^(1/2.2);
vor_wout_filt_rgb(isnan(vor_wout_filt_rgb) | isinf(vor_wout_filt_rgb)) = 0.0;
vor_wout_filt_rgb_gc_lin(:, 1) = vor_wout_filt_rgb(:, 1).^gR;
vor_wout_filt_rgb_gc_lin(:, 2) = vor_wout_filt_rgb(:, 2).^gG;
vor_wout_filt_rgb_gc_lin(:, 3) = vor_wout_filt_rgb(:, 3).^gB;
vor_wout_filter_dkl = real(rgb2dkl(RGB2DKL_T, vor_wout_filt_rgb_gc_lin'))';

% lets look at what happens when we change voronoi map to dkl and apply conv map transform
% from glaven
n = length(find(filter_idxs));

vor_glaven_filtered_dkl = vor_wout_filter_dkl;
vor_glaven_filtered_dkl(filter_idxs, :) = (fmin_M_glaven*vor_wout_filter_dkl(filter_idxs, :)')' + repmat(fmin_t_glaven, n, 1);
vor_glaven_filtered_rgb = real(dkl2rgb(DKL2RGB, vor_glaven_filtered_dkl'))';

% and if we compare that also to the application of the conv map for the best matching filter
vor_filter_filtered_dkl = vor_wout_filter_dkl;
vor_filter_filtered_dkl(filter_idxs, :) = (fmin_M_filt*vor_wout_filter_dkl(filter_idxs, :)')' + repmat(fmin_t_filt, n, 1);
vor_filter_filtered_rgb = real(dkl2rgb(DKL2RGB, vor_filter_filtered_dkl'))';

% let's plot all of them in a figure and see
figure(1);
clf;
subplot(1, 3, 1);
title('Pop Avg Match');
imshow(reshape(vor_rgb, 256, 256, 3));

subplot(1, 3, 2);
vor_glaven_filtered_rgb_gc(:, 1) =  vor_glaven_filtered_rgb(:, 1).^(1/gR);
vor_glaven_filtered_rgb_gc(:, 2) =  vor_glaven_filtered_rgb(:, 2).^(1/gG);
vor_glaven_filtered_rgb_gc(:, 3) =  vor_glaven_filtered_rgb(:, 3).^(1/gB);
title('Glaven Affine Map');
imshow(reshape(vor_glaven_filtered_rgb_gc, 256, 256, 3));

subplot(1, 3, 3);
vor_filter_filtered_rgb_gc(:, 1) =  vor_filter_filtered_rgb(:, 1).^(1/gR);
vor_filter_filtered_rgb_gc(:, 2) =  vor_filter_filtered_rgb(:, 2).^(1/gG);
vor_filter_filtered_rgb_gc(:, 3) =  vor_filter_filtered_rgb(:, 3).^(1/gB);
title('Filter Match Affine Map');
imshow(reshape(vor_filter_filtered_rgb_gc, 256, 256, 3));

export_fig('../figures/conv_map_fig_filter_examples.png');

generated_diff_vecs = vor_filter_filtered_dkl - vor_wout_filter_dkl;
user_match_diff_vecs = vor_filter_dkl - vor_wout_filter_dkl;

fis = find(filter_idxs);
idxs = randi([0, length(fis)], 10);

h = figure(2);
fullfig(h);
clf;
subplot(1, 2, 1);
hold on;
quiver(vor_wout_filter_dkl(fis(idxs, :), 2), vor_wout_filter_dkl(fis(idxs, :), 3), generated_diff_vecs(fis(idxs, :), 2), generated_diff_vecs(fis(idxs, :), 3), 0, 'r');
quiver(vor_wout_filter_dkl(fis(idxs, :), 2), vor_wout_filter_dkl(fis(idxs, :), 3), user_match_diff_vecs(fis(idxs, :), 2), user_match_diff_vecs(fis(idxs, :), 3), 0, 'k');
axis square;

subplot(1, 2, 2);
hold on;
quiver(vor_wout_filter_dkl(fis(idxs, :), 3), vor_wout_filter_dkl(fis(idxs, :), 1), generated_diff_vecs(fis(idxs, :), 3), generated_diff_vecs(fis(idxs, :), 1), 0, 'r');
quiver(vor_wout_filter_dkl(fis(idxs, :), 3), vor_wout_filter_dkl(fis(idxs, :), 1), user_match_diff_vecs(fis(idxs, :), 3), user_match_diff_vecs(fis(idxs, :), 1), 0, 'k');
axis square;

export_fig('../figures/conv_map_fig_quiver_plots.png');

figure(3);
clf;
subplot(1, 2, 1);
hist(generated_diff_vecs(filter_idxs(:), :));
subplot(1, 2, 2);
hist(user_match_diff_vecs(filter_idxs(:), :));

export_fig('../figures/conv_map_fig_histos.png');

% compute best fitting affine map for flat filter and for filter generated from best fitting
clear the_orig_data;
global the_orig_data;
the_orig_data = vor_wout_filter_dkl;
clear the_filter_data;
global the_filter_data;
the_filter_data = vor_filter_filtered_dkl;
x0 = [1.0, 1.0, 1.0;
      1.0, 1.0, 1.0;
      1.0, 1.0, 1.0;
      0.0, 0.0, 0.0];
n = size(vor_wout_filter_dkl, 1);
x_filter_filtered = fminsearch(@(x) conv_transform_to_min(x), x0, optimset('MaxFunEvals', 5000));

clear the_orig_data;
global the_orig_data;
the_orig_data = vor_wout_filter_dkl;
clear the_filter_data;
global the_filter_data;
the_filter_data = vor_filter_dkl;
x0 = [1.0, 1.0, 1.0;
      1.0, 1.0, 1.0;
      1.0, 1.0, 1.0;
      0.0, 0.0, 0.0];
n = size(vor_wout_filter_dkl, 1);
x_user_match = fminsearch(@(x) conv_transform_to_min(x), x0, optimset('MaxFunEvals', 5000));

small_img = [-1.0000   -0.0030         0;
    -0.6017   -0.0020    0.0977;
    -0.4097   -0.0006    0.1199;
    -0.2197   -0.0000    0.1358;
    -0.0310    0.0013    0.1472]
(x_user_match(1:3,:)*small_img')' + repmat(x_user_match(4,:), 5, 1)
after_trans = [-0.9411    0.0092    0.0109;
   -0.6350   -0.0066    0.0881;
   -0.4806   -0.0048    0.1111;
   -0.3261   -0.0016    0.1303;
   -0.1715    0.0036    0.1469];
vecs = after_trans - small_img;

figure;
quiver3(small_img(:,1), small_img(:,2), small_img(:,3), vecs(:,1),vecs(:,2),vecs(:,3));

idxs = [17526, 17534, 17777, 18051];
before_filter = [-0.6017   -0.0030    0.0977;
   -0.4097   -0.0020    0.1199;
   -0.2197   -0.0006    0.1358;
   -0.0310    0.0013    0.1472];
after_filter = [   -0.9337   -0.0246    0.0423;
   -0.9018   -0.0360    0.0573;
   -0.8702   -0.0472    0.0709;
   -0.8389   -0.0582    0.0834];
vecs = after_filter - before_filter;

figure;
quiver3(before_filter(:,1), before_filter(:,2), before_filter(:,3), vecs(:,1),vecs(:,2),vecs(:,3));

clear the_orig_data;
global the_orig_data;
the_orig_data = before_filter;
clear the_filter_data;
global the_filter_data;
the_filter_data = after_filter;
x0 = [1.0, 1.0, 1.0;
      1.0, 1.0, 1.0;
      1.0, 1.0, 1.0;
      0.0, 0.0, 0.0];
n = size(before_filter, 1);
x_just_few_pxs = fminsearch(@(x) conv_transform_to_min(x), x0, optimset('MaxFunEvals', 10000));
% from maroon
x_just_few_pxs = [0.6540   -0.9136   -5.5430;
   -0.0147   -3.8008   -0.1802;
    0.0929   -0.5689   -0.1904;
    0.0132   -0.0277    0.1160];
(x_just_few_pxs(1:3,:)*before_filter')' + repmat(x_just_few_pxs(4,:), 4, 1)
after_trans = [-0.9191   -0.0251    0.0432;
   -0.9175   -0.0357    0.0562;
   -0.8827   -0.0467    0.0701;
   -0.8242   -0.0587    0.0844];