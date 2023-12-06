clear all; close all; clc;

gcRGB = csvread('../calibration/eizo_gamma.csv');

DKL2RGB = dkl2rgbFromCalib('../calibration/eizo_chroma.csv');
RGB2DKL_T = inv(DKL2RGB);

glaven_for_plot = im2double(imread(['../images/patterned/mitsuba_caustics_rgb_trans_illum_blue_refl_munsell_red_lighter.png']));
empty_scene_for_plot = im2double(imread('../images/empty_scene_blue_illum.png'));

img_with_filter = gammaCorr(im2double(imread(['../images/patterned/mitsuba_caustics_rgb_trans_illum_blue_refl_munsell_red_lighter.png'])), gcRGB);
img_wout_filter = gammaCorr(im2double(imread('../images/empty_scene_blue_illum.png')), gcRGB);

mask = imread('../images/masks/obj_mask_correct.png');
mask = squeeze(mask(:, :, 1));
mask = logical(mask);

rwof = img_wout_filter(:, :, 1); rwof = rwof(:); r_wout_filter = rwof(mask(:));
gwof = img_wout_filter(:, :, 2); gwof = gwof(:); g_wout_filter = gwof(mask(:));
bwof = img_wout_filter(:, :, 3); bwof = bwof(:); b_wout_filter = bwof(mask(:));

dkl_wout_filter = real(rgb2dkl(RGB2DKL_T, [r_wout_filter(:), g_wout_filter(:), b_wout_filter(:)]'))';

r = img_with_filter(:, :, 1); r = r(:); r_img = r(mask(:));
g = img_with_filter(:, :, 2); g = g(:); g_img = g(mask(:));
b = img_with_filter(:, :, 3); b = b(:); b_img = b(mask(:));

dkl_img = real(rgb2dkl(RGB2DKL_T, [r_img(:), g_img(:), b_img(:)]'))';

[ldgr_d1, ldgr_d3, rggr_d1, rggr_d2, yvgr_d2, yvgr_d3, ld_interp_d1, ld_interp_d3, rg_interp_d1, rg_interp_d2, yv_interp_d2, yv_interp_d3] = divergence_map(dkl_wout_filter, dkl_img);

%%%%%%%%%%%%
% compute best affine map for full DKL data
global Bfit_affinemap;
compute_forms_of_convergence(dkl_wout_filter, dkl_img, false);

h = figure(1);
set(gcf, 'color', 'w');
fullfig(h);
clf;
subplot(2, 3, 1);
imshow(empty_scene_for_plot);
title('A');
axis square;

subplot(2, 3, 2);
imshow(glaven_for_plot);
title('B');
axis square;

% example result for best fitting affine map
subplot(2, 3, 3);
hold off;
hist3([dkl_img(:, 2), dkl_img(:, 3)], 'FaceColor', 'b', 'FaceAlpha', 0.65, 'Ctrs', ...
      {linspace(-0.2, 1.5, 20), linspace(-0.8, 1.8, 20)});
hold on;
hist3([Bfit_affinemap(:, 2), Bfit_affinemap(:, 3)], 'FaceColor', 'r', 'FaceAlpha', 0.65, 'Ctrs', ...
      {linspace(-0.2, 1.5, 20), linspace(-0.8, 1.8, 20)});
xlabel('L-M');
ylabel('S-(L+M)');
zlabel('Number of pixels');
title('C');
axis square;

subplot(2, 3, 4);
quiver(rggr_d2, yvgr_d2, rg_interp_d2, yv_interp_d2, 10);
axis([-0.5, 1.2, -1, 2]);
xlabel('L-M');
ylabel('S-(L+M)');
title('D');
axis square;

subplot(2, 3, 5);
hold on;
streamslice(rggr_d2, yvgr_d2, rg_interp_d2, yv_interp_d2);
plot(0.05, 0, 'rs', 'MarkerSize', 8, 'LineWidth', 3);
axis([-0.5, 1.2, -1, 1.6]);
xlabel('L-M');
ylabel('S-(L+M)');
title('E');
axis square;

subplot(2, 3, 6);
hold on;
streamslice(yvgr_d3, ldgr_d3, yv_interp_d3, ld_interp_d3, 8, 'method', 'nearest');
plot(0, -0.92, 'rs', 'MarkerSize', 8, 'LineWidth', 3);
axis([-1, 2, -1.1, 0.6]);
xlabel('S-(L+M)');
ylabel('L+M+S');
title('F');
axis square;

export_fig(['../figures/fig1.png'], '-m2', '-r500');

%% make figure for methods section
h2 = figure(2);
set(gcf, 'color', 'w');
fullfig(h2);
clf;
subplot(2, 2, 1);
imshow(empty_scene_for_plot);
title('A');
axis square;

subplot(2, 2, 2);
imshow(glaven_for_plot);
title('B');
axis square;

subplot(2, 2, 3);
hold on;
vecs = dkl_img - dkl_wout_filter;
quiver(dkl_wout_filter(:, 2), dkl_wout_filter(:, 3), vecs(:, 2), vecs(:, 3), 10);
plot(0.465, 0.375, 'rs', 'MarkerSize', 8, 'LineWidth', 3);
axis([-0.5, 1.2, -1, 2]);
xlabel('L-M');
ylabel('S-(L+M)');
title('C');
axis square;

subplot(2, 2, 4);
quiver(dkl_wout_filter(:, 2), dkl_wout_filter(:, 3), vecs(:, 2), vecs(:, 3), 5, 'MaxHeadSize', 0.05);
axis([0.45, 0.48, 0.36, 0.39]);
xlabel('L-M');
ylabel('S-(L+M)');
title('D');
axis square;

export_fig(['../figures/figA1.png'], '-m2', '-r500');
