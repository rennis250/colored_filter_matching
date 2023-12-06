% try stream slice on new divergence format
clear all; close all; clc;

gcRGB = csvread('../calibration/eizo_gamma.csv');

DKL2RGB = dkl2rgbFromCalib('../calibration/eizo_chroma.csv');
RGB2DKL_T = inv(DKL2RGB);

blue_img_wout_filter = gammaCorr(im2double(imread('../images/empty_scene_blue_illum.png')), gcRGB);

masks = {'../images/masks/obj_mask.png'};

img_with_filter = gammaCorr(im2double(imread(['../images/patterned/mitsuba_caustics_rgb_trans_illum_blue_refl_munsell_red_lighter.png'])), gcRGB);
img_wout_filter = blue_img_wout_filter;

mc = 1;
mask = imread(masks{mc}); mask = logical(mask);

r = img_with_filter(:,:,1); r = r(:); r_filter = r(mask(:));
g = img_with_filter(:,:,2); g = g(:); g_filter = g(mask(:));
b = img_with_filter(:,:,3); b = b(:); b_filter = b(mask(:));

dkl_filter = real(rgb2dkl(RGB2DKL_T, [r_filter(:), g_filter(:), b_filter(:)]'))';

rwof = img_wout_filter(:,:,1); rwof = rwof(:); r_wout_filter = rwof(mask(:));
gwof = img_wout_filter(:,:,2); gwof = gwof(:); g_wout_filter = gwof(mask(:));
bwof = img_wout_filter(:,:,3); bwof = bwof(:); b_wout_filter = bwof(mask(:));

dkl_wout_filter = real(rgb2dkl(RGB2DKL_T, [r_wout_filter(:), g_wout_filter(:), b_wout_filter(:)]'))';

[xgr, xgr_d3, ygr_d1, ygr_d2, zgr, zgr_d3, ld_interp, ld_interp_d3, rg_interp_d1, rg_interp_d2, yv_interp, yv_interp_d3] = divergence_map(dkl_wout_filter, dkl_filter);
figure(1);
clf;
subplot(1,3,1);
streamslice(xgr, ygr_d1, ld_interp, rg_interp_d1, 3);
axis square;
subplot(1,3,2);
streamslice(ygr_d2, zgr, rg_interp_d2, yv_interp);
axis square;
subplot(1,3,3);
streamslice(xgr_d3, zgr_d3, ld_interp_d3, yv_interp_d3, 3);
axis square;

fh = figure(2);
clf;
subplot(1,3,1);
quiver(ygr_d1, xgr, rg_interp_d1, ld_interp, 5);
axis([-1.1, 1, 0, 1]);
axis square;
subplot(1,3,2);
quiver(ygr_d2, zgr, rg_interp_d2, yv_interp, 5);
axis square;
subplot(1,3,3);
quiver(xgr_d3, zgr_d3, ld_interp_d3, yv_interp_d3, 5);
axis square;
axis([-1.1, 1, 0, 1]);
fh.Position(3) = fh.Position(4);
set(gcf, 'units', 'inches', 'position', 1.5*[0, 2.3, 8, 2.3]);
export_fig example_different_planes_of_filter_vec_field.tiff

div = divergence(ygr_d2, zgr, rg_interp_d2, yv_interp);
figure(3);
pcolor(ygr_d2, zgr, div);
shading interp;
colorbar;
colorposneg;
axis square;
