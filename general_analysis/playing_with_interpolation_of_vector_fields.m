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

rwof = img_wout_filter(:, :, 1); rwof = rwof(:); r_wout_filter = rwof(mask(:));
gwof = img_wout_filter(:, :, 2); gwof = gwof(:); g_wout_filter = gwof(mask(:));
bwof = img_wout_filter(:, :, 3); bwof = bwof(:); b_wout_filter = bwof(mask(:));

dkl_wout_filter = real(rgb2dkl(RGB2DKL_T, [r_wout_filter(:), g_wout_filter(:), b_wout_filter(:)]'))';

rwof = img_with_filter(:, :, 1); rwof = rwof(:); r_filter = rwof(mask(:));
gwof = img_with_filter(:, :, 2); gwof = gwof(:); g_filter = gwof(mask(:));
bwof = img_with_filter(:, :, 3); bwof = bwof(:); b_filter = bwof(mask(:));

dkl_filter = real(rgb2dkl(RGB2DKL_T, [r_filter(:), g_filter(:), b_filter(:)]'))';

vecs = dkl_filter - dkl_wout_filter;

ymin = min(dkl_wout_filter(:, 2)) - 0.01;
zmin = min(dkl_wout_filter(:, 3)) - 0.01;

ymax = max(dkl_wout_filter(:, 2)) + 0.01;
zmax = max(dkl_wout_filter(:, 3)) + 0.01;

steps = 300;
ys = linspace(ymin, ymax, steps);
zs = linspace(zmin, zmax, steps);
[ygr, zgr] = meshgrid(ys, zs);
rgs = squeeze(dkl_wout_filter(:, 2));
yvs = squeeze(dkl_wout_filter(:, 3));
vec_rgs = squeeze(vecs(:, 2));
vec_yvs = squeeze(vecs(:, 3));

rg_interp = scatteredInterpolant(rgs, yvs, vec_rgs);
yv_interp = scatteredInterpolant(rgs, yvs, vec_yvs);
yvgrid = yv_interp(ygr, zgr);
rggrid = rg_interp(ygr, zgr);

div = divergence(ygr, zgr, rggrid, yvgrid);

figure;
pcolor(ygr, zgr, div);
colorposneg;
colorbar;
shading interp;
export_fig divergence_of_interpolated_filter_vec_field.tiff

figure;
quiver(ygr, zgr, rggrid, yvgrid, 100);
axis square;
export_fig interpolated_filter_vec_field.tiff