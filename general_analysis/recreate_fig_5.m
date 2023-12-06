clear all; close all; clc;

gcRGB = csvread('../calibration/eizo_gamma.csv');
monxyY = csvread('../calibration/eizo_chroma.csv');
DKL2RGB = dkl2rgbFromCalib('../calibration/eizo_chroma.csv');
RGB2LMS = rgb2lmsFromCalib('../calibration/eizo_mon_spectra.csv');
monxyz = xyY2XYZ(monxyY);
RGB2DKL_T = inv(DKL2RGB);
LMS2RGB = inv(RGB2LMS);

mean3 = @(x)(mean(x, 1));

load('../data/obsstats.mat', 'data_table');

% mitsuba_caustics_rgb_trans_lightest_illum_blue_refl_munsell_yellow_lightest.png

obs_data = grpstats(data_table, {'exp_name', 'obs_name', 'mask_name', 'illum_glaven', 'body_glaven', 'hilo_glaven'}, 'mean');
obs_data = removevars(obs_data, {'obs_name'});
pop_data = grpstats(obs_data, {'exp_name', 'mask_name', 'illum_glaven', 'body_glaven', 'hilo_glaven'}, 'mean');

idxs = pop_data.exp_name == 'patch_dye' & ...
    pop_data.mask_name == 'obj_mask' & ...
    pop_data.illum_glaven == 'blue' & ...
    pop_data.body_glaven == 'yellow' & ...
    pop_data.hilo_glaven == '0';
patch_match = pop_data(idxs, :);
ld = patch_match.mean_mean_gwld_match;
rg = patch_match.mean_mean_gwrg_match;
yv = patch_match.mean_mean_gwyv_match;
round(255.*dkl2rgb(DKL2RGB, ld, rg, yv)'.^(1.0./2.2))

vor_texture = csvread('vor.csv');
lms_absorp_curves = csvread('linss2_10e_1.csv');
lms_absorp_curves = lms_absorp_curves(1:14:end, 2:end);

redtrans_file = 'munsell_red_better.spd';
greentrans_file = 'munsell_green_better.spd';
bluetrans_file = 'munsell_blue_better.spd';
yellowtrans_file = 'munsell_yellow_better.spd';

monitor_spectra = csvread('../calibration/eizo_mon_spectra.csv');

[vor_filter_map, filter_idxs, redtrans, bluetrans, greentrans, yellowtrans] = ...
        gen_voronoi_filter_map(monitor_spectra, redtrans_file, bluetrans_file, greentrans_file, yellowtrans_file);

% all values below are calculated in obsstats.m
ld = 0.6209;
rg_mix = 0.5889;
by_mix = 0.0831;

filter = ld.*(((rg_mix.*redtrans(:,2) + (1 - rg_mix).*greentrans(:,2)) + (by_mix.*bluetrans(:,2) + (1 - by_mix).*yellowtrans(:,2))));
vor_col2 = vor_filter_map;
vor_col2(filter_idxs, :) = bsxfun(@times, vor_filter_map(filter_idxs, :), filter'.^2.2);
vor_col2(isnan(vor_col2) | isinf(vor_col2)) = 0.0;
vor_rgb = (vor_col2*lms_absorp_curves*LMS2RGB').^(1.0/2.2);
% vor_rgb = vor_rgb.^(1/2.2);
imwrite(reshape(vor_rgb, 256, 256, 3), '../images/fig5_filter_match_example.png');

% and we will recreate fig 15 in here

% first bottom left
ld = 0.7053;
rg_mix = 0.5498;
by_mix = 0.2676;

filter = ld.*(((rg_mix.*redtrans(:,2) + (1 - rg_mix).*greentrans(:,2)) + (by_mix.*bluetrans(:,2) + (1 - by_mix).*yellowtrans(:,2))));
vor_col2 = vor_filter_map;
vor_col2(filter_idxs, :) = bsxfun(@times, vor_filter_map(filter_idxs, :), filter'.^2.2);
vor_col2(isnan(vor_col2) | isinf(vor_col2)) = 0.0;
vor_rgb = (vor_col2*lms_absorp_curves*LMS2RGB').^(1.0/2.2);
% vor_rgb = vor_rgb.^(1/2.2);
imwrite(reshape(vor_rgb, 256, 256, 3), '../images/fig15_bottom_left_filter_match_example.png');


% then bottom right
ld = 0.4476;
rg_mix = 0.5050;
by_mix = 0.2699;

filter = ld.*(((rg_mix.*redtrans(:,2) + (1 - rg_mix).*greentrans(:,2)) + (by_mix.*bluetrans(:,2) + (1 - by_mix).*yellowtrans(:,2))));
vor_col2 = vor_filter_map;
vor_col2(filter_idxs, :) = bsxfun(@times, vor_filter_map(filter_idxs, :), filter'.^2.2);
vor_col2(isnan(vor_col2) | isinf(vor_col2)) = 0.0;
vor_rgb = (vor_col2*lms_absorp_curves*LMS2RGB').^(1.0/2.2);
% vor_rgb = vor_rgb.^(1/2.2);
imwrite(reshape(vor_rgb, 256, 256, 3), '../images/fig15_bottom_right_filter_match_example.png');

