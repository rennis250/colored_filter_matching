clear all; close all; clc;

img1 = imread('../images/exp5/mitsuba_comp_rmc_rsd_illum_3_rg_1_by_4_ld_2_bkgd_voronoi_red_brighter.png');
refl = dlmread('../exp5_rmc_rsd_comp/stimuli/rmc_rsd_comp/refls/rg_-1_by_0.9999998_ld_1.spd');
illum = dlmread('../exp5_rmc_rsd_comp/stimuli/rmc_rsd_comp/illums/red.spd');

figure(1);
clf;
subplot(1, 2, 1);
imshow(img1);
subplot(1, 2, 2);
hold on;
plot(refl(:, 1), refl(:, 2), 'r');
plot(illum(:, 1), illum(:, 2), 'k');
axis square
export_fig('../figures/purple_example.png');

img2 = imread('../images/exp5/mitsuba_comp_rmc_rsd_illum_3_rg_1_by_2_ld_2_bkgd_voronoi_green_brighter.png');
refl = dlmread('../exp5_rmc_rsd_comp/stimuli/rmc_rsd_comp/refls/rg_-1_by_-0.3333334_ld_1.spd');
illum = dlmread('../exp5_rmc_rsd_comp/stimuli/rmc_rsd_comp/illums/red.spd');

figure(2);
clf;
subplot(1, 2, 1);
imshow(img2);
subplot(1, 2, 2);
hold on;
plot(refl(:, 1), refl(:, 2), 'r');
plot(illum(:, 1), illum(:, 2), 'k');
axis square
export_fig('../figures/green_example.png');