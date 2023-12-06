clear all; close all; clc;

monxyY = [0.6788 0.3075 32.5020
        0.2062 0.6941 76.9090
        0.1535 0.0522 7.1828];
monxyz = xyY2XYZ(monxyY);

pm = imread('../../../../images/masks/obj_mask.png');
pm = logical(pm);
gisp = find(pm);
pm_wout_high = pm;

bm = imread('../../../../images/masks/high_mask.png');
bm = logical(bm);
gisb = find(bm);

pm_wout_high(gisb) = logical(0.0); % don't include highlight in estimate of color
gisp_wout_high = find(pm_wout_high);

sm = logical(ones(size(pm,1),size(pm,2)));

xc = 121; yc = 129;
backx = 212; backy = 61;
shadx = 600; shady = 281;

gR = 2.2; gG = 2.2; gB = 2.2;

test_img = im2double(imread('mitsuba_caustics_rgb_trans_illum_blue_refl_munsell_blue_lighter.png'));

[imp, imh] = maskAndLinearTransparent(test_img, pm, bm, sm, xc, yc, backx, backy, shadx, shady, gR, gG, gB);
[imp_wout_high, ~] = maskAndLinearTransparent(test_img, pm_wout_high, bm, sm, xc, yc, backx, backy, shadx, shady, gR, gG, gB);

rp = imp(:,:,1); rp = rp(:);
gp = imp(:,:,2); gp = gp(:);
bp = imp(:,:,3); bp = bp(:);
obj_rgb = [rp(gisp) gp(gisp) bp(gisp)];
obj_lab = rgb2labRob(squeeze(obj_rgb), monxyz);
mean_obj_rgb = mean(squeeze(obj_rgb), 1);
mean_obj_lab = mean(squeeze(obj_lab), 1);

rp = imp_wout_high(:,:,1); rp = rp(:);
gp = imp_wout_high(:,:,2); gp = gp(:);
bp = imp_wout_high(:,:,3); bp = bp(:);
obj_rgb_wout_high = [rp(gisp_wout_high) gp(gisp_wout_high) bp(gisp_wout_high)];
obj_lab_wout_high = rgb2labRob(squeeze(obj_rgb_wout_high), monxyz);
mean_obj_rgb_wout_high = mean(squeeze(obj_rgb_wout_high), 1);
mean_obj_lab_wout_high = mean(squeeze(obj_lab_wout_high), 1);

rh = imh(:,:,1); rh = rh(:);
gh = imh(:,:,2); gh = gh(:);
bh = imh(:,:,3); bh = bh(:);
high_rgb = [rh(gisb) gh(gisb) bh(gisb)];
high_lab = rgb2labRob(squeeze(high_rgb), monxyz);
mean_high_rgb = mean(squeeze(high_rgb), 1);
mean_high_lab = mean(squeeze(high_lab), 1);

rust_stats = readtable('mitsuba_caustics_rgb_trans_illum_blue_refl_munsell_blue_lighter.png_masked.csv');
rust_obj_stats = rust_stats(1, :);
rust_high_stats = rust_stats(2, :);
rust_obj_wout_high_stats = rust_stats(4, :);

rgb = mean_high_rgb;

rgbScaled = 2.0.*(rgb - 0.5);

T= [1  1  0.1029
    1 -0.4357  -0.1447
    1  0.0285  1];

Idkl = inv(T)*rgbScaled';
