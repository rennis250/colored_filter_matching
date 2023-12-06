clear all; close all; clc;

expns = {'11b'};
white_walls_filter_exp = 1;
filter_exps = [white_walls_filter_exp];

% masks = {'../images/masks/obj_mask_correct.png'};
masks = {'../images/masks/white_walls_mask.png'};
mask_names = {'obj_mask'};
stats = [struct()]; % exps x imgs x masks

expc = 1;
expn = expns{expc};
imgns = dir(['../images/exp' expn '/*.png']);

gcRGB = csvread('../calibration/oled_gamma.csv');
RGB2LMS = rgb2lmsFromCalib('../calibration/oled_mon_spectra.csv');
LMS2RGB = inv(RGB2LMS);

imgc = 1;
imgn = imgns(imgc);
img = gammaCorr(im2double(imread(['../images/exp' expn '/' imgn.name])), gcRGB);

maskc = 1;
mask = imread(masks{maskc});
mask_r = mask(:, :, 1);
mask_g = mask(:, :, 2);
mask_b = mask(:, :, 3);
bkgd_pxs = mask_r == 0 & mask_g == 0 & mask_b == 0;
bkgd_pxs = logical(bkgd_pxs(:));

r = img(:, :, 1); rs = r(:);
g = img(:, :, 2); gs = g(:);
b = img(:, :, 3); bs = b(:);

rgb_gc_lin = [rs(bkgd_pxs), gs(bkgd_pxs), bs(bkgd_pxs)];
lms = real(rgb2lms(RGB2LMS, rgb_gc_lin'))';

vor_img = gammaCorr(im2double(imread('voronoi_filter_ld_1.3190000_rg_0.6110000_by_0.5880000.png')), gcRGB);
vor_rgb = reshape(vor_img, size(vor_img, 1)*size(vor_img, 2), 3);
vor_lms = real(rgb2lms(RGB2LMS, vor_rgb'))';

disp('lms');
disp(mean(lms, 1));
disp(mean(vor_lms, 1));

disp('rgb');
disp(mean(rgb_gc_lin, 1));
disp(mean(vor_rgb, 1));

figure(1);
clf;
subplot(1, 2, 1);
imshow(sqrt(img));
subplot(1, 2, 2);
imshow(reshape(sqrt(vor_rgb), 256, 256, 3));

disp(mean(mean(lms, 1)./mean(vor_lms, 1)));
