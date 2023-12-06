% clear all; close all; clc;

%% looking at the cone excitation ratios the correct way

rgb2lms = rgb2lmsFromCalib('mon_spectra.csv');
lms2rgb = inv(rgb2lms);

monitor_spectra = csvread('mon_spectra.csv');

lms = csvread('linss2_10e_1.csv');
wlns = lms(:,1);
lms = lms(1:14:end,2:end);

vor_texture = csvread('vor.csv');

whiteillum = max(max(monitor_spectra))*ones(length(wlns(1:14:end)),1);

bluetrans = dlmread('../spectra/munsell_blue_better.spd');
bluetranstmp = [];
bluetranstmp(:,1) = wlns(1:14:end);
bluetranstmp(:,2) = interp1(bluetrans(:,1), bluetrans(:,2), wlns(1:14:end));
bluetrans = bluetranstmp;

greentrans = dlmread('../spectra/munsell_green_better.spd');
greentranstmp = [];
greentranstmp(:,1) = wlns(1:14:end);
greentranstmp(:,2) = interp1(greentrans(:,1), greentrans(:,2), wlns(1:14:end));
greentrans = greentranstmp;

redtrans = dlmread('../spectra/munsell_red_better.spd');
redtranstmp = [];
redtranstmp(:,1) = wlns(1:14:end);
redtranstmp(:,2) = interp1(redtrans(:,1), redtrans(:,2), wlns(1:14:end));
redtrans = redtranstmp;

yellowtrans = dlmread('../spectra/munsell_yellow_better.spd');
yellowtranstmp = [];
yellowtranstmp(:,1) = wlns(1:14:end);
yellowtranstmp(:,2) = interp1(yellowtrans(:,1), yellowtrans(:,2), wlns(1:14:end));
yellowtrans = yellowtranstmp;

vor_col = vor_texture.*3;
vor_col = bsxfun(@times, vor_col, whiteillum');
vor_col2 = vor_col;

rg_mix = 0.2;
by_mix = 0.8;
ld = 0.5;

filter = ld.*(((rg_mix.*redtrans(:,2) + (1 - rg_mix).*greentrans(:,2)) + (by_mix.*bluetrans(:,2) + (1 - by_mix).*yellowtrans(:,2))));
vor_col3 = bsxfun(@times, vor_col, filter'.^2.2);
vor_col3(isnan(vor_col3) | isinf(vor_col3)) = 0.0;

figure();
clf;
subplot(1,3,1);
plot(vor_col2(:, 1), vor_col3(:, 1), 'ko');
subplot(1,3,2);
plot(vor_col2(:, 2), vor_col3(:, 2), 'ko');
subplot(1,3,3);
plot(vor_col2(:, 3), vor_col3(:, 3), 'ko');

l_filt_mean = mean(vor_col3(:,1));
m_filt_mean = mean(vor_col3(:,2));
s_filt_mean = mean(vor_col3(:,3));

l_vor_mean = mean(vor_col2(:,1));
m_vor_mean = mean(vor_col2(:,2));
s_vor_mean = mean(vor_col2(:,3));

%% let's look at the cone excitation ratios for one of the images

pm = imread('../images/objMask.png');
pm = pm(:,:,1); pm = logical(pm);
gisp = find(pm);

img = im2double(imread('../images/patterned/mitsuba_caustics_rgb_trans_lightest_illum_blue_refl_munsell_green_lightest.png'));
img = reshape(img, size(img,1)*size(img,2), 3);
lms_img = img*rgb2lms;

lms_obj = lms_img(gisp, :);
lms_bkg = lms_img(~pm, :);

l_obj_mean = mean(lms_obj(:,1));
m_obj_mean = mean(lms_obj(:,2));
s_obj_mean = mean(lms_obj(:,3));

l_bkg_mean = mean(lms_bkg(:,1));
m_bkg_mean = mean(lms_bkg(:,2));
s_bkg_mean = mean(lms_bkg(:,3));