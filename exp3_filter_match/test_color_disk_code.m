clear all; close all; clc;

monxyY = [.6804 .3073 30.94
    .2029 .6968 74.22
    .1527 .0508 6.74];

rgb2lms = rgb2lmsFromCalib('~/ownCloud/chromaticSceneStatistics/experiment/laptopTest/analysis/Labor313EXP_Eizo_spectra_05-Dec-2016230530_Konica.csv');
lms2rgb = inv(rgb2lms);

monitor_spectra = csvread('~/ownCloud/chromaticSceneStatistics/experiment/laptopTest/analysis/Labor313EXP_Eizo_spectra_05-Dec-2016230530_Konica.csv');

lms = csvread('linss2_10e_1.csv');
wlns = lms(:,1);
lms = lms(1:14:end,2:end);

redrefl = dlmread('spectra/munsell_red_lightest.spd');
redrefltmp = [];
redrefltmp(:,1) = wlns(1:14:end);
redrefltmp(:,2) = interp1(redrefl(:,1), redrefl(:,2), wlns(1:14:end));
redrefl = redrefltmp;

grayrefl = [];
grayrefl(:,1) = wlns(1:14:end);
grayrefl(:,2) = ones(length(wlns(1:14:end)),1);

greenspect = dlmread('spectra/illum_greenOrtho.spd');
greenspecttmp = [];
greenspecttmp(:,1) = wlns(1:14:end);
greenspecttmp(:,2) = interp1(greenspect(:,1), greenspect(:,2), wlns(1:14:end));
greenspect = greenspecttmp;

scaler = max(max(monitor_spectra))/max(greenspect(:,2));
greenspect = greenspect.*scaler;

bluetrans = dlmread('spectra/munsell_blue_lighter.spd');
bluetranstmp = [];
bluetranstmp(:,1) = wlns(1:14:end);
bluetranstmp(:,2) = interp1(bluetrans(:,1), bluetrans(:,2), wlns(1:14:end));
bluetrans = bluetranstmp;

transmittedspect = [];
transmittedspect(:,1) = wlns(1:14:end);
transmittedspect(:,2) = grayrefl(:,2).*greenspect(:,2).*bluetrans(:,2).^2.2;
transmittedspect(isnan(transmittedspect(:,2)) | isinf(transmittedspect(:,2)), 2) = 0.0;

lms_coords = lms'*transmittedspect(:,2);
rgb = lms2rgb*lms_coords;

vor_texture = csvread('../voronoi_for_transparency/vor.csv');

circ_mask = zeros(256,256);
for r = 1:60
    for theta = 0.001:0.001:2*pi;
        y = r*sin(theta);
        x = r*cos(theta);
        if sqrt(x^2 + y^2) < 256
            y = ceil(y) + 256/2;
            x = ceil(x) + 256/2;
            circ_mask(y, x) = 1.0;
        end
    end
end
idxs = find(circ_mask == 1.0);
[y, x] = ind2sub(size(circ_mask), find(circ_mask == 1.0));

whiteillum = max(max(monitor_spectra))*ones(length(wlns(1:14:end)),1);

vor_col = vor_texture;
vor_col = bsxfun(@times, vor_col, whiteillum');
vor_col(idxs, :) = bsxfun(@times, vor_col(idxs, :), bluetrans(:,2)'.^2.2);
vor_col(isnan(vor_col) | isinf(vor_col)) = 0.0;
vor_rgb = reshape(vor_col*lms*lms2rgb', 256, 256, 3);
imshow(imgaussfilt(vor_rgb.^(1/2.2).*2, 0.55));
