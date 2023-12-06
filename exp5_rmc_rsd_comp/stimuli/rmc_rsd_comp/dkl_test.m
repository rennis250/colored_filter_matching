clear all; close all; clc;

RGB2LMS = rgb2lmsFromCalib('../../calibration/eizo_mon_spectra.csv');

in_mask = logical(imread('../../images/masks/obj_mask.png')); obj_mask = in_mask(:);
out_mask = ~in_mask; bkgd_mask = out_mask(:);

imgs = dir('scenes/*png');
rmc = zeros(size(imgs, 1), 3);
rsd = zeros(size(imgs, 1), 3);

for ic = 1:length(imgs)
    img = im2double(imread(['scenes/' imgs(ic).name]));
    imggc = img.^(2.2);
    rgb = reshape(imggc, size(img,1)*size(img,2), size(img,3));
    lms = rgb2lms(rgb', RGB2LMS);

    obj_lms = lms(:, obj_mask);
    bkgd_lms = lms(:, bkgd_mask);

    rmc(ic, :) = mean(obj_lms, 2) ./ mean(bkgd_lms, 2);
    rsd(ic, :) = std(obj_lms, 0, 2) ./ std(bkgd_lms, 0, 2);

    ic
end

figure(1);
clf;
subplot(1, 3, 1);
hold on;
plot(rmc(:, 1), rsd(:, 1), 'ko');
axis square;
subplot(1, 3, 2);
hold on;
plot(rmc(:, 2), rsd(:, 2), 'ko');
axis square;
subplot(1, 3, 3);
hold on;
plot(rmc(:, 3), rsd(:, 3), 'ko');
axis square;


