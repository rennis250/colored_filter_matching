clear all; close all; clc;

% first, do color filter stimulus

monxyY = csvread('../calibration/eizo_chroma.csv');
monxyz = xyY2XYZ(monxyY);

wd = 400;
he = 400;

rgb2lms = rgb2lmsFromCalib('../calibration/eizo_mon_spectra.csv');
lms2rgb = inv(rgb2lms);

monitor_spectra = csvread('../calibration/eizo_mon_spectra.csv');

lms = csvread('../aux/linss2_10e_1.csv');
wlns = lms(:,1);
lms = lms(1:14:end,2:end);

vor_texture = csvread('../exp3_filter_match/vor.csv');

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
circ_idxs = find(circ_mask == 1.0);
[y, x] = ind2sub(size(circ_mask), find(circ_mask == 1.0));

% by_mult = 1;
% bluetrans = dlmread('spectra/munsell_blue_better.spd');
% bluetranstmp = [];
% bluetranstmp(:,1) = wlns(1:14:end);
% bluetranstmp(:,2) = interp1(bluetrans(:,1), bluetrans(:,2), wlns(1:14:end));
% bluetrans = bluetranstmp;

% greentrans = dlmread('spectra/munsell_green_better.spd');
% greentranstmp = [];
% greentranstmp(:,1) = wlns(1:14:end);
% greentranstmp(:,2) = interp1(greentrans(:,1), greentrans(:,2), wlns(1:14:end));
% greentrans = greentranstmp;

% redtrans = dlmread('spectra/munsell_red_better.spd');
% redtranstmp = [];
% redtranstmp(:,1) = wlns(1:14:end);
% redtranstmp(:,2) = interp1(redtrans(:,1), redtrans(:,2), wlns(1:14:end));
% redtrans = redtranstmp;

% yellowtrans = dlmread('spectra/munsell_yellow_better.spd');
% yellowtranstmp = [];
% yellowtranstmp(:,1) = wlns(1:14:end);
% yellowtranstmp(:,2) = interp1(yellowtrans(:,1), yellowtrans(:,2), wlns(1:14:end));
% yellowtrans = yellowtranstmp;

by_mult = 3;
bluetrans = dlmread('../base_stimuli/spectra/refl_blueBall.spd');
bluetranstmp = [];
bluetranstmp(:,1) = wlns(1:14:end);
bluetranstmp(:,2) = interp1(bluetrans(:,1), bluetrans(:,2), wlns(1:14:end));
bluetrans = bluetranstmp;

greentrans = dlmread('../base_stimuli/spectra/refl_greenBall.spd');
greentranstmp = [];
greentranstmp(:,1) = wlns(1:14:end);
greentranstmp(:,2) = interp1(greentrans(:,1), greentrans(:,2), wlns(1:14:end));
greentrans = greentranstmp;

redtrans = dlmread('../base_stimuli/spectra/refl_redBall.spd');
redtranstmp = [];
redtranstmp(:,1) = wlns(1:14:end);
redtranstmp(:,2) = interp1(redtrans(:,1), redtrans(:,2), wlns(1:14:end));
redtrans = redtranstmp;

yellowtrans = dlmread('../base_stimuli/spectra/refl_yellowBall.spd');
yellowtranstmp = [];
yellowtranstmp(:,1) = wlns(1:14:end);
yellowtranstmp(:,2) = interp1(yellowtrans(:,1), yellowtrans(:,2), wlns(1:14:end));
yellowtrans = yellowtranstmp;

whiteillum = max(max(monitor_spectra))*ones(length(wlns(1:14:end)),1);

vor_col = vor_texture.*3;
vor_col = bsxfun(@times, vor_col, whiteillum');
vor_col2 = vor_col;

rg_mix = 0.6126;
by_mix = 0.0077;
ld = 0.5273;

filter = ld.*(((rg_mix.*redtrans(:,2) + (1 - rg_mix).*greentrans(:,2)) + (by_mult.*by_mix.*bluetrans(:,2) + (1 - by_mix).*yellowtrans(:,2))));
vor_col2(circ_idxs, :) = bsxfun(@times, vor_col(circ_idxs, :), filter'.^2.2);
vor_col2(isnan(vor_col2) | isinf(vor_col2)) = 0.0;
vor_image = real(vor_col2*lms*lms2rgb');
vor_image = reshape(vor_image.^(1/2.2), 256, 256, 3);
imshow(vor_image)

% next, do square patch stimulus (extract from match data)

r = 0.8574; % 219
g = 0.8106; % 207
b = 0.5747; % 147
