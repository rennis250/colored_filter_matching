clear all; close all; clc;

DKL2RGB = dkl2rgbFromCalib('../calibration/eizo_chroma.csv');

rgb_red = real(dkl2rgb(DKL2RGB, 0, 1, 0'))';
rgb_green = real(dkl2rgb(DKL2RGB, 0, -1, 0'))';

img = zeros(400, 400, 3);

% mix with 65% red coverage overlay
a = 0.5;
for x = 1:3
  % green bkgd
  img(:, :, x) = rgb_green(x);
end
for x = 1:3
  % red overlay
  img(150:250, 150:250, x) = a*rgb_red(x) + (1 - a)*rgb_green(x);
end

figure(1);
clf;
imshow(sqrt(img));
