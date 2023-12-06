clear all; close all; clc;

img = im2double(imread('../images/patterned/mitsuba_caustics_rgb_trans_lightest_illum_yellow_refl_munsell_blue_lightest.png'));

obj_mask = imread('../images/masks/obj_mask.png'); obj_mask = logical(obj_mask);
high_mask = imread('../images/masks/high_mask.png'); high_mask = logical(high_mask);
refl_mask = imread('../images/masks/refl_mask.png'); refl_mask = logical(refl_mask);
obj_wout_high_mask = imread('../images/masks/obj_wout_high_mask.png'); obj_wout_high_mask = logical(obj_wout_high_mask);

r = img(:,:,1); r(obj_mask == 0) = 0.8;
g = img(:,:,2); g(obj_mask == 0) = 0.8;
b = img(:,:,3); b(obj_mask == 0) = 0.8;
img_obj = zeros(size(img));
img_obj(:,:,1) = r;
img_obj(:,:,2) = g;
img_obj(:,:,3) = b;
imwrite(img_obj, '../figures/obj.png');

r = img(:,:,1); r(high_mask == 0) = 0.8;
g = img(:,:,2); g(high_mask == 0) = 0.8;
b = img(:,:,3); b(high_mask == 0) = 0.8;
img_high = zeros(size(img));
img_high(:,:,1) = r;
img_high(:,:,2) = g;
img_high(:,:,3) = b;
imwrite(img_high, '../figures/high.png');

r = img(:,:,1); r(refl_mask == 0) = 0.8;
g = img(:,:,2); g(refl_mask == 0) = 0.8;
b = img(:,:,3); b(refl_mask == 0) = 0.8;
img_refl = zeros(size(img));
img_refl(:,:,1) = r;
img_refl(:,:,2) = g;
img_refl(:,:,3) = b;
imwrite(img_refl, '../figures/refl.png');

r = img(:,:,1); r(obj_wout_high_mask == 0) = 0.8;
g = img(:,:,2); g(obj_wout_high_mask == 0) = 0.8;
b = img(:,:,3); b(obj_wout_high_mask == 0) = 0.8;
img_obj_wout_high = zeros(size(img));
img_obj_wout_high(:,:,1) = r;
img_obj_wout_high(:,:,2) = g;
img_obj_wout_high(:,:,3) = b;
imwrite(img_obj_wout_high, '../figures/img_obj_wout_high.png');
