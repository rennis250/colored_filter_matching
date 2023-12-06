clear all; close all; clc;

monxyY = [.6804 .3073 30.94
    .2029 .6968 74.22
    .1527 .0508 6.74];
monxyz = xyY2XYZ(monxyY);

% extract from patch of plane

% xs = 288; xe = 328;
% ys = 166; ye = 206;

% blue_wall = im2double(imread('images/gray_wall_illum_extraction__illum_blue_rotated_plane.png')).^(2.2);
% red_wall = im2double(imread('images/gray_wall_illum_extraction__illum_red_rotated_plane.png')).^(2.2);
% % green_wall = im2double(imread('images/gray_wall_illum_extraction__illum_green_rotated_plane.png')).^(2.2);
% yellow_wall = im2double(imread('images/gray_wall_illum_extraction__illum_yellow_rotated_plane.png')).^(2.2);

% blue_patch = blue_wall(ys:ye, xs:xe, :); blue_patch = reshape(blue_patch, size(blue_patch, 1)*size(blue_patch, 2), 3);
% red_patch = red_wall(ys:ye, xs:xe, :); red_patch = reshape(red_patch, size(red_patch, 1)*size(red_patch, 2), 3);
% % green_patch = green_wall(ys:ye, xs:xe, :); green_patch = reshape(green_patch, size(green_patch, 1)*size(green_patch, 2), 3);
% yellow_patch = yellow_wall(ys:ye, xs:xe, :); yellow_patch = reshape(yellow_patch, size(yellow_patch, 1)*size(yellow_patch, 2), 3);

% blue_lab = rgb2labRob(blue_patch, monxyz);
% % green_lab = rgb2labRob(green_patch, monxyz);
% red_lab = rgb2labRob(red_patch, monxyz);
% yellow_lab = rgb2labRob(yellow_patch, monxyz);

% blue_illum = mean(blue_lab, 1);
% red_illum = mean(red_lab, 1);
% % green_illum = mean(green_lab, 1);
% yellow_illum = mean(yellow_lab, 1);

% extract from illuminant in scene

% xs = 227; xe = 248;
% ys = 198; ye = 219;

% blue_wall = im2double(imread('images/illum_center_extraction__illum_blue_rotated_plane.png')).^(2.2);
% red_wall = im2double(imread('images/illum_center_extraction__illum_red_rotated_plane.png')).^(2.2);
% % green_wall = im2double(imread('images/illum_center_extraction__illum_green_rotated_plane.png')).^(2.2);
% yellow_wall = im2double(imread('images/illum_center_extraction__illum_yellow_rotated_plane.png')).^(2.2);

% blue_patch = blue_wall(ys:ye, xs:xe, :); blue_patch = reshape(blue_patch, size(blue_patch, 1)*size(blue_patch, 2), 3);
% red_patch = red_wall(ys:ye, xs:xe, :); red_patch = reshape(red_patch, size(red_patch, 1)*size(red_patch, 2), 3);
% % green_patch = green_wall(ys:ye, xs:xe, :); green_patch = reshape(green_patch, size(green_patch, 1)*size(green_patch, 2), 3);
% yellow_patch = yellow_wall(ys:ye, xs:xe, :); yellow_patch = reshape(yellow_patch, size(yellow_patch, 1)*size(yellow_patch, 2), 3);

% blue_lab = rgb2labRob(blue_patch, monxyz);
% % green_lab = rgb2labRob(green_patch, monxyz);
% red_lab = rgb2labRob(red_patch, monxyz);
% yellow_lab = rgb2labRob(yellow_patch, monxyz);

% blue_illum = mean(blue_lab, 1);
% red_illum = mean(red_lab, 1);
% % green_illum = mean(green_lab, 1);
% yellow_illum = mean(yellow_lab, 1);

% save('illums_lab.mat', 'blue_illum', 'red_illum', 'green_illum', 'yellow_illum');

% extract from patch of plane in voronoi scene

xs = 288; xe = 328;
ys = 166; ye = 206;

blue_wall = im2double(imread('images/vor_wall_illum_extraction__illum_blue_rotated_plane.png')).^(2.2);
red_wall = im2double(imread('images/vor_wall_illum_extraction__illum_red_rotated_plane.png')).^(2.2);
% green_wall = im2double(imread('images/vor_wall_illum_extraction__illum_green_rotated_plane.png')).^(2.2);
yellow_wall = im2double(imread('images/vor_wall_illum_extraction__illum_yellow_rotated_plane.png')).^(2.2);

blue_patch = blue_wall(ys:ye, xs:xe, :); blue_patch = reshape(blue_patch, size(blue_patch, 1)*size(blue_patch, 2), 3);
red_patch = red_wall(ys:ye, xs:xe, :); red_patch = reshape(red_patch, size(red_patch, 1)*size(red_patch, 2), 3);
% green_patch = green_wall(ys:ye, xs:xe, :); green_patch = reshape(green_patch, size(green_patch, 1)*size(green_patch, 2), 3);
yellow_patch = yellow_wall(ys:ye, xs:xe, :); yellow_patch = reshape(yellow_patch, size(yellow_patch, 1)*size(yellow_patch, 2), 3);

blue_lab = rgb2labRob(blue_patch, monxyz);
% green_lab = rgb2labRob(green_patch, monxyz);
red_lab = rgb2labRob(red_patch, monxyz);
yellow_lab = rgb2labRob(yellow_patch, monxyz);

blue_illum = mean(blue_lab, 1);
red_illum = mean(red_lab, 1);
% green_illum = mean(green_lab, 1);
yellow_illum = mean(yellow_lab, 1);

save('illums_lab.mat', 'blue_illum', 'red_illum', 'yellow_illum');
