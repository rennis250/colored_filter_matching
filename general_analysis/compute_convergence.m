clear all; close all; clc;

img = imread('../images/empty_scene.png');
img_w_filter = imread('../images/patterned/mitsuba_caustics_rgb_trans_lightest_illum_blue_refl_munsell_blue_lightest.png');

lab_vals = csvread('../images/patterned/mitsuba_caustics_rgb_trans_lightest_illum_blue_refl_munsell_blue_lightest_lab_coords.csv');

lab_filter = lab_vals(:,1:3);
lab_wout_filter = lab_vals(:,4:6);

lab_coords = lab_wout_filter;
lab_vecs = lab_filter - lab_wout_filter;

div = divergence(lab_coords(:,2), lab_coords(:,3), lab_vecs(:,2), lab_vecs(:,3));
