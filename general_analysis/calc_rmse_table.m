% clear all; close all; clc;

% load('output/divergence/rmse_vals.mat');

% blc, ic, bc, mc, darkerc, labc, chromc
% t = squeeze(rmse.affinemap(1,:,:,:,:,:,:));
% ic, bc, mc, darkerc, labc, chromc
% m = squeeze(mean(t, [1, 2, 4, 6]));
% s = squeeze(std(t, 0, [1, 2, 4, 6]));
% mc, labc

% now do it for no or little refraction

clear rmse;

load('output/divergence/rmse_vals.mat');

% blc, ic, bc, mc, darkerc, labc, chromc, refc
t = squeeze(rmse.affinemap(1,:,:,:,:,:,:,:));
% ic, bc, mc, labc, chromc, refc
m = squeeze(mean(t, [1, 2, 5]));
s = squeeze(std(t, 0, [1, 2, 5]));
% mc, labc, refc
