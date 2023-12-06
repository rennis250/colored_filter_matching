clear all; close all; clc;

load('lms_prep.mat');
t = [data(:,1).spectralData(:) data(:,2).spectralData(:) data(:,3).spectralData(:)];
csvwrite('mon_spectra.csv', t);