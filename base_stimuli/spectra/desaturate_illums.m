clear all; close all; clc;

red = dlmread('spectra/brighter_illum_redOrtho.spd');
green = dlmread('spectra/brighter_illum_greenOrtho.spd');
blue = dlmread('spectra/brighter_illum_blueOrtho.spd');
yellow = dlmread('spectra/brighter_illum_yellowOrtho.spd');

dered = zeros(size(red));
degreen = zeros(size(red));
deblue = zeros(size(red));
deyellow = zeros(size(red));

dered(:,1) = red(:,1);
degreen(:,1) = green(:,1);
deblue(:,1) = blue(:,1);
deyellow(:,1) = yellow(:,1);

dered(:,2) = 0.30.*green(:,2) + 0.70.*red(:,2);
degreen(:,2) = 0.70.*green(:,2) + 0.30.*red(:,2);
deblue(:,2) = 0.70.*blue(:,2) + 0.30.*yellow(:,2);
deyellow(:,2) = 0.30.*blue(:,2) + 0.70.*yellow(:,2);

dlmwrite('spectra/de_brighter_illum_redOrtho.spd', dered, ' ');
dlmwrite('spectra/de_brighter_illum_greenOrtho.spd', degreen, ' ');
dlmwrite('spectra/de_brighter_illum_blueOrtho.spd', deblue, ' ');
dlmwrite('spectra/de_brighter_illum_yellowOrtho.spd', deyellow, ' ');