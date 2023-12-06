clear all; close all; clc;

load munsell380_800_1.mat

% CIE1931 XYZ CMFs
cmfData = csvread('ciexyz31_1.csv');
wavelength_cmf = cmfData(:,1);
x_bar = cmfData(:,2);
y_bar = cmfData(:,3);
z_bar = cmfData(:,4);

wlns = 380:800;
cmf = interp1(wavelength_cmf(20:422), [x_bar(20:422) y_bar(20:422) z_bar(20:422)], wlns, 'spline');
cmf(cmf < 0) = 0;
cmf(isna(cmf)) = 0;

XYZ_CMFs = [cmf(:,1) cmf(:,2) cmf(:,3)];
xyz = XYZ_CMFs' * munsell;
xyY = [xyz(1,:)./sum(xyz,1); xyz(2,:)./sum(xyz,1); xyz(2,:)];

rgb = xyz2rgbOLED(xyz);

% most 'transparent' points from munsell
% find green lightest point
green_lightestidx = find(xyY(1,:) > 0.3172 & xyY(1,:) < 0.3174 & xyY(2,:) > 0.3483 & xyY(2,:) < 0.3485)
disp('green munsell refl')
S(green_lightestidx,:) % 5G 2.5/2

% find red lightest point
red_lightestidx = find(xyY(1,:) > 0.3510 & xyY(1,:) < 0.3512 & xyY(2,:) > 0.3332 & xyY(2,:) < 0.3334)
disp('red munsell refl')
S(red_lightestidx,:) % 2.5R 7/2

% find blue lightest point
blue_lightestidx = find(xyY(1,:) > 0.3164 & xyY(1,:) < 0.3166 & xyY(2,:) > 0.3271 & xyY(2,:) < 0.3273)
disp('blue munsell refl')
S(blue_lightestidx,:) % 5B 4/1

% find yellow lightest point
yellow_lightestidx = find(xyY(1,:) > 0.3480 & xyY(1,:) < 0.3482 & xyY(2,:) > 0.3542 & xyY(2,:) < 0.3544)
disp('green munsell refl')
S(yellow_lightestidx,:) % 10Y 5/1

figure(1);
clf;
hold on;
axis([0 1 0 1]);
scatter(xyY(1,:), xyY(2,:), 10, rgb');

lasso = csvread('xyYcoords.csv');
plot(lasso(:,1), lasso(:,2), 'k-');

plank = csvread('planck.csv');
plot(plank(:,1), plank(:,2), 'r-', 'LineWidth', 2);

line([xyY(1,green_lightestidx) xyY(1,red_lightestidx)], [xyY(2,green_lightestidx) xyY(2,red_lightestidx)], 'Color', [0.0 1.0 0.0]);

line([xyY(1,blue_lightestidx) xyY(1,yellow_lightestidx)], [xyY(2,blue_lightestidx) xyY(2,yellow_lightestidx)], 'Color', [1.0 0.0 1.0]);

% less "transparent" points for munsell
% find green point
greenidx = find(xyY(1,:) > 0.2802 & xyY(1,:) < 0.2804 & xyY(2,:) > 0.3694 & xyY(2,:) < 0.3696)

% find red point
redidx = find(xyY(1,:) > 0.3931 & xyY(1,:) < 0.3933 & xyY(2,:) > 0.3194 & xyY(2,:) < 0.3196)

% find blue point
blueidx = find(xyY(1,:) > 0.2898 & xyY(1,:) < 0.2900 & xyY(2,:) > 0.3058 & xyY(2,:) < 0.3060)

% find yellow point
yellowidx = find(xyY(1,:) > 0.3692 & xyY(1,:) < 0.3694 & xyY(2,:) > 0.3764 & xyY(2,:) < 0.3766)

figure(2);
clf;
hold on;
plot(munsell(:,green_lightestidx)./max(munsell(:,green_lightestidx)), 'g-');
plot(munsell(:,red_lightestidx)./max(munsell(:,red_lightestidx)), 'r-');
plot(munsell(:,blue_lightestidx)./max(munsell(:,blue_lightestidx)), 'b-');
plot(munsell(:,yellow_lightestidx)./max(munsell(:,yellow_lightestidx)), 'y-');

% write out for mitsuba
green_lightest = [wlns(:) munsell(:,green_lightestidx)./max(munsell(:,green_lightestidx))];
dlmwrite("munsell_green_lightest.spd", green_lightest, 'delimiter', ' ')
red_lightest = [wlns(:) munsell(:,red_lightestidx)./max(munsell(:,red_lightestidx))];
dlmwrite("munsell_red_lightest.spd", red_lightest, 'delimiter', ' ')
blue_lightest = [wlns(:)  munsell(:,blue_lightestidx)./max(munsell(:,blue_lightestidx))];
dlmwrite("munsell_blue_lightest.spd", blue_lightest, 'delimiter', ' ')
yellow_lightest = [wlns(:) munsell(:,yellow_lightestidx)./max(munsell(:,yellow_lightestidx))];
dlmwrite("munsell_yellow_lightest.spd", yellow_lightest, 'delimiter', ' ')

% write out scaled versions for mitsuba
green_lighter = [wlns(:) 0.68.*munsell(:,green_lightestidx)./max(munsell(:,green_lightestidx))];
dlmwrite("munsell_green_lighter.spd", green_lighter, 'delimiter', ' ')
red_lighter = [wlns(:) 0.68.*munsell(:,red_lightestidx)./max(munsell(:,red_lightestidx))];
dlmwrite("munsell_red_lighter.spd", red_lighter, 'delimiter', ' ')
blue_lighter = [wlns(:) 0.68.* munsell(:,blue_lightestidx)./max(munsell(:,blue_lightestidx))];
dlmwrite("munsell_blue_lighter.spd", blue_lighter, 'delimiter', ' ')
yellow_lighter = [wlns(:) 0.68.*munsell(:,yellow_lightestidx)./max(munsell(:,yellow_lightestidx))];
dlmwrite("munsell_yellow_lighter.spd", yellow_lighter, 'delimiter', ' ')

figure(3);
clf;
hold on;
plot(green_lighter(:,2), 'g-');
plot(red_lighter(:,2), 'r-');
plot(blue_lighter(:,2), 'b-');
plot(yellow_lighter(:,2), 'y-');
axis([0 421 0 1])

% new blue
blueidx = find(xyY(1,:) > 0.2836 & xyY(1,:) < 0.2838 & xyY(2,:) > 0.2982 & xyY(2,:) < 0.2984)
disp('blue filter match chip')
S(blueidx, :) % 7.5B 2.5/2
blue_better = [wlns(:)  munsell(:,blueidx)./max(munsell(:,blueidx))];
dlmwrite("spectra/munsell_blue_better.spd", blue_better, 'delimiter', ' ')

redidx = find(xyY(1,:) > 0.4103 & xyY(1,:) < 0.4105 & xyY(2,:) > 0.3104 & xyY(2,:) < 0.3106)
disp('red filter match chip')
S(redidx, :) % 7.5RP 6/8
red_better = [wlns(:)  munsell(:,redidx)./max(munsell(:,redidx))];
dlmwrite("spectra/munsell_red_better.spd", red_better, 'delimiter', ' ')

greenidx = find(xyY(1,:) > 0.2653 & xyY(1,:) < 0.2655 & xyY(2,:) > 0.3762 & xyY(2,:) < 0.3764)
disp('green filter match chip')
S(greenidx, :) % 10G 7/8
green_better = [wlns(:)  munsell(:,greenidx)./max(munsell(:,greenidx))];
dlmwrite("spectra/munsell_green_better.spd", green_better, 'delimiter', ' ')

yellowidx = find(xyY(1,:) > 0.3653 & xyY(1,:) < 0.3655 & xyY(2,:) > 0.3734 & xyY(2,:) < 0.3736)
disp('yellow filter match chip')
S(yellowidx, :) % 7.5Y 5/2
yellow_better = [wlns(:)  munsell(:,yellowidx)./max(munsell(:,yellowidx))];
dlmwrite("spectra/munsell_yellow_better.spd", yellow_better, 'delimiter', ' ')

figure(5);
clf;
hold on;
plot(blue_better(:,2), 'b-');
plot(red_better(:,2), 'r-');
plot(green_better(:,2), 'g-');
plot(yellow_better(:,2), 'k-');
axis([0 421 0 1])

figure(1);
hold on;
plot(0.2837, 0.2983, 'ks', 'MarkerSize', 12);
plot(0.4104, 0.3105, 'ks', 'MarkerSize', 12);
plot(0.2654, 0.3763, 'ks', 'MarkerSize', 12);
plot(0.3654, 0.3735, 'ks', 'MarkerSize', 12);

blueextremidx = find(xyY(1,:) > 0.2171 & xyY(1,:) < 0.2173 & xyY(2,:) > 0.2120 & xyY(2,:) < 0.2122)
blue_ext = [wlns(:) munsell(:,blueextremidx)./max(munsell(:,blueextremidx))];
dlmwrite("spectra/munsell_blue_extrem.spd", blue_ext, 'delimiter', ' ')

greenextremidx = find(xyY(1,:) > 0.2793 & xyY(1,:) < 0.2795 & xyY(2,:) > 0.4588 & xyY(2,:) < 0.4590)
green_ext = [wlns(:) munsell(:,greenextremidx)./max(munsell(:,greenextremidx))];
dlmwrite("spectra/munsell_green_extrem.spd", green_ext, 'delimiter', ' ')

yellowextremidx = find(xyY(1,:) > 0.4361 & xyY(1,:) < 0.4363 & xyY(2,:) > 0.4255 & xyY(2,:) < 0.4257)
yellow_ext = [wlns(:) munsell(:,yellowextremidx)./max(munsell(:,yellowextremidx))];
dlmwrite("spectra/munsell_yellow_extrem.spd", yellow_ext, 'delimiter', ' ')

redextremidx = find(xyY(1,:) > 0.4416 & xyY(1,:) < 0.4418 & xyY(2,:) > 0.2656 & xyY(2,:) < 0.2658)
red_ext = [wlns(:) munsell(:,redextremidx)./max(munsell(:,redextremidx))];
dlmwrite("spectra/munsell_red_extrem.spd", red_ext, 'delimiter', ' ')
