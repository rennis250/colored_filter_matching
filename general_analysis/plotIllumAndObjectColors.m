clear all;
close all;
clc;

% Change default axes fonts.
set(0,'DefaultAxesFontSize', 15);

% Change default text fonts.
set(0,'DefaultTextFontSize', 15);

object_colors = zeros(8, 3);

% load distributions for illum and object colors
lightredrefl = dlmread('/Users/rje/my_docs/simpler_pathmarcher/spectra/munsell_red_lightest.spd');
lightbluerefl = dlmread('/Users/rje/my_docs/simpler_pathmarcher/spectra/munsell_blue_lightest.spd');
lightgreenrefl = dlmread('/Users/rje/my_docs/simpler_pathmarcher/spectra/munsell_green_lightest.spd');
lightyellowrefl = dlmread('/Users/rje/my_docs/simpler_pathmarcher/spectra/munsell_yellow_lightest.spd');

[lightredrgb,lightredx] = colormatchPlanck(lightredrefl(:,2), lightredrefl(:,1));
lightredrgb = lightredrgb./650;
object_colors(1, :) = lightredrgb
lightredxyy = [lightredx(1)/sum(lightredx) lightredx(2)/sum(lightredx) lightredx(2)];

[lightgreenrgb,lightgreenx] = colormatchPlanck(lightgreenrefl(:,2), lightgreenrefl(:,1));
lightgreenrgb = lightgreenrgb./650;
object_colors(2, :) = lightgreenrgb;
lightgreenxyy = [lightgreenx(1)/sum(lightgreenx) lightgreenx(2)/sum(lightgreenx) lightgreenx(2)];

[lightbluergb,lightbluex] = colormatchPlanck(lightbluerefl(:,2), lightbluerefl(:,1));
lightbluergb = lightbluergb./650;
object_colors(3, :) = lightbluergb;
lightbluexyy = [lightbluex(1)/sum(lightbluex) lightbluex(2)/sum(lightbluex) lightbluex(2)];

[lightyellowrgb,lightyellowx] = colormatchPlanck(lightyellowrefl(:,2), lightyellowrefl(:,1));
lightyellowrgb = lightyellowrgb./650;
object_colors(4, :) = lightyellowrgb;
lightyellowxyy = [lightyellowx(1)/sum(lightyellowx) lightyellowx(2)/sum(lightyellowx) lightyellowx(2)];

% plot light transmissions - one at a time

figure(3);
clf;
hold on;
xlim([lightredrefl(1,1), lightredrefl(end,1)]);
ylim([0,100]);
ylabel('Transmission (%)');
xlabel('Wavelength (nm)');
plot(lightredrefl(:,1), 100.0*lightredrefl(:,2), 'r-', 'LineWidth', 2, 'Color', lightredrgb);
axis square;
export_fig('light_red_tranmission.pdf', '-dpdf');

figure(4);
clf;
hold on;
xlim([lightredrefl(1,1), lightredrefl(end,1)]);
ylim([0,100]);
ylabel('Transmission (%)');
xlabel('Wavelength (nm)');
plot(lightgreenrefl(:,1), 100.0*lightgreenrefl(:,2), 'g-', 'LineWidth', 2, 'Color', lightgreenrgb);
axis square;
export_fig('light_green_tranmission.pdf', '-dpdf');

figure(5);
clf;
hold on;
xlim([lightredrefl(1,1), lightredrefl(end,1)]);
ylim([0,100]);
ylabel('Transmission (%)');
xlabel('Wavelength (nm)');
plot(lightbluerefl(:,1), 100.0*lightbluerefl(:,2), 'b-', 'LineWidth', 2, 'Color', lightbluergb);
axis square;
export_fig('light_blue_tranmission.pdf', '-dpdf');

figure(6);
clf;
hold on;
xlim([lightredrefl(1,1), lightredrefl(end,1)]);
ylim([0,100]);
ylabel('Transmission (%)');
xlabel('Wavelength (nm)');
plot(lightyellowrefl(:,1), 100.0*lightyellowrefl(:,2), 'y-', 'LineWidth', 2, 'Color', lightyellowrgb);
axis square;
export_fig('light_yellow_tranmission.pdf', '-dpdf');

darkredrefl = dlmread('/Users/rje/my_docs/simpler_pathmarcher/spectra/munsell_red_lighter.spd');
darkbluerefl = dlmread('/Users/rje/my_docs/simpler_pathmarcher/spectra/munsell_blue_lighter.spd');
darkgreenrefl = dlmread('/Users/rje/my_docs/simpler_pathmarcher/spectra/munsell_green_lighter.spd');
darkyellowrefl = dlmread('/Users/rje/my_docs/simpler_pathmarcher/spectra/munsell_yellow_lighter.spd');

[darkredrgb,darkredx] = colormatchPlanck(darkredrefl(:,2), darkredrefl(:,1));
darkredrgb = darkredrgb./650;
object_colors(5, :) = darkredrgb;
darkredxyy = [darkredx(1)/sum(darkredx) darkredx(2)/sum(darkredx) darkredx(2)];

[darkgreenrgb,darkgreenx] = colormatchPlanck(darkgreenrefl(:,2), darkgreenrefl(:,1));
darkgreenrgb = darkgreenrgb./650;
object_colors(6, :) = darkgreenrgb;
darkgreenxyy = [darkgreenx(1)/sum(darkgreenx) darkgreenx(2)/sum(darkgreenx) darkgreenx(2)];

[darkbluergb,darkbluex] = colormatchPlanck(darkbluerefl(:,2), darkbluerefl(:,1));
darkbluergb = darkbluergb./650;
object_colors(7, :) = darkbluergb;
darkbluexyy = [darkbluex(1)/sum(darkbluex) darkbluex(2)/sum(darkbluex) darkbluex(2)];

[darkyellowrgb,darkyellowx] = colormatchPlanck(darkyellowrefl(:,2), darkyellowrefl(:,1));
darkyellowrgb = darkyellowrgb./650;
object_colors(8, :) = darkyellowrgb;
darkyellowxyy = [darkyellowx(1)/sum(darkyellowx) darkyellowx(2)/sum(darkyellowx) darkyellowx(2)];

% plot dark transmissions - one at a time

figure(3);
clf;
hold on;
axis square;
xlim([lightredrefl(1,1), lightredrefl(end,1)]);
ylim([0,100]);
ylabel('Transmission (%)');
xlabel('Wavelength (nm)');
plot(darkredrefl(:,1), 100.0*darkredrefl(:,2), 'r-', 'LineWidth', 2, 'Color', darkredrgb);
axis square;
export_fig('dark_red_tranmission.pdf', '-dpdf');

figure(4);
clf;
hold on;
axis square;
xlim([lightredrefl(1,1), lightredrefl(end,1)]);
ylim([0,100]);
ylabel('Transmission (%)');
xlabel('Wavelength (nm)');
plot(darkgreenrefl(:,1), 100.0*darkgreenrefl(:,2), 'g-', 'LineWidth', 2, 'Color', darkgreenrgb);
axis square;
export_fig('dark_green_tranmission.pdf', '-dpdf');

figure(5);
clf;
hold on;
axis square;
xlim([lightredrefl(1,1), lightredrefl(end,1)]);
ylim([0,100]);
ylabel('Transmission (%)');
xlabel('Wavelength (nm)');
plot(darkbluerefl(:,1), 100.0*darkbluerefl(:,2), 'b-', 'LineWidth', 2, 'Color', darkbluergb);
axis square;
export_fig('dark_blue_tranmission.pdf', '-dpdf');

figure(6);
clf;
hold on;
axis square;
xlim([lightredrefl(1,1), lightredrefl(end,1)]);
ylim([0,100]);
ylabel('Transmission (%)');
xlabel('Wavelength (nm)');
plot(darkyellowrefl(:,1), 100.0*darkyellowrefl(:,2), 'y-', 'LineWidth', 2, 'Color', darkyellowrgb);
axis square;
export_fig('dark_yellow_tranmission.pdf', '-dpdf');

redspectra = dlmread('/Users/rje/my_docs/simpler_pathmarcher/spectra/de_brighter_illum_red.spd');
bluespectra = dlmread('/Users/rje/my_docs/simpler_pathmarcher/spectra/de_brighter_illum_blue.spd');
greenspectra = dlmread('/Users/rje/my_docs/simpler_pathmarcher/spectra/de_brighter_illum_green.spd');
yellowspectra = dlmread('/Users/rje/my_docs/simpler_pathmarcher/spectra/de_brighter_illum_yellow.spd');

illum_colors = zeros(2, 3);

[redsrgb,redsx] = colormatchPlanck(redspectra(:,2), redspectra(:,1));
redsrgb = redsrgb./2500;
redsxyy = [redsx(1)/sum(redsx) redsx(2)/sum(redsx) redsx(2)];
[greensrgb,greensx] = colormatchPlanck(greenspectra(:,2), greenspectra(:,1));
greensrgb = greensrgb./2500;
greensxyy = [greensx(1)/sum(greensx) greensx(2)/sum(greensx) greensx(2)];

[bluesrgb,bluesx] = colormatchPlanck(bluespectra(:,2), bluespectra(:,1));
bluesrgb = bluesrgb./2500;
illum_colors(1, :) = bluesrgb;
bluesxyy = [bluesx(1)/sum(bluesx) bluesx(2)/sum(bluesx) bluesx(2)];

[yellowsrgb,yellowsx] = colormatchPlanck(yellowspectra(:,2), yellowspectra(:,1));
yellowsrgb = yellowsrgb./2500;
illum_colors(2, :) = yellowsrgb;
yellowsxyy = [yellowsx(1)/sum(yellowsx) yellowsx(2)/sum(yellowsx) yellowsx(2)];

save('object_illum_colors_for_plot.mat', 'object_colors', 'illum_colors');

pause;

% plot illum spectra - one at a time

% note that somehow, someway, the green and red spectra were accidentally switched into
% different files. as usual, i have no idea how...

% green spectrum
figure(3);
clf;
hold on;
axis square;
ylim([0,14]);
xlabel('Wavelength (nm)');
ylabel('Radiance (W/(sr * m^{2}))');
plot(redspectra(:,1), redspectra(:,2), 'r-', 'LineWidth', 2, 'Color', redsrgb);
export_fig('illum_green_spectra.pdf', '-dpdf');

% red spectrum
figure(4);
clf;
hold on;
axis square;
ylim([0,14]);
xlabel('Wavelength (nm)');
ylabel('Radiance (W/(sr * m^{2}))');
plot(greenspectra(:,1), greenspectra(:,2), 'g-', 'LineWidth', 2, 'Color', greensrgb);
export_fig('illum_red_spectra.pdf', '-dpdf');

figure(5);
clf;
hold on;
axis square;
ylim([0,14]);
xlabel('Wavelength (nm)');
ylabel('Radiance (W/(sr * m^{2}))');
plot(bluespectra(:,1), bluespectra(:,2), 'b-', 'LineWidth', 2, 'Color', bluesrgb);
export_fig('illum_blue_spectra.pdf', '-dpdf');

figure(6);
clf;
hold on;
axis square;
ylim([0,14]);
xlabel('Wavelength (nm)');
ylabel('Radiance (W/(sr * m^{2}))');
plot(yellowspectra(:,1), yellowspectra(:,2), 'y-', 'LineWidth', 2, 'Color', yellowsrgb);
export_fig('illum_yellow_spectra.pdf', '-dpdf');

% grab pre-computed xyY coords of CIE lasso and planckian locus for plotting
lassocs = csvread('xyYcoords.csv');
planck = csvread('planck.csv');

% first plot showing transmission colors (only need to plot lighter or darker, since they have
% same chroma, but different luminance)
figure(1);
clf;
hold on;
plot(lassocs(:,1), lassocs(:,2),'k-');
plot(planck(:,1), planck(:,2),'k-');
plot(lightredxyy(1), lightredxyy(2), 'ro', 'MarkerFaceColor', lightredrgb, 'MarkerEdgeColor', lightredrgb);
plot(lightgreenxyy(1), lightgreenxyy(2), 'go', 'MarkerFaceColor', lightgreenrgb, 'MarkerEdgeColor', lightgreenrgb);
plot(lightbluexyy(1), lightbluexyy(2), 'bo', 'MarkerFaceColor', lightbluergb, 'MarkerEdgeColor', lightbluergb);
plot(lightyellowxyy(1), lightyellowxyy(2), 'yo', 'MarkerFaceColor', lightyellowrgb, 'MarkerEdgeColor', lightyellowrgb);
axis([0 1 0 1]);
axis square;
xlabel('x');
ylabel('y');
title('"Chromaticity coordinates" of transmission distributions');
axis square

export_fig('transmission_colors_xyy.pdf', '-dpdf');

% second plot showing blue-yellow illuminant colors
figure(2);
clf;
hold on;
plot(lassocs(:,1), lassocs(:,2),'k-');
plot(planck(:,1), planck(:,2),'k-');
plot(bluesxyy(1), bluesxyy(2), 'bo', 'MarkerFaceColor', bluesrgb, 'MarkerEdgeColor', bluesrgb);
plot(yellowsxyy(1), yellowsxyy(2), 'yo', 'MarkerFaceColor', yellowsrgb, 'MarkerEdgeColor', yellowsrgb);
axis([0 1 0 1]);
axis square;
xlabel('x');
ylabel('y');
title('Chromaticity coordinates of "blue-yellow" illuminants');
axis square

export_fig('blue_yellow_illum_xyy.pdf', '-dpdf');

% third plot showing red-green illuminant colors
figure(2);
clf;
hold on;
plot(lassocs(:,1), lassocs(:,2),'k-');
plot(planck(:,1), planck(:,2),'k-');
% red and green were accidentally swapped
plot(redsxyy(1), redsxyy(2), 'ro', 'MarkerFaceColor', redsrgb, 'MarkerEdgeColor', redsrgb);
plot(greensxyy(1), greensxyy(2), 'go', 'MarkerFaceColor', greensrgb, 'MarkerEdgeColor', greensrgb);
axis([0 1 0 1]);
axis square;
xlabel('x');
ylabel('y');
title('Chromaticity coordinates of "red-green" illuminants');
axis square

export_fig('red_green_illum_xyy.pdf', '-dpdf');
