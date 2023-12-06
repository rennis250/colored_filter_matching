clear all; close all; clc;

%% basic stuff for plotting

% Change default axes fonts.
set(0,'DefaultAxesFontSize', 15);

% Change default text fonts.
set(0,'DefaultTextFontSize', 15);

load('object_illum_colors_for_plot.mat');
load('illums_lab.mat');

trans_color_for_plot(1, :) = [0 0 1];
trans_color_for_plot(2, :) = [1.0 0.5 0.0];
trans_color_for_plot(3, :) = [1 0 0];
trans_color_for_plot(4, :) = [0 1 0];
trans_color_for_plot(5, :) = [0 0 0.5];
trans_color_for_plot(6, :) = [0.5, 0.25, 0.0];
trans_color_for_plot(7, :) = [0.5 0 0];
trans_color_for_plot(8, :) = [0 0.5 0];

illum_color_for_plot(1, :) = object_colors(1, :);
illum_color_for_plot(2, :) = object_colors(2, :);

%% plot exps 2 and 3 first

load('../data/exp2.mat', 'mean_obj_lab', 'popavgdata_lab', 'colors_for_scatter_plot', 'which_illum');
mean_obj_lab_exp2 = mean_obj_lab;
popavgdata_lab_exp2 = popavgdata_lab;

load('../data/exp3.mat', 'mean_obj_lab', 'popavgdata_lab');
mean_obj_lab_exp3 = mean_obj_lab;
popavgdata_lab_exp3 = popavgdata_lab;

%% exps 2 and 3

figure(1);
clf;
set(gcf,'color','w');
subplot(1,3,1);
hold on;
axis square;
xlabel("L* mean object");
ylabel("L* mean match");
l = line([-100 100], [-100 100]);
set(l, 'LineWidth', 2);
set(l, 'LineStyle', '--');
axis([0 100 0 100]);
axis square;

coeffs = polyfit(mean_obj_lab_exp2(:, 1), popavgdata_lab_exp2(:, 1), 1);
fittedX = linspace(0, 100, 200);
fittedY = polyval(coeffs, fittedX);
l1 = plot(fittedX, fittedY, 'k-', 'LineWidth', 3);

p1 = scatter(mean_obj_lab_exp2(which_illum == 1, 1), popavgdata_lab_exp2(which_illum == 1, 1), 100, colors_for_scatter_plot(which_illum == 1, :), 'o', 'filled');
p2 = scatter(mean_obj_lab_exp2(which_illum == 2, 1), popavgdata_lab_exp2(which_illum == 2, 1), 100, colors_for_scatter_plot(which_illum == 2, :), 'v');

coeffs = polyfit(mean_obj_lab_exp3(:, 1), popavgdata_lab_exp3(:, 1), 1);
fittedX = linspace(0, 100, 200);
fittedY = polyval(coeffs, fittedX);
l2 = plot(fittedX, fittedY, 'r-', 'LineWidth', 3, 'Color', [0.5 0.5 0.5]);

scatter(mean_obj_lab_exp3(which_illum == 1, 1), popavgdata_lab_exp3(which_illum == 1, 1), 100, colors_for_scatter_plot(which_illum == 1, :), 'o', 'filled');
scatter(mean_obj_lab_exp3(which_illum == 2, 1), popavgdata_lab_exp3(which_illum == 2, 1), 100, colors_for_scatter_plot(which_illum == 2, :), 'v');

l = line([-100 100], [-100 100]);
set(l, 'LineWidth', 2);
set(l, 'LineStyle', '--');

legend([p1, p2, l1, l2], 'Blue Illuminant', 'Yellow Illuminant', 'Experiment 1', 'Experiment 2', 'Location','SouthEast');

subplot(1,3,2);
hold on;
axis square;
xlabel("a* mean object");
ylabel("a* mean match");

coeffs = polyfit(mean_obj_lab_exp2(:, 2), popavgdata_lab_exp2(:, 2), 1);
fittedX = linspace(-100, 100, 200);
fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'k-', 'LineWidth', 3);

p1 = scatter(mean_obj_lab_exp2(which_illum == 1, 2), popavgdata_lab_exp2(which_illum == 1, 2), 100, colors_for_scatter_plot(which_illum == 1, :), 'o', 'filled');
p2 = scatter(mean_obj_lab_exp2(which_illum == 2, 2), popavgdata_lab_exp2(which_illum == 2, 2), 100, colors_for_scatter_plot(which_illum == 2, :), 'v');

coeffs = polyfit(mean_obj_lab_exp3(:, 2), popavgdata_lab_exp3(:, 2), 1);
fittedX = linspace(-100, 100, 200);
fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'r-', 'LineWidth', 3, 'Color', [0.5 0.5 0.5]);

scatter(mean_obj_lab_exp3(which_illum == 1, 2), popavgdata_lab_exp3(which_illum == 1, 2), 100, colors_for_scatter_plot(which_illum == 1, :), 'o', 'filled');
scatter(mean_obj_lab_exp3(which_illum == 2, 2), popavgdata_lab_exp3(which_illum == 2, 2), 100, colors_for_scatter_plot(which_illum == 2, :), 'v');

l = line([-100 100], [-100 100]);
set(l, 'LineWidth', 2);
set(l, 'LineStyle', '--');
axis square;

subplot(1,3,3);
hold on;
axis square;
xlabel("b* mean object");
ylabel("b* mean match");

coeffs = polyfit(mean_obj_lab_exp2(:, 3), popavgdata_lab_exp2(:, 3), 1);
fittedX = linspace(-100, 100, 200);
fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'k-', 'LineWidth', 3);

p1 = scatter(mean_obj_lab_exp2(which_illum == 1, 3), popavgdata_lab_exp2(which_illum == 1, 3), 100, colors_for_scatter_plot(which_illum == 1, :), 'o', 'filled');
p2 = scatter(mean_obj_lab_exp2(which_illum == 2, 3), popavgdata_lab_exp2(which_illum == 2, 3), 100, colors_for_scatter_plot(which_illum == 2, :), 'v');

coeffs = polyfit(mean_obj_lab_exp3(:, 3), popavgdata_lab_exp3(:, 3), 1);
fittedX = linspace(-100, 100, 200);
fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'r-', 'LineWidth', 3, 'Color', [0.5 0.5 0.5]);

scatter(mean_obj_lab_exp3(which_illum == 1, 3), popavgdata_lab_exp3(which_illum == 1, 3), 100, colors_for_scatter_plot(which_illum == 1, :), 'o', 'filled');
scatter(mean_obj_lab_exp3(which_illum == 2, 3), popavgdata_lab_exp3(which_illum == 2, 3), 100, colors_for_scatter_plot(which_illum == 2, :), 'v');

l = line([-100 100], [-100 100]);
set(l, 'LineWidth', 2);
set(l, 'LineStyle', '--');
axis square;
export_fig('../figures/exps_2_and_3_together.pdf', '-pdf');

%% exps 4 and 5

load('../data/achromatic_exp.mat', 'mean_obj_lab', 'popavgdata_lab', 'colors_for_scatter_plot');
mean_obj_lab_exp4 = mean_obj_lab;
popavgdata_lab_exp4 = popavgdata_lab;

load('../data/achromatic_w_patch.mat', 'mean_obj_lab', 'popavgdata_lab');
mean_obj_lab_exp5 = mean_obj_lab;
popavgdata_lab_exp5 = popavgdata_lab;

figure(1);
clf;
set(gcf,'color','w');
subplot(1,3,1);
hold on;
axis square;
xlabel("L* mean object");
ylabel("L* mean match");
l = line([-100 100], [-100 100]);
set(l, 'LineWidth', 2);
set(l, 'LineStyle', '--');
axis([0 100 0 100]);
axis square;

coeffs = polyfit(mean_obj_lab_exp4(:, 1), popavgdata_lab_exp4(:, 1), 1);
fittedX = linspace(0, 100, 200);
fittedY = polyval(coeffs, fittedX);
l1 = plot(fittedX, fittedY, 'r-', 'LineWidth', 3, 'Color', [0.5 0.5 0.5]);

p1 = scatter(mean_obj_lab_exp4(:, 1), popavgdata_lab_exp4(:, 1), 100, colors_for_scatter_plot(:, :), 's', 'filled');

coeffs = polyfit(mean_obj_lab_exp5(:, 1), popavgdata_lab_exp5(:, 1), 1);
fittedX = linspace(0, 100, 200);
fittedY = polyval(coeffs, fittedX);
l2 = plot(fittedX, fittedY, 'k-', 'LineWidth', 3);

scatter(mean_obj_lab_exp5(:, 1), popavgdata_lab_exp5(:, 1), 100, colors_for_scatter_plot(:, :), 's', 'filled');

l = line([-100 100], [-100 100]);
set(l, 'LineWidth', 2);
set(l, 'LineStyle', '--');

legend([p1, l2, l1], 'Achromatic Illuminant', 'Experiment 3', 'Experiment 4', 'Location','SouthEast');

subplot(1,3,2);
hold on;
axis square;
xlabel("a* mean object");
ylabel("a* mean match");

coeffs = polyfit(mean_obj_lab_exp4(:, 2), popavgdata_lab_exp4(:, 2), 1);
fittedX = linspace(-100, 100, 200);
fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'r-', 'LineWidth', 3, 'Color', [0.5 0.5 0.5]);

p1 = scatter(mean_obj_lab_exp4(:, 2), popavgdata_lab_exp4(:, 2), 100, colors_for_scatter_plot(:, :), 's', 'filled');

coeffs = polyfit(mean_obj_lab_exp5(:, 2), popavgdata_lab_exp5(:, 2), 1);
fittedX = linspace(-100, 100, 200);
fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'k-', 'LineWidth', 3);

scatter(mean_obj_lab_exp5(:, 2), popavgdata_lab_exp5(:, 2), 100, colors_for_scatter_plot(:, :), 's', 'filled');

l = line([-100 100], [-100 100]);
set(l, 'LineWidth', 2);
set(l, 'LineStyle', '--');
axis square;

subplot(1,3,3);
hold on;
axis square;
xlabel("b* mean object");
ylabel("b* mean match");

coeffs = polyfit(mean_obj_lab_exp4(:, 3), popavgdata_lab_exp4(:, 3), 1);
fittedX = linspace(-100, 100, 200);
fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'r-', 'LineWidth', 3, 'Color', [0.5 0.5 0.5]);

p1 = scatter(mean_obj_lab_exp4(:, 3), popavgdata_lab_exp4(:, 3), 100, colors_for_scatter_plot(:, :), 's', 'filled');

coeffs = polyfit(mean_obj_lab_exp5(:, 3), popavgdata_lab_exp5(:, 3), 1);
fittedX = linspace(-100, 100, 200);
fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'k-', 'LineWidth', 3);

scatter(mean_obj_lab_exp5(:, 3), popavgdata_lab_exp5(:, 3), 100, colors_for_scatter_plot(:, :), 's', 'filled');

l = line([-100 100], [-100 100]);
set(l, 'LineWidth', 2);
set(l, 'LineStyle', '--');
axis square;
export_fig('../figures/exps_4_and_5_together.pdf', '-pdf');
