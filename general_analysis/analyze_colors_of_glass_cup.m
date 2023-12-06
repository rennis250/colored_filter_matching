clear all; close all; clc;

% Change default axes fonts.
set(0,'DefaultAxesFontSize', 15);

% Change default text fonts.
set(0,'DefaultTextFontSize', 15);

monxyY = [.6804 .3073 30.94
    .2029 .6968 74.22
    .1527 .0508 6.74]; % assume that image is shown on same monitor from experiment
monxyz = xyY2XYZ(monxyY);

red_book = im2double(imread('glass_red_book.jpg'));
green_book = im2double(imread('glass_green_book.jpg'));
blue_book = im2double(imread('glass_blue_book.jpg'));

red_mask = im2double(imread('glass_red_bookMask.png'));
green_mask = im2double(imread('glass_green_bookMask.png'));
blue_mask = im2double(imread('glass_blue_bookMask.png'));

imgs = {};

imgs{1} = red_book;
imgs{2} = green_book;
imgs{3} = blue_book;

rgbs = {};
labs = {};

[imgp, ~, gis, ~] = betterMaskingFunction(red_mask);
[rgb, lab, mean_lab] = turnRegionsToLABandRGB(imgp, monxyz, gis);

% only use a random sub-set of pixels, since images are huge
rand_idxs = randperm(length(gis));
rand_idxs = rand_idxs(1:10000);

rgbs{1} = rgb(rand_idxs, :);
labs{1} = lab(rand_idxs, :);
mean_rgbs(1, :) = mean(rgb, 1);
mean_labs(1, :) = mean_lab;

[imgp, ~, gis, ~] = betterMaskingFunction(green_mask);
[rgb, lab, mean_lab] = turnRegionsToLABandRGB(imgp, monxyz, gis);

rgbs{2} = rgb(rand_idxs, :);
labs{2} = lab(rand_idxs, :);
mean_rgbs(2, :) = mean(rgb, 1);
mean_labs(2, :) = mean_lab;

[imgp, ~, gis, ~] = betterMaskingFunction(blue_mask);
[rgb, lab, mean_lab] = turnRegionsToLABandRGB(imgp, monxyz, gis);

rgbs{3} = rgb(rand_idxs, :);
labs{3} = lab(rand_idxs, :);
mean_rgbs(3, :) = mean(rgb, 1);
mean_labs(3, :) = mean_lab;

% now make big plot showing color distributions from regions of interest and mean colors.
figure(1);
clf;

for x = 1:3
    subplot(2,3,x);
    hold on;
    xlabel('a*');
    ylabel('b*');
    axis([-150 150 -150 150]);
    axis square;
    grid;
    lstp = labs{x};
    rstp = rgbs{x};
    p_obj_pts = scatter(lstp(:, 2), lstp(:, 3), [], rstp.^(1/2.2), 'filled');
    p_obj = plot(mean_labs(x, 2), mean_labs(x, 3), 's', 'MarkerEdgeColor', [1 1 1], 'MarkerFaceColor', mean_rgbs(x, :).^(1/2.2), 'MarkerSize', 15);

    h_leg = legend([p_obj_pts, p_obj], {'Object', 'Mean Object'});

    subplot(2,3,x+3);
    imshow(imgs{x});
end

export_fig('glassCups_Book_LAB_colors.pdf', '-pdf');
