clear all; close all; clc;

% Change default axes fonts.
set(0,'DefaultAxesFontSize', 15);

% Change default text fonts.
set(0,'DefaultTextFontSize', 15);

monxyY = [.6804 .3073 30.94
    .2029 .6968 74.22
    .1527 .0508 6.74]; % assume that image is shown on same monitor from experiment
monxyz = xyY2XYZ(monxyY);

de_dist = 15;

% prepare masks to extract object and highlight colors
pm = imread('glassMask.png');
pm = pm(:,:,1); pm = logical(pm);
gisp = find(pm);
pm_wout_high = pm;

bm = imread('glassHighlight.png');
br = bm(:,:,1); bg = bm(:,:,2); bb = bm(:,:,3);
idxs = br ~= 0 | bg ~= 0 | bb ~= 0;
hr2 = br; hg2 = bg; hb2 = bb;
hr2(idxs) = 255; hg2(idxs) = 255; hb2(idxs) = 255;
highlight2(:,:,1) = hr2;
highlight2(:,:,2) = hg2;
highlight2(:,:,3) = hb2;
bm = highlight2(:,:,1); bm = logical(bm);
gisb = find(bm);

pm_wout_high(gisb) = logical(0.0); % don't include highlight in estimate of color
gisp_wout_high = find(pm_wout_high);

sm = logical(ones(size(pm,1),size(pm,2)));

xc = 121; yc = 129;
backx = 212; backy = 61;
shadx = 600; shady = 281;

gR = 2.2; gG = 2.2; gB = 2.2;

curr_img = im2double(imread('glassbottle.png'));

obj_lab = [];
high_lab = [];
mean_obj_lab = [];
mean_high_lab = [];
lab_lum_splits = {};
lum_split_diffs = zeros(1,3);

% first let's compute mean colors of regions of interest
[imp, imh] = maskAndLinearTransparent(curr_img, pm, bm, sm, xc, yc, backx, backy, shadx, shady, gR, gG, gB);
[imp_wout_high, ~] = maskAndLinearTransparent(curr_img, pm_wout_high, bm, sm, xc, yc, backx, backy, shadx, shady, gR, gG, gB);

rp = imp(:,:,1); rp = rp(:);
gp = imp(:,:,2); gp = gp(:);
bp = imp(:,:,3); bp = bp(:);
obj_rgb = [rp(gisp) gp(gisp) bp(gisp)];
obj_lab = rgb2labRob(obj_rgb, monxyz);
mean_obj_lab = mean(squeeze(obj_lab(:, :)), 1);

rp = imp_wout_high(:,:,1); rp = rp(:);
gp = imp_wout_high(:,:,2); gp = gp(:);
bp = imp_wout_high(:,:,3); bp = bp(:);
obj_rgb_wout_high = [rp(gisp_wout_high) gp(gisp_wout_high) bp(gisp_wout_high)];
obj_lab_wout_high = rgb2labRob(obj_rgb_wout_high, monxyz);
mean_obj_rgb_wout_high = mean(obj_rgb_wout_high, 1);
mean_obj_lab_wout_high = mean(squeeze(obj_lab_wout_high(:, :)), 1);

rh = imh(:,:,1); rh = rh(:);
gh = imh(:,:,2); gh = gh(:);
bh = imh(:,:,3); bh = bh(:);
high_rgb = [rh(gisb) gh(gisb) bh(gisb)];
high_lab = rgb2labRob(high_rgb, monxyz);
mean_high_rgb = mean(high_rgb, 1);
mean_high_lab = mean(squeeze(high_lab), 1);

lum_regions = linspace(0, squeeze(max(max(obj_lab(:, 1)))), 4);

% next, let's try to get regions of different luminance and see what the mean colors are
c = 1;
for lc = 1:length(lum_regions)-1
    ls = squeeze(obj_lab_wout_high(:, 1));
    lab_lum_splits{c} = squeeze(obj_lab_wout_high(ls >= lum_regions(lc) & ls < lum_regions(lc+1), :));
    mean_lab_lum_splits(c, :) = mean(lab_lum_splits{c}, 1);
    c = c + 1;
end

% and we can find the closest point in the luminance bands to the matched color
for lc = 1:c-1
    lum_split_diffs(lc) = norm(squeeze(mean_lab_lum_splits(lc, :)) - mean_obj_lab_wout_high');
end
best_lum_split = find(lum_split_diffs(:) == min(lum_split_diffs(:)));

% let's compute DE units for difference between mean of obj color and mean of highlight color
de_units = deltaE2000(repmat(mean_obj_lab_wout_high, size(obj_lab, 1), 1), squeeze(obj_lab));

% output a grayscale version with pixels that are 'de_dist' DE units colored in,
% as well as a version where the different luminance bands are highlighted
curr_img_r = curr_img(:,:,1); curr_img_r(~pm) = 0.0;
curr_img_g = curr_img(:,:,2); curr_img_g(~pm) = 0.0;
curr_img_b = curr_img(:,:,3); curr_img_b(~pm) = 0.0;

curr_img2(:,:,1) = curr_img_r;
curr_img2(:,:,2) = curr_img_g;
curr_img2(:,:,3) = curr_img_b;

curr_img_gray_tmp = rgb2gray(curr_img2);
curr_img_gray = zeros(size(curr_img2));
curr_img_gray(:,:,1) = curr_img_gray_tmp;
curr_img_gray(:,:,2) = curr_img_gray_tmp;
curr_img_gray(:,:,3) = curr_img_gray_tmp;

idxs = gisp(find(de_units <= de_dist));
curr_img_col = reshape(curr_img2, size(curr_img2,1)*size(curr_img2,2), 3);
curr_img_gray_col = reshape(curr_img_gray, size(curr_img2,1)*size(curr_img2,2), 3);

curr_img_gray_col(idxs, :) = curr_img_col(idxs, :);

imwrite(reshape(curr_img_gray_col, size(curr_img2,1), size(curr_img2,2), size(curr_img2,3)), 'glassbottle_de_units.png');

lum_regions = linspace(0, squeeze(max(max(obj_lab(:, 1)))), 4);

for lc = 1:length(lum_regions)-1
    ls = squeeze(obj_lab(:, 1));

    idxs = gisp(find(ls >= lum_regions(lc) & ls < lum_regions(lc+1)));
    curr_img_col = reshape(curr_img2, size(curr_img2,1)*size(curr_img2,2), 3);
    curr_img_gray_col = reshape(curr_img_gray, size(curr_img2,1)*size(curr_img2,2), 3);

    curr_img_gray_col(idxs, :) = curr_img_col(idxs, :);

    imwrite(reshape(curr_img_gray_col, size(curr_img2,1), size(curr_img2,2), size(curr_img2,3)), ['glassbottle_lum_region_' num2str(lc) '.png']);
end

% now make big plot showing color distributions from regions of interest and mean colors.
% make multiple so that they can build up over slides

bnd_idxs = boundary(squeeze(obj_lab(:, 2:3)));
bnd_idxs_wout_high = boundary(squeeze(obj_lab_wout_high(:, 2:3)));
bnd_idxs_high = boundary(squeeze(high_lab(:, 2:3)));

obj_rgb_mean = mean_obj_rgb_wout_high(:).^(1/2.2);
high_rgb_mean = mean_high_rgb(:).^(1/2.2);

% first plot - show LAB distribution of object only
figure(2);
clf;
subplot(1,2,1);
hold on;
xlabel('a*');
ylabel('b*');
axis([-150 150 -150 150]);
axis square;
grid;
title(['Mean Luminance - ' num2str(mean_obj_lab_wout_high(1))]);
p_obj_pts = scatter(obj_lab(:, 2), obj_lab(:, 3), [], obj_rgb.^(1/2.2), 'filled');
p_obj = plot(mean_obj_lab(2), mean_obj_lab(3), 's', 'MarkerEdgeColor', [1 1 1], 'MarkerFaceColor', obj_rgb_mean, 'MarkerSize', 15);

h_leg = legend([p_obj_pts, p_obj], {'Object', 'Mean Object'});

subplot(1,2,2);
imshow(curr_img);
export_fig('glassbottle_lab_obj_only.pdf', '-pdf');

% second plot - show highlight also
figure(2);
clf;
subplot(1,2,1);
hold on;
xlabel('a*');
ylabel('b*');
axis([-150 150 -150 150]);
axis square;
grid;
title(['Mean Luminance - ' num2str(mean_obj_lab_wout_high(1))]);
p_obj_pts = plot(obj_lab_wout_high(bnd_idxs_wout_high, 2), obj_lab_wout_high(bnd_idxs_wout_high, 3), 'k-');
p_high_pts = scatter(high_lab(:, 2), high_lab(:, 3), [], high_rgb.^(1/2.2), 'filled');
p_obj = plot(mean_obj_lab_wout_high(2), mean_obj_lab_wout_high(3), 's', 'MarkerFaceColor', obj_rgb_mean, 'MarkerSize', 15);
p_high = plot(mean_high_lab(2), mean_high_lab(3), 'o', 'MarkerFaceColor', high_rgb_mean, 'MarkerSize', 15);

h_leg = legend([p_obj_pts, p_obj, p_high_pts, p_high], {'Object', 'Mean Object', 'Highlight', 'Mean Highlight'});

subplot(1,2,2);
imshow(curr_img)
export_fig('glassbottle_lab_obj_high.pdf', '-pdf');

%%%% save examples of the masks

imtom = curr_img;

% show highlight
immr = imtom(:,:,1); immg = imtom(:,:,2); immb = imtom(:,:,3);
immr(~bm) = 0.5; immg(~bm) = 0.5; immb(~bm) = 0.5;
imtos = zeros(size(imtom));
imtos(:,:,1) = immr; imtos(:,:,2) = immg; imtos(:,:,3) = immb;
figure(10);
clf;
hold on;
imshow(imtos);
imwrite(imtos, 'glassbottle_masked_high.png');

% show object
immr = imtom(:,:,1); immg = imtom(:,:,2); immb = imtom(:,:,3);
immr(~pm) = 0.5; immg(~pm) = 0.5; immb(~pm) = 0.5;
imtos = zeros(size(imtom));
imtos(:,:,1) = immr; imtos(:,:,2) = immg; imtos(:,:,3) = immb;
figure(10);
clf;
hold on;
imshow(imtos);
imwrite(imtos, 'glassbottle_masked_obj.png');
