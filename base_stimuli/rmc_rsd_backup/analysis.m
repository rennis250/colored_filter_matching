clear all; close all; clc;

bkg_mask = imread('../../../images/masks/bkg_mask.png');
bkg_mask = logical(bkg_mask(:,:,1) == 1);
bkg_mask = squeeze(reshape(bkg_mask, size(bkg_mask,1)*size(bkg_mask,2), 1));

obj_mask = imread('../../../images/masks/obj_mask.png');
obj_mask = logical(obj_mask(:,:,1) == 1);
obj_mask = squeeze(reshape(obj_mask, size(obj_mask,1)*size(obj_mask,2), 1));

mon_xyY = csvread('/home/me/my_docs/science_projects/color_filter_matching/new_structure/color_filter_matching/exp3_filter_match/calibration/Labor313EXP_Eizo_chroma_05-Dec-2016230530_Konica.csv');
DKL2RGB = dkl2rgbFromCalib('/home/me/my_docs/science_projects/color_filter_matching/new_structure/color_filter_matching/exp3_filter_match/calibration/Labor313EXP_Eizo_chroma_05-Dec-2016230530_Konica.csv');
RGB2DKL = inv(DKL2RGB);

DKL2RGB_cw = dkl2rgbFromCalib_cw(mon_xyY);
RGB2DKL_cw = inv(DKL2RGB_cw);

RGB2LMS = rgb2lmsFromCalib('/home/me/my_docs/science_projects/color_filter_matching/new_structure/color_filter_matching/exp3_filter_match/calibration/Labor313EXP_Eizo_spectra_05-Dec-2016230530_Konica.csv');

mon_XYZ_cw = xyY2XYZ_chris(mon_xyY);
wp_XYZ_cw = sum(mon_XYZ_cw, 1);
wp_xyY_cw = XYZ2xyY(wp_XYZ_cw);

fns = dir('*.png');

rmc_cw = zeros(length(fns)-2,3);
rmc_rob = zeros(length(fns)-2,3);

rsd_cw = zeros(length(fns)-2,3);
rsd_rob = zeros(length(fns)-2,3);

imc = 1;

for fnc = 1:length(fns)
	fn = fns(fnc).name;
	if strcmp(fn, 'bkg_mask.png') || strcmp(fn, 'obj_mask.png')
		continue;
	end

	img = im2double(imread(fn));
	rgb = img.^(2.2);
	dkl_img = rgb2dklIMG(RGB2DKL, rgb);

	colrs0 = reshape(rgb, size(rgb, 1)*size(rgb, 2), size(rgb, 3));
	lms_rob = rgb2lms(colrs0', RGB2LMS);

	dkl_cw = rgb2dkl_cw(RGB2DKL_cw, colrs0(:,1), colrs0(:,2), colrs0(:,3));

	colrs_from_cw_xyz_func = colourconverter(dkl_cw', 'dkl', 2, wp_xyY_cw, mon_xyY);

	bkg_lms_cw = colrs_from_cw_xyz_func.lms(bkg_mask,:);
	obj_lms_cw = colrs_from_cw_xyz_func.lms(obj_mask,:);

	bkg_lms_rob = lms_rob(:,bkg_mask)';
	obj_lms_rob = lms_rob(:,obj_mask)';

	rmc_cw(imc, :) = mean(obj_lms_cw,1)./mean(bkg_lms_cw,1);
	rmc_rob(imc, :) = mean(obj_lms_rob,1)./mean(bkg_lms_rob,1);

	rsd_cw(imc, :) = std(obj_lms_cw,1)./std(bkg_lms_cw,1);
	rsd_rob(imc, :) = std(obj_lms_rob,1)./std(bkg_lms_rob,1);

	imc = imc + 1;

	imc
end

% compare results of cw's functions and mine

figure(1);
clf;
subplot(3,2,1);
axis square;
hold on;
plot(rmc_cw(:,1), rmc_rob(:,1), 'ko');

subplot(3,2,2);
axis square;
hold on;
plot(rmc_cw(:,2), rmc_rob(:,2), 'ko');

subplot(3,2,3);
axis square;
hold on;
plot(rmc_cw(:,3), rmc_rob(:,3), 'ko');

subplot(3,2,4);
axis square;
hold on;
plot(rsd_cw(:,1), rsd_rob(:,1), 'ko');

subplot(3,2,5);
axis square;
hold on;
plot(rsd_cw(:,2), rsd_rob(:,2), 'ko');

subplot(3,2,6);
axis square;
hold on;
plot(rsd_cw(:,3), rsd_rob(:,3), 'ko');

% what is actual difference

mean(rmc_cw(:,1)./rmc_rob(:,1))
mean(rmc_cw(:,2)./rmc_rob(:,2))
mean(rmc_cw(:,3)./rmc_rob(:,3))

mean(rsd_cw(:,1)./rsd_rob(:,1))
mean(rsd_cw(:,2)./rsd_rob(:,2))
mean(rsd_cw(:,3)./rsd_rob(:,3))

std(rmc_cw(:,1)./rmc_rob(:,1))
std(rmc_cw(:,2)./rmc_rob(:,2))
std(rmc_cw(:,3)./rmc_rob(:,3))

std(rsd_cw(:,1)./rsd_rob(:,1))
std(rsd_cw(:,2)./rsd_rob(:,2))
std(rsd_cw(:,3)./rsd_rob(:,3))

% plot actual rmc vs rsd for each cone class, looking at cw's output and mine and the rust output

rust_rmc_rsd = csv2cell('test.csv');
rust_mean_lms = [rust_rmc_rsd{2:end,13}; rust_rmc_rsd{2:end,14}; rust_rmc_rsd{2:end,15}];
rust_sd_lms = [rust_rmc_rsd{2:end,16}; rust_rmc_rsd{2:end,17}; rust_rmc_rsd{2:end,18}];
bkg_or_obj_str = {rust_rmc_rsd{2:end,19}};

bkg_or_obj = logical(zeros(size(bkg_or_obj_str,2),1));
for n = 1:size(bkg_or_obj_str,2)
    bkg_or_obj(n) = strcmp(bkg_or_obj_str{n}, 'bkg_mask.png');
end
bkg = bkg_or_obj == 1;
obj = bkg_or_obj == 0;

rust_mean_obj_lms = rust_mean_lms(:,obj);
rust_mean_bkg_lms = rust_mean_lms(:,bkg);
rust_sd_obj_lms = rust_sd_lms(:,obj);
rust_sd_bkg_lms = rust_sd_lms(:,bkg);

rmc_rust = rust_mean_obj_lms./rust_mean_bkg_lms; rmc_rust = rmc_rust';
rsd_rust = rust_sd_obj_lms./rust_sd_bkg_lms; rsd_rust = rsd_rust';

figure(2);
clf;
subplot(3,3,1);
axis square;
hold on;
plot(rmc_cw(:,1), rsd_cw(:,1), 'ko');

subplot(3,3,2);
axis square;
hold on;
plot(rmc_cw(:,2), rsd_cw(:,2), 'ko');

subplot(3,3,3);
axis square;
hold on;
plot(rmc_cw(:,3), rsd_cw(:,3), 'ko');

subplot(3,3,4);
axis square;
hold on;
plot(rmc_rob(:,1), rsd_rob(:,1), 'ko');

subplot(3,3,5);
axis square;
hold on;
plot(rmc_rob(:,2), rsd_rob(:,2), 'ko');

subplot(3,3,6);
axis square;
hold on;
plot(rmc_rob(:,3), rsd_rob(:,3), 'ko');

subplot(3,3,7);
axis square;
axis([0.3 1 0.5 1]);
hold on;
plot(rmc_rust(:,1), rsd_rust(:,1), 'ko');

subplot(3,3,8);
axis square;
axis([0.5 1 0.7 1.3]);
hold on;
plot(rmc_rust(:,2), rsd_rust(:,2), 'ko');

subplot(3,3,9);
axis square;
axis([0 1 0 1]);
hold on;
plot(rmc_rust(:,3), rsd_rust(:,3), 'ko');