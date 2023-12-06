clear all; close all; clc;

gcRGB = csvread('../calibration/eizo_gamma.csv');

DKL2RGB = dkl2rgbFromCalib('../calibration/eizo_chroma.csv');
RGB2DKL_T = inv(DKL2RGB);

img_with_filter = gammaCorr(im2double(imread(['../images/patterned/mitsuba_caustics_rgb_trans_illum_blue_refl_munsell_red_lighter.png'])), gcRGB);
img_with_plastic = gammaCorr(im2double(imread(['../images/convergence_diffuse_comparision_blender/mitsuba_plastic_blue_red_lighter.png'])), gcRGB);
img_with_lightest_glass = gammaCorr(im2double(imread(['../images/patterned/mitsuba_caustics_rgb_trans_lightest_illum_blue_refl_munsell_red_lightest.png'])), gcRGB);
img_with_dark_glass = gammaCorr(im2double(imread(['../images/convergence_diffuse_comparision_blender/mitsuba_dark_glass_blue_red_lighter.png'])), gcRGB);
img_with_dark_plastic = gammaCorr(im2double(imread(['../images/convergence_diffuse_comparision_blender/mitsuba_darker_plastic_blue_red_lighter.png'])), gcRGB);
img_wout_filter = gammaCorr(im2double(imread('../images/empty_scene_blue_illum.png')), gcRGB);

images = {img_with_filter, img_with_plastic, img_with_dark_glass, img_with_dark_plastic, img_with_lightest_glass};
names = {'normal filter', 'plastic', 'dark filter', 'dark plastic', 'lightest glass'};
fnames = {'normal', 'plastic', 'dark_filter', 'dark_plastic', 'lightest_plastic'};

mask = imread('../images/masks/obj_mask.png'); mask = logical(mask);

rwof = img_wout_filter(:, :, 1); rwof = rwof(:); r_wout_filter = rwof(mask(:));
gwof = img_wout_filter(:, :, 2); gwof = gwof(:); g_wout_filter = gwof(mask(:));
bwof = img_wout_filter(:, :, 3); bwof = bwof(:); b_wout_filter = bwof(mask(:));

dkl_wout_filter = real(rgb2dkl(RGB2DKL_T, [r_wout_filter(:), g_wout_filter(:), b_wout_filter(:)]'))';

% for ic = 1:length(images)
for ic = 1
	img = images{ic};

	r = img(:, :, 1); r = r(:); r_img = r(mask(:));
	g = img(:, :, 2); g = g(:); g_img = g(mask(:));
	b = img(:, :, 3); b = b(:); b_img = b(mask(:));

	dkl_img = real(rgb2dkl(RGB2DKL_T, [r_img(:), g_img(:), b_img(:)]'))';

	[ldgr_d1, ldgr_d3, rggr_d1, rggr_d2, yvgr_d2, yvgr_d3, ld_interp_d1, ld_interp_d3, rg_interp_d1, rg_interp_d2, yv_interp_d2, yv_interp_d3] = divergence_map(dkl_wout_filter, dkl_img);
    
	figure;
% 	subplot(1,3,1);
	streamslice(rggr_d1, ldgr_d1, rg_interp_d1, ld_interp_d1, 3);
	xlabel('RG');
	ylabel('LD');
	axis square;
	export_fig(['streamslice_' fnames{ic} '_ld_rg.pdf']);
    
    figure;
% 	subplot(1,3,2);
	streamslice(rggr_d2, yvgr_d2, rg_interp_d2, yv_interp_d2);
	xlabel('RG');
	ylabel('YV');
	axis square;
	export_fig(['streamslice_' fnames{ic} '_rg_yv.pdf']);
    
    figure;
% 	subplot(1,3,3);
	streamslice(yvgr_d3, ldgr_d3, yv_interp_d3, ld_interp_d3, 8, 'method', 'nearest');
	xlabel('YV');
	ylabel('LD');
	axis square;
	set(gcf, 'units', 'inches', 'position', 1.5*[3, 0, 9, 3]);
    set(gcf, 'color', 'w');
	export_fig(['streamslice_' fnames{ic} '_ld_yv.pdf']);

	figure;
% 	subplot(3,1,1);
% 	quiver(rggr_d1, ldgr_d1, rg_interp_d1, ld_interp_d1, 5);
% 	xlabel('RG');
% 	ylabel('LD');
% 	axis square;
% 	subplot(3,1,2);
    quiver(rggr_d2, yvgr_d2, rg_interp_d2, yv_interp_d2, 8);
	xlabel('L-M');
	ylabel('S-(L+M)');
    axis square;
% 	subplot(3,1,3);
%     quiver(yvgr_d3, ldgr_d3, yv_interp_d3, ld_interp_d3, 5);
% 	xlabel('YV');
% 	ylabel('LD');
%     axis square;
% 	set(gcf, 'units', 'inches', 'position', 1.5*[3, 0, 3, 9]);
%     set(gcf, 'color', 'w');
	export_fig(['vector_field_interp_' fnames{ic} '.pdf']);

	div_rg_ld = divergence(rggr_d1, ldgr_d1, rg_interp_d1, ld_interp_d1);
	div_rg_yv = divergence(rggr_d2, yvgr_d2, rg_interp_d2, yv_interp_d2);
	div_yv_ld = divergence(yvgr_d3, ldgr_d3, yv_interp_d3, ld_interp_d3);
    
	figure;
	subplot(3,1,1);
	pcolor(rggr_d1, ldgr_d1, div_rg_ld);
	xlabel('RG');
	ylabel('LD');
	shading interp;
	colorbar;
	colorposneg;
	axis square;
	subplot(3,1,2);
	pcolor(rggr_d2, yvgr_d2, div_rg_yv);
	xlabel('RG');
	ylabel('YV');
	shading interp;
	colorbar;
	colorposneg;
	axis square;
	subplot(3,1,3);
	pcolor(yvgr_d3, ldgr_d3, div_yv_ld);
	xlabel('YV');
	ylabel('LD');
	shading interp;
	colorbar;
	colorposneg;
	axis square;
	set(gcf, 'units', 'inches', 'position', 1.5*[3, 0, 3, 9]);
    set(gcf, 'color', 'w');
	export_fig(['divergence_' fnames{ic} '.tiff']);
end
