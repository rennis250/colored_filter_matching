clear all; close all; clc;

expns = {'2'};
exp_names = {'patch_dye'};
patch_dye_exp = 1;
eizo_exps = [patch_dye_exp]; % which exp in the iteration loop is it?
gamma_fiasco_exps = [patch_dye_exp];

avg_bkgd_color = [0.8, 0.8, 0.8];

% exp1, exp2, exp11a
gRs = [2.1102, 2.1102, 2.1102];
gGs = [2.1243, 2.1243, 2.1243];
gBs = [2.0170, 2.0170, 2.0170];

% exp1, exp2, exp11a
gR_acts = [2.1442, 2.1442, 2.1101];
gG_acts = [2.1514, 2.1514, 2.1244];
gB_acts = [1.9483, 1.9483, 2.0172];

gRGB_filter_exp = csvread('../calibration/eizo_gamma.csv');

expc = 1;
expn = expns{expc};

monxyY = csvread('../calibration/eizo_chroma.csv');
DKL2RGB = dkl2rgbFromCalib('../calibration/eizo_chroma.csv');
RGB2LMS = rgb2lmsFromCalib('../calibration/eizo_mon_spectra.csv');
LMS2RGB = inv(RGB2LMS);
monxyz = xyY2XYZ(monxyY);
RGB2DKL_T = inv(DKL2RGB);

is_it_a_gamma_corr_fiasco_exp = expc == gamma_fiasco_exps;
if any(is_it_a_gamma_corr_fiasco_exp)
  disp(['doing a gamma corr fiasco exp: ', expn]);

  which_exp_is_it = find(is_it_a_gamma_corr_fiasco_exp);
  if length(which_exp_is_it) > 1
  	error('error in selecting gamma corr fiasco exp: too many possibilities?!');
  end
  gR = gRs(which_exp_is_it);
  gG = gGs(which_exp_is_it);
  gB = gBs(which_exp_is_it);

  gR_act = gR_acts(which_exp_is_it);
  gG_act = gG_acts(which_exp_is_it);
  gB_act = gB_acts(which_exp_is_it);
end

fns = dir(['../data/exp' expn '/*.mat']);
imgns = dir(['../images/exp' expn '/*.png']);

% exps x imgs x masks
load('../data/imgstats.mat');

stats_for_curr_exp = squeeze(stats(expc, :, :));

rgbs = [];
goodics = [];
illums = {};
hilos = [];
bodies = [];
ic = 1;
for imgc = 1:length(imgns)
  imgn = imgns(imgc).name;
  fparts = strsplit(imgn, '_');
  if strcmp(fparts{end}, 'lighter.png')
    if strcmp(fparts{6}, 'red') || strcmp(fparts{6}, 'green')
      continue;
    end
    illum = fparts{6};
    body = fparts{9};
    hilo = 1;
  else
    if strcmp(fparts{7}, 'red') || strcmp(fparts{7}, 'green')
      continue;
    end
    illum = fparts{7};
    body = fparts{10};
    hilo = 0;
  end

  oc = 1;
  for obsc = 1:length(fns)
    fn = fns(obsc).name;
    load(['../data/exp' expn '/' fn], 'data', 'subID', 'trialOrder');

    if expc == patch_dye_exp
      if strcmp(subID, '05') || strcmp(subID, 'go') || strcmp(subID, 'rr') || strcmp(subID, 'vc')
        continue;
      end
    end

    idxs = find(trialOrder == imgc);
    trialc = 1;
    for idx = idxs
      % reproduce obs match
      rgbs(ic, oc, trialc, :) = squeeze(data(idx, 1, :))';
      trialc = trialc + 1;
    end
    oc = oc + 1;
  end

  goodics = [goodics; imgc];
  illums{ic} = illum;
  hilos(ic) = hilo;
  bodies{ic} = body;
  ic = ic + 1;
end

rgb = squeeze(mean(rgbs, 3)); % obs mean
rgbm = squeeze(mean(rgb, 2)); % pop mean

the_covms_ld_yv = cell(1, size(rgbs, 1));
the_covms_rg_yv = cell(1, size(rgbs, 1));
for ic = 1:length(goodics)
	temp_dkl = real(rgb2dkl(RGB2DKL_T, squeeze(rgb(ic, :, :))')');
	temp_ld_yv = temp_dkl(:, [3 1]);
	the_covms_ld_yv{ic} = cov(temp_ld_yv);
	temp_rg_yv = temp_dkl(:, [2 3]);
	the_covms_rg_yv{ic} = cov(temp_rg_yv);
end

rgb_sent_to_monitor = real(rgbm.^[gR, gG, gB]);
rgb_gc_lin = real(rgb_sent_to_monitor.^[gR_act, gG_act, gB_act]);
dkls = real(rgb2dkl(RGB2DKL_T, rgb_gc_lin'))';

%% then load the corresponding glaven and calculate convergence for it

mask = imread('../images/masks/obj_mask_correct.png');
mask_r = mask(:, :, 1);
mask_g = mask(:, :, 2);
mask_b = mask(:, :, 3);
mask = mask_r;
reject_pxs = mask_r == 127 & mask_g == 127 & mask_b == 127;
reject_pxs = logical(reject_pxs(:));
bkgd_pxs = mask_r == 0 & mask_g == 0 & mask_b == 0;
bkgd_pxs = logical(bkgd_pxs(:));
obj_pxs = ~(reject_pxs | bkgd_pxs);
obj_pxs = obj_pxs(:);
mask(obj_pxs) = 1;
mask = logical(mask(:));

h = figure(1);
fullfig(h);
for ic = 1:length(the_covms_ld_yv)
	imgn = imgns(goodics(ic)).name;
        orig_img_with_filter = im2double(imread(['../images/exp' expn '/' imgn]));
	img_with_filter = gammaCorr(orig_img_with_filter, gRGB_filter_exp);
	img_wout_filter = gammaCorr(im2double(imread(['../images/empty_scene_' illums{ic} '_illum.png'])), gRGB_filter_exp);

	rwof = img_wout_filter(:, :, 1); rwof = rwof(:); r_wout_filter = rwof(obj_pxs(:));
	gwof = img_wout_filter(:, :, 2); gwof = gwof(:); g_wout_filter = gwof(obj_pxs(:));
	bwof = img_wout_filter(:, :, 3); bwof = bwof(:); b_wout_filter = bwof(obj_pxs(:));

	dkl_wout_filter = real(rgb2dkl(RGB2DKL_T, [r_wout_filter(:), g_wout_filter(:), b_wout_filter(:)]'))';

	r = img_with_filter(:, :, 1); r = r(:); r_img = r(obj_pxs(:));
	g = img_with_filter(:, :, 2); g = g(:); g_img = g(obj_pxs(:));
	b = img_with_filter(:, :, 3); b = b(:); b_img = b(obj_pxs(:));

	dkl_img = real(rgb2dkl(RGB2DKL_T, [r_img(:), g_img(:), b_img(:)]'))';

	[ldgr_d1, ldgr_d3, rggr_d1, rggr_d2, yvgr_d2, yvgr_d3, ld_interp_d1, ld_interp_d3, rg_interp_d1, rg_interp_d2, yv_interp_d2, yv_interp_d3] = divergence_map(dkl_wout_filter, dkl_img);

	dkl = dkls(ic, :);

	clf;
        set(gcf, 'color', 'w');

        subplot(2, 2, 1);
        hold on;
        imshow(orig_img_with_filter);

        subplot(2, 2, 2);
        hold on;
        pmatch = 0.8*ones(400, 400, 3);
        pmatch(end/2-50:end/2+50, end/2-50:end/2+50, 1) = rgb_sent_to_monitor(ic, 1);
        pmatch(end/2-50:end/2+50, end/2-50:end/2+50, 2) = rgb_sent_to_monitor(ic, 2);
        pmatch(end/2-50:end/2+50, end/2-50:end/2+50, 3) = rgb_sent_to_monitor(ic, 3);
        imshow(pmatch);

	subplot(2, 2, 3);
	hold on;
	h1 = streamslice(yvgr_d3, ldgr_d3, yv_interp_d3, ld_interp_d3);
	set(h1, 'Color', [0, 0, 0]);
	plot(dkl(3), dkl(1), 'rs');
	error_ellipse(the_covms_ld_yv{ic}, [dkl(3), dkl(1)], 0.3934, 'style', 'r-');
	if hilos(ic) == 1
		for x = 1:length(the_covms_ld_yv)
			if hilos(x) == 0 && strcmp(illums{ic}, illums{x}) && strcmp(bodies{ic}, bodies{x})
				break;
			end
		end
		companion_dkl_to_plot = dkls(x, :);
		plot(companion_dkl_to_plot(3), companion_dkl_to_plot(1), 'bs');
		error_ellipse(the_covms_ld_yv{x}, [companion_dkl_to_plot(3), companion_dkl_to_plot(1)], 0.3934, 'style', 'b-');
	end
	ylabel('L+M+S');
	xlabel('S-(L+M)');
	axis square;
	
	subplot(2, 2, 4);
	hold on;
	h1 = streamslice(rggr_d2, yvgr_d2, rg_interp_d2, yv_interp_d2);
	set(h1, 'Color', [0, 0, 0]);
	plot(dkl(2), dkl(3), 'rs');
	error_ellipse(the_covms_rg_yv{ic}, [dkl(2), dkl(3)], 0.3934, 'style', 'r-');
	if hilos(ic) == 1
		plot(companion_dkl_to_plot(2), companion_dkl_to_plot(3), 'bs');
		error_ellipse(the_covms_rg_yv{x}, [companion_dkl_to_plot(2), companion_dkl_to_plot(3)], 0.3934, 'style', 'b-');
	end
	ylabel('S-(L+M)');
	xlabel('L-M');
	axis square;

	export_fig(['../figures/glaven_patch_conv_comp_' num2str(imgn) '.tif'], '-m2', '-r500');
end

