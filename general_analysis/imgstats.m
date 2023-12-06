clear all; close all; clc;

expns = {'1', '2', '3', '5', '11a', '11b', '6'};
patch_app_exp = 1;
patch_dye_exp = 2;
filter_exp = 3;
rmc_rsd_comp_exp = 4;
white_walls_patch_exp = 5;
white_walls_filter_exp = 6;
robust_rsd_exp = 7;
eizo_exps = [patch_app_exp, patch_dye_exp, filter_exp, rmc_rsd_comp_exp, robust_rsd_exp]; % which exp in the iteration loop is it?
filter_exps = [filter_exp, rmc_rsd_comp_exp, white_walls_filter_exp, robust_rsd_exp];
gamma_fiasco_exps = [patch_app_exp, patch_dye_exp, white_walls_patch_exp];

masks = {'../images/masks/obj_mask_correct.png', '../images/masks/refl_mask_correct.png', '../images/masks/obj_wout_high_mask_correct.png'};
mask_names = {'obj_mask', 'refl_mask', 'obj_wout_high_mask'};
stats = [struct()]; % exps x imgs x masks

% the dark_thresh is always 5% of the max luminance in the Glaven
% (a more strict threshold; it excludes A LOT of pixels in the dark Glavens),
% but we don't want it to be biased by the highlight, so load that mask into
% a separate array for later
wout_high_mask = imread('../images/masks/obj_wout_high_mask_correct.png');
wout_high_mask_r = wout_high_mask(:, :, 1);
wout_high_mask_g = wout_high_mask(:, :, 2);
wout_high_mask_b = wout_high_mask(:, :, 3);
wout_high_mask = wout_high_mask_r;
reject_pxs = wout_high_mask_r == 127 & wout_high_mask_g == 127 & wout_high_mask_b == 127;
reject_pxs = logical(reject_pxs(:));
bkgd_pxs = wout_high_mask_r == 0 & wout_high_mask_g == 0 & wout_high_mask_b == 0;
bkgd_pxs = logical(bkgd_pxs(:));
obj_pxs = ~(reject_pxs | bkgd_pxs);
obj_pxs = logical(obj_pxs(:));
wout_high_mask(obj_pxs) = 1;
wout_high_mask = logical(wout_high_mask(:));

tic;
for expc = 1:length(expns)
    expn = expns{expc};
    imgns = dir(['../images/exp' expn '/*.png']);
    
    if any(expc == eizo_exps)
        gcRGB = csvread('../calibration/eizo_gamma.csv');
        monxyY = csvread('../calibration/eizo_chroma.csv');
        DKL2RGB = dkl2rgbFromCalib('../calibration/eizo_chroma.csv');
        RGB2LMS = rgb2lmsFromCalib('../calibration/eizo_mon_spectra.csv');
    else
        gcRGB = csvread('../calibration/oled_gamma.csv');
        monxyY = csvread('../calibration/oled_chroma.csv');
        DKL2RGB = dkl2rgbFromCalib('../calibration/oled_chroma.csv');
        RGB2LMS = rgb2lmsFromCalib('../calibration/oled_mon_spectra.csv');
    end
    monxyz = xyY2XYZ(monxyY);
    RGB2DKL_T = inv(DKL2RGB);
    
    for imgc = 1:length(imgns)
        imgn = imgns(imgc);
        img = gammaCorr(im2double(imread(['../images/exp' expn '/' imgn.name])), gcRGB);
        
        for maskc = 1:length(masks)
            mask = imread(masks{maskc});
            mask_r = mask(:, :, 1);
            mask_g = mask(:, :, 2);
            mask_b = mask(:, :, 3);
            mask = mask_r;
            reject_pxs = mask_r == 127 & mask_g == 127 & mask_b == 127;
            reject_pxs = logical(reject_pxs(:));
            bkgd_pxs = mask_r == 0 & mask_g == 0 & mask_b == 0;
            bkgd_pxs = logical(bkgd_pxs(:));
            obj_pxs = ~(reject_pxs | bkgd_pxs);
            obj_pxs = logical(obj_pxs(:));
            mask(obj_pxs) = 1;
            mask = logical(mask(:));
            
            r = img(:, :, 1); rs = r(:);
            g = img(:, :, 2); gs = g(:);
            b = img(:, :, 3); bs = b(:);
            
            rgb_gc_lin = [rs, gs, bs];
            
            rgb_gc_masked = [rs(obj_pxs(:)), gs(obj_pxs(:)), bs(obj_pxs(:))];
            mean_rgb_gc_masked = mean(rgb_gc_masked, 1);
            gwlab_from_mean_rgb = real(rgb2labRob(monxyz, mean_rgb_gc_masked));
            
            lab = real(rgb2labRob(monxyz, rgb_gc_lin));
            dkl = real(rgb2dkl(RGB2DKL_T, rgb_gc_lin'))';
            lms = real(rgb2lms(RGB2LMS, rgb_gc_lin'))';
            xyz = real(rgb2xyzRob(monxyz, rgb_gc_lin'))';
            
            [cr, ~] = lrc(lms, bkgd_pxs, obj_pxs);
            wpdkl = wp(dkl, bkgd_pxs, obj_pxs);
            wplab = wp(lab, bkgd_pxs, obj_pxs);
            wplab_top5 = wp_top_5(lab, bkgd_pxs, obj_pxs, xyz);
            msatdkl = msat_dkl(dkl, bkgd_pxs, obj_pxs);
            msatlab = msat_lab(lab, bkgd_pxs, obj_pxs);
            if maskc == 3
                mfreq_lab = mfreq_col(lab, bkgd_pxs, obj_pxs);
            else
                mfreq_lab = [NaN, NaN, NaN];
            end
            [gwdkl, ~] = gw_sd(dkl, bkgd_pxs, obj_pxs);
            [gwlab, sdlab] = gw_sd(lab, bkgd_pxs, obj_pxs);
            hv = huevar(dkl, bkgd_pxs, obj_pxs);
            [lmsm, lmssd] = gw_sd(lms, bkgd_pxs, obj_pxs);
            [rmc, rsd] = rmc_rsd(lms, bkgd_pxs, obj_pxs);
            [tau_simpler, ~] = robust_rsd(lms, bkgd_pxs, obj_pxs, 'simpler');
            [tau_general, ~] = robust_rsd(lms, bkgd_pxs, obj_pxs, 'general');
            
            stats(expc, imgc, maskc).mask_name = mask_names{maskc};
            stats(expc, imgc, maskc).imgn = imgn;
            stats(expc, imgc, maskc).lrc = cr;
            stats(expc, imgc, maskc).wpdkl = wpdkl;
            stats(expc, imgc, maskc).wplab = wplab;
            stats(expc, imgc, maskc).wplab_top5 = wplab_top5;
            stats(expc, imgc, maskc).gwdkl = gwdkl;
            stats(expc, imgc, maskc).gwlab = gwlab;
            stats(expc, imgc, maskc).sdlab = sdlab;
            stats(expc, imgc, maskc).gwlab_from_mean_rgb = gwlab_from_mean_rgb;
            stats(expc, imgc, maskc).msatdkl = msatdkl;
            stats(expc, imgc, maskc).msatlab = msatlab;
            stats(expc, imgc, maskc).mfreq_lab = mfreq_lab;
            stats(expc, imgc, maskc).hv = hv;
            stats(expc, imgc, maskc).lmsm = lmsm;
            stats(expc, imgc, maskc).lmssd = lmssd;
            stats(expc, imgc, maskc).rmc = rmc;
            stats(expc, imgc, maskc).rsd = rsd;
            stats(expc, imgc, maskc).tau_simpler = tau_simpler;
            stats(expc, imgc, maskc).tau_general = tau_general;

            % dark excluded stats
            dark_thresh = 0.05*max(squeeze(xyz(wout_high_mask, 2)));
            % dark_thresh = 3;
            dark_pxs = squeeze(xyz(:, 2)) < dark_thresh;
            reject_pxs(dark_pxs) = true;
            obj_pxs(dark_pxs) = false;
            bkgd_pxs(dark_pxs) = false;
%             no_dark_mask = mask & not_dark_pxs;
%             rs2 = rs; rs2(~no_dark_mask) = 1;
%             gs2 = gs; gs2(~no_dark_mask) = 1;
%             bs2 = bs; bs2(~no_dark_mask) = 1;
%             imshow(reshape([rs2, gs2, bs2].^(1/2), 400, 400, 3));
%             pause;
            
            [cr, ~] = lrc(lms, bkgd_pxs, obj_pxs);
            wpdkl = wp(dkl, bkgd_pxs, obj_pxs);
            wplab = wp(lab, bkgd_pxs, obj_pxs);
            wplab_top5 = wp_top_5(lab, bkgd_pxs, obj_pxs, xyz);
            msatdkl = msat_dkl(dkl, bkgd_pxs, obj_pxs);
            msatlab = msat_lab(lab, bkgd_pxs, obj_pxs);
            if maskc == 3
                mfreq_lab = mfreq_col(lab, bkgd_pxs, obj_pxs);
            else
                mfreq_lab = [NaN, NaN, NaN];
            end
            [gwdkl, ~] = gw_sd(dkl, bkgd_pxs, obj_pxs);
            [gwlab, sdlab] = gw_sd(lab, bkgd_pxs, obj_pxs);
            hv = huevar(dkl, bkgd_pxs, obj_pxs);
            [lmsm, lmssd] = gw_sd(lms, bkgd_pxs, obj_pxs);
            [rmc, rsd] = rmc_rsd(lms, bkgd_pxs, obj_pxs);
            [tau_simpler, ~] = robust_rsd(lms, bkgd_pxs, obj_pxs, 'simpler');
            [tau_general, ~] = robust_rsd(lms, bkgd_pxs, obj_pxs, 'general');
            
            stats(expc, imgc, maskc).dark_exc_lrc = cr;
            stats(expc, imgc, maskc).dark_exc_wpdkl = wpdkl;
            stats(expc, imgc, maskc).dark_exc_wplab = wplab;
            stats(expc, imgc, maskc).dark_exc_wplab_top5 = wplab_top5;
            stats(expc, imgc, maskc).dark_exc_gwdkl = gwdkl;
            stats(expc, imgc, maskc).dark_exc_gwlab = gwlab;
            stats(expc, imgc, maskc).dark_exc_sdlab = sdlab;
            stats(expc, imgc, maskc).dark_exc_gwlab_from_mean_rgb = gwlab_from_mean_rgb;
            stats(expc, imgc, maskc).dark_exc_msatdkl = msatdkl;
            stats(expc, imgc, maskc).dark_exc_msatlab = msatlab;
            stats(expc, imgc, maskc).dark_exc_mfreq_lab = mfreq_lab;
            stats(expc, imgc, maskc).dark_exc_hv = hv;
            stats(expc, imgc, maskc).dark_exc_lmsm = lmsm;
            stats(expc, imgc, maskc).dark_exc_lmssd = lmssd;
            stats(expc, imgc, maskc).dark_exc_rmc = rmc;
            stats(expc, imgc, maskc).dark_exc_rsd = rsd;
            stats(expc, imgc, maskc).dark_exc_tau_simpler = tau_simpler;
            stats(expc, imgc, maskc).dark_exc_tau_general = tau_general;
        end
    end
end
toc;

save('../data/imgstats.mat', 'stats');