clear all; close all; clc;

% a separate file for rmc_rsd_stats, since those images are in a separate
% directory. they were also not tested in an experiment, so it is good to
% keep an organisation that reflects that they are conceptually distinct.

masks = {'../images/masks/obj_mask_correct.png', '../images/masks/obj_wout_high_mask_correct.png'};
mask_names = {'obj_mask', 'obj_wout_high_mask'};

wout_high_mask = imread('../images/masks/obj_mask_correct.png');
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

rmc_rsd_stats = [struct()]; % imgs

imgns = dir('../exp5_rmc_rsd_comp/stimuli/rmc_rsd_comp/scenes/*.png');

gcRGB = csvread('../calibration/eizo_gamma.csv');
monxyY = csvread('../calibration/eizo_chroma.csv');
DKL2RGB = dkl2rgbFromCalib('../calibration/eizo_chroma.csv');
RGB2LMS = rgb2lmsFromCalib('../calibration/eizo_mon_spectra.csv');
monxyz = xyY2XYZ(monxyY);
RGB2DKL_T = inv(DKL2RGB);

tic;
for imgc = 1:length(imgns)
    imgn = imgns(imgc);
    img = gammaCorr(im2double(imread(['../exp5_rmc_rsd_comp/stimuli/rmc_rsd_comp/scenes/' imgn.name])), gcRGB);
    
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
        
        r = img(:,:,1);
        g = img(:,:,2);
        b = img(:,:,3);
        
        rgb_gc_lin = [r(:), g(:), b(:)];
        
        rgb_gc_masked = [r(obj_pxs(:)), g(obj_pxs(:)), b(obj_pxs(:))];
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
        %         mfreq_lab = mfreq_col(lab, bkgd_pxs, obj_pxs);
        [gwdkl, ~] = gw_sd(dkl, bkgd_pxs, obj_pxs);
        [gwlab, sdlab] = gw_sd(lab, bkgd_pxs, obj_pxs);
        hv = huevar(dkl, bkgd_pxs, obj_pxs);
        [lmsm, lmssd] = gw_sd(lms, bkgd_pxs, obj_pxs);
        [rmc, rsd] = rmc_rsd(lms, bkgd_pxs, obj_pxs);
        [tau_simpler, ~] = robust_rsd(lms, bkgd_pxs, obj_pxs, 'simpler');
        [tau_general, ~] = robust_rsd(lms, bkgd_pxs, obj_pxs, 'general');
        
        rmc_rsd_stats(imgc, maskc).mask_name = mask_names{maskc};
        rmc_rsd_stats(imgc, maskc).imgn = imgn;
        rmc_rsd_stats(imgc, maskc).lrc = cr;
        rmc_rsd_stats(imgc, maskc).wpdkl = wpdkl;
        rmc_rsd_stats(imgc, maskc).wplab = wplab;
        rmc_rsd_stats(imgc, maskc).wplab_top5 = wplab_top5;
        rmc_rsd_stats(imgc, maskc).gwdkl = gwdkl;
        rmc_rsd_stats(imgc, maskc).gwlab = gwlab;
        rmc_rsd_stats(imgc, maskc).sdlab = sdlab;
        rmc_rsd_stats(imgc, maskc).gwlab_from_mean_rgb = gwlab_from_mean_rgb;
        rmc_rsd_stats(imgc, maskc).msatdkl = msatdkl;
        rmc_rsd_stats(imgc, maskc).msatlab = msatlab;
        %         rmc_rsd_stats(imgc, maskc).mfreq_lab = mfreq_lab;
        rmc_rsd_stats(imgc, maskc).hv = hv;
        rmc_rsd_stats(imgc, maskc).lmsm = lmsm;
        rmc_rsd_stats(imgc, maskc).lmssd = lmssd;
        rmc_rsd_stats(imgc, maskc).rmc = rmc;
        rmc_rsd_stats(imgc, maskc).rsd = rsd;
        rmc_rsd_stats(imgc, maskc).tau_simpler = tau_simpler;
        rmc_rsd_stats(imgc, maskc).tau_general = tau_general;
        
        % dark excluded stats
        dark_thresh = 0.05*max(squeeze(xyz(wout_high_mask & xyz(:, 2) < 50, 2)));
        dark_pxs = squeeze(xyz(:, 2)) < dark_thresh;
        reject_pxs(dark_pxs) = true;
        obj_pxs(dark_pxs) = false;
        bkgd_pxs(dark_pxs) = false;
        
        [cr, ~] = lrc(lms, bkgd_pxs, obj_pxs);
        wpdkl = wp(dkl, bkgd_pxs, obj_pxs);
        wplab = wp(lab, bkgd_pxs, obj_pxs);
        wplab_top5 = wp_top_5(lab, bkgd_pxs, obj_pxs, xyz);
        msatdkl = msat_dkl(dkl, bkgd_pxs, obj_pxs);
        msatlab = msat_lab(lab, bkgd_pxs, obj_pxs);
        %         mfreq_lab = mfreq_col(lab, bkgd_pxs, obj_pxs);
        [gwdkl, ~] = gw_sd(dkl, bkgd_pxs, obj_pxs);
        [gwlab, sdlab] = gw_sd(lab, bkgd_pxs, obj_pxs);
        hv = huevar(dkl, bkgd_pxs, obj_pxs);
        [lmsm, lmssd] = gw_sd(lms, bkgd_pxs, obj_pxs);
        [rmc, rsd] = rmc_rsd(lms, bkgd_pxs, obj_pxs);
        [tau_simpler, ~] = robust_rsd(lms, bkgd_pxs, obj_pxs, 'simpler');
        [tau_general, ~] = robust_rsd(lms, bkgd_pxs, obj_pxs, 'general');
        
        rmc_rsd_stats(imgc, maskc).dark_exc_lrc = cr;
        rmc_rsd_stats(imgc, maskc).dark_exc_wpdkl = wpdkl;
        rmc_rsd_stats(imgc, maskc).dark_exc_wplab = wplab;
        rmc_rsd_stats(imgc, maskc).dark_exc_wplab_top5 = wplab_top5;
        rmc_rsd_stats(imgc, maskc).dark_exc_gwdkl = gwdkl;
        rmc_rsd_stats(imgc, maskc).dark_exc_gwlab = gwlab;
        rmc_rsd_stats(imgc, maskc).dark_exc_sdlab = sdlab;
        rmc_rsd_stats(imgc, maskc).dark_exc_gwlab_from_mean_rgb = gwlab_from_mean_rgb;
        rmc_rsd_stats(imgc, maskc).dark_exc_msatdkl = msatdkl;
        rmc_rsd_stats(imgc, maskc).dark_exc_msatlab = msatlab;
        %         rmc_rsd_stats(imgc, maskc).dark_exc_mfreq_lab = mfreq_lab;
        rmc_rsd_stats(imgc, maskc).dark_exc_hv = hv;
        rmc_rsd_stats(imgc, maskc).dark_exc_lmsm = lmsm;
        rmc_rsd_stats(imgc, maskc).dark_exc_lmssd = lmssd;
        rmc_rsd_stats(imgc, maskc).dark_exc_rmc = rmc;
        rmc_rsd_stats(imgc, maskc).dark_exc_rsd = rsd;
        rmc_rsd_stats(imgc, maskc).dark_exc_tau_simpler = tau_simpler;
        rmc_rsd_stats(imgc, maskc).dark_exc_tau_general = tau_general;
    end
end
toc;

save('../data/rmc_rsd_stats.mat', 'rmc_rsd_stats');

%% here are four good images to choose for testing people again
% format = [rmc, rsd]
% first four are from normal, without highlight
% last two are from dark excluded, without highlight
imgs_to_test_coords = [0.622412075242242, 0.360768191612012;
                       0.344654558176497, 0.327026662259461;
                       0.633516643649885, 0.991066538187762;
                       0.338605212100397, 0.978995998293276;
                
                       0.780935828570723, 0.290205873072069;
                       0.785629115965031, 0.975065535980824];

imgs_to_test = zeros(size(imgs_to_test_coords, 1), 1);
for ic = 1:size(imgs_to_test_coords, 1)
    min_rmc_diff = Inf;
    for rc = 1:size(rmc_rsd_stats, 1)
        if ic < 5
            rmc_m = rmc_rsd_stats(rc, 2).rmc(2);
            rsd_m = rmc_rsd_stats(rc, 2).tau_general(2);
        else
            rmc_m = rmc_rsd_stats(rc, 2).dark_exc_rmc(2);
            rsd_m = rmc_rsd_stats(rc, 2).dark_exc_tau_general(2);
        end
        
        curr_rmc_m = imgs_to_test_coords(ic, 1);
        curr_rsd_m = imgs_to_test_coords(ic, 2);
        
        rmc_diff = sqrt((rmc_m - curr_rmc_m)^2 + (rsd_m - curr_rsd_m)^2);
        if rmc_diff < min_rmc_diff
            imgs_to_test(ic) = rc;
            min_rmc_diff = rmc_diff;
        end
    end
end

%% let's quickly plot rmc vs robust rsd for all of these images

% first without the highlight
figure(1);
clf;
hold on;
rmc = [rmc_rsd_stats(:, 2).rmc];
rmc = [rmc(1:3:end); rmc(2:3:end); rmc(3:3:end)];
robust_rsd = [rmc_rsd_stats(:, 2).tau_general];
robust_rsd = [robust_rsd(1:3:end); robust_rsd(2:3:end); robust_rsd(3:3:end)];
subplot(3, 3, 1);
plot(rmc(1,:), robust_rsd(1,:), 'ko');
axis square;
subplot(3, 3, 2);
plot(rmc(2,:), robust_rsd(2,:), 'ko');
axis square;
subplot(3, 3, 3);
plot(rmc(3,:), robust_rsd(3,:), 'ko');
axis square;

dark_exc_rmc = [rmc_rsd_stats(:, 2).dark_exc_rmc];
dark_exc_rmc = [dark_exc_rmc(1:3:end); dark_exc_rmc(2:3:end); dark_exc_rmc(3:3:end)];
dark_exc_robust_rsd = [rmc_rsd_stats(:, 2).dark_exc_tau_general];
dark_exc_robust_rsd = [dark_exc_robust_rsd(1:3:end); dark_exc_robust_rsd(2:3:end); dark_exc_robust_rsd(3:3:end)];
subplot(3, 3, 4);
plot(dark_exc_rmc(1,:), dark_exc_robust_rsd(1,:), 'ko');
axis square;
subplot(3, 3, 5);
plot(dark_exc_rmc(2,:), dark_exc_robust_rsd(2,:), 'ko');
axis square;
subplot(3, 3, 6);
plot(dark_exc_rmc(3,:), dark_exc_robust_rsd(3,:), 'ko');
axis square;

rmc = [rmc_rsd_stats(:, 2).rmc];
rmc = [rmc(1:3:end); rmc(2:3:end); rmc(3:3:end)];
simpler_robust_rsd = [rmc_rsd_stats(:, 2).tau_simpler];
simpler_robust_rsd = [simpler_robust_rsd(1:3:end); simpler_robust_rsd(2:3:end); simpler_robust_rsd(3:3:end)];
subplot(3, 3, 7);
plot(rmc(1,:), simpler_robust_rsd(1,:), 'ko');
axis square;
subplot(3, 3, 8);
plot(rmc(2,:), simpler_robust_rsd(2,:), 'ko');
axis square;
subplot(3, 3, 9);
plot(rmc(3,:), simpler_robust_rsd(3,:), 'ko');
axis square;

% and then, with the highlight
figure(2);
clf;
hold on;
rmc = [rmc_rsd_stats(:, 1).rmc];
rmc = [rmc(1:3:end); rmc(2:3:end); rmc(3:3:end)];
robust_rsd = [rmc_rsd_stats(:, 1).tau_general];
robust_rsd = [robust_rsd(1:3:end); robust_rsd(2:3:end); robust_rsd(3:3:end)];
subplot(3, 3, 1);
plot(rmc(1,:), robust_rsd(1,:), 'ko');
axis square;
subplot(3, 3, 2);
plot(rmc(2,:), robust_rsd(2,:), 'ko');
axis square;
subplot(3, 3, 3);
plot(rmc(3,:), robust_rsd(3,:), 'ko');
axis square;

dark_exc_rmc = [rmc_rsd_stats(:, 1).dark_exc_rmc];
dark_exc_rmc = [dark_exc_rmc(1:3:end); dark_exc_rmc(2:3:end); dark_exc_rmc(3:3:end)];
dark_exc_robust_rsd = [rmc_rsd_stats(:, 1).dark_exc_tau_general];
dark_exc_robust_rsd = [dark_exc_robust_rsd(1:3:end); dark_exc_robust_rsd(2:3:end); dark_exc_robust_rsd(3:3:end)];
subplot(3, 3, 4);
plot(dark_exc_rmc(1,:), dark_exc_robust_rsd(1,:), 'ko');
axis square;
subplot(3, 3, 5);
plot(dark_exc_rmc(2,:), dark_exc_robust_rsd(2,:), 'ko');
axis square;
subplot(3, 3, 6);
plot(dark_exc_rmc(3,:), dark_exc_robust_rsd(3,:), 'ko');
axis square;

rmc = [rmc_rsd_stats(:, 1).rmc];
rmc = [rmc(1:3:end); rmc(2:3:end); rmc(3:3:end)];
simpler_robust_rsd = [rmc_rsd_stats(:, 1).tau_simpler];
simpler_robust_rsd = [simpler_robust_rsd(1:3:end); simpler_robust_rsd(2:3:end); simpler_robust_rsd(3:3:end)];
subplot(3, 3, 7);
plot(rmc(1,:), simpler_robust_rsd(1,:), 'ko');
axis square;
subplot(3, 3, 8);
plot(rmc(2,:), simpler_robust_rsd(2,:), 'ko');
axis square;
subplot(3, 3, 9);
plot(rmc(3,:), simpler_robust_rsd(3,:), 'ko');
axis square;