clear all; close all; clc;

expns = {'11b'};
white_walls_filter_exp = 1;
filter_exps = [white_walls_filter_exp];

masks = {'../images/masks/white_walls_mask.png'};
mask_names = {'white_walls_mask'};
stats = [struct()]; % exps x imgs x masks

tic;
for expc = 1:length(expns)
    expn = expns{expc};
    imgns = dir(['../images/exp' expn '/*.png']);
    
        gcRGB = csvread('../calibration/oled_gamma.csv');
        monxyY = csvread('../calibration/oled_chroma.csv');
        DKL2RGB = dkl2rgbFromCalib('../calibration/oled_chroma.csv');
        RGB2LMS = rgb2lmsFromCalib('../calibration/oled_mon_spectra.csv');
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
            reject_pxs = mask_r == 255 & mask_g == 0 & mask_b == 0;
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
        end
    end
end
toc;

save('../data/white_walls_imgstats.mat', 'stats');