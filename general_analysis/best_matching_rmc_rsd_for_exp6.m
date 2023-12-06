clear all; close all; clc;

% exp2 - eizo
gR = 2.1102;
gG = 2.1243;
gB = 2.0170;

% exp2 - eizo
gR_act = 2.1442;
gG_act = 2.1514;
gB_act = 1.9483;

gRGB_rmc_rsd_comp_exp = csvread('../calibration/eizo_gamma.csv');
gR = gRGB_rmc_rsd_comp_exp(1);
gG = gRGB_rmc_rsd_comp_exp(2);
gB = gRGB_rmc_rsd_comp_exp(3);

vor_texture = csvread('vor.csv');
lms_absorp_curves = csvread('linss2_10e_1.csv');
lms_absorp_curves = lms_absorp_curves(1:14:end, 2:end);

monxyY = csvread('../calibration/eizo_chroma.csv');
DKL2RGB = dkl2rgbFromCalib('../calibration/eizo_chroma.csv');
RGB2LMS = rgb2lmsFromCalib('../calibration/eizo_mon_spectra.csv');

LMS2RGB = inv(RGB2LMS);
monxyz = xyY2XYZ(monxyY);
RGB2DKL_T = inv(DKL2RGB);

redtrans_file = '../base_stimuli/spectra/munsell_red_EXTREEEMMMEE.spd';
greentrans_file = '../base_stimuli/spectra/munsell_green_EXTREEEMMMEE.spd';
bluetrans_file = '../base_stimuli/spectra/munsell_blue_EXTREEEMMMEE.spd';
yellowtrans_file = '../base_stimuli/spectra/munsell_yellow_EXTREEEMMMEE.spd';

monitor_spectra = csvread('../calibration/eizo_mon_spectra.csv');

[vor_filter_map, filter_idxs, redtrans, bluetrans, greentrans, yellowtrans] = ...
    gen_voronoi_filter_map(monitor_spectra, redtrans_file, bluetrans_file, greentrans_file, yellowtrans_file, false);

%% let's try generating "all" possible filters for our experiment and see which RSD and RMC combos we get out

rmc_sphere = [];
rsd_sphere = [];
ld_rg_bys = [];
for ld = linspace(0, 1, 20)
    for r = linspace(0, 1, 20)
        for theta = linspace(0, 2*pi, 20)
            rg_mix = r*cos(theta);
            by_mix = r*sin(theta);
            
            filter = ld.*(((rg_mix.*redtrans(:,2) + (1 - rg_mix).*greentrans(:,2)) + (by_mix.*bluetrans(:,2) + (1 - by_mix).*yellowtrans(:,2))));
            vor_col2 = vor_filter_map;
            vor_col2(filter_idxs(:), :) = bsxfun(@times, vor_filter_map(filter_idxs(:), :), filter'.^2.2);
            vor_col2(isnan(vor_col2) | isinf(vor_col2)) = 0.0;
            vor_rgb = (vor_col2*lms_absorp_curves*LMS2RGB').^(1.0/2.2);
            
            rs = vor_rgb(:, 1);
            gs = vor_rgb(:, 2);
            bs = vor_rgb(:, 3);
            
            rgbs_sent_to_monitor = [rs(:), gs(:), bs(:)];
            
            if any(rgbs_sent_to_monitor > 1 | rgbs_sent_to_monitor < 0)
                continue;
            end
            
            % simulate the monitor's action of clamping RGB values to [0, 1]
            rgbs_sent_to_monitor(rgbs_sent_to_monitor > 1) = 1;
            rgbs_sent_to_monitor(rgbs_sent_to_monitor < 0) = 0;
            
            clear rgbs_gc_lin;
            rgbs_gc_lin(:, 1) = rgbs_sent_to_monitor(:, 1).^gR;
            rgbs_gc_lin(:, 2) = rgbs_sent_to_monitor(:, 2).^gG;
            rgbs_gc_lin(:, 3) = rgbs_sent_to_monitor(:, 3).^gB;
            
            lms = real(rgb2lms(RGB2LMS, rgbs_gc_lin'))';
            [rmc, ~] = rmc_rsd(lms, ~filter_idxs, filter_idxs);
            [tau_general, ~] = robust_rsd(lms, bkgd_pxs, obj_pxs, 'general');
            
            ld_rg_bys = [ld_rg_bys; ld, rg_mix, by_mix];
            rmc_sphere = [rmc_sphere; rmc];
            rsd_sphere = [rsd_sphere; tau_general];
        end
    end
end

%% okay, now let's see what happens, when we find the filters that are most similar to 
% the rmc_rsd_comp glavens in the paper

rmc_rsd_comp_ins = dir('../exp6_rmc_robust_rsd_comp/stimuli/images/*png');
rmc_rsd_comp_images = zeros(9, 400, 400, 3);
for x = 1:10    
    rmc_rsd_comp_images = im2double(imread(['../exp6_rmc_robust_rsd_comp/stimuli/images/' rmc_rsd_comp_ins(x).name]));
end

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

figure;
for gc = 1:size(rmc_rsd_comp_images, 1)
	img = gammaCorr(im2double(imread(rmc_rsd_comp_images(gc).name), [gR, gG, gB]);
        
    r = img(:,:,1);
    g = img(:,:,2);
    b = img(:,:,3);
        
    rgb_gc_lin = [r(:), g(:), b(:)];

	lms = real(rgb2lms(RGB2LMS, rgb_gc_lin'))';

    [rmc, ~] = rmc_rsd(lms, bkgd_pxs, obj_pxs);
    [tau_general, ~] = robust_rsd(lms, bkgd_pxs, obj_pxs, 'general');

    rmc_eqdiff = sqrt(sum((rmc - rmc_sphere).^2, 2));
    [~, rmc_idx] = min(rmc_eqdiff);
    
    rsd_eqdiff = sqrt(sum((tau_general - rsd_sphere).^2, 2));
    [~, rsd_idx] = min(rsd_eqdiff);
    
    ld = ld_rg_bys(rsd_idx, 1);
    rg_mix = ld_rg_bys(rsd_idx, 2);
    by_mix = ld_rg_bys(rsd_idx, 3);
    filter = ld.*(((rg_mix.*redtrans(:,2) + (1 - rg_mix).*greentrans(:,2)) + (by_mix.*bluetrans(:,2) + (1 - by_mix).*yellowtrans(:,2))));
    vor_col2 = vor_filter_map;
    vor_col2(filter_idxs(:), :) = bsxfun(@times, vor_filter_map(filter_idxs(:), :), filter'.^2.2);
    vor_col2(isnan(vor_col2) | isinf(vor_col2)) = 0.0;
    vor_rgb_best_rsd = (vor_col2*lms_absorp_curves*LMS2RGB').^(1.0/2.2);
    
    ld = ld_rg_bys(rmc_idx, 1);
    rg_mix = ld_rg_bys(rmc_idx, 2);
    by_mix = ld_rg_bys(rmc_idx, 3);
    filter = ld.*(((rg_mix.*redtrans(:,2) + (1 - rg_mix).*greentrans(:,2)) + (by_mix.*bluetrans(:,2) + (1 - by_mix).*yellowtrans(:,2))));
    vor_col2 = vor_filter_map;
    vor_col2(filter_idxs(:), :) = bsxfun(@times, vor_filter_map(filter_idxs(:), :), filter'.^2.2);
    vor_col2(isnan(vor_col2) | isinf(vor_col2)) = 0.0;
    vor_rgb_best_rmc = (vor_col2*lms_absorp_curves*LMS2RGB').^(1.0/2.2);

    clf;
    
    subplot(2, 3, 1);
    imshow(squeeze(rmc_rsd_comp_images(good_idxs(gc), :, :, :)));
    subplot(2, 3, 4);
    imshow(reshape(vor_rgb_best_rmc, 256, 256, 3));
    title('rmc');
    subplot(2, 3, 6);
    imshow(reshape(vor_rgb_best_rsd, 256, 256, 3));
    title('rsd');
    
    export_fig(['../figures/best_rmc_robust_rsd_filter_' good_idxs(gc) '.png']);
    
    pause;
end