clear all; close all; clc;

expns = {'1', '2', '3', '5', '11a', '11b'};
exp_names = {'patch_app', 'patch_dye', 'filter', 'rmc_rsd_comp', 'white_walls_patch', 'white_walls_filter'};
patch_app_exp = 1;
patch_dye_exp = 2;
filter_exp = 3;
rmc_rsd_comp_exp = 4;
white_walls_patch_exp = 5;
white_walls_filter_exp = 6;
eizo_exps = [patch_app_exp, patch_dye_exp, filter_exp, rmc_rsd_comp_exp]; % which exp in the iteration loop is it?
filter_exps = [filter_exp, rmc_rsd_comp_exp, white_walls_filter_exp];
gamma_fiasco_exps = [patch_app_exp, patch_dye_exp, white_walls_patch_exp];

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
gRGB_rmc_rsd_comp_exp = csvread('../calibration/eizo_gamma.csv');
gRGB_white_walls_filter_exp = csvread('../calibration/oled_gamma.csv');

vor_texture = csvread('vor.csv');
lms_absorp_curves = csvread('linss2_10e_1.csv');
lms_absorp_curves = lms_absorp_curves(1:14:end, 2:end);

global_counter = 1;

total_data = zeros(length(expns), 1);
skipped_data = zeros(length(expns), 1);

% exps x imgs x masks
load('../data/imgstats.mat');

tic;
try
    for expc = 1:length(expns)
        expn = expns{expc};
        
        if any(expc == eizo_exps)
            monxyY = csvread('../calibration/eizo_chroma.csv');
            DKL2RGB = dkl2rgbFromCalib('../calibration/eizo_chroma.csv');
            RGB2LMS = rgb2lmsFromCalib('../calibration/eizo_mon_spectra.csv');
        else
            monxyY = csvread('../calibration/oled_chroma.csv');
            DKL2RGB = dkl2rgbFromCalib('../calibration/oled_chroma.csv');
            RGB2LMS = rgb2lmsFromCalib('../calibration/oled_mon_spectra.csv');
        end
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
        
        stats_for_curr_exp = squeeze(stats(expc, :, :));
        num_masks = size(stats_for_curr_exp, 2);
        
        for obsc = 1:length(fns)
            fn = fns(obsc).name;
            load(['../data/exp' expn '/' fn], 'data', 'subID', 'trialOrder');
            
            if expc == filter_exp || expc == white_walls_filter_exp
                if strcmp(subID, 'll')
                    continue;
                end
            end
            
            redtrans_file = 'munsell_red_better.spd';
            greentrans_file = 'munsell_green_better.spd';
            bluetrans_file = 'munsell_blue_better.spd';
            yellowtrans_file = 'munsell_yellow_better.spd';
            if any(expc == eizo_exps)
                if expc == rmc_rsd_comp_exp
                    if strcmp(subID, 're')
                        continue
                    end
                    
                    if ~strcmp(subID, '02') || ~strcmp(subID, '03')
                        redtrans_file = 'munsell_red_EXTREEEMMMEE.spd';
                        greentrans_file = 'munsell_green_EXTREEEMMMEE.spd';
                        bluetrans_file = 'munsell_blue_EXTREEEMMMEE.spd';
                        yellowtrans_file = 'munsell_yellow_EXTREEEMMMEE.spd';
                    end
                end
                
                monitor_spectra = csvread('../calibration/eizo_mon_spectra.csv');
            else
                monitor_spectra = csvread('../calibration/oled_mon_spectra.csv');
            end
            [vor_filter_map, filter_idxs, redtrans, bluetrans, greentrans, yellowtrans] = ...
                gen_voronoi_filter_map(monitor_spectra, redtrans_file, bluetrans_file, greentrans_file, yellowtrans_file);
            
            for imgc = 1:size(imgns, 1)
                imgn = imgns(imgc).name;
                fparts = strsplit(imgn, '_');
                if any(expc == [patch_app_exp, patch_dye_exp, filter_exp])
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
                elseif expc == white_walls_patch_exp || expc == white_walls_filter_exp
                    if strcmp(fparts{end}, 'lighter.png')
                        hilo = 1;
                    else
                        hilo = 0;
                    end
                    illum = fparts{5};
                    body = fparts{9};
                elseif expc == rmc_rsd_comp_exp
                    % illum 3 = red
                    % illum 4 = green
                    % ld 1 = darker
                    % ld 2 = lighter
                    
                    if strcmp(fparts{6}, '3')
                        illum = 'red';
                    else
                        illum = 'green';
                    end
                    rg = fparts{8};
                    by = fparts{10};
                    if strcmp(fparts{12}, '1')
                        hilo = 1;
                        % lighter_darker = 'darker';
                    else
                        hilo = 0;
                        % lighter_darker = 'lighter';
                    end
                else
                    error('analyzing an experiment that does not have any info in the image file names?!');
                end
                
                for maskc = 1:num_masks
                    img_idx = 1;
                    for img_idx = 1:length(imgns)
                        if strcmp(imgn, stats_for_curr_exp(img_idx, maskc).imgn.name)
                            break
                        end
                    end
                    
                    idxs = find(trialOrder == imgc);
                    for idx = idxs
                        % reproduce obs match
                        if any(expc == filter_exps)
                            rg_mix = data(idx, 1);
                            by_mix = data(idx, 2);
                            ld = data(idx, 3);
                            
                            if expc == rmc_rsd_comp_exp
                                qualc = data(idx, 4);
                                qualityc(global_counter, 1) = qualc;
                            else
                                qualityc(global_counter, 1) = NaN;
                            end
                            
                            total_data(expc) = total_data(expc) + 1;
                            
                            filter = ld.*(((rg_mix.*redtrans(:,2) + (1 - rg_mix).*greentrans(:,2)) + (by_mix.*bluetrans(:,2) + (1 - by_mix).*yellowtrans(:,2))));
                            vor_col2 = vor_filter_map;
                            vor_col2(filter_idxs, :) = bsxfun(@times, vor_filter_map(filter_idxs, :), filter'.^2.2);
                            vor_col2(isnan(vor_col2) | isinf(vor_col2)) = 0.0;
                            vor_rgb = (vor_col2*lms_absorp_curves*LMS2RGB').^(1.0/2.2);
                            
                            rs = vor_rgb(:,1);
                            gs = vor_rgb(:,2);
                            bs = vor_rgb(:,3);
                            
                            rgbs_sent_to_monitor = [rs(:), gs(:), bs(:)];
                            
                            if any(rgbs_sent_to_monitor(:) > 1) || any(rgbs_sent_to_monitor(:) < 0)
                                %                     subID
                                %                     squeeze(data(x, 1, :))
                                skipped_data(expc) = skipped_data(expc) + 1;
                                continue
                            end
                            
                            if expc == filter_exp
                                gR = gRGB_filter_exp(1);
                                gG = gRGB_filter_exp(2);
                                gB = gRGB_filter_exp(3);
                            elseif expc == rmc_rsd_comp_exp
                                gR = gRGB_rmc_rsd_comp_exp(1);
                                gG = gRGB_rmc_rsd_comp_exp(2);
                                gB = gRGB_rmc_rsd_comp_exp(3);
                            elseif expc == white_walls_filter_exp
                                gR = gRGB_white_walls_filter_exp(1);
                                gG = gRGB_white_walls_filter_exp(2);
                                gB = gRGB_white_walls_filter_exp(3);
                            else
                                error('somehow making voronoi filters for an experiment that was not done with filters?!');
                            end
                            
                            clear rgbs_gc_lin;
                            rgbs_gc_lin(:, 1) = rgbs_sent_to_monitor(:, 1).^gR;
                            rgbs_gc_lin(:, 2) = rgbs_sent_to_monitor(:, 2).^gG;
                            rgbs_gc_lin(:, 3) = rgbs_sent_to_monitor(:, 3).^gB;
                            
                            rgb_gc_masked = [rgbs_gc_lin(filter_idxs, 1), rgbs_gc_lin(filter_idxs, 2), rgbs_gc_lin(filter_idxs, 3)];
                            mean_rgb_gc_masked = mean(rgb_gc_masked, 1);
                            lab_from_mean_rgb = real(rgb2labRob(monxyz, mean_rgb_gc_masked));
                            
                            lab = real(rgb2labRob(monxyz, rgbs_gc_lin'))';
                            dkl = real(rgb2dkl(RGB2DKL_T, rgbs_gc_lin'))';
                            lms = real(rgb2lms(RGB2LMS, rgbs_gc_lin'))';
                            xyz = real(rgb2xyzRob(monxyz, rgbs_gc_lin'))';
                            
                            % dark excluded stats
                            dark_thresh = 0.05*max(squeeze(xyz(filter_idxs, 2)));
                            not_dark_pxs = squeeze(xyz(:, 2)) > dark_thresh;
                            if isempty(find(not_dark_pxs & filter_idxs))
                                % some filter settings are so dark that
                                % excluding all pixels would make no sense.
                                no_dark_mask = filter_idxs;
                            else
                                no_dark_mask = filter_idxs & not_dark_pxs;
                            end
                            
                            wplab = wp(lab, no_dark_mask);
                            wplab_top5 = wp_top_5(lab, no_dark_mask);
                            [gwdkl, ~] = gw_sd(dkl, no_dark_mask);
                            [gwlab, ~] = gw_sd(lab, no_dark_mask);
                            msatdkl = msat_dkl(dkl, no_dark_mask);
                            msatlab = msat_lab(lab, no_dark_mask);
                            mfreq_lab = mfreq_col(lab, no_dark_mask);
                            [lmsm, ~] = gw_sd(lms, no_dark_mask);
                            [rmc, rsd] = rmc_rsd(lms, no_dark_mask);
                            
                            qualityc(global_counter, 1) = NaN;
                            
                            gwld_match(global_counter, :) = gwdkl(1);
                            gwrg_match(global_counter, :) = gwdkl(2);
                            gwyv_match(global_counter, :) = gwdkl(3);
                            
                            gwL_match(global_counter, :) = gwlab(1);
                            gwa_match(global_counter, :) = gwlab(2);
                            gwb_match(global_counter, :) = gwlab(3);
                            
                            gwL_from_mean_rgb_match(global_counter, :) = lab_from_mean_rgb(1);
                            gwa_from_mean_rgb_match(global_counter, :) = lab_from_mean_rgb(2);
                            gwb_from_mean_rgb_match(global_counter, :) = lab_from_mean_rgb(3);
                            
                            lm_match(global_counter, :) = lmsm(1);
                            mm_match(global_counter, :) = lmsm(2);
                            sm_match(global_counter, :) = lmsm(3);
                            
                            rmc_l_match(global_counter, :) = rmc(1);
                            rmc_m_match(global_counter, :) = rmc(2);
                            rmc_s_match(global_counter, :) = rmc(3);
                            
                            rsd_l_filter(global_counter, :) = rsd(1);
                            rsd_m_filter(global_counter, :) = rsd(2);
                            rsd_s_filter(global_counter, :) = rsd(3);
                            
                            wpL_filter(global_counter, :) = wplab(1);
                            wpa_filter(global_counter, :) = wplab(2);
                            wpb_filter(global_counter, :) = wplab(3);
                            
                            wp_top5L_filter(global_counter, :) = wplab_top5(1);
                            wp_top5a_filter(global_counter, :) = wplab_top5(2);
                            wp_top5b_filter(global_counter, :) = wplab_top5(3);
                            
                            msat_ld_filter(global_counter, :) = msatdkl(1);
                            msat_rg_filter(global_counter, :) = msatdkl(2);
                            msat_yv_filter(global_counter, :) = msatdkl(3);
                            
                            msat_L_filter(global_counter, :) = msatlab(1);
                            msat_a_filter(global_counter, :) = msatlab(2);
                            msat_b_filter(global_counter, :) = msatlab(3);
                            
                            mfreq_L_filter(global_counter, :) = mfreq_lab(1);
                            mfreq_a_filter(global_counter, :) = mfreq_lab(2);
                            mfreq_b_filter(global_counter, :) = mfreq_lab(3);
                        else
                            rgb = squeeze(data(idx, 1, :))';
                            rgb_sent_to_monitor = real(rgb.^[gR, gG, gB]);
                            rgb_gc_lin = real(rgb_sent_to_monitor.^[gR_act, gG_act, gB_act]);
                            
                            total_data(expc) = total_data(expc) + 1;
                            
                            if any(rgb_sent_to_monitor(:) > 1) || any(rgb_sent_to_monitor(:) < 0)
                                %                     subID
                                %                     squeeze(data(x, 1, :))
                                skipped_data(expc) = skipped_data(expc) + 1;
                                continue
                            end
                            
                            dkl = real(rgb2dkl(RGB2DKL_T, rgb_gc_lin'))';
                            lms = real(rgb2lms(RGB2LMS, rgb_gc_lin'))';
                            lab = real(rgb2labRob(monxyz, rgb_gc_lin));
                            
                            avg_bkgd_gc_lin = avg_bkgd_color.^[gR_act, gG_act, gB_act];
                            lms_bkgd = real(rgb2lms(RGB2LMS, avg_bkgd_gc_lin'))';
                            rmc = lms./lms_bkgd;
                            
                            lm_match(global_counter, :) = lms(1);
                            mm_match(global_counter, :) = lms(2);
                            sm_match(global_counter, :) = lms(3);
                            
                            gwld_match(global_counter, :) = dkl(1);
                            gwrg_match(global_counter, :) = dkl(2);
                            gwyv_match(global_counter, :) = dkl(3);
                            
                            gwL_match(global_counter, :) = lab(1);
                            gwa_match(global_counter, :) = lab(2);
                            gwb_match(global_counter, :) = lab(3);
                            
                            % the mean of a uniform patch is the same as
                            % any single pixel
                            gwL_from_mean_rgb_match(global_counter, :) = lab(1);
                            gwa_from_mean_rgb_match(global_counter, :) = lab(2);
                            gwb_from_mean_rgb_match(global_counter, :) = lab(3);
                            
                            qualityc(global_counter, 1) = NaN;
                            
                            rmc_l_match(global_counter, :) = rmc(1);
                            rmc_m_match(global_counter, :) = rmc(2);
                            rmc_s_match(global_counter, :) = rmc(3);
                            
                            rsd_l_filter(global_counter, :) = NaN;
                            rsd_m_filter(global_counter, :) = NaN;
                            rsd_s_filter(global_counter, :) = NaN;
                            
                            wpL_filter(global_counter, :) = NaN;
                            wpa_filter(global_counter, :) = NaN;
                            wpb_filter(global_counter, :) = NaN;
                            
                            wp_top5L_filter(global_counter, :) = NaN;
                            wp_top5a_filter(global_counter, :) = NaN;
                            wp_top5b_filter(global_counter, :) = NaN;
                            
                            msat_ld_filter(global_counter, :) = NaN;
                            msat_rg_filter(global_counter, :) = NaN;
                            msat_yv_filter(global_counter, :) = NaN;
                            
                            msat_L_filter(global_counter, :) = NaN;
                            msat_a_filter(global_counter, :) = NaN;
                            msat_b_filter(global_counter, :) = NaN;
                            
                            mfreq_L_filter(global_counter, :) = NaN;
                            mfreq_a_filter(global_counter, :) = NaN;
                            mfreq_b_filter(global_counter, :) = NaN;
                        end
                        
                        exp_name{global_counter, 1} = exp_names{expc};
                        obs_name{global_counter, 1} = subID;
                        illum_glaven{global_counter, 1} = illum;
                        body_glaven{global_counter, 1} = body;
                        hilo_glaven(global_counter, 1) = hilo;
                        mask_name{global_counter, 1} = stats_for_curr_exp(img_idx, maskc).mask_name;
                        
                        lrc_glaven(global_counter, 1) = stats_for_curr_exp(img_idx, maskc).dark_exc_lrc;
                        
                        wpld_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).dark_exc_wpdkl(1);
                        wprg_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).dark_exc_wpdkl(2);
                        wpyv_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).dark_exc_wpdkl(3);
                        
                        wpL_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).dark_exc_wplab(1);
                        wpa_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).dark_exc_wplab(2);
                        wpb_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).dark_exc_wplab(3);
                        
                        wp_top5L_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).dark_exc_wplab_top5(1);
                        wp_top5a_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).dark_exc_wplab_top5(2);
                        wp_top5b_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).dark_exc_wplab_top5(3);
                        
                        gwld_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).dark_exc_gwdkl(1);
                        gwrg_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).dark_exc_gwdkl(2);
                        gwyv_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).dark_exc_gwdkl(3);
                        
                        gwL_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).dark_exc_gwlab(1);
                        gwa_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).dark_exc_gwlab(2);
                        gwb_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).dark_exc_gwlab(3);
                        
                        gwL_from_mean_rgb_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).dark_exc_gwlab_from_mean_rgb(1);
                        gwa_from_mean_rgb_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).dark_exc_gwlab_from_mean_rgb(2);
                        gwb_from_mean_rgb_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).dark_exc_gwlab_from_mean_rgb(3);
                        
                        msat_ld_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).dark_exc_msatdkl(1);
                        msat_rg_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).dark_exc_msatdkl(2);
                        msat_yv_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).dark_exc_msatdkl(3);
                        
                        msat_L_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).dark_exc_msatlab(1);
                        msat_a_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).dark_exc_msatlab(2);
                        msat_b_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).dark_exc_msatlab(3);
                        
                        mfreq_L_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).dark_exc_mfreq_lab(1);
                        mfreq_a_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).dark_exc_mfreq_lab(2);
                        mfreq_b_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).dark_exc_mfreq_lab(3);
                        
                        hv_glaven(global_counter, 1) = stats_for_curr_exp(img_idx, maskc).dark_exc_hv;
                        
                        lm_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).dark_exc_lmsm(1);
                        mm_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).dark_exc_lmsm(2);
                        sm_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).dark_exc_lmsm(3);
                        
                        lsd_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).dark_exc_lmssd(1);
                        msd_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).dark_exc_lmssd(2);
                        ssd_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).dark_exc_lmssd(3);
                        
                        rmc_l_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).dark_exc_rmc(1);
                        rmc_m_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).dark_exc_rmc(2);
                        rmc_s_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).dark_exc_rmc(3);
                        
                        rsd_l_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).dark_exc_rsd(1);
                        rsd_m_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).dark_exc_rsd(2);
                        rsd_s_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).dark_exc_rsd(3);
                        
                        global_counter = global_counter + 1;
                    end
                end
            end
        end
    end
    
    disp('skipped data');
    skipped_data
    total_data
    skipped_data./total_data
catch e
    e
end
toc;

exp_name = categorical(exp_name);
obs_name = categorical(obs_name);
mask_name = categorical(mask_name);
illum_glaven = categorical(illum_glaven);
body_glaven = categorical(body_glaven);
hilo_glaven = categorical(hilo_glaven);

data_table = table(exp_name, ...
    obs_name, ...
    mask_name, ...
    illum_glaven, ...
    body_glaven, ...
    hilo_glaven, ...
    lrc_glaven, ...
    wpld_glaven, ...
    wprg_glaven, ...
    wpyv_glaven, ...
    wpL_glaven, ...
    wpa_glaven, ...
    wpb_glaven, ...
    wp_top5L_glaven, ...
    wp_top5a_glaven, ...
    wp_top5b_glaven, ...
    gwld_glaven, ...
    gwrg_glaven, ...
    gwyv_glaven, ...
    gwL_glaven, ...
    gwa_glaven, ...
    gwb_glaven, ...
    gwL_from_mean_rgb_glaven, ...
    gwa_from_mean_rgb_glaven, ...
    gwb_from_mean_rgb_glaven, ...
    msat_ld_glaven, ...
    msat_rg_glaven, ...
    msat_yv_glaven, ...
    msat_L_glaven, ...
    msat_a_glaven, ...
    msat_b_glaven, ...
    mfreq_L_glaven, ...
    mfreq_a_glaven, ...
    mfreq_b_glaven, ...
    hv_glaven, ...
    lm_glaven, ...
    mm_glaven, ...
    sm_glaven, ...
    lsd_glaven, ...
    msd_glaven, ...
    ssd_glaven, ...
    rmc_l_glaven, ...
    rmc_m_glaven, ...
    rmc_s_glaven, ...
    rsd_l_glaven, ...
    rsd_m_glaven, ...
    rsd_s_glaven, ...
    lm_match, ...
    mm_match, ...
    sm_match, ...
    gwld_match, ...
    gwrg_match, ...
    gwyv_match, ...
    gwL_match, ...
    gwa_match, ...
    gwb_match, ...
    gwL_from_mean_rgb_match, ...
    gwa_from_mean_rgb_match, ...
    gwb_from_mean_rgb_match, ...
    qualityc, ...
    rmc_l_match, ...
    rmc_m_match, ...
    rmc_s_match, ...
    rsd_l_filter, ...
    rsd_m_filter, ...
    rsd_s_filter, ...
    wpL_filter, ...
    wpa_filter, ...
    wpb_filter, ...
    wp_top5L_filter, ...
    wp_top5a_filter, ...
    wp_top5b_filter, ...
    msat_ld_filter, ...
    msat_rg_filter, ...
    msat_yv_filter, ...
    msat_L_filter, ...
    msat_a_filter, ...
    msat_b_filter, ...
    mfreq_L_filter, ...
    mfreq_a_filter, ...
    mfreq_b_filter);

save('../data/dark_exc_obsstats.mat', 'data_table');
