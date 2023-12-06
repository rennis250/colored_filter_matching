clear all; close all; clc;

expns = {'1', '2', '3', '5', '11a', '11b', '6'};
exp_names = {'patch_app', 'patch_dye', 'filter', 'rmc_rsd_comp', 'white_walls_patch', 'white_walls_filter', 'robust_rsd'};
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
gRGB_robust_rsd_comp_exp = csvread('../calibration/eizo_gamma.csv');
gRGB_white_walls_filter_exp = csvread('../calibration/oled_gamma.csv');

vor_texture = csvread('vor.csv');
lms_absorp_curves = csvread('linss2_10e_1.csv');
lms_absorp_curves = lms_absorp_curves(1:14:end, 2:end);

global_counter = 1;
filter_match_counter = 1;
filter_match_counter2 = 1;
filter_match_counter3 = 1;

rmc_to_examine_for_fig_13 = [];
rmc_counter = [1, 1];

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

            if expc == patch_dye_exp
                if strcmp(subID, '05') || strcmp(subID, 'go') || strcmp(subID, 'rr') || strcmp(subID, 'vc')
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

                if expc == robust_rsd_exp
                    redtrans_file = 'munsell_red_EXTREEEMMMEE.spd';
                    greentrans_file = 'munsell_green_EXTREEEMMMEE.spd';
                    bluetrans_file = 'munsell_blue_EXTREEEMMMEE.spd';
                    yellowtrans_file = 'munsell_yellow_EXTREEEMMMEE.spd';
                end

                monitor_spectra = csvread('../calibration/eizo_mon_spectra.csv');
            else
                monitor_spectra = csvread('../calibration/oled_mon_spectra.csv');
            end
            [vor_filter_map, filter_idxs, redtrans, bluetrans, greentrans, yellowtrans] = ...
                gen_voronoi_filter_map(monitor_spectra, redtrans_file, bluetrans_file, greentrans_file, yellowtrans_file, false);

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
                        if expc == white_walls_patch_exp
                            continue;
                        else
                            hilo = 1;
                        end
                    else
                        hilo = 0;
                    end
                    illum = fparts{5};
                    body = fparts{9};
                elseif expc == rmc_rsd_comp_exp || expc == robust_rsd_exp
                    % illum 1 = blue
                    % illum 2 = yellow
                    % illum 3 = red
                    % illum 4 = green
                    % ld 1 = darker
                    % ld 2 = lighter

                    if strcmp(fparts{6}, '1')
                        illum = 'blue';
                    elseif strcmp(fparts{6}, '2')
                        illum = 'yellow';
                    elseif strcmp(fparts{6}, '3')
                        illum = 'red';
                    else
                        illum = 'green';
                    end
                    rg = fparts{8};
                    by = fparts{10};
                    if strcmp(fparts{12}, '1')
                        hilo = 0;
                        % lighter_darker = 'darker';
                    else
                        hilo = 1;
                        % lighter_darker = 'lighter';
                    end
                    if strcmp(fparts{end}, 'dist.png')
                        bkgd = 'diff_dist.png';
                    elseif strcmp(fparts{end}, 'brighter.png')
                        if strcmp(fparts{end-1}, 'dist')
                            bkgd = 'diff_dist_brighter.png';
                        else
                            bkgd = [fparts{end-1}, '_brighter.png'];
                        end
                    else
                        bkgd = fparts{end};
                    end
                else
                    error('analyzing an experiment that does not have any info in the image file names?!');
                end

                % needed for stupid change in number of stimuli tested
                % between observers in the white walls patch experiment
                if expc == white_walls_patch_exp
                    if ~strcmp(subID, '01')
                        imgc = imgc/2;
                    end
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

                            if maskc == 1 % otherwise, we count the same data point 3 times, once for each mask
                                total_data(expc) = total_data(expc) + 1;
                            end

                            filter = ld.*(((rg_mix.*redtrans(:,2) + (1 - rg_mix).*greentrans(:,2)) + (by_mix.*bluetrans(:,2) + (1 - by_mix).*yellowtrans(:,2))));
                            vor_col2 = vor_filter_map;
                            vor_col2(filter_idxs(:), :) = bsxfun(@times, vor_filter_map(filter_idxs(:), :), filter'.^2.2);
                            vor_col2(isnan(vor_col2) | isinf(vor_col2)) = 0.0;
                            vor_rgb = (vor_col2*lms_absorp_curves*LMS2RGB').^(1.0/2.2);

                            rs = vor_rgb(:,1);
                            gs = vor_rgb(:,2);
                            bs = vor_rgb(:,3);

                            rgbs_sent_to_monitor = [rs(:), gs(:), bs(:)];

                            if any(rgbs_sent_to_monitor(:) > 1) || any(rgbs_sent_to_monitor(:) < 0)
                                %                     subID
                                %                     squeeze(data(x, 1, :))
                                if maskc == 1 % otherwise, we count the same data point 3 times, once for each mask
                                    skipped_data(expc) = skipped_data(expc) + 1;
                                end
                                continue
                            end

                            if expc == filter_exp && strcmp(illum, 'blue') && strcmp(body, 'yellow') && hilo == 0 && maskc == 1
                                ld_obs(filter_match_counter) = ld;
                                rg_mix_obs(filter_match_counter) = rg_mix;
                                by_mix_obs(filter_match_counter) = by_mix;
                                filter_match_counter = filter_match_counter + 1;
                                disp('did it!');
                            end

                            if expc == rmc_rsd_comp_exp || expc == robust_rsd_exp
                                disp({expc, rg, by, illum, hilo, bkgd, maskc});
                            end

                            save_that_rmc = 0;

                            if expc == rmc_rsd_comp_exp && strcmp(rg, '3') && strcmp(by, '1') && strcmp(illum, 'green') && hilo == 1 && strcmp(bkgd, 'diff_dist.png') && maskc == 1
                                ld_rmc_darker_obs(filter_match_counter2) = ld;
                                rg_mix_rmc_darker_obs(filter_match_counter2) = rg_mix;
                                by_mix_rmc_darker_obs(filter_match_counter2) = by_mix;
                                filter_match_counter2 = filter_match_counter2 + 1;
                                disp('did it! 2');
                                save_that_rmc = 1;
                            end

                            if expc == rmc_rsd_comp_exp && strcmp(rg, '3') && strcmp(by, '1') && strcmp(illum, 'green') && hilo == 1 && strcmp(bkgd, 'diff_dist_brighter.png') && maskc == 1
                                ld_rmc_lighter_obs(filter_match_counter3) = ld;
                                rg_mix_rmc_lighter_obs(filter_match_counter3) = rg_mix;
                                by_mix_rmc_lighter_obs(filter_match_counter3) = by_mix;
                                filter_match_counter3 = filter_match_counter3 + 1;
                                disp('did it! 3');
                                save_that_rmc = 2;
                            end

                            if expc == filter_exp
                                gR = gRGB_filter_exp(1);
                                gG = gRGB_filter_exp(2);
                                gB = gRGB_filter_exp(3);
                            elseif expc == rmc_rsd_comp_exp
                                gR = gRGB_rmc_rsd_comp_exp(1);
                                gG = gRGB_rmc_rsd_comp_exp(2);
                                gB = gRGB_rmc_rsd_comp_exp(3);
                            elseif expc == robust_rsd_exp
                                gR = gRGB_robust_rsd_comp_exp(1);
                                gG = gRGB_robust_rsd_comp_exp(2);
                                gB = gRGB_robust_rsd_comp_exp(3);
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

                            rgb_gc_masked = [rgbs_gc_lin(filter_idxs(:), 1), rgbs_gc_lin(filter_idxs(:), 2), rgbs_gc_lin(filter_idxs(:), 3)];
                            mean_rgb_gc_masked = mean(rgb_gc_masked, 1);
                            lab_from_mean_rgb = real(rgb2labRob(monxyz, mean_rgb_gc_masked));

                            lab = real(rgb2labRob(monxyz, rgbs_gc_lin));
                            dkl = real(rgb2dkl(RGB2DKL_T, rgbs_gc_lin'))';
                            lms = real(rgb2lms(RGB2LMS, rgbs_gc_lin'))';
                            xyz = real(rgb2xyzRob(monxyz, rgbs_gc_lin'))';

                            wplab = wp(lab, ~filter_idxs, filter_idxs);
                            wplab_top5 = wp_top_5(lab, ~filter_idxs, filter_idxs, xyz);
                            [gwdkl, ~] = gw_sd(dkl, ~filter_idxs, filter_idxs);
                            [gwlab, sdlab] = gw_sd(lab, ~filter_idxs, filter_idxs);
                            msatdkl = msat_dkl(dkl, ~filter_idxs, filter_idxs);
                            msatlab = msat_lab(lab, ~filter_idxs, filter_idxs);
                            mfreq_lab = mfreq_col(lab, ~filter_idxs, filter_idxs);
                            [lmsm, ~] = gw_sd(lms, ~filter_idxs, filter_idxs);
                            [rmc, rsd] = rmc_rsd(lms, ~filter_idxs, filter_idxs);
                            [tau_simpler, ~] = robust_rsd(lms, ~filter_idxs, filter_idxs, 'simpler');
                            [tau_general, ~] = robust_rsd(lms, ~filter_idxs, filter_idxs, 'general');

                            if save_that_rmc ~= 0
                                rmc_to_examine_for_fig_13(rmc_counter(save_that_rmc), save_that_rmc, :) = rmc;
                                rmc_counter(save_that_rmc) = rmc_counter(save_that_rmc) + 1;
                            end

                            ld_or_r(global_counter, 1) = ld;
                            rg_or_g(global_counter, 1) = rg_mix;
                            by_or_b(global_counter, 1) = by_mix;

                            if expc == rmc_rsd_comp_exp || expc == robust_rsd_exp
                                qualityc(global_counter, 1) = data(idx, 4);
                            else
                                qualityc(global_counter, 1) = NaN;
                            end

                            gwld_match(global_counter, :) = gwdkl(1);
                            gwrg_match(global_counter, :) = gwdkl(2);
                            gwyv_match(global_counter, :) = gwdkl(3);

                            gwL_match(global_counter, :) = gwlab(1);
                            gwa_match(global_counter, :) = gwlab(2);
                            gwb_match(global_counter, :) = gwlab(3);

                            sdL_match(global_counter, :) = sdlab(1);
                            sda_match(global_counter, :) = sdlab(2);
                            sdb_match(global_counter, :) = sdlab(3);

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

                            tau_simpler_l_filter(global_counter, :) = tau_simpler(1);
                            tau_simpler_m_filter(global_counter, :) = tau_simpler(2);
                            tau_simpler_s_filter(global_counter, :) = tau_simpler(3);

                            tau_general_l_filter(global_counter, :) = tau_general(1);
                            tau_general_m_filter(global_counter, :) = tau_general(2);
                            tau_general_s_filter(global_counter, :) = tau_general(3);

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
                        else % exp1, exp2, exp11a - that's why it's okay that i do not check for
                            % whether or not to apply gamma corr here
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

                            ld_or_r(global_counter, 1) = rgb_gc_lin(1);
                            rg_or_g(global_counter, 1) = rgb_gc_lin(2);
                            by_or_b(global_counter, 1) = rgb_gc_lin(3);

                            % there were no quality ratings given for
                            % uniform patch experiments
                            qualityc(global_counter, 1) = NaN;

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

                            sdL_match(global_counter, :) = NaN;
                            sda_match(global_counter, :) = NaN;
                            sdb_match(global_counter, :) = NaN;

                            rsd_l_filter(global_counter, :) = NaN;
                            rsd_m_filter(global_counter, :) = NaN;
                            rsd_s_filter(global_counter, :) = NaN;

                            tau_simpler_l_filter(global_counter, :) = NaN;
                            tau_simpler_m_filter(global_counter, :) = NaN;
                            tau_simpler_s_filter(global_counter, :) = NaN;

                            tau_general_l_filter(global_counter, :) = NaN;
                            tau_general_m_filter(global_counter, :) = NaN;
                            tau_general_s_filter(global_counter, :) = NaN;

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
                        img_name{global_counter, 1} =  imgn;
                        illum_glaven{global_counter, 1} = illum;
                        if expc == rmc_rsd_comp_exp || expc == robust_rsd_exp
                            body_glaven{global_counter, 1} = 'no_body_name';
                            rg_glaven(global_counter, 1) = str2double(rg);
                            by_glaven(global_counter, 1) = str2double(by);
                            bkgd_glaven{global_counter, 1} = bkgd;
                        else
                            body_glaven{global_counter, 1} = body;
                            rg_glaven(global_counter, 1) = -1; % indicates that this is not relevant for non-rmc_rsd_comp exps
                            by_glaven(global_counter, 1) = -1;
                            bkgd_glaven{global_counter, 1} = 'not_a_rmc_rsd_comp_image';
                        end
                        hilo_glaven(global_counter, 1) = hilo;
                        mask_name{global_counter, 1} = stats_for_curr_exp(img_idx, maskc).mask_name;

                        lrc_glaven(global_counter, 1) = stats_for_curr_exp(img_idx, maskc).lrc;

                        wpld_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).wpdkl(1);
                        wprg_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).wpdkl(2);
                        wpyv_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).wpdkl(3);

                        wpL_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).wplab(1);
                        wpa_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).wplab(2);
                        wpb_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).wplab(3);

                        wp_top5L_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).wplab_top5(1);
                        wp_top5a_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).wplab_top5(2);
                        wp_top5b_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).wplab_top5(3);

                        gwld_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).gwdkl(1);
                        gwrg_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).gwdkl(2);
                        gwyv_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).gwdkl(3);

                        gwL_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).gwlab(1);
                        gwa_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).gwlab(2);
                        gwb_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).gwlab(3);

                        sdL_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).sdlab(1);
                        sda_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).sdlab(2);
                        sdb_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).sdlab(3);

                        gwL_from_mean_rgb_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).gwlab_from_mean_rgb(1);
                        gwa_from_mean_rgb_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).gwlab_from_mean_rgb(2);
                        gwb_from_mean_rgb_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).gwlab_from_mean_rgb(3);

                        msat_ld_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).msatdkl(1);
                        msat_rg_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).msatdkl(2);
                        msat_yv_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).msatdkl(3);

                        msat_L_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).msatlab(1);
                        msat_a_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).msatlab(2);
                        msat_b_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).msatlab(3);

                        mfreq_L_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).mfreq_lab(1);
                        mfreq_a_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).mfreq_lab(2);
                        mfreq_b_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).mfreq_lab(3);

                        hv_glaven(global_counter, 1) = stats_for_curr_exp(img_idx, maskc).hv;

                        lm_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).lmsm(1);
                        mm_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).lmsm(2);
                        sm_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).lmsm(3);

                        lsd_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).lmssd(1);
                        msd_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).lmssd(2);
                        ssd_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).lmssd(3);

                        rmc_l_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).rmc(1);
                        rmc_m_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).rmc(2);
                        rmc_s_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).rmc(3);

                        rsd_l_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).rsd(1);
                        rsd_m_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).rsd(2);
                        rsd_s_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).rsd(3);

                        tau_simpler_l_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).tau_simpler(1);
                        tau_simpler_m_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).tau_simpler(2);
                        tau_simpler_s_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).tau_simpler(3);

                        tau_general_l_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).tau_general(1);
                        tau_general_m_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).tau_general(2);
                        tau_general_s_glaven(global_counter, :) = stats_for_curr_exp(img_idx, maskc).tau_general(3);

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
bkgd_glaven = categorical(bkgd_glaven);
img_name = categorical(img_name);

data_table = table(exp_name, ...
    obs_name, ...
    img_name, ...
    mask_name, ...
    illum_glaven, ...
    body_glaven, ...
    hilo_glaven, ...
    rg_glaven, ...
    by_glaven, ...
    bkgd_glaven, ...
    ld_or_r, ...
    rg_or_g, ...
    by_or_b, ...
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
    tau_simpler_l_glaven, ...
    tau_simpler_m_glaven, ...
    tau_simpler_s_glaven, ...
    tau_general_l_glaven, ...
    tau_general_m_glaven, ...
    tau_general_s_glaven, ...
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
    tau_simpler_l_filter, ...
    tau_simpler_m_filter, ...
    tau_simpler_s_filter, ...
    tau_general_l_filter, ...
    tau_general_m_filter, ...
    tau_general_s_filter, ...
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

save('../data/obsstats.mat', 'data_table');

disp('needed for re-creating filter in fig 5');
[mean(ld_obs), mean(rg_mix_obs), mean(by_mix_obs)]

% there's one response that is clearly an outlier as it has a large
% influence on the result
disp('needed for re-creating filter in fig 15 (bottom left)');
[mean(ld_rmc_darker_obs(ld_rmc_darker_obs > 0)), mean(rg_mix_rmc_darker_obs(ld_rmc_darker_obs > 0)), mean(by_mix_rmc_darker_obs(ld_rmc_darker_obs > 0))]

% there's one response that is clearly an outlier as it has a large
% influence on the result
disp('needed for re-creating filter in fig 15 (bottom right)');
[mean(ld_rmc_lighter_obs(ld_rmc_lighter_obs > 0)), mean(rg_mix_rmc_lighter_obs(ld_rmc_lighter_obs > 0)), mean(by_mix_rmc_lighter_obs(ld_rmc_lighter_obs > 0))]