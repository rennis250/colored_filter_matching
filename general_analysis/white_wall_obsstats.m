clear all; close all; clc;

expns = {'11b'};
exp_names = {'white_walls_filter'};
white_walls_filter_exp = 1;
filter_exps = [white_walls_filter_exp];

avg_bkgd_color = [0.8, 0.8, 0.8];

gRGB_white_walls_filter_exp = csvread('../calibration/oled_gamma.csv');

vor_texture = csvread('vor.csv');
lms_absorp_curves = csvread('linss2_10e_1.csv');
lms_absorp_curves = lms_absorp_curves(1:14:end, 2:end);

global_counter = 1;

total_data = zeros(length(expns), 1);
skipped_data = zeros(length(expns), 1);

% exps x imgs x masks
load('../data/white_walls_imgstats.mat');

tic;
try
    for expc = 1:length(expns)
        expn = expns{expc};

            monxyY = csvread('../calibration/oled_chroma.csv');
            DKL2RGB = dkl2rgbFromCalib('../calibration/oled_chroma.csv');
            RGB2LMS = rgb2lmsFromCalib('../calibration/oled_mon_spectra.csv');
        LMS2RGB = inv(RGB2LMS);
        monxyz = xyY2XYZ(monxyY);
        RGB2DKL_T = inv(DKL2RGB);

        fns = dir(['../data/exp' expn '/*.mat']);
        imgns = dir(['../images/exp' expn '/*.png']);

        stats_for_curr_exp = squeeze(stats(expc, :, :));
        % num_masks = size(stats_for_curr_exp, 2);
        num_masks = 1;

        for obsc = 1:length(fns)
            fn = fns(obsc).name;
            load(['../data/exp' expn '/' fn], 'data', 'subID', 'trialOrder');

            if expc == white_walls_filter_exp
                if strcmp(subID, 'll')
                    continue;
                end
            end

            redtrans_file = 'munsell_red_better.spd';
            greentrans_file = 'munsell_green_better.spd';
            bluetrans_file = 'munsell_blue_better.spd';
            yellowtrans_file = 'munsell_yellow_better.spd';
            monitor_spectra = csvread('../calibration/oled_mon_spectra.csv');
            [vor_filter_map, filter_idxs, redtrans, bluetrans, greentrans, yellowtrans] = ...
                gen_voronoi_filter_map(monitor_spectra, redtrans_file, bluetrans_file, greentrans_file, yellowtrans_file, false);

            for imgc = 1:size(imgns, 1)
                imgn = imgns(imgc).name;
                fparts = strsplit(imgn, '_');
                if expc == white_walls_filter_exp
                    if strcmp(fparts{end}, 'lighter.png')
                            hilo = 1;
                    else
                        hilo = 0;
                    end
                    illum = fparts{5};
                    body = fparts{9};
                else
                    error('analyzing an experiment that does not have any info in the image file names?!');
                end

                for maskc = 1:num_masks
                    img_idx = 1;
                    for img_idx = 1:length(imgns)
                        if strcmp(imgn, stats_for_curr_exp(1, img_idx).imgn.name)
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

                                qualityc(global_counter, 1) = NaN;

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

                            if expc == white_walls_filter_exp
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

                            gwld_match(global_counter, :) = gwdkl(1);
                            gwrg_match(global_counter, :) = gwdkl(2);
                            gwyv_match(global_counter, :) = gwdkl(3);

                            gwL_match(global_counter, :) = gwlab(1);
                            gwa_match(global_counter, :) = gwlab(2);
                            gwb_match(global_counter, :) = gwlab(3);

                            sdL_match(global_counter, :) = sdlab(1);
                            sda_match(global_counter, :) = sdlab(2);
                            sdb_match(global_counter, :) = sdlab(3);

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
                        end

                        exp_name{global_counter, 1} = exp_names{expc};
                        obs_name{global_counter, 1} = subID;
                        img_name{global_counter, 1} =  imgn;
                        illum_glaven{global_counter, 1} = illum;
                            body_glaven{global_counter, 1} = body;
                            rg_glaven(global_counter, 1) = -1; % indicates that this is not relevant for non-rmc_rsd_comp exps
                            by_glaven(global_counter, 1) = -1;
                            bkgd_glaven{global_counter, 1} = 'not_a_rmc_rsd_comp_image';
                        hilo_glaven(global_counter, 1) = hilo;
                        mask_name{global_counter, 1} = stats_for_curr_exp(1, img_idx).mask_name;

                        lrc_glaven(global_counter, 1) = stats_for_curr_exp(1, img_idx).lrc;

                        wpld_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).wpdkl(1);
                        wprg_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).wpdkl(2);
                        wpyv_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).wpdkl(3);

                        wpL_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).wplab(1);
                        wpa_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).wplab(2);
                        wpb_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).wplab(3);

                        wp_top5L_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).wplab_top5(1);
                        wp_top5a_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).wplab_top5(2);
                        wp_top5b_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).wplab_top5(3);

                        gwld_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).gwdkl(1);
                        gwrg_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).gwdkl(2);
                        gwyv_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).gwdkl(3);

                        gwL_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).gwlab(1);
                        gwa_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).gwlab(2);
                        gwb_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).gwlab(3);

                        sdL_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).sdlab(1);
                        sda_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).sdlab(2);
                        sdb_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).sdlab(3);

                        msat_ld_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).msatdkl(1);
                        msat_rg_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).msatdkl(2);
                        msat_yv_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).msatdkl(3);

                        msat_L_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).msatlab(1);
                        msat_a_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).msatlab(2);
                        msat_b_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).msatlab(3);

                        mfreq_L_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).mfreq_lab(1);
                        mfreq_a_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).mfreq_lab(2);
                        mfreq_b_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).mfreq_lab(3);

                        hv_glaven(global_counter, 1) = stats_for_curr_exp(1, img_idx).hv;

                        lm_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).lmsm(1);
                        mm_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).lmsm(2);
                        sm_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).lmsm(3);

                        lsd_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).lmssd(1);
                        msd_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).lmssd(2);
                        ssd_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).lmssd(3);

                        rmc_l_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).rmc(1);
                        rmc_m_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).rmc(2);
                        rmc_s_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).rmc(3);

                        rsd_l_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).rsd(1);
                        rsd_m_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).rsd(2);
                        rsd_s_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).rsd(3);

                        tau_simpler_l_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).tau_simpler(1);
                        tau_simpler_m_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).tau_simpler(2);
                        tau_simpler_s_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).tau_simpler(3);

                        tau_general_l_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).tau_general(1);
                        tau_general_m_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).tau_general(2);
                        tau_general_s_glaven(global_counter, :) = stats_for_curr_exp(1, img_idx).tau_general(3);

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

save('../data/white_walls_obsstats.mat', 'data_table');