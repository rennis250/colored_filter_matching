clear all; clc;

pcolors = {'k', 'b', 'r', 'g'};
psymbols = {'o', 'v', 'o', 'v'};

%% organize the data!

mean3 = @(x)(mean(x, 1));

load('../data/obsstats.mat', 'data_table');

obs_data = grpstats(data_table, {'exp_name', 'obs_name', 'img_name', 'mask_name', 'illum_glaven', 'rg_glaven', 'by_glaven', 'bkgd_glaven', 'body_glaven', 'hilo_glaven'}, {'mean', 'std'});
obs_data = removevars(obs_data, {'obs_name'});
pop_data = grpstats(obs_data, {'exp_name', 'img_name', 'mask_name', 'rg_glaven', 'by_glaven', 'bkgd_glaven', 'illum_glaven', 'body_glaven', 'hilo_glaven'}, {'mean', 'std'});

%% scatter plots - LAB
figure(1);
clf;
exp_names = {categorical({'patch_app'}), categorical({'patch_dye'}), categorical({'filter'}), categorical({'filter'})};
subplot_offset = 0;
for expc = 1:length(exp_names)
    if expc == 3
        subplot_offset = 3;
    end
    expn = exp_names{expc};
    idxs = pop_data.exp_name == expn & pop_data.mask_name == categorical({'obj_wout_high_mask'});
    curr_exp_data = pop_data(idxs, :);
    
    xs = [curr_exp_data.mean_mean_gwL_glaven, curr_exp_data.mean_mean_gwa_glaven, curr_exp_data.mean_mean_gwb_glaven];
    if expc ~= 4
        ys = [curr_exp_data.mean_mean_gwL_match, curr_exp_data.mean_mean_gwa_match, curr_exp_data.mean_mean_gwb_match];
    else
        ys = [curr_exp_data.mean_mean_wpL_filter, curr_exp_data.mean_mean_wpa_filter, curr_exp_data.mean_mean_wpb_filter];
    end
    plotLABStat(xs, ys, 2, 3, subplot_offset, pcolors, psymbols, expc, 'Image Stat - Matching Stimulus', 'Image Stat - Glaven', false);
end

%% scatter plots - LAB from mean RGB
figure(2);
clf;
exp_names = {categorical({'patch_app'}), categorical({'patch_dye'}), categorical({'filter'}), categorical({'filter'})};
subplot_offset = 0;
for expc = 1:length(exp_names)
    if expc == 3
        subplot_offset = 3;
    end
    expn = exp_names{expc};
    idxs = pop_data.exp_name == expn & pop_data.mask_name == categorical({'obj_wout_high_mask'});
    curr_exp_data = pop_data(idxs, :);
    
    xs = [curr_exp_data.mean_mean_gwL_from_mean_rgb_glaven, curr_exp_data.mean_mean_gwa_from_mean_rgb_glaven, curr_exp_data.mean_mean_gwb_from_mean_rgb_glaven];
    if expc ~= 4
        ys = [curr_exp_data.mean_mean_gwL_from_mean_rgb_match, curr_exp_data.mean_mean_gwa_from_mean_rgb_match, curr_exp_data.mean_mean_gwb_from_mean_rgb_match];
    else
        ys = [curr_exp_data.mean_mean_wpL_filter, curr_exp_data.mean_mean_wpa_filter, curr_exp_data.mean_mean_wpb_filter];
    end
    plotLABStat(xs, ys, 2, 3, subplot_offset, pcolors, psymbols, expc, 'Image Stat - Matching Stimulus', 'Image Stat - Glaven', false);
end

%% scatter plots - DKL
figure(3);
clf;
exp_names = {categorical({'patch_app'}), categorical({'patch_dye'}), categorical({'filter'}), categorical({'filter'})};
subplot_offset = 0;
for expc = 1:length(exp_names)
    if expc == 3
        subplot_offset = 3;
    end
    expn = exp_names{expc};
    idxs = pop_data.exp_name == expn & pop_data.mask_name == categorical({'obj_wout_high_mask'});
    curr_exp_data = pop_data(idxs, :);
    
    xs = [curr_exp_data.mean_mean_gwld_glaven, curr_exp_data.mean_mean_gwrg_glaven, curr_exp_data.mean_mean_gwyv_glaven];
    ys = [curr_exp_data.mean_mean_gwld_match, curr_exp_data.mean_mean_gwrg_match, curr_exp_data.mean_mean_gwyv_match];
    plotDKLStat(xs, ys, 2, 3, subplot_offset, pcolors, psymbols, expc, 'Image Stat - Matching Stimulus', 'Image Stat - Glaven', false);
end

%% scatter plots - LMS
figure(4);
clf;
exp_names = {categorical({'patch_app'}), categorical({'patch_dye'}), categorical({'filter'}), categorical({'filter'})};
subplot_offset = 0;
for expc = 1:length(exp_names)
    if expc == 3
        subplot_offset = 3;
    end
    expn = exp_names{expc};
    idxs = pop_data.exp_name == expn & pop_data.mask_name == categorical({'obj_wout_high_mask'});
    curr_exp_data = pop_data(idxs, :);
    
    xs = [curr_exp_data.mean_mean_lm_glaven, curr_exp_data.mean_mean_mm_glaven, curr_exp_data.mean_mean_sm_glaven];
    ys = [curr_exp_data.mean_mean_lm_match, curr_exp_data.mean_mean_mm_match, curr_exp_data.mean_mean_sm_match];
    plotLMSStat(xs, ys, 2, 3, subplot_offset, pcolors, psymbols, expc, 'Image Stat - Matching Stimulus', 'Image Stat - Glaven', false);
end

%% 2-d chromaticity plots
figure(5);
clf;
exp_names = {categorical({'patch_app'}), categorical({'patch_dye'}), categorical({'filter'}), categorical({'filter'})};
for expc = 1:length(exp_names)
    expn = exp_names{expc};
    idxs = pop_data.exp_name == expn & pop_data.mask_name == categorical({'obj_wout_high_mask'});
    curr_exp_data = pop_data(idxs, :);
    
    xg = curr_exp_data.mean_mean_gwa_glaven;
    xm = curr_exp_data.mean_mean_gwa_match;
    yg = curr_exp_data.mean_mean_gwb_glaven;
    ym = curr_exp_data.mean_mean_gwb_match;
    
    subplot(2, 2, expc);
    hold on;
    plot(xm, ym, 'ko');
    plot(xg, yg, 'rs');
    quiver(xg, yg, xm - xg, ym - yg, 0);
    
    xlim([-100, 100]);
    ylim([-100, 100]);
    axis square;
end

%% white point for uniform patch and filter
figure(6)
clf;
subplot_offset = 0;
exp_names = {categorical({'patch_app'}), categorical({'patch_dye'}), categorical({'filter'})};
for expc = 1:length(exp_names)
    if expc == 3
        subplot_offset = 3;
    end
    expn = exp_names{expc};
    idxs = pop_data.exp_name == expn & pop_data.mask_name == categorical({'obj_wout_high_mask'});
    curr_exp_data = pop_data(idxs, :);
    
    xs = [curr_exp_data.mean_mean_wpL_glaven, curr_exp_data.mean_mean_wpa_glaven, curr_exp_data.mean_mean_wpb_glaven];
    if expc ~= 3
        ys = [curr_exp_data.mean_mean_gwL_match, curr_exp_data.mean_mean_gwa_match, curr_exp_data.mean_mean_gwb_match];
    else
        ys = [curr_exp_data.mean_mean_wpL_filter, curr_exp_data.mean_mean_wpa_filter, curr_exp_data.mean_mean_wpb_filter];
    end
    plotLABStat(xs, ys, 2, 3, subplot_offset, pcolors, psymbols, expc, 'Image Stat - Matching Stimulus', 'Image Stat - Glaven', false);
end

%% RMC/RSD for uniform patch
figure(7)
clf;
subplot_offset = 0;
exp_names = {categorical({'patch_app'}), categorical({'patch_dye'}), categorical({'filter'})};
for expc = 1:length(exp_names)
    if expc == 3
        subplot_offset = 3;
    end
    expn = exp_names{expc};
    idxs = pop_data.exp_name == expn & pop_data.mask_name == categorical({'obj_wout_high_mask'});
    curr_exp_data = pop_data(idxs, :);
    
    subplot(2, 3, subplot_offset+1);
    hold on;
    plot(curr_exp_data.mean_mean_rmc_l_glaven, curr_exp_data.mean_mean_rmc_l_match, [pcolors{expc} psymbols{expc}]);
    line([0, 1], [0, 1]);
    xlim([0, 1]);
    ylim([0, 1]);
    if expc == 1
        title('RMC for uniform patches');
    elseif expc == 2
        legend('Appearance instructions', 'Dye match instructions', 'Location', 'SouthEast');
    elseif expc == 3
        title('RMC for flat filters');
    end
    xlabel('RMC for Glaven');
    ylabel('RMC for matching stimulus');
    axis square;
    
    subplot(2, 3, subplot_offset+2);
    hold on;
    plot(curr_exp_data.mean_mean_rmc_m_glaven, curr_exp_data.mean_mean_rmc_m_match, [pcolors{expc} psymbols{expc}]);
    line([0, 1], [0, 1]);
    xlim([0, 1]);
    ylim([0, 1]);
    xlabel('RMC for Glaven');
    axis square;
    
    subplot(2, 3, subplot_offset+3);
    hold on;
    plot(curr_exp_data.mean_mean_rmc_s_glaven, curr_exp_data.mean_mean_rmc_s_match, [pcolors{expc} psymbols{expc}]);
    line([0, 1], [0, 1]);
    xlim([0, 1]);
    ylim([0, 1]);
    xlabel('RMC for Glaven');
    axis square;
end

%% most saturated
figure(8)
clf;
subplot_offset = 0;
exp_names = {categorical({'patch_app'}), categorical({'patch_dye'}), categorical({'filter'})};
for expc = 1:length(exp_names)
    if expc == 3
        subplot_offset = 3;
    end
    expn = exp_names{expc};
    idxs = pop_data.exp_name == expn & pop_data.mask_name == categorical({'obj_wout_high_mask'});
    curr_exp_data = pop_data(idxs, :);
    
    xs = [curr_exp_data.mean_mean_msat_L_glaven, curr_exp_data.mean_mean_msat_a_glaven, curr_exp_data.mean_mean_msat_b_glaven];
    ys = [curr_exp_data.mean_mean_gwL_match, curr_exp_data.mean_mean_gwa_match, curr_exp_data.mean_mean_gwb_match];
    plotLABStat(xs, ys, 2, 3, subplot_offset, pcolors, psymbols, expc, 'Image Stat - Matching Stimulus', 'Image Stat - Glaven', false);
end

%% most frequent for uniform patch and flat filter
figure(9)
clf;
subplot_offset = 0;
exp_names = {categorical({'patch_app'}), categorical({'patch_dye'}), categorical({'filter'})};
for expc = 1:length(exp_names)
    if expc == 3
        subplot_offset = 3;
    end
    expn = exp_names{expc};
    idxs = pop_data.exp_name == expn & pop_data.mask_name == categorical({'obj_wout_high_mask'});
    curr_exp_data = pop_data(idxs, :);
    
    xs = [curr_exp_data.mean_mean_mfreq_L_glaven, curr_exp_data.mean_mean_mfreq_a_glaven, curr_exp_data.mean_mean_mfreq_b_glaven];
    if expc ~= 3
        ys = [curr_exp_data.mean_mean_gwL_match, curr_exp_data.mean_mean_gwa_match, curr_exp_data.mean_mean_gwb_match];
    else
        ys = [curr_exp_data.mean_mean_mfreq_L_filter, curr_exp_data.mean_mean_mfreq_a_filter, curr_exp_data.mean_mean_mfreq_b_filter];
    end
    plotLABStat(xs, ys, 2, 3, subplot_offset, pcolors, psymbols, expc, 'Image Stat - Matching Stimulus', 'Image Stat - Glaven', false);
end

%% scatter plots - RMC/RSD
figure(10);
clf;
subplot_offset = 0;
exp_names = {categorical({'filter'}), categorical({'filter'})};
for expc = 1:length(exp_names)
    if expc == 2
        subplot_offset = 3;
    end
    expn = exp_names{expc};
    idxs = pop_data.exp_name == expn & pop_data.mask_name == categorical({'obj_wout_high_mask'});
    curr_exp_data = pop_data(idxs, :);
    
    if expc == 1
        xs = [curr_exp_data.mean_mean_rmc_l_glaven, curr_exp_data.mean_mean_rmc_m_glaven, curr_exp_data.mean_mean_rmc_s_glaven];
        ys = [curr_exp_data.mean_mean_rmc_l_match, curr_exp_data.mean_mean_rmc_m_match, curr_exp_data.mean_mean_rmc_s_match];
    else
        xs = [curr_exp_data.mean_mean_rsd_l_glaven, curr_exp_data.mean_mean_rsd_m_glaven, curr_exp_data.mean_mean_rsd_s_glaven];
        ys = [curr_exp_data.mean_mean_rsd_l_filter, curr_exp_data.mean_mean_rsd_m_filter, curr_exp_data.mean_mean_rsd_s_filter];
    end
    plotLMSStat(xs, ys, 2, 3, subplot_offset, pcolors, psymbols, expc, 'Image Stat - Matching Stimulus', 'Image Stat - Glaven', false);
end

%% RMC/RSD - white walls - filter
figure(11)
clf;
subplot_offset = 0;
psymbols = {'o', 'v', 'd', 's'};
exp_names = {categorical({'white_walls_filter'}), categorical({'white_walls_filter'})};
for expc = 1:length(exp_names)
    if expc == 2
        subplot_offset = 3;
    end
    expn = exp_names{expc};
    idxs = pop_data.exp_name == expn & pop_data.mask_name == categorical({'obj_wout_high_mask'});
    curr_exp_data = pop_data(idxs, :);
    
    if expc == 1
        xs = [curr_exp_data.mean_mean_rmc_l_glaven, curr_exp_data.mean_mean_rmc_m_glaven, curr_exp_data.mean_mean_rmc_s_glaven];
        ys = [curr_exp_data.mean_mean_rmc_l_match, curr_exp_data.mean_mean_rmc_m_match, curr_exp_data.mean_mean_rmc_s_match];
    else
        xs = [curr_exp_data.mean_mean_rsd_l_glaven, curr_exp_data.mean_mean_rsd_m_glaven, curr_exp_data.mean_mean_rsd_s_glaven];
        ys = [curr_exp_data.mean_mean_rsd_l_filter, curr_exp_data.mean_mean_rsd_m_filter, curr_exp_data.mean_mean_rsd_s_filter];
    end
    plotLMSStat(xs, ys, 2, 3, subplot_offset, pcolors, psymbols, expc, 'Image Stat - Glaven', 'Image Stat - Matching Stimulus', true);
end

%% white point (top 5%) for uniform patch and filter
figure(12)
clf;
subplot_offset = 0;
exp_names = {categorical({'patch_app'}), categorical({'patch_dye'}), categorical({'filter'})};
for expc = 1:length(exp_names)
    if expc == 3
        subplot_offset = 3;
    end
    expn = exp_names{expc};
    idxs = pop_data.exp_name == expn & pop_data.mask_name == categorical({'obj_wout_high_mask'});
    curr_exp_data = pop_data(idxs, :);
    
    xs = [curr_exp_data.mean_mean_wp_top5L_glaven, curr_exp_data.mean_mean_wp_top5a_glaven, curr_exp_data.mean_mean_wp_top5b_glaven];
    if expc ~= 3
        ys = [curr_exp_data.mean_mean_gwL_match, curr_exp_data.mean_mean_gwa_match, curr_exp_data.mean_mean_gwb_match];
    else
        ys = [curr_exp_data.mean_mean_wp_top5L_filter, curr_exp_data.mean_mean_wp_top5a_filter, curr_exp_data.mean_mean_wp_top5b_filter];
    end
    plotLABStat(xs, ys, 2, 3, subplot_offset, pcolors, psymbols, expc, 'Image Stat - Matching Stimulus', 'Image Stat - Glaven', false);
end

%% scatter plots - robust RSD
figure(13);
clf;
subplot_offset = 0;
exp_names = {categorical({'filter'}), categorical({'filter'})};
for expc = 1:length(exp_names)
    if expc == 2
        subplot_offset = 3;
    end
    expn = exp_names{expc};
    idxs = pop_data.exp_name == expn & pop_data.mask_name == categorical({'obj_wout_high_mask'});
    curr_exp_data = pop_data(idxs, :);
    
    if expc == 1
        xs = [curr_exp_data.mean_mean_tau_simpler_l_glaven, curr_exp_data.mean_mean_tau_simpler_m_glaven, curr_exp_data.mean_mean_tau_simpler_s_glaven];
        ys = [curr_exp_data.mean_mean_tau_simpler_l_filter, curr_exp_data.mean_mean_tau_simpler_m_filter, curr_exp_data.mean_mean_tau_simpler_s_filter];
    else
        xs = [curr_exp_data.mean_mean_rsd_l_glaven, curr_exp_data.mean_mean_rsd_m_glaven, curr_exp_data.mean_mean_rsd_s_glaven];
        ys = [curr_exp_data.mean_mean_rsd_l_filter, curr_exp_data.mean_mean_rsd_m_filter, curr_exp_data.mean_mean_rsd_s_filter];
        % xs = [curr_exp_data.mean_mean_tau_general_l_glaven, curr_exp_data.mean_mean_tau_general_m_glaven, curr_exp_data.mean_mean_tau_general_s_glaven];
        % ys = [curr_exp_data.mean_mean_tau_general_l_filter, curr_exp_data.mean_mean_tau_general_m_filter, curr_exp_data.mean_mean_tau_general_s_filter];
    end
    plotLMSStat(xs, ys, 2, 3, subplot_offset, pcolors, psymbols, expc, 'Image Stat - Matching Stimulus', 'Image Stat - Glaven', false);
end

%% scatter plots - robust RSD
figure(14);
clf;
subplot_offset = 0;
exp_names = {categorical({'robust_rsd'}), categorical({'robust_rsd'})};
for expc = 1:length(exp_names)
    if expc == 2
        subplot_offset = 3;
    end
    expn = exp_names{expc};
    idxs = pop_data.exp_name == expn & pop_data.mask_name == categorical({'obj_wout_high_mask'});
    curr_exp_data = pop_data(idxs, :);
    
    if expc == 1
        xs = [curr_exp_data.mean_mean_tau_general_l_glaven, curr_exp_data.mean_mean_tau_general_m_glaven, curr_exp_data.mean_mean_tau_general_s_glaven];
        ys = [curr_exp_data.mean_mean_tau_general_l_filter, curr_exp_data.mean_mean_tau_general_m_filter, curr_exp_data.mean_mean_tau_general_s_filter];
    else
        xs = [curr_exp_data.mean_mean_rmc_l_glaven, curr_exp_data.mean_mean_rmc_m_glaven, curr_exp_data.mean_mean_rmc_s_glaven];
        ys = [curr_exp_data.mean_mean_rmc_l_match, curr_exp_data.mean_mean_rmc_m_match, curr_exp_data.mean_mean_rmc_s_match];
        % xs = [curr_exp_data.mean_mean_tau_general_l_glaven, curr_exp_data.mean_mean_tau_general_m_glaven, curr_exp_data.mean_mean_tau_general_s_glaven];
        % ys = [curr_exp_data.mean_mean_tau_general_l_filter, curr_exp_data.mean_mean_tau_general_m_filter, curr_exp_data.mean_mean_tau_general_s_filter];
    end
    plotLMSStat(xs, ys, 2, 3, subplot_offset, pcolors, psymbols, expc, 'Image Stat - Matching Stimulus', 'Image Stat - Glaven', false);
end

%% utility functions for plotting
function plotLABStat(xs, ys, nrows, ncols, subplot_offset, pcolors, psymbols, expc, xlab, ylab, with_rsq)
for c = 1:3
    subplot_er(nrows, ncols, subplot_offset+c);
    hold on;
    plot(xs(:, c), ys(:, c), [pcolors{expc} psymbols{expc}]);
    if c == 1
        line([0, 100], [0, 100], 'Color', 'k');
        xlim([0, 100]);
        ylim([0, 100]);
    else
        line([-100, 100], [-100, 100], 'Color', 'k');
        xlim([-100, 100]);
        ylim([-100, 100]);
    end
    xlabel(xlab);
    if c == 1
        ylabel(ylab);
    end
    if with_rsq
        title(['Adjusted R^2 = ', num2str(rsq(xs(:, c), ys(:, c), 1))]);
    end
    axis square;
end
end

function plotDKLStat(xs, ys, nrows, ncols, subplot_offset, pcolors, psymbols, expc, xlab, ylab, with_rsq)
for c = 1:3
    subplot_er(nrows, ncols, subplot_offset+c);
    hold on;
    plot(xs(:, c), ys(:, c), [pcolors{expc} psymbols{expc}]);
    line([-1, 1], [-1, 1], 'Color', 'k');
    xlim([-1, 1]);
    ylim([-1, 1]);
    xlabel(xlab);
    if c == 1
        ylabel(ylab);
    end
    if with_rsq
        title(['Adjusted R^2 = ', num2str(rsq(xs(:, c), ys(:, c), 1))]);
    end
    axis square;
end
end

function plotLMSStat(xs, ys, nrows, ncols, subplot_offset, pcolors, psymbols, expc, xlab, ylab, with_rsq)
coneNs = {'L', 'M', 'S'};
for c = 1:3
    subplot_er(nrows, ncols, subplot_offset+c);
    hold on;
    for cc = 1:size(xs, 1)
        plot(xs(cc, c), ys(cc, c), [pcolors{expc} psymbols{expc}]);
    end
    line([0, 2], [0, 2], 'Color', 'k');
    p = polyfit(xs(:, c), ys(:, c), 1);
    yfit = polyval(p, [0, 2]);
    plot([0, 2], yfit, 'k--');
    xlim([0, 2]);
    ylim([0, 2]);
    xlabel(xlab);
    if c == 1
        ylabel(ylab);
    end
    if with_rsq
        title([coneNs{c}, ' - Adjusted R^2 = ', num2str(rsq(xs(:, c), ys(:, c), 1))]);
    else
        title(coneNs{c});
    end
    axis square;
end
end