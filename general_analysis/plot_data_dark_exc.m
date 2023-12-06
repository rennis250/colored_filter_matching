clear all; clc;

%% organize the data!

mean3 = @(x)(mean(x, 1));

load('../data/dark_exc_obsstats.mat', 'data_table');

obs_data = grpstats(data_table, {'exp_name', 'obs_name', 'mask_name', 'illum_glaven', 'body_glaven', 'hilo_glaven'}, 'mean');
obs_data = removevars(obs_data, {'obs_name'});
pop_data = grpstats(obs_data, {'exp_name', 'mask_name', 'illum_glaven', 'body_glaven', 'hilo_glaven'}, 'mean');

%% scatter plots - LAB
figure(1);
clf;
pcolors = {'k', 'b', 'r', 'g'};
psymbols = {'o', 'v', 'o', 'v'};
exp_names = {categorical({'patch_app'}), categorical({'patch_dye'}), categorical({'filter'}), categorical({'filter'})};
subplot_offset = 0;
for expc = 1:length(exp_names)
    if expc == 3
        subplot_offset = 3;
    end
    expn = exp_names{expc};
    
    idxs = pop_data.exp_name == expn & pop_data.mask_name == categorical({'obj_mask'});
    curr_exp_data = pop_data(idxs, :);
    
    if expc ~= 4
        subplot(2, 3, subplot_offset+1);
        hold on;
        plot(curr_exp_data.mean_mean_gwL_glaven, curr_exp_data.mean_mean_gwL_match, [pcolors{expc} psymbols{expc}]);
        line([0, 100], [0, 100]);
        xlim([0, 100]);
        ylim([0, 100]);
        axis square;
        
        subplot(2, 3, subplot_offset+2);
        hold on;
        plot(curr_exp_data.mean_mean_gwa_glaven, curr_exp_data.mean_mean_gwa_match, [pcolors{expc} psymbols{expc}]);
        line([-100, 100], [-100, 100]);
        xlim([-100, 100]);
        ylim([-100, 100]);
        axis square;
        
        subplot(2, 3, subplot_offset+3);
        hold on;
        plot(curr_exp_data.mean_mean_gwb_glaven, curr_exp_data.mean_mean_gwb_match, [pcolors{expc} psymbols{expc}]);
        line([-100, 100], [-100, 100]);
        xlim([-100, 100]);
        ylim([-100, 100]);
        axis square;
    else
        subplot(2, 3, subplot_offset+1);
        hold on;
        plot(curr_exp_data.mean_mean_gwL_glaven, curr_exp_data.mean_mean_wpL_filter, [pcolors{expc} psymbols{expc}]);
        line([0, 100], [0, 100]);
        xlim([0, 100]);
        ylim([0, 100]);
        axis square;
        
        subplot(2, 3, subplot_offset+2);
        hold on;
        plot(curr_exp_data.mean_mean_gwa_glaven, curr_exp_data.mean_mean_wpa_filter, [pcolors{expc} psymbols{expc}]);
        line([-100, 100], [-100, 100]);
        xlim([-100, 100]);
        ylim([-100, 100]);
        axis square;
        
        subplot(2, 3, subplot_offset+3);
        hold on;
        plot(curr_exp_data.mean_mean_gwb_glaven, curr_exp_data.mean_mean_wpb_filter, [pcolors{expc} psymbols{expc}]);
        line([-100, 100], [-100, 100]);
        xlim([-100, 100]);
        ylim([-100, 100]);
        axis square;
    end
end

%% scatter plots - LAB from mean RGB
figure(2);
clf;
pcolors = {'k', 'b', 'r', 'g'};
psymbols = {'o', 'v', 'o', 'v'};
exp_names = {categorical({'patch_app'}), categorical({'patch_dye'}), categorical({'filter'}), categorical({'filter'})};
subplot_offset = 0;
for expc = 1:length(exp_names)
    if expc == 3
        subplot_offset = 3;
    end
    expn = exp_names{expc};
    
    idxs = pop_data.exp_name == expn & pop_data.mask_name == categorical({'obj_mask'});
    curr_exp_data = pop_data(idxs, :);
    
    if expc ~= 4
        subplot(2, 3, subplot_offset+1);
        hold on;
        plot(curr_exp_data.mean_mean_gwL_from_mean_rgb_glaven, curr_exp_data.mean_mean_gwL_from_mean_rgb_match, [pcolors{expc} psymbols{expc}]);
        line([0, 100], [0, 100]);
        xlim([0, 100]);
        ylim([0, 100]);
        axis square;
        
        subplot(2, 3, subplot_offset+2);
        hold on;
        plot(curr_exp_data.mean_mean_gwa_from_mean_rgb_glaven, curr_exp_data.mean_mean_gwa_from_mean_rgb_match, [pcolors{expc} psymbols{expc}]);
        line([-100, 100], [-100, 100]);
        xlim([-100, 100]);
        ylim([-100, 100]);
        axis square;
        
        subplot(2, 3, subplot_offset+3);
        hold on;
        plot(curr_exp_data.mean_mean_gwb_from_mean_rgb_glaven, curr_exp_data.mean_mean_gwb_from_mean_rgb_match, [pcolors{expc} psymbols{expc}]);
        line([-100, 100], [-100, 100]);
        xlim([-100, 100]);
        ylim([-100, 100]);
        axis square;
    else
        subplot(2, 3, subplot_offset+1);
        hold on;
        plot(curr_exp_data.mean_mean_gwL_from_mean_rgb_glaven, curr_exp_data.mean_mean_wpL_filter, [pcolors{expc} psymbols{expc}]);
        line([0, 100], [0, 100]);
        xlim([0, 100]);
        ylim([0, 100]);
        axis square;
        
        subplot(2, 3, subplot_offset+2);
        hold on;
        plot(curr_exp_data.mean_mean_gwa_from_mean_rgb_glaven, curr_exp_data.mean_mean_wpa_filter, [pcolors{expc} psymbols{expc}]);
        line([-100, 100], [-100, 100]);
        xlim([-100, 100]);
        ylim([-100, 100]);
        axis square;
        
        subplot(2, 3, subplot_offset+3);
        hold on;
        plot(curr_exp_data.mean_mean_gwb_from_mean_rgb_glaven, curr_exp_data.mean_mean_wpb_filter, [pcolors{expc} psymbols{expc}]);
        line([-100, 100], [-100, 100]);
        xlim([-100, 100]);
        ylim([-100, 100]);
        axis square;
    end
end

%% scatter plots - DKL
figure(3);
clf;
pcolors = {'k', 'b', 'r', 'g'};
psymbols = {'o', 'v', 'o', 'v'};
exp_names = {categorical({'patch_app'}), categorical({'patch_dye'}), categorical({'filter'}), categorical({'filter'})};
subplot_offset = 0;
for expc = 1:length(exp_names)
    if expc == 3
        subplot_offset = 3;
    end
    expn = exp_names{expc};
    
    idxs = pop_data.exp_name == expn & pop_data.mask_name == categorical({'obj_mask'});
    curr_exp_data = pop_data(idxs, :);
    
    subplot(2, 3, subplot_offset+1);
    hold on;
    plot(curr_exp_data.mean_mean_gwld_glaven, curr_exp_data.mean_mean_gwld_match, [pcolors{expc} psymbols{expc}]);
    line([-1, 1], [-1, 1]);
    xlim([-1, 1]);
    ylim([-1, 1]);
    axis square;
    
    subplot(2, 3, subplot_offset+2);
    hold on;
    plot(curr_exp_data.mean_mean_gwrg_glaven, curr_exp_data.mean_mean_gwrg_match, [pcolors{expc} psymbols{expc}]);
    line([-1, 1], [-1, 1]);
    xlim([-1, 1]);
    ylim([-1, 1]);
    axis square;
    
    subplot(2, 3, subplot_offset+3);
    hold on;
    plot(curr_exp_data.mean_mean_gwyv_glaven, curr_exp_data.mean_mean_gwyv_match, [pcolors{expc} psymbols{expc}]);
    line([-1, 1], [-1, 1]);
    xlim([-1, 1]);
    ylim([-1, 1]);
    axis square;
end

%% scatter plots - LMS
figure(4);
clf;
pcolors = {'k', 'b', 'r', 'g'};
psymbols = {'o', 'v', 'o', 'v'};
exp_names = {categorical({'patch_app'}), categorical({'patch_dye'}), categorical({'filter'}), categorical({'filter'})};
subplot_offset = 0;
for expc = 1:length(exp_names)
    if expc == 3
        subplot_offset = 3;
    end
    expn = exp_names{expc};
    
    idxs = pop_data.exp_name == expn & pop_data.mask_name == categorical({'obj_mask'});
    curr_exp_data = pop_data(idxs, :);
    
    subplot(2, 3, subplot_offset+1);
    hold on;
    plot(curr_exp_data.mean_mean_lm_glaven, curr_exp_data.mean_mean_lm_match, [pcolors{expc} psymbols{expc}]);
    line([0, 1], [0, 1]);
    xlim([0, 0.1]);
    ylim([0, 0.1]);
    axis square;
    
    subplot(2, 3, subplot_offset+2);
    hold on;
    plot(curr_exp_data.mean_mean_mm_glaven, curr_exp_data.mean_mean_mm_match, [pcolors{expc} psymbols{expc}]);
    line([0, 1], [0, 1]);
    xlim([0, 0.1]);
    ylim([0, 0.1]);
    axis square;
    
    subplot(2, 3, subplot_offset+3);
    hold on;
    plot(curr_exp_data.mean_mean_sm_glaven, curr_exp_data.mean_mean_sm_match, [pcolors{expc} psymbols{expc}]);
    line([0, 1], [0, 1]);
    xlim([0, 0.1]);
    ylim([0, 0.1]);
    axis square;
end

%% 2-d chromaticity plots
figure(5);
clf;
exp_names = {categorical({'patch_app'}), categorical({'patch_dye'}), categorical({'filter'}), categorical({'filter'})};
for expc = 1:length(exp_names)
    expn = exp_names{expc};
    idxs = pop_data.exp_name == expn & pop_data.mask_name == categorical({'obj_mask'});
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
pcolors = {'k', 'b', 'r', 'g'};
psymbols = {'o', 'v', 'o', 'v'};
exp_names = {categorical({'patch_app'}), categorical({'patch_dye'}), categorical({'filter'})};
for expc = 1:length(exp_names)
    if expc == 3
        subplot_offset = 3;
    end
    expn = exp_names{expc};
    idxs = pop_data.exp_name == expn & pop_data.mask_name == categorical({'obj_mask'});
    curr_exp_data = pop_data(idxs, :);
    
    if expc ~= 3
        subplot(2, 3, subplot_offset+1);
        hold on;
        plot(curr_exp_data.mean_mean_wpL_glaven, curr_exp_data.mean_mean_gwL_match, [pcolors{expc} psymbols{expc}]);
        line([0, 100], [0, 100]);
        xlim([0, 100]);
        ylim([0, 100]);
        axis square;
        
        subplot(2, 3, subplot_offset+2);
        hold on;
        plot(curr_exp_data.mean_mean_wpa_glaven, curr_exp_data.mean_mean_gwa_match, [pcolors{expc} psymbols{expc}]);
        line([-100, 100], [-100, 100]);
        xlim([-100, 100]);
        ylim([-100, 100]);
        axis square;
        
        subplot(2, 3, subplot_offset+3);
        hold on;
        plot(curr_exp_data.mean_mean_wpb_glaven, curr_exp_data.mean_mean_gwb_match, [pcolors{expc} psymbols{expc}]);
        line([-100, 100], [-100, 100]);
        xlim([-100, 100]);
        ylim([-100, 100]);
        axis square;
    else
        subplot(2, 3, subplot_offset+1);
        hold on;
        plot(curr_exp_data.mean_mean_wpL_glaven, curr_exp_data.mean_mean_wpL_filter, [pcolors{expc} psymbols{expc}]);
        line([0, 100], [0, 100]);
        xlim([0, 100]);
        ylim([0, 100]);
        axis square;
        
        subplot(2, 3, subplot_offset+2);
        hold on;
        plot(curr_exp_data.mean_mean_wpa_glaven, curr_exp_data.mean_mean_wpa_filter, [pcolors{expc} psymbols{expc}]);
        line([-100, 100], [-100, 100]);
        xlim([-100, 100]);
        ylim([-100, 100]);
        axis square;
        
        subplot(2, 3, subplot_offset+3);
        hold on;
        plot(curr_exp_data.mean_mean_wpb_glaven, curr_exp_data.mean_mean_wpb_filter, [pcolors{expc} psymbols{expc}]);
        line([-100, 100], [-100, 100]);
        xlim([-100, 100]);
        ylim([-100, 100]);
        axis square;
    end
end

%% RMC/RSD for uniform patch
figure(7)
clf;
subplot_offset = 0;
pcolors = {'k', 'b', 'r', 'g'};
psymbols = {'o', 'v', 'o', 'v'};
exp_names = {categorical({'patch_app'}), categorical({'patch_dye'}), categorical({'filter'})};
for expc = 1:length(exp_names)
    if expc == 3
        subplot_offset = 3;
    end
    expn = exp_names{expc};
    idxs = pop_data.exp_name == expn & pop_data.mask_name == categorical({'obj_mask'});
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

%% most saturated for uniform patch
figure(8)
clf;
subplot_offset = 0;
pcolors = {'k', 'b', 'r', 'g'};
psymbols = {'o', 'v', 'o', 'v'};
exp_names = {categorical({'patch_app'}), categorical({'patch_dye'}), categorical({'filter'})};
for expc = 1:length(exp_names)
    if expc == 3
        subplot_offset = 3;
    end
    expn = exp_names{expc};
    idxs = pop_data.exp_name == expn & pop_data.mask_name == categorical({'obj_mask'});
    curr_exp_data = pop_data(idxs, :);
    
    subplot(2, 3, subplot_offset+1);
    hold on;
    plot(curr_exp_data.mean_mean_msat_L_glaven, curr_exp_data.mean_mean_gwL_match, [pcolors{expc} psymbols{expc}]);
    line([0, 100], [0, 100]);
    xlim([0, 100]);
    ylim([0, 100]);
    axis square;
    
    subplot(2, 3, subplot_offset+2);
    hold on;
    plot(curr_exp_data.mean_mean_msat_a_glaven, curr_exp_data.mean_mean_gwa_match, [pcolors{expc} psymbols{expc}]);
    line([-100, 100], [-100, 100]);
    xlim([-100, 100]);
    ylim([-100, 100]);
    axis square;
    
    subplot(2, 3, subplot_offset+3);
    hold on;
    plot(curr_exp_data.mean_mean_msat_b_glaven, curr_exp_data.mean_mean_gwb_match, [pcolors{expc} psymbols{expc}]);
    line([-100, 100], [-100, 100]);
    xlim([-100, 100]);
    ylim([-100, 100]);
    axis square;
end

%% most frequent for uniform patch and flat filter
figure(9)
clf;
subplot_offset = 0;
pcolors = {'k', 'b', 'r', 'g'};
psymbols = {'o', 'v', 'o', 'v'};
exp_names = {categorical({'patch_app'}), categorical({'patch_dye'}), categorical({'filter'})};
for expc = 1:length(exp_names)
    if expc == 3
        subplot_offset = 3;
    end
    expn = exp_names{expc};
    idxs = pop_data.exp_name == expn & pop_data.mask_name == categorical({'obj_wout_high_mask'});
    curr_exp_data = pop_data(idxs, :);
    
    if expc ~= 3
        subplot(2, 3, subplot_offset+1);
        hold on;
        plot(curr_exp_data.mean_mean_mfreq_L_glaven, curr_exp_data.mean_mean_gwL_match, [pcolors{expc} psymbols{expc}]);
        line([0, 100], [0, 100]);
        xlim([0, 100]);
        ylim([0, 100]);
        axis square;
        
        subplot(2, 3, subplot_offset+2);
        hold on;
        plot(curr_exp_data.mean_mean_mfreq_a_glaven, curr_exp_data.mean_mean_gwa_match, [pcolors{expc} psymbols{expc}]);
        line([-100, 100], [-100, 100]);
        xlim([-100, 100]);
        ylim([-100, 100]);
        axis square;
        
        subplot(2, 3, subplot_offset+3);
        hold on;
        plot(curr_exp_data.mean_mean_mfreq_b_glaven, curr_exp_data.mean_mean_gwb_match, [pcolors{expc} psymbols{expc}]);
        line([-100, 100], [-100, 100]);
        xlim([-100, 100]);
        ylim([-100, 100]);
        axis square;
    else
        subplot(2, 3, subplot_offset+1);
        hold on;
        plot(curr_exp_data.mean_mean_mfreq_L_glaven, curr_exp_data.mean_mean_mfreq_L_filter, [pcolors{expc} psymbols{expc}]);
        line([0, 100], [0, 100]);
        xlim([0, 100]);
        ylim([0, 100]);
        axis square;
        
        subplot(2, 3, subplot_offset+2);
        hold on;
        plot(curr_exp_data.mean_mean_mfreq_a_glaven, curr_exp_data.mean_mean_mfreq_a_filter, [pcolors{expc} psymbols{expc}]);
        line([-100, 100], [-100, 100]);
        xlim([-100, 100]);
        ylim([-100, 100]);
        axis square;
        
        subplot(2, 3, subplot_offset+3);
        hold on;
        plot(curr_exp_data.mean_mean_mfreq_b_glaven, curr_exp_data.mean_mean_mfreq_b_filter, [pcolors{expc} psymbols{expc}]);
        line([-100, 100], [-100, 100]);
        xlim([-100, 100]);
        ylim([-100, 100]);
        axis square;
    end
end

%% scatter plots - RMC/RSD
figure(10);
clf;
pcolors = {'k', 'b', 'r', 'g'};
psymbols = {'o', 'v', 'o', 'v'};
exp_names = {categorical({'filter'}), categorical({'filter'})};
for expc = 1:length(exp_names)
    expn = exp_names{expc};
    
    idxs = pop_data.exp_name == expn & pop_data.mask_name == categorical({'obj_mask'});
    curr_exp_data = pop_data(idxs, :);
    
    if expc == 1
        subplot(1, 3, 1);
        hold on;
        plot(curr_exp_data.mean_mean_rmc_l_glaven, curr_exp_data.mean_mean_rmc_l_match, [pcolors{expc} psymbols{expc}]);
        line([0, 1.5], [0, 1.5]);
        xlim([0, 1.5]);
        ylim([0, 1.5]);
        axis square;
        
        subplot(1, 3, 2);
        hold on;
        plot(curr_exp_data.mean_mean_rmc_m_glaven, curr_exp_data.mean_mean_rmc_m_match, [pcolors{expc} psymbols{expc}]);
        line([0, 1.5], [0, 1.5]);
        xlim([0, 1.5]);
        ylim([0, 1.5]);
        axis square;
        
        subplot(1, 3, 3);
        hold on;
        plot(curr_exp_data.mean_mean_rmc_s_glaven, curr_exp_data.mean_mean_rmc_s_match, [pcolors{expc} psymbols{expc}]);
        line([0, 1.5], [0, 1.5]);
        xlim([0, 1.5]);
        ylim([0, 1.5]);
        axis square;
    else
        subplot(1, 3, 1);
        hold on;
        plot(curr_exp_data.mean_mean_rsd_l_glaven, curr_exp_data.mean_mean_rsd_l_filter, [pcolors{expc} psymbols{expc}]);
        line([0, 1.5], [0, 1.5]);
        xlim([0, 1.5]);
        ylim([0, 1.5]);
        axis square;
        
        subplot(1, 3, 2);
        hold on;
        plot(curr_exp_data.mean_mean_rsd_m_glaven, curr_exp_data.mean_mean_rsd_m_filter, [pcolors{expc} psymbols{expc}]);
        line([0, 1.5], [0, 1.5]);
        xlim([0, 1.5]);
        ylim([0, 1.5]);
        axis square;
        
        subplot(1, 3, 3);
        hold on;
        plot(curr_exp_data.mean_mean_rsd_s_glaven, curr_exp_data.mean_mean_rsd_s_filter, [pcolors{expc} psymbols{expc}]);
        line([0, 1.5], [0, 1.5]);
        xlim([0, 1.5]);
        ylim([0, 1.5]);
        axis square;
    end
end

%% RMC/RSD for uniform patch - white walls
figure(11)
clf;
subplot_offset = 0;
pcolors = {'k', 'b', 'r', 'g'};
psymbols = {'o', 'v', 'o', 'v'};
exp_names = {categorical({'white_walls_patch'}), categorical({'white_walls_filter'})};
for expc = 1:length(exp_names)
    expn = exp_names{expc};
    idxs = pop_data.exp_name == expn & pop_data.mask_name == categorical({'obj_mask'});
    curr_exp_data = pop_data(idxs, :);
    
    subplot(1, 3, 1);
    hold on;
    plot(curr_exp_data.mean_mean_rmc_l_glaven, curr_exp_data.mean_mean_rmc_l_match, [pcolors{expc} psymbols{expc}]);
    line([0, 1], [0, 1]);
    xlim([0, 1]);
    ylim([0, 1]);
    xlabel('RMC for Glaven');
    ylabel('RMC for matching stimulus');
    axis square;
    
    subplot(1, 3, 2);
    hold on;
    plot(curr_exp_data.mean_mean_rmc_m_glaven, curr_exp_data.mean_mean_rmc_m_match, [pcolors{expc} psymbols{expc}]);
    line([0, 1], [0, 1]);
    xlim([0, 1]);
    ylim([0, 1]);
    xlabel('RMC for Glaven');
    axis square;
    
    subplot(1, 3, 3);
    hold on;
    plot(curr_exp_data.mean_mean_rmc_s_glaven, curr_exp_data.mean_mean_rmc_s_match, [pcolors{expc} psymbols{expc}]);
    line([0, 1], [0, 1]);
    xlim([0, 1]);
    ylim([0, 1]);
    xlabel('RMC for Glaven');
    axis square;
end

