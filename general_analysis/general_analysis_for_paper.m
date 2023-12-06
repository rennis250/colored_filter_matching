clear all; close all; clc;

load('../data/obsstats.mat');

% obs_data = grpstats(data_table, {'exp_name', 'obs_name', 'mask_name', 'illum_glaven', 'body_glaven', 'hilo_glaven'}, 'mean');
% obs_data = removevars(obs_data, {'obs_name'});
% pop_data = grpstats(obs_data, {'exp_name', 'mask_name', 'illum_glaven', 'body_glaven', 'hilo_glaven'}, 'mean');

DKL2RGB = dkl2rgbFromCalib('../calibration/eizo_chroma.csv');
RGB2DKL_T = inv(DKL2RGB);

pop_mean_data = {};
obs_mean_data = {};
obs_covms_mean = {};
obs_pas_mean = {};
obs_covms_max = {};
obs_pas_max = {};

%% collect the data together
exps = unique(data_table.exp_name);
for cc = 1:length(exps)
    exp = exps(cc);
    d = data_table(data_table.exp_name == exp, :);

    obs_data = {};
    sub_ids = unique(d.obs_name);
    illums = unique(d.illum_glaven);
    bodies = unique(d.body_glaven);
    highlows = unique(d.hilo_glaven);
    for sc = 1:length(sub_ids)
        sub_id = sub_ids(sc);
        for ic = 1:length(illums)
            illum = illums(ic);
            for bc = 1:length(bodies)
                body = bodies(bc);
                for hc = 1:length(highlows)
                    highlow = highlows(hc);
                    sub_data = d(d.mask_name == 'obj_wout_high_mask' & d.obs_name == sub_id & d.illum_glaven == illum & d.body_glaven == body & d.hilo_glaven == highlow, :);
                    if isempty(sub_data)
                        continue;
                    end

                    % from image
                    t = {exp, mean(sub_data.gwL_glaven), std(sub_data.gwL_glaven), ...
                              mean(sub_data.gwa_glaven), std(sub_data.gwa_glaven), ...
                              mean(sub_data.gwb_glaven), std(sub_data.gwb_glaven)};
                    t = [t, {illum, body, highlow}];

                    % from obs
                    t = [t, {mean(sub_data.gwL_match), std(sub_data.gwL_match), ...
                             mean(sub_data.gwa_match), std(sub_data.gwa_match), ...
                             mean(sub_data.gwb_match), std(sub_data.gwb_match)}];
                         
                    if exp == 'filter' || exp == 'rmc_rsd_comp' || exp == 'white_walls_filter'
                        t = [t, {mean(sub_data.wp_top5L_filter), std(sub_data.wp_top5L_filter), ...
                                 mean(sub_data.wp_top5a_filter), std(sub_data.wp_top5a_filter), ...
                                 mean(sub_data.wp_top5b_filter), std(sub_data.wp_top5b_filter)}];
                    else
                        t = [t, {NaN, NaN, NaN, NaN, NaN, NaN}];
                    end

                    obs_data = [obs_data; t];
                end
            end
        end
    end
    obs_data = cell2table(obs_data, 'VariableNames', {'exp_name', 'mean_gwL_glaven', 'sd_gwL_glaven', ...
                                                                  'mean_gwa_glaven', 'sd_gwa_glaven', ...
                                                                  'mean_gwb_glaven', 'sd_gwb_glaven', 'illum', 'body', 'hilo', ...
                                                                  ...
                                                                  'mean_gwL_match', 'sd_gwL_match', ...
                                                                  'mean_gwa_match', 'sd_gwa_match', ...
                                                                  'mean_gwb_match', 'sd_gwb_match', ...
                                                                  ...
                                                                  'mean_wp_top5L_filter', 'sd_wp_top5L_filter', ...
                                                                  'mean_wp_top5a_filter', 'sd_wp_top5a_filter', ...
                                                                  'mean_wp_top5b_filter', 'sd_wp_top5b_filter'});

    ct = 1;
    covms_mean = {};
    pas_mean = {};
    covms_max = {};
    pas_max = {};
    tmp_data = {};
    illums = unique(obs_data.illum);
    bodies = unique(obs_data.body);
    highlows = unique(obs_data.hilo);
    for ic = 1:length(illums)
        illum = illums(ic);
        for bc = 1:length(bodies)
            body = bodies(bc);
            for hc = 1:length(highlows)
                highlow = highlows(hc);
                sub_data = obs_data(obs_data.illum == illum & obs_data.body == body & obs_data.hilo == highlow, :);
                if isempty(sub_data)
                    continue;
                end

                % from image
                t = {exp, mean(sub_data.mean_gwL_glaven), std(sub_data.mean_gwL_glaven), ...
                          mean(sub_data.mean_gwa_glaven), std(sub_data.mean_gwa_glaven), ...
                          mean(sub_data.mean_gwb_glaven), std(sub_data.mean_gwb_glaven)};
                t = [t, {illum, body, highlow}];

                % from obs
                t = [t, {mean(sub_data.mean_gwL_match), std(sub_data.mean_gwL_match), ...
                         mean(sub_data.mean_gwa_match), std(sub_data.mean_gwa_match), ...
                         mean(sub_data.mean_gwb_match), std(sub_data.mean_gwb_match)}];
                     
                covms_mean{ct} = cov([sub_data.mean_gwa_match(:), sub_data.mean_gwb_match(:)]);
                pas_mean{ct} = pca([sub_data.mean_gwa_match(:), sub_data.mean_gwb_match(:)]);
                if exp == 'filter' || exp == 'rmc_rsd_comp' || exp == 'white_walls_filter'
                    t = [t, {mean(sub_data.mean_wp_top5L_filter), std(sub_data.mean_wp_top5L_filter), ...
                             mean(sub_data.mean_wp_top5a_filter), std(sub_data.mean_wp_top5a_filter), ...
                             mean(sub_data.mean_wp_top5b_filter), std(sub_data.mean_wp_top5b_filter)}];
                         
                    covms_max{ct} = cov([sub_data.mean_wp_top5a_filter(:), sub_data.mean_wp_top5b_filter(:)]);
                    pas_max{ct} = pca([sub_data.mean_wp_top5a_filter(:), sub_data.mean_wp_top5b_filter(:)]);
                else
                    t = [t, {NaN, NaN, NaN, NaN, NaN, NaN}];
                end

                tmp_data = [tmp_data; t];

                ct = ct + 1;
            end
        end
    end
    tmp_data = cell2table(tmp_data, 'VariableNames', {'exp_name', 'mean_mean_gwL_glaven', 'se_gwL_glaven', ...
                                                                  'mean_mean_gwa_glaven', 'se_gwa_glaven', ...
                                                                  'mean_mean_gwb_glaven', 'se_gwb_glaven', 'illum', 'body', 'hilo', ...
                                                                  ...
                                                                  'mean_mean_gwL_match', 'se_gwL_match', ...
                                                                  'mean_mean_gwa_match', 'se_gwa_match', ...
                                                                  'mean_mean_gwb_match', 'se_gwb_match', ...
                                                                  ...
                                                                  'mean_mean_wp_top5L_filter', 'se_wp_top5L_filter', ...
                                                                  'mean_mean_wp_top5a_filter', 'se_wp_top5a_filter', ...
                                                                  'mean_mean_wp_top5b_filter', 'se_wp_top5b_filter'});

    color = zeros(size(tmp_data, 1), 3);
    colors = {'red', 'green', 'blue', 'yellow'};
    hilos = {categorical(0), categorical(1)};
    dkls = {dkl2rgb(RGB2DKL_T, [0, 1, 0]')', ...
        dkl2rgb(RGB2DKL_T, [0, -1, 0]')', ...
        dkl2rgb(RGB2DKL_T, [0, 0, 1]')', ...
        dkl2rgb(RGB2DKL_T, [0, 0, -1]')', ...

        dkl2rgb(RGB2DKL_T, [-0.6, 1, 0]')', ...
        dkl2rgb(RGB2DKL_T, [-0.6, -1, 0]')', ...
        dkl2rgb(RGB2DKL_T, [-0.6, 0, 1]')', ...
        dkl2rgb(RGB2DKL_T, [-0.6, 0, -1]')'};
    ct = 1;
    for c = colors
        for hl = hilos
            idxs = find(tmp_data.body == c{1} & tmp_data.hilo == hl{1});
            color(idxs, :) = repmat(dkls{ct}, length(idxs), 1);
            ct = ct + 1;
        end
    end
    tmp_data.color = color;

    pch = cell(size(tmp_data, 1), 1);
    if cc == 4 || cc == 5 || cc == 6
        for x = 1:size(pch, 1)
            pch{x} = 's';
        end
    else
        blueidxs = find(tmp_data.illum == 'blue');
        yellowidxs = find(tmp_data.illum == 'yellow');
        for bc = 1:length(blueidxs)
            pch{blueidxs(bc)} = 'o';
        end
        for yc = 1:length(yellowidxs)
            pch{yellowidxs(yc)} = 'v';
        end
    end
    tmp_data.pch = pch;

    pop_mean_data{cc} = tmp_data;
    obs_mean_data{cc} = obs_data;
    obs_covms_mean{cc} = covms_mean;
    obs_pas_mean{cc} = pas_mean;
    obs_covms_max{cc} = covms_max;
    obs_pas_max{cc} = pas_max;
end

%% make a plot showing the color matches in a chromaticity diagram
figure(1)
clf;
exps = {'Patch proximal match', 'Patch dye match', 'Flat filter - mean color', 'Flat filter - White Point'};
for i = 1:4
    subplot(2, 2, i);
    hold on;
    axis square;
    xlim([-100, 100]);
    ylim([-100, 100]);
    xlabel('a*');
    ylabel('b*');
    title(exps{i});

    if i == 1
        pdata = pop_mean_data{2};
        odata = obs_mean_data{2};
        the_covm = obs_covms_mean{2};
    elseif i == 2
        pdata = pop_mean_data{3};
        odata = obs_mean_data{3};
        the_covm = obs_covms_mean{3};
    elseif i == 3
        pdata = pop_mean_data{1};
        odata = obs_mean_data{1};
        the_covm = obs_covms_mean{1};
    elseif i == 4
        pdata = pop_mean_data{1};
        odata = obs_mean_data{1};
        the_covm = obs_covms_max{1};
    end

    if i < 3
        img_a = pdata.mean_mean_gwa_glaven;
        img_b = pdata.mean_mean_gwb_glaven;

        pop_a = pdata.mean_mean_gwa_match;
        pop_b = pdata.mean_mean_gwb_match;
    else
        img_a = pdata.mean_mean_gwa_glaven;
        img_b = pdata.mean_mean_gwb_glaven;

        pop_a = pdata.mean_mean_wp_top5a_filter;
        pop_b = pdata.mean_mean_wp_top5b_filter;
    end

    for x = 1:length(the_covm)
        error_ellipse(the_covm{x}, [pop_a(x), pop_b(x)], 0.3934, 'style', 'k-');
        plot(img_a(x), img_b(x), 'ko', 'MarkerEdgeColor', lab2rgb([50, img_a(x)*0.7, img_b(x)*0.7]), 'MarkerFaceColor', [1, 1, 1]);
        plot(pop_a(x), pop_b(x), 'ko', 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', lab2rgb([50, img_a(x)*0.7, img_b(x)*0.7]));
        quiver(img_a(x), img_b(x), pop_a(x) - img_a(x), pop_b(x) - img_b(x), 0, 'Color', [0, 0, 0]);
    end

    % axis(1, at = seq(-50, 50, 50));
    % axis(1, at = 100, labels = "100", col.axis = 4);
    % axis(1, at = -100, labels = "-100", col.axis = 5);

    % axis(2, at = seq(-50, 50, 50));
    % axis(2, at = 100, labels = "100", col.axis = 7);
    % axis(2, at = -100, labels = "-100", col.axis = 6);
end

%% make a separate chromaticity diagram plot for the white walls stimuli
figure(2)
clf;
hold on;
axis square;
xlim([-100, 100]);
ylim([-100, 100]);
xlabel('a*');
ylabel('b*');
title('White walls');

pdata = pop_mean_data{5};
odata = obs_mean_data{5};
the_covm = obs_covms_max{5};

img_a = pdata.mean_mean_gwa_glaven;
img_b = pdata.mean_mean_gwb_glaven;

pop_a = pdata.mean_mean_gwa_match;
pop_b = pdata.mean_mean_gwb_match;

for x = 1:length(the_covm)
    error_ellipse(the_covm{x}, [pop_a(x), pop_b(x)], 0.3934, 'style', 'k-');
    plot(img_a(x), img_b(x), 'ko', 'MarkerEdgeColor', clamp(lab2rgb([50, img_a(x)*0.7, img_b(x)*0.7])), 'MarkerFaceColor', [1, 1, 1]);
    plot(pop_a(x), pop_b(x), 'ko', 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', clamp(lab2rgb([50, img_a(x)*0.7, img_b(x)*0.7])));
    quiver(img_a(x), img_b(x), pop_a(x) - img_a(x), pop_b(x) - img_b(x), 0, 'Color', [0, 0, 0]);
end

% axis(1, at = seq(-50, 50, 50));
% axis(1, at = 100, labels = "100", col.axis = 4);
% axis(1, at = -100, labels = "-100", col.axis = 5);

% axis(2, at = seq(-50, 50, 50));
% axis(2, at = 100, labels = "100", col.axis = 7);
% axis(2, at = -100, labels = "-100", col.axis = 6);

%% plotting results for whole obj in terms of mean LAB
figure(3);
clf;
cols = {'k', 'b', 'r', 'g'};
dc = 1;
for d = {'l', 'a', 'b'}
    for i = 1:2
        subplot(2, 3, dc);
        hold on;
        axis square;

        if i == 1
            pdata = pop_mean_data{2};
            odata = obs_mean_data{2};
        elseif i == 2
            pdata = pop_mean_data{3};
            odata = obs_mean_data{3};
        end

        if strcmp(d, 'l')
            pop_x = pdata.mean_mean_gwL_glaven;
            pop_y = pdata.mean_mean_gwL_match;

            pop_sd = pdata.se_gwL_match;

            obs_x = odata.mean_gwL_glaven;
            obs_y = odata.mean_gwL_match;
        elseif strcmp(d, 'a')
            pop_x = pdata.mean_mean_gwa_glaven;
            pop_y = pdata.mean_mean_gwa_match;

            pop_sd = pdata.se_gwa_match;

            obs_x = odata.mean_gwa_glaven;
            obs_y = odata.mean_gwa_match;
        elseif strcmp(d, 'b')
            pop_x = pdata.mean_mean_gwb_glaven;
            pop_y = pdata.mean_mean_gwb_match;

            pop_sd = pdata.se_gwb_match;

            obs_x = odata.mean_gwb_glaven;
            obs_y = odata.mean_gwb_match;
        end

        plotch = pdata.pch;

        if strcmp(d, 'l')
            minlim = 0;
            maxlim = 100;
        else
            minlim = -100;
            maxlim = 100;
        end
        xlim([minlim, maxlim]);
        ylim([minlim, maxlim]);

        blueidxs = find(pdata.illum == 'blue');
        yellowidxs = find(pdata.illum == 'yellow');

        blue_x = pop_x(blueidxs);
        blue_y = pop_y(blueidxs);
        blue_sd = pop_sd(blueidxs);
        errorbar(blue_x, blue_y, blue_sd, [cols{i}, 'o'], 'MarkerFaceColor', cols{i}, 'LineWidth', 2);

        yellow_x = pop_x(yellowidxs);
        yellow_y = pop_y(yellowidxs);
        yellow_sd = pop_sd(yellowidxs);
        errorbar(yellow_x, yellow_y, yellow_sd, [cols{i}, 'v'], 'LineWidth', 2);

        p = polyfit(obs_x, obs_y, 1);
        f1 = polyval(p, [minlim, maxlim]);
        plot([minlim, maxlim], f1, [cols{i}, '--'], 'LineWidth', 2);
    end
    rl = refline(1, 0);
    rl.Color = 'k';
    rl.LineWidth = 1.5;
    if strcmp(d, 'l')
        xlabel(['Mean Image ', upper(d{1}), '*']);
        ylabel(['Observer Match ', upper(d{1}), '*']);
    else
        xlabel(['Mean Image ', d{1}, '*']);
        ylabel(['Observer Match ', d{1}, '*']);
    end

    if strcmp(d, 'l')
        %         axis(1, at = seq(0, 100, 20));
        %         axis(2, at = seq(0, 100, 20));
        legend('bottomright', {'Blue Illuminant', 'White Illuminant', 'Mean Proximal Match', 'Mean Dye Match'});
    elseif strcmp(d, 'a')
        %         axis(1, at = seq(-50, 50, 50));
        %         axis(2, at = seq(-50, 50, 50));
        %         axis(1, at = -100, labels = '-100', col.axis = 5);
        %         axis(1, at = 100, labels = '100', col.axis = 4);
        %         axis(2, at = -100, labels = '-100', col.axis = 5);
        %         axis(2, at = 100, labels = '100', col.axis = 4);
    elseif strcmp(d, 'b')
        %         axis(1, at = seq(-50, 50, 50));
        %         axis(2, at = seq(-50, 50, 50));
        %         axis(1, at = -100, labels = '-100', col.axis = 7);
        %         axis(1, at = 100, labels = '100', col.axis = 6);
        %         axis(2, at = -100, labels = '-100', col.axis = 7);
        %         axis(2, at = 100, labels = '100', col.axis = 6);
    end

    dc = dc + 1;
end

for d = {'l', 'a', 'b'}
    for i = 3:4
        subplot(2, 3, dc);
        hold on;
        axis square;

        if i == 3
            pdata = pop_mean_data{1};
            odata = obs_mean_data{1};
        elseif i == 4
            pdata = pop_mean_data{1};
            odata = obs_mean_data{1};
        end

        if strcmp(d, 'l')
            pop_x = pdata.mean_mean_gwL_glaven;
            if i ~= 3
                pop_y = pdata.mean_mean_gwL_match;
            else
                pop_y = pdata.mean_mean_wp_top5L_filter;
            end

            if i ~= 3
                pop_sd = pdata.se_gwL_match;
            else
                pop_sd = pdata.se_wp_top5L_filter;
            end

            obs_x = odata.mean_gwL_glaven;
            if i ~= 3
                obs_y = odata.mean_gwL_match;
            else
                obs_y = odata.mean_wp_top5L_filter;
            end
        elseif strcmp(d, 'a')
            pop_x = pdata.mean_mean_gwa_glaven;
            if i ~= 3
                pop_y = pdata.mean_mean_gwa_match;
            else
                pop_y = pdata.mean_mean_wp_top5a_filter;
            end

            if i ~= 3
                pop_sd = pdata.se_gwa_match;
            else
                pop_sd = pdata.se_wp_top5a_filter;
            end

            obs_x = odata.mean_gwa_glaven;
            if i ~= 3
                obs_y = odata.mean_gwa_match;
            else
                obs_y = odata.mean_wp_top5a_filter;
            end
        elseif strcmp(d, 'b')
            pop_x = pdata.mean_mean_gwb_glaven;
            if i ~= 3
                pop_y = pdata.mean_mean_gwb_match;
            else
                pop_y = pdata.mean_mean_wp_top5b_filter;
            end

            if i ~= 3
                pop_sd = pdata.se_gwb_match;
            else
                pop_sd = pdata.se_wp_top5b_filter;
            end

            obs_x = odata.mean_gwb_glaven;
            if i ~= 3
                obs_y = odata.mean_gwb_match;
            else
                obs_y = odata.mean_wp_top5b_filter;
            end
        end

        plotch = pdata.pch;

        if strcmp(d, 'l')
            minlim = 0;
            maxlim = 100;
        else
            minlim = -100;
            maxlim = 100;
        end
        xlim([minlim, maxlim]);
        ylim([minlim, maxlim]);

        blueidxs = find(pdata.illum == 'blue');
        yellowidxs = find(pdata.illum == 'yellow');

        blue_x = pop_x(blueidxs);
        blue_y = pop_y(blueidxs);
        blue_sd = pop_sd(blueidxs);
        errorbar(blue_x, blue_y, blue_sd, [cols{i}, 'o'], 'MarkerFaceColor', cols{i}, 'LineWidth', 2);

        yellow_x = pop_x(yellowidxs);
        yellow_y = pop_y(yellowidxs);
        yellow_sd = pop_sd(yellowidxs);
        errorbar(yellow_x, yellow_y, yellow_sd, [cols{i}, 'v'], 'LineWidth', 2);

        p = polyfit(obs_x, obs_y, 1);
        f1 = polyval(p, [minlim, maxlim]);
        plot([minlim, maxlim], f1, [cols{i}, '--'], 'LineWidth', 2);
    end
    rl = refline(1, 0);
    rl.Color = 'k';
    rl.LineWidth = 1.5;
    if strcmp(d, 'l')
        xlabel(['Mean Image ', upper(d{1}), '*']);
        ylabel(['Observer Match ', upper(d{1}), '*']);
    else
        xlabel(['Mean Image ', d{1}, '*']);
        ylabel(['Observer Match ', d{1}, '*']);
    end

    if strcmp(d, 'l')
        legend('bottomright', {'Blue Illuminant', 'White Illuminant', 'Mean Filter Match', 'White Point Filter Match'});
    elseif strcmp(d, 'a')
        ticklabels = get(gca,'YTickLabel');
        ticklabels_new = cell(size(ticklabels));
        ticklabels_new{1} = ['\color{green} ' ticklabels{1}];
        for i = 2:length(ticklabels)-1
            ticklabels_new{i} = ticklabels{i};
        end
        ticklabels_new{end} = ['\color{red} ' ticklabels{end}];
        set(gca, 'YTickLabel', ticklabels_new);

        ticklabels = get(gca,'XTickLabel');
        ticklabels_new = cell(size(ticklabels));
        ticklabels_new{1} = ['\color{green} ' ticklabels{1}];
        for i = 2:length(ticklabels)-1
            ticklabels_new{i} = ticklabels{i};
        end
        ticklabels_new{end} = ['\color{red} ' ticklabels{end}];
        set(gca, 'XTickLabel', ticklabels_new);
    elseif strcmp(d, 'b')
        ticklabels = get(gca,'YTickLabel');
        ticklabels_new = cell(size(ticklabels));
        ticklabels_new{1} = ['\color{yellow} ' ticklabels{1}];
        for i = 2:length(ticklabels)-1
            ticklabels_new{i} = ticklabels{i};
        end
        ticklabels_new{end} = ['\color{blue} ' ticklabels{end}];
        set(gca, 'YTickLabel', ticklabels_new);

        ticklabels = get(gca,'XTickLabel');
        ticklabels_new = cell(size(ticklabels));
        ticklabels_new{1} = ['\color{yellow} ' ticklabels{1}];
        for i = 2:length(ticklabels)-1
            ticklabels_new{i} = ticklabels{i};
        end
        ticklabels_new{end} = ['\color{blue} ' ticklabels{end}];
        set(gca, 'XTickLabel', ticklabels_new);
    end

    dc = dc + 1;
end

%% supporting functions
function [y] = clamp(x)
    y = x;
    y(x < 0) = 0;
    y(x > 1) = 1;
    return;
end