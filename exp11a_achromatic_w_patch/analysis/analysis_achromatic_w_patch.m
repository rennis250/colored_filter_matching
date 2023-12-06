clear all; close all; clc;

% Change default axes fonts.
set(0,'DefaultAxesFontSize', 15);

% Change default text fonts.
set(0,'DefaultTextFontSize', 15);

load('object_illum_colors_for_plot.mat');
load('illums_lab.mat');

monxyY = [.6804 .3073 30.94
    .2029 .6968 74.22
    .1527 .0508 6.74];
monxyz = xyY2XYZ(monxyY);

which_exp = 'achromatic_w_patch';
conditions = {'p'};

wd = 400;
he = 400;

de_dist = 25;
num_lum_bands = 5;

trans_color_for_plot(1, :) = [0 0 1];
trans_color_for_plot(2, :) = [1.0 0.5 0.0];
trans_color_for_plot(3, :) = [1 0 0];
trans_color_for_plot(4, :) = [0 1 0];
trans_color_for_plot(5, :) = [0 0 0.5];
trans_color_for_plot(6, :) = [0.5, 0.25, 0.0];
trans_color_for_plot(7, :) = [0.5 0 0];
trans_color_for_plot(8, :) = [0 0.5 0];

illum_color_for_plot(1, :) = object_colors(1, :);
illum_color_for_plot(2, :) = object_colors(2, :);

% prepare masks to extract object, highlight, and caustic colors
pm = imread('images/objMask.png');
pm = pm(:,:,1); pm = logical(pm);
gisp = find(pm);
pm_wout_high = pm;

bm = imread('images/highMask.png');
br = bm(:,:,1); bg = bm(:,:,2); bb = bm(:,:,3);
idxs = br < 255 | bg < 255 | bb < 255;
hr2 = br; hg2 = bg; hb2 = bb;
hr2(idxs) = 0; hg2(idxs) = 0; hb2(idxs) = 0;
highlight2(:,:,1) = hr2;
highlight2(:,:,2) = hg2;
highlight2(:,:,3) = hb2;
bm = highlight2(:,:,1); bm = logical(bm);
gisb = find(bm);

pm_wout_high(gisb) = logical(0.0); % don't include highlight in estimate of color
gisp_wout_high = find(pm_wout_high);

rm = imread('images/reflMask.png');
rr = rm(:,:,1); rg = rm(:,:,2); rb = rm(:,:,3);
gidxs = rr == 255 & rg == 0 & rb == 0;
hr2 = rr; hg2 = rg; hb2 = rb;
hr2(gidxs) = 255; hg2(gidxs) = 255; hb2(gidxs) = 255;
bidxs = rr < 255 | rg < 255 | rb < 255;
hr2(bidxs) = 255; hg2(bidxs) = 255; hb2(bidxs) = 255;
hr2(gidxs) = 0; hg2(gidxs) = 0; hb2(gidxs) = 0;
refl2(:,:,1) = hr2;
refl2(:,:,2) = hg2;
refl2(:,:,3) = hb2;
rm = refl2(:,:,1); rm = ~logical(rm);
gisp_wout_refls = find(rm);

pm_only_refls = imread('images/reflMask.png');
rr = pm_only_refls(:,:,1); rg = pm_only_refls(:,:,2); rb = pm_only_refls(:,:,3);
bidxs = rr == 255 & rg == 0 & rb == 0;
hr2 = rr; hg2 = rg; hb2 = rb;
hr2(bidxs) = 255; hg2(bidxs) = 255; hb2(bidxs) = 255;
gidxs = hr2 < 255 | hg2 < 255 | hb2 < 255;
hr2(gidxs) = 0; hg2(gidxs) = 0; hb2(gidxs) = 0;
refl2(:,:,1) = hr2;
refl2(:,:,2) = hg2;
refl2(:,:,3) = hb2;
pm_only_refls = refl2(:,:,1); pm_only_refls = ~logical(pm_only_refls);
gisp_refls = find(pm_only_refls);

sm = logical(ones(size(pm,1),size(pm,2)));

xc = 121; yc = 129;
backx = 212; backy = 61;
shadx = 600; shady = 281;

gR = 2.2; gG = 2.2; gB = 2.2;

for x = 1:length(conditions)
    condition = conditions{x};
    fns = dir(['data/' which_exp '/*_' condition '.mat']);

    if strcmp(condition, 'u')
        ims = dir('images/uniform/mitsuba_caus*.png');
    elseif strcmp(condition, 'p')
        ims = dir('images/achromatic_exp/*.png');
    end
    imgs = zeros(1,wd,he,3);
    ims2 = {};
    which_trans = [];
    xc = 1;
    for x = 1:length(ims)
        % we restrict analysis to blue-yellow illums only now, so only focus on that
        fparts = strsplit(ims(x).name, '_');
        if strcmp(fparts{end}, 'lighter.png')
            refln = fparts{9};
            ld = 'darker';
        else
            refln = fparts{9};
            ld = 'lighter';
        end
        imgs(xc,:,:,:) = im2double(imread(['images/achromatic_exp/' ims(x).name]));
        ims2{xc} = ims(x).name;

        if strcmp(refln, 'blue')
            if strcmp(ld, 'lighter')
                colors_for_scatter_plot(xc, :) = [0 0 1];
                which_trans(xc) = 1;
            else
                colors_for_scatter_plot(xc, :) = [0 0 0.5];
                which_trans(xc) = 5;
            end
        elseif strcmp(refln, 'yellow')
            if strcmp(ld, 'lighter')
                colors_for_scatter_plot(xc, :) = [1 0.5 0];
                which_trans(xc) = 2;
            else
                colors_for_scatter_plot(xc, :) = [0.5 0.25 0];
                which_trans(xc) = 6;
            end
        elseif strcmp(refln, 'red')
            if strcmp(ld, 'lighter')
                colors_for_scatter_plot(xc, :) = [1 0 0];
                which_trans(xc) = 3;
            else
                colors_for_scatter_plot(xc, :) = [0.5 0 0];
                which_trans(xc) = 7;
            end
        elseif strcmp(refln, 'green')
            if strcmp(ld, 'lighter')
                colors_for_scatter_plot(xc, :) = [0 1 0];
                which_trans(xc) = 4;
            else
                colors_for_scatter_plot(xc, :) = [0 0.5 0];
                which_trans(xc) = 8;
            end
        end
        xc = xc + 1;
    end

    avgdata = [];
    avgdata_lab = [];
    data_lab = cell(length(fns), size(imgs, 1));
    data_rgb = cell(length(fns), size(imgs, 1));

    for fc = 1:length(fns)
        clear data;
        clear subID;
        clear trialOrder;
        clear ntrials;
        load(['../data/' which_exp '/' fns(fc).name], 'data', 'subID', 'trialOrder', 'ntrials');

        c = 1;
        for tc = 1:size(ims, 1)
            idxs = find(trialOrder == tc);
            avgdata(fc, c, :) = real(mean(squeeze(data(idxs, 1, :)), 1));
            data_lab{fc, c} = [];
            data_rgb{fc, c} = [];
            for zc = 1:length(idxs)
                if any(data(idxs(zc), 1, :) < 0 | data(idxs(zc), 1, :) > 1)
                    continue;
                else
                    tl = real(rgb2labRob(squeeze(data(idxs(zc), 1, :))', monxyz));
                    tr = real(squeeze(data(idxs(zc), 1, :))');
                    data_lab{fc, c} = [data_lab{fc, c}; tl];
                    data_rgb{fc, c} = [data_rgb{fc, c}; tr];
                end
            end
            avgdata_lab(fc, c, :) = real(mean(data_lab{fc, c}, 1));
            stddata_lab(fc, c, :) = real(std(data_lab{fc, c}, 0, 1));
            avgdata_rgb(fc, c, :) = real(mean(data_rgb{fc, c}, 1));

            c = c + 1;
        end
    end

    % average and SEM of all observer data
    popavgdata_rgb = [];
    popavgdata_lab = [];
    popsemdata_lab = [];
    % temporary b/c only one obs right now
    % for tc = 1:size(imgs, 1)
    %    popavgdata_rgb(tc, :) = real(mean(squeeze(avgdata(:, tc, :)), 1));
    %    popavgdata_lab(tc, :) = real(mean(squeeze(avgdata_lab(:, tc, :)), 1));
    %    popsemdata_lab(tc, :) = real(std(squeeze(avgdata_lab(:, tc, :)), 0, 1));
    % end
    for tc = 1:size(imgs, 1)
       popavgdata_rgb(tc, :) = real(avgdata(:, tc, :));
       popavgdata_lab(tc, :) = real(avgdata_lab(:, tc, :));
       popsemdata_lab(tc, :) = real(avgdata_lab(:, tc, :));
    end

    obj_lab = [];
    high_lab = [];
    mean_obj_lab = [];
    mean_high_lab = [];
    lab_lum_splits = {};
    lum_split_diffs = zeros(1,3);
    lum_split_diffs_wout_refls = zeros(1,3);

    for ic = 1:size(imgs, 1)
        curr_img = squeeze(imgs(ic, :, :, :));

        mean_pop_match = squeeze(popavgdata_lab(ic, :));

        % first let's compute mean colors of regions of interest
        [imp, imh] = maskAndLinearTransparent(curr_img, pm, bm, sm, xc, yc, backx, backy, shadx, shady, gR, gG, gB);
        [imp_wout_high, ~] = maskAndLinearTransparent(curr_img, pm_wout_high, bm, sm, xc, yc, backx, backy, shadx, shady, gR, gG, gB);
        [imp_wout_refls, ~] = maskAndLinearTransparent(curr_img, rm, bm, sm, xc, yc, backx, backy, shadx, shady, gR, gG, gB);
        [imp_refls, ~] = maskAndLinearTransparent(curr_img, pm_only_refls, bm, sm, xc, yc, backx, backy, shadx, shady, gR, gG, gB);

        rp = imp(:,:,1); rp = rp(:);
        gp = imp(:,:,2); gp = gp(:);
        bp = imp(:,:,3); bp = bp(:);
        obj_rgb(ic, :, :) = [rp(gisp) gp(gisp) bp(gisp)];
        obj_lab(ic, :, :) = rgb2labRob(squeeze(obj_rgb(ic, :, :)), monxyz);
        obj_lms(ic, :, :) = rgb2lms*squeeze(obj_rgb(ic, :, :))';
        mean_obj_rgb(ic, :) = mean(squeeze(obj_rgb(ic, :, :)), 1);
        mean_obj_lab(ic, :) = mean(squeeze(obj_lab(ic, :, :)), 1);

        rb = imb(:,:,1); rb = rb(:);
        gb = imb(:,:,2); gb = gb(:);
        bb = imb(:,:,3); bb = bb(:);
        bkg_rgb(ic, :, :) = [rb(gisb) gb(gisb) bb(gisb)];
        bkg_lab(ic, :, :) = rgb2labRob(squeeze(bkg_rgb(ic, :, :)), monxyz);
        bkg_lms(ic, :, :) = rgb2lms*squeeze(bkg_rgb(ic, :, :))';
        mean_bkg_rgb(ic, :) = mean(squeeze(bkg_rgb(ic, :, :)), 1);
        mean_bkg_lab(ic, :) = mean(squeeze(bkg_lab(ic, :, :)), 1);

        rp = imp_wout_high(:,:,1); rp = rp(:);
        gp = imp_wout_high(:,:,2); gp = gp(:);
        bp = imp_wout_high(:,:,3); bp = bp(:);
        obj_rgb_wout_high(ic, :, :) = [rp(gisp_wout_high) gp(gisp_wout_high) bp(gisp_wout_high)];
        obj_lab_wout_high(ic, :, :) = rgb2labRob(squeeze(obj_rgb_wout_high(ic, :, :)), monxyz);
        obj_lms_wout_high(ic, :, :) = rgb2lms*squeeze(obj_rgb_wout_high(ic, :, :))';
        mean_obj_rgb_wout_high(ic, :) = mean(squeeze(obj_rgb_wout_high(ic, :, :)), 1);
        mean_obj_lab_wout_high(ic, :) = mean(squeeze(obj_lab_wout_high(ic, :, :)), 1);

        rp = imp_wout_refls(:,:,1); rp = rp(:);
        gp = imp_wout_refls(:,:,2); gp = gp(:);
        bp = imp_wout_refls(:,:,3); bp = bp(:);
        obj_rgb_wout_refls(ic, :, :) = [rp(gisp_wout_refls) gp(gisp_wout_refls) bp(gisp_wout_refls)];
        obj_lab_wout_refls(ic, :, :) = rgb2labRob(squeeze(obj_rgb_wout_refls(ic, :, :)), monxyz);
        obj_lms_wout_refls(ic, :, :) = rgb2lms*squeeze(obj_rgb_wout_refls(ic, :, :))';
        mean_obj_rgb_wout_refls(ic, :) = mean(squeeze(obj_lab_wout_refls(ic, :, :)), 1);
        mean_obj_lab_wout_refls(ic, :) = mean(squeeze(obj_lab_wout_refls(ic, :, :)), 1);

        rp = imp_refls(:,:,1); rp = rp(:);
        gp = imp_refls(:,:,2); gp = gp(:);
        bp = imp_refls(:,:,3); bp = bp(:);
        refls_rgb(ic, :, :) = [rp(gisp_refls) gp(gisp_refls) bp(gisp_refls)];
        refls_lab(ic, :, :) = rgb2labRob(squeeze(refls_rgb(ic, :, :)), monxyz);
        refls_lms(ic, :, :) = rgb2lms*squeeze(refls_lms(ic, :, :))';
        mean_refls_rgb(ic, :) = mean(squeeze(refls_rgb(ic, :, :)), 1);
        mean_refls_lab(ic, :) = mean(squeeze(refls_lab(ic, :, :)), 1);

        rh = imh(:,:,1); rh = rh(:);
        gh = imh(:,:,2); gh = gh(:);
        bh = imh(:,:,3); bh = bh(:);
        high_rgb(ic, :, :) = [rh(gisb) gh(gisb) bh(gisb)];
        high_lab(ic, :, :) = rgb2labRob(squeeze(high_rgb(ic, :, :)), monxyz);
        high_lms(ic, :, :) = rgb2lms*squeeze(high_lms(ic, :, :))';
        mean_high_rgb(ic, :) = mean(squeeze(high_rgb(ic, :, :)), 1);
        mean_high_lab(ic, :) = mean(squeeze(high_lab(ic, :, :)), 1);

        lum_regions = linspace(0, squeeze(max(max(obj_lab(ic, :, 1)))), num_lum_bands);
        lum_regions_wout_refls = linspace(0, squeeze(max(max(obj_lab_wout_refls(ic, :, 1)))), num_lum_bands);

        % next, let's try to get regions of different luminance and see what the mean colors are
        c = 1;
        for lc = 1:length(lum_regions)-1
            ls = squeeze(obj_lab(ic, :, 1));
            lab_lum_splits{ic, c} = squeeze(obj_lab(ic, ls >= lum_regions(lc) & ls < lum_regions(lc+1), :));
            mean_lab_lum_splits(ic, c, :) = mean(lab_lum_splits{ic, c}, 1);

            ls_wout_refls = squeeze(obj_lab_wout_refls(ic, :, 1));
            lab_lum_splits_wout_refls{ic, c} = squeeze(obj_lab_wout_refls(ic, ls_wout_refls >= lum_regions_wout_refls(lc) & ls_wout_refls < lum_regions_wout_refls(lc+1), :));
            mean_lab_lum_splits_wout_refls(ic, c, :) = mean(lab_lum_splits_wout_refls{ic, c}, 1);

            c = c + 1;
        end

        % and we can find the closest point in the luminance bands to the matched color
        for sc = 1:length(fns)
            for lc = 1:size(mean_lab_lum_splits, 2)
                lum_split_diffs(sc, ic, lc) = norm(squeeze(mean_lab_lum_splits(ic, lc, :)) - squeeze(avgdata_lab(sc, ic, :)));
            end
            best_lum_split(sc, ic) = find(lum_split_diffs(sc, ic, :) == min(lum_split_diffs(sc, ic, :)));

            % check the same for obj without reflections included
           for lc = 1:c-1
                lum_split_diffs_wout_refls(sc, ic, lc) = norm(squeeze(mean_lab_lum_splits_wout_refls(ic, lc, :)) - squeeze(avgdata_lab(sc, ic, :)));
            end
            best_lum_split_wout_refls(sc, ic) = find(lum_split_diffs_wout_refls(sc, ic, :) == min(lum_split_diffs_wout_refls(sc, ic, :)));

            for lc = 1:c-1
                lum_split_diffs_wout_lum(sc, ic, lc) = norm(squeeze(mean_lab_lum_splits(ic, lc, 2:3)) - squeeze(avgdata_lab(sc, ic, 2:3)));
            end
            best_lum_split_wout_lum(sc, ic) = find(lum_split_diffs_wout_lum(sc, ic, :) == min(lum_split_diffs_wout_lum(sc, ic, :)));
        end

        % now, let's try finding if matches were closer to mean of obj or mean of highlight
        % compute for chromaticity only, since observers are overestimating luminance by
        % a huge degree
        for sc = 1:length(fns)
            dist_from_obj(sc, ic) = norm(squeeze(avgdata_lab(sc, ic, 2:3))' - mean_obj_lab_wout_high(ic, 2:3));
            dist_from_obj_wout_refls(sc, ic) = norm(squeeze(avgdata_lab(sc, ic, 2:3))' - mean_obj_lab_wout_refls(ic, 2:3));
            dist_from_high(sc, ic) = norm(squeeze(avgdata_lab(sc, ic, 2:3))' - mean_high_lab(ic, 2:3));
        end

        % let's grab pixels of object that were 5 DE units away from color match
        de_units(ic, :) = deltaE2000(repmat(mean_pop_match, size(obj_lab, 2), 1), squeeze(obj_lab(ic, :, :)));
        de_units_wout_refls(ic, :) = deltaE2000(repmat(mean_pop_match, size(obj_lab_wout_refls, 2), 1), squeeze(obj_lab_wout_refls(ic, :, :)));
    end
end

% save data to plot experiments together
save('exp5_data.mat');

do_crazy_amount_of_plots = false;

if do_crazy_amount_of_plots
    %%%% summary plots of all images at once

    %%%%%%%%%%%%%%%%%%% correlation plots

    % plot showing correlation between mean colors of objects and mean matches
    figure(4);
    clf;
    subplot(1,3,1);
    hold on;
    axis square;
    xlabel("L* mean object");
    ylabel("L* mean match");
    l = line([-100 100], [-100 100]);
    set(l, 'LineWidth', 2);
    set(l, 'LineStyle', '--');
    axis([0 100 0 100]);
    axis square;

    rs_ps = [];

    [r, p] = corr(mean_obj_lab(:, 1), popavgdata_lab(:, 1));

    rs_ps(1, 1) = r;
    rs_ps(1, 2) = p;

    p1 = scatter(mean_obj_lab(:, 1), popavgdata_lab(:, 1), 100, colors_for_scatter_plot(:, :), 'o', 'filled');

    title(['r = ' num2str(r, 3) '*, p < 0.0001']);

    subplot(1,3,2);
    hold on;
    axis square;
    xlabel("a* mean object");
    ylabel("a* mean match");

    [r, p] = corr(mean_obj_lab(:, 2), popavgdata_lab(:, 2));

    rs_ps(2, 1) = r;
    rs_ps(2, 2) = p;

    p1 = scatter(mean_obj_lab(:, 2), popavgdata_lab(:, 2), 100, colors_for_scatter_plot(:, :), 'o', 'filled');

    title(['r = ' num2str(r, 3) '*, p < 0.0001']);

    l = line([-100 100], [-100 100]);
    set(l, 'LineWidth', 2);
    set(l, 'LineStyle', '--');
    axis square;

    subplot(1,3,3);
    hold on;
    axis square;
    xlabel("b* mean object");
    ylabel("b* mean match");

    [r, p] = corr(mean_obj_lab(:, 3), popavgdata_lab(:, 3));

    rs_ps(3, 1) = r;
    rs_ps(3, 2) = p;

    p1 = scatter(mean_obj_lab(:, 3), popavgdata_lab(:, 3), 100, colors_for_scatter_plot(:, :), 'o', 'filled');

    title(['r = ' num2str(r, 3) '*, p < 0.009']);

    l = line([-100 100], [-100 100]);
    set(l, 'LineWidth', 2);
    set(l, 'LineStyle', '--');
    axis square;
    export_fig(['~/Dropbox/transparent_color_skin_of_things_nov_2018/' which_exp '/plots/mean_match_obj_correlation.pdf'], '-pdf');

    dlmwrite(['~/Dropbox/transparent_color_skin_of_things_nov_2018/' which_exp '/correlations.csv'], rs_ps);

    % just write out correlations of matches with mean of highlight color

    % first for all illums together
    rs_ps = [];

    [r, p] = corr(mean_high_lab(:, 1), popavgdata_lab(:, 1));

    rs_ps(1, 1) = r;
    rs_ps(1, 2) = p;

    [r, p] = corr(mean_high_lab(:, 2), popavgdata_lab(:, 2));

    rs_ps(2, 1) = r;
    rs_ps(2, 2) = p;

    [r, p] = corr(mean_high_lab(:, 3), popavgdata_lab(:, 3));

    rs_ps(3, 1) = r;
    rs_ps(3, 2) = p;

    dlmwrite(['~/Dropbox/transparent_color_skin_of_things_nov_2018/' which_exp '/correlations_high_only.csv'], rs_ps);

    %%%%%%%%%%%%%%%%% end correlation plots

    %%%%% luminance region plots

    % illums together
    figure(5);
    for x = 1:2 % with luminance in distance computation and without luminance
        if x == 1
            lum_split_to_plot = best_lum_split;
        else
            lum_split_to_plot = best_lum_split_wout_lum;
        end
        % plot showing lum regions with best match
        clf;
        hold on;
        for sc = 1:size(lum_split_to_plot, 1)
            bs = [length(find(lum_split_to_plot(sc, :) == 1)) length(find(lum_split_to_plot(sc, :) == 2)) length(find(lum_split_to_plot(sc, :) == 3)) length(find(lum_split_to_plot(sc, :) == 4))];
            tot_best_lum_split(sc, :) = bs;
        end
        avg_best_lum_split = mean(tot_best_lum_split, 1);
        sd_best_lum_split = std(tot_best_lum_split, 0, 1);
        b = bar(1:4, avg_best_lum_split, 'k');
        b.FaceColor = [1 1 1];
        eb = errorbar(1:4, avg_best_lum_split, sd_best_lum_split, '.');
        eb.Color = [0 0 0];
        eb.LineWidth = 2;

        % bs = [length(find(best_lum_split_wout_refls == 1)) length(find(best_lum_split_wout_refls == 2)) length(find(best_lum_split_wout_refls == 3))];
        % p = bar(1:3, bs);

        axis([0.5 4.5 0 15]);
        axis square;
        xlabel('Luminance (L*) band');
        ylabel('Frequency');
        title('Freq. of avg. color of L* band being closest to avg. color matches');
        xticks([1 2 3 4]);
        xticklabels({'1st quartile', '2nd quartile', '3rd quartile', 'Brightest quartile'});
        if x == 1
            export_fig(['~/Dropbox/transparent_color_skin_of_things_nov_2018/' which_exp '/plots/hist_lum_regions_w_lum.pdf'], '-pdf');
        else
            export_fig(['~/Dropbox/transparent_color_skin_of_things_nov_2018/' which_exp '/plots/hist_lum_regions_wout_lum.pdf'], '-pdf');
        end
    end

    %%%%%%%%%%%%% end luminance region plots

    % plot showing histogram of distances from obj and high to get an idea of how often matches are
    % pulled towards the highlight color
    figure(6);
    clf;
    hold on;
    [h,p,~,st] = ttest(dist_from_obj, dist_from_high);
    hist(dist_from_obj_wout_refls(:));
    hist(dist_from_high(:));
    hist(dist_from_obj(:));
    hs = findobj(gca,'Type','patch');
    h1 = hs(1);
    h1.FaceColor = [0 0 0];
    h1.EdgeColor = 'w';
    h2 = hs(2);
    h2.FaceColor = [0 0.5 0.5];
    h2.EdgeColor = 'w';
    h3 = hs(3);
    h3.FaceColor = [0.5 0.5 0.5];
    h3.EdgeColor = 'w';
    axis([0 165 0 45]);
    axis square;
    xlabel('Euclidean Distance (CIELAB units)');
    ylabel('Frequency');
    title('Distance from avg. color matches');
    legend([hs(3) hs(1) hs(2)], 'Dist. from avg. object color', 'Dist. from obj w/out spec. refls.', 'Dist. from avg. highlight color');
    export_fig(['~/Dropbox/transparent_color_skin_of_things_nov_2018/' which_exp '/plots/hist_dist_matches.pdf'], '-pdf');

    % loop through each image and output a version with pixels that are 5 DE units from match in white
    % as well as a version where the different luminance bands are highlighted in white, all computed
    % for the object
    mean_de_25 = [];
    for ic = 1:size(imgs, 1)
        fn = ims2{ic};
        curr_img = squeeze(imgs(ic, :, :, :));

        curr_img_r = curr_img(:,:,1); curr_img_r(~pm) = 0.0;
        curr_img_g = curr_img(:,:,2); curr_img_g(~pm) = 0.0;
        curr_img_b = curr_img(:,:,3); curr_img_b(~pm) = 0.0;

        curr_img(:,:,1) = curr_img_r;
        curr_img(:,:,2) = curr_img_g;
        curr_img(:,:,3) = curr_img_b;

        curr_img_gray_tmp = rgb2gray(curr_img);
        curr_img_gray = zeros(size(curr_img));
        curr_img_gray(:,:,1) = curr_img_gray_tmp;
        curr_img_gray(:,:,2) = curr_img_gray_tmp;
        curr_img_gray(:,:,3) = curr_img_gray_tmp;

        idxs = gisp(find(de_units(ic, :) <= de_dist));
        curr_img_col = reshape(curr_img, size(curr_img,1)*size(curr_img,2), 3);
        curr_img_gray_col = reshape(curr_img_gray, size(curr_img,1)*size(curr_img,2), 3);

        curr_img_gray_col(idxs, :) = curr_img_col(idxs, :);

        % compute mean color of 25 DE unit pixels

        mean_de_25(ic, :) = mean(squeeze(obj_lab(ic, de_units(ic, :) <= de_dist, :)), 1);

        % imwrite(reshape(curr_img_gray_col, size(curr_img,1), size(curr_img,2), size(curr_img,3)), ['~/Dropbox/transparent_color_skin_of_things_nov_2018/' which_exp '/images/de_25_units/' fn(1:end-4) '_de_25.png']);

        lum_regions = linspace(0, squeeze(max(max(obj_lab(ic, :, 1)))), num_lum_bands);

        for lc = 1:length(lum_regions)-1
            ls = squeeze(obj_lab(ic, :, 1));

            idxs = gisp(find(ls >= lum_regions(lc) & ls < lum_regions(lc+1)));
            curr_img_col = reshape(curr_img, size(curr_img,1)*size(curr_img,2), 3);
            curr_img_gray_col = reshape(curr_img_gray, size(curr_img,1)*size(curr_img,2), 3);

            curr_img_gray_col(idxs, :) = curr_img_col(idxs, :);

            % imwrite(reshape(curr_img_gray_col, size(curr_img,1), size(curr_img,2), size(curr_img,3)), ['~/Dropbox/transparent_color_skin_of_things_nov_2018/' which_exp '/images/lum_regions/' fn(1:end-4) '_lum_region_' num2str(lc) '.png']);
        end
    end

    %%%% compute correlation of mean of DE 25 unit pixels with mean of observer match

    % first for all illums together
    rs_ps = [];

    [r, p] = corr(mean_de_25(:, 1), popavgdata_lab(:, 1));

    rs_ps(1, 1) = r;
    rs_ps(1, 2) = p;

    [r, p] = corr(mean_de_25(:, 2), popavgdata_lab(:, 2));

    rs_ps(2, 1) = r;
    rs_ps(2, 2) = p;

    [r, p] = corr(mean_de_25(:, 3), popavgdata_lab(:, 3));

    rs_ps(3, 1) = r;
    rs_ps(3, 2) = p;

    dlmwrite(['~/Dropbox/transparent_color_skin_of_things_nov_2018/' which_exp '/correlations_DE_25_only.csv'], rs_ps);

    % make a plot showing the correlations for DE 25 pixels

    figure(4);
    clf;
    subplot(1,3,1);
    hold on;
    axis square;
    xlabel("L* mean object");
    ylabel("L* mean match");
    l = line([-100 100], [-100 100]);
    set(l, 'LineWidth', 2);
    set(l, 'LineStyle', '--');
    axis([0 100 0 100]);
    axis square;

    [r, p] = corr(mean_de_25(:, 1), popavgdata_lab(:, 1));

    p1 = scatter(mean_de_25(:, 1), popavgdata_lab(:, 1), 100, colors_for_scatter_plot(:, :), 'o', 'filled');

    title(['r = ' num2str(r, 3) '*, p < 0.0001']);

    subplot(1,3,2);
    hold on;
    axis square;
    xlabel("a* mean object");
    ylabel("a* mean match");

    [r, p] = corr(mean_de_25(:, 2), popavgdata_lab(:, 2));

    p1 = scatter(mean_de_25(:, 2), popavgdata_lab(:, 2), 100, colors_for_scatter_plot(:, :), 'o', 'filled');

    title(['r = ' num2str(r, 3) '*, p < 0.0001']);

    l = line([-100 100], [-100 100]);
    set(l, 'LineWidth', 2);
    set(l, 'LineStyle', '--');
    axis square;

    subplot(1,3,3);
    hold on;
    axis square;
    xlabel("b* mean object");
    ylabel("b* mean match");

    [r, p] = corr(mean_de_25(:, 3), popavgdata_lab(:, 3));

    p1 = scatter(mean_de_25(:, 3), popavgdata_lab(:, 3), 100, colors_for_scatter_plot(:, :), 'o', 'filled');

    title(['r = ' num2str(r, 3) '*, p < 0.009']);

    l = line([-100 100], [-100 100]);
    set(l, 'LineWidth', 2);
    set(l, 'LineStyle', '--');
    axis square;
    export_fig(['~/Dropbox/transparent_color_skin_of_things_nov_2018/' which_exp '/plots/mean_match_DE_25_correlation.pdf'], '-pdf');

    % loop through each image and output a version with pixels that are 5 DE units from match marked
    % but this time without reflections included in computation
    for ic = 1:size(imgs, 1)
        fn = ims2{ic};
        curr_img = squeeze(imgs(ic, :, :, :));

        curr_img_r = curr_img(:,:,1); curr_img_r(~pm) = 0.0;
        curr_img_g = curr_img(:,:,2); curr_img_g(~pm) = 0.0;
        curr_img_b = curr_img(:,:,3); curr_img_b(~pm) = 0.0;

        curr_img(:,:,1) = curr_img_r;
        curr_img(:,:,2) = curr_img_g;
        curr_img(:,:,3) = curr_img_b;

        curr_img_gray_tmp = rgb2gray(curr_img);
        curr_img_gray = zeros(size(curr_img));
        curr_img_gray(:,:,1) = curr_img_gray_tmp;
        curr_img_gray(:,:,2) = curr_img_gray_tmp;
        curr_img_gray(:,:,3) = curr_img_gray_tmp;

        idxs = gisp_wout_refls(find(de_units_wout_refls(ic, :) <= de_dist));
        curr_img_col = reshape(curr_img, size(curr_img,1)*size(curr_img,2), 3);
        curr_img_gray_col = reshape(curr_img_gray, size(curr_img,1)*size(curr_img,2), 3);

        curr_img_gray_col(idxs, :) = curr_img_col(idxs, :);

        imwrite(reshape(curr_img_gray_col, size(curr_img,1), size(curr_img,2), size(curr_img,3)), ['~/Dropbox/transparent_color_skin_of_things_nov_2018/' which_exp '/images/de_25_units_wout_refls/' fn(1:end-4) '_de_25.png']);

        lum_regions_wout_refls = linspace(0, squeeze(max(max(obj_lab_wout_refls(ic, :, 1)))), num_lum_bands);

        for lc = 1:length(lum_regions_wout_refls)-1
            ls = squeeze(obj_lab_wout_refls(ic, :, 1));

            idxs = gisp_wout_refls(find(ls >= lum_regions_wout_refls(lc) & ls < lum_regions_wout_refls(lc+1)));
            curr_img_col = reshape(curr_img, size(curr_img,1)*size(curr_img,2), 3);
            curr_img_gray_col = reshape(curr_img_gray, size(curr_img,1)*size(curr_img,2), 3);

            curr_img_gray_col(idxs, :) = curr_img_col(idxs, :);

            imwrite(reshape(curr_img_gray_col, size(curr_img,1), size(curr_img,2), size(curr_img,3)), ['~/Dropbox/transparent_color_skin_of_things_nov_2018/' which_exp '/images/lum_regions_wout_refls/' fn(1:end-4) '_lum_region_' num2str(lc) '.png']);
        end
    end

    do_big_plot = true;

    % now make big per image plots showing average settings. make multiple so that they
    % can build up over slides

    figure(2);
    for ic = 1:size(imgs, 1)
        % now make plots for each image
        if do_big_plot
            bnd_idxs = boundary(squeeze(obj_lab(ic, :, 2:3)));
            bnd_idxs_refls = boundary(squeeze(refls_lab(ic, :, 2:3)));
            bnd_idxs_high = boundary(squeeze(high_lab(ic, :, 2:3)));

            mean_match_rgb = popavgdata_rgb(ic, :).^(1/2.2);

            % first plot - show LAB distribution of object only
            clf;
            subplot(1,2,1);
            hold on;
            xlabel('a*');
            ylabel('b*');
            axis([-150 150 -150 150]);
            axis square;
            grid;
            title(['Mean Lum - ' num2str(mean_obj_lab(ic, 1))]);
            p_obj_pts = scatter(obj_lab(ic, :, 2), obj_lab(ic, :, 3), [], squeeze(obj_rgb(ic, :, :).^(1/2.2)), 'filled');
            p_obj = plot(mean_obj_lab(ic, 2), mean_obj_lab(ic, 3), 's', 'MarkerFaceColor', mean_obj_rgb(ic, :).^(1/2.2), 'MarkerSize', 15);

            h_leg = legend([p_obj_pts, p_obj], {'Object', 'Mean Object'}, 'Location', 'NorthEastOutside');

            subplot(1,2,2);
            imshow(squeeze(imgs(ic, :, :, :)));

            fnparts = strsplit(ims2{ic}, '_');
            if strcmp(fnparts{end}, 'lightest.png')
               transmission = [fnparts{9} '_less_sat'];
            else
               transmission = [fnparts{9} '_more_sat'];
            end
            export_fig(['~/Dropbox/transparent_color_skin_of_things_nov_2018/' which_exp '/plots/achromatic_transmission_' transmission '_obj.pdf'], '-pdf');

            % third plot - show highlights also
            clf;
            subplot(1,2,1);
            hold on;
            xlabel('a*');
            ylabel('b*');
            axis([-150 150 -150 150]);
            axis square;
            grid;
            title(['Mean Lum - ' num2str(mean_obj_lab(ic, 1))]);
            p_obj_pts = plot(obj_lab(ic, bnd_idxs, 2), obj_lab(ic, bnd_idxs, 3), 'k-');
            % p_refl_pts = plot(refls_lab(ic, bnd_idxs_refls, 2), refls_lab(ic, bnd_idxs_refls, 3), 'b-');
            p_high_pts = scatter(high_lab(ic, :, 2), high_lab(ic, :, 3), [], squeeze(high_rgb(ic, :, :).^(1/2.2)), 'filled');
            p_obj = plot(mean_obj_lab(ic, 2), mean_obj_lab(ic, 3), 's', 'MarkerFaceColor', mean_obj_rgb(ic, :).^(1/2.2), 'MarkerSize', 15);
            % p_refl = plot(mean_refls_lab(ic, 2), mean_refls_lab(ic, 3), 'o', 'MarkerFaceColor', mean_refls_rgb(ic, :).^(1/2.2), 'MarkerSize', 15);
            p_high = plot(mean_high_lab(ic, 2), mean_high_lab(ic, 3), 'o', 'MarkerFaceColor', mean_high_rgb(ic, :).^(1/2.2), 'MarkerSize', 15);

            h_leg = legend([p_obj_pts, p_obj, p_high_pts, p_high], {'Object', 'Mean Object', 'Highlight', 'Mean Highlight'}, 'Location', 'NorthEastOutside');

            subplot(1,2,2);
            imshow(squeeze(imgs(ic, :, :, :)))

            fnparts = strsplit(ims2{ic}, '_');
            if strcmp(fnparts{end}, 'lightest.png')
               transmission = [fnparts{9} '_less_sat'];
            else
               transmission = [fnparts{9} '_more_sat'];
            end
            export_fig(['~/Dropbox/transparent_color_skin_of_things_nov_2018/' which_exp '/plots/achromatic_transmission_' transmission '_obj_high.pdf'], '-pdf');

            % fourth plot - show obs matches
            clf;
            subplot(1,2,1);
            hold on;
            xlabel('a*');
            ylabel('b*');
            axis([-150 150 -150 150]);
            axis square;
            grid;
            title(['Mean Lum - ' num2str(mean_obj_lab(ic, 1)) ' - Match Lum - ' num2str(popavgdata_lab(ic, 1))]);
            p_obj_pts = plot(obj_lab(ic, bnd_idxs, 2), obj_lab(ic, bnd_idxs, 3), 'k-');
            % p_refl_pts = plot(refls_lab(ic, bnd_idxs_refls, 2), refls_lab(ic, bnd_idxs_refls, 3), 'b-');
            % p_high_pts = plot(high_lab(ic, bnd_idxs_high, 2), high_lab(ic, bnd_idxs_high, 3), 'r-', 'Color', [0.5 0 0]);
            p_obj = plot(mean_obj_lab(ic, 2), mean_obj_lab(ic, 3), 's', 'MarkerFaceColor', mean_obj_rgb(ic, :).^(1/2.2), 'MarkerSize', 15);
            % p_refl = plot(mean_refls_lab(ic, 2), mean_refls_lab(ic, 3), 'o', 'MarkerFaceColor', mean_refls_rgb(ic, :).^(1/2.2), 'MarkerSize', 15);
            p_high = plot(mean_high_lab(ic, 2), mean_high_lab(ic, 3), 'o', 'MarkerFaceColor', mean_high_rgb(ic, :).^(1/2.2), 'MarkerSize', 15);

            for sc = 1:size(avgdata_lab, 1)
                p_matches = scatter(avgdata_lab(sc, ic, 2), avgdata_lab(sc, ic, 3), [], squeeze(avgdata_rgb(sc, ic, :).^(1/2.2))', 'filled');
            end

            h_leg = legend([p_obj_pts, p_obj, p_high, p_matches], {'Object', 'Mean Object', 'Mean Highlight', 'Obs. matches'}, 'Location', 'NorthEastOutside');

            subplot(1,2,2);
            imshow(squeeze(imgs(ic, :, :, :)))

            fnparts = strsplit(ims2{ic}, '_');
            if strcmp(fnparts{end}, 'lightest.png')
               transmission = [fnparts{9} '_less_sat'];
            else
               transmission = [fnparts{9} '_more_sat'];
            end
            export_fig(['~/Dropbox/transparent_color_skin_of_things_nov_2018/' which_exp '/plots/achromatic_transmission_' transmission '_obj_high_obs_matches.pdf'], '-pdf');

            % fifth plot - show mean of obs matches
            clf;
            subplot(1,2,1);
            hold on;
            xlabel('a*');
            ylabel('b*');
            axis([-150 150 -150 150]);
            axis square;
            grid;
            title(['Mean Lum - ' num2str(mean_obj_lab(ic, 1)) ' - Match Lum - ' num2str(popavgdata_lab(ic, 1))]);
            p_obj_pts = plot(obj_lab(ic, bnd_idxs, 2), obj_lab(ic, bnd_idxs, 3), 'k-');
            % p_refl_pts = plot(refls_lab(ic, bnd_idxs_refls, 2), refls_lab(ic, bnd_idxs_refls, 3), 'b-');
            % p_high_pts = plot(high_lab(ic, bnd_idxs_high, 2), high_lab(ic, bnd_idxs_high, 3), 'r-', 'Color', [0.5 0 0]);
            p_obj = plot(mean_obj_lab(ic, 2), mean_obj_lab(ic, 3), 's', 'MarkerFaceColor', mean_obj_rgb(ic, :).^(1/2.2), 'MarkerSize', 15);
            % p_refl = plot(mean_refls_lab(ic, 2), mean_refls_lab(ic, 3), 'o', 'MarkerFaceColor', mean_refls_rgb(ic, :).^(1/2.2), 'MarkerSize', 15);
            p_high = plot(mean_high_lab(ic, 2), mean_high_lab(ic, 3), 'o', 'MarkerFaceColor', mean_high_rgb(ic, :).^(1/2.2), 'MarkerSize', 15);

            for sc = 1:size(avgdata_lab, 1)
                p_matches = scatter(avgdata_lab(sc, ic, 2), avgdata_lab(sc, ic, 3), [], squeeze(avgdata_rgb(sc, ic, :).^(1/2.2))', 'filled');
            end

            p_match = plot(popavgdata_lab(ic, 2), popavgdata_lab(ic, 3), 'rd', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', mean_match_rgb, 'MarkerSize', 15);

            h_leg = legend([p_obj_pts, p_obj, p_high, p_matches, p_match], {'Object', 'Mean Object', 'Mean Highlight', 'Obs. matches', 'Grand mean matches'}, 'Location', 'NorthEastOutside');

            subplot(1,2,2);
            imshow(squeeze(imgs(ic, :, :, :)))

            fnparts = strsplit(ims2{ic}, '_');
            if strcmp(fnparts{end}, 'lightest.png')
               transmission = [fnparts{9} '_less_sat'];
            else
               transmission = [fnparts{9} '_more_sat'];
            end
            export_fig(['~/Dropbox/transparent_color_skin_of_things_nov_2018/' which_exp '/plots/achromatic_transmission_' transmission '_obj_high_obs_matches_grand_mean.pdf'], '-pdf');
        end
    end

    % now make big plots showing distributions of a+b and L+a, as well as observer settings in the
    % cloud

    figure(12);
    for ic = 1:size(imgs, 1)
        % now make plots for each image
        if do_big_plot
            bnd_idxs_a_b = boundary(squeeze(obj_lab(ic, :, 2:3)));
            bnd_idxs_l_a = boundary(squeeze(obj_lab(ic, :, 1:2)));

            mean_match_rgb = popavgdata_rgb(ic, :).^(1/2.2);

            clf;
            subplot(2,2,1);
            hold on;
            xlabel('a*');
            ylabel('b*');
            axis([-150 150 -150 150]);
            axis square;
            grid;
            p_obj_pts = scatter(obj_lab(ic, :, 2), obj_lab(ic, :, 3), [], squeeze(obj_rgb(ic, :, :).^(1/2.2)), 'filled');
            p_obj = plot(mean_obj_lab(ic, 2), mean_obj_lab(ic, 3), 's', 'MarkerFaceColor', mean_obj_rgb(ic, :).^(1/2.2), 'MarkerSize', 15);

            h_leg = legend([p_obj_pts, p_obj], {'Object', 'Mean Object'}, 'Location', 'NorthEastOutside');

            subplot(2,2,2);
            hold on;
            xlabel('a*');
            ylabel('L*');
            axis([-150 150 0 100]);
            axis square;
            grid;
            p_obj_pts = scatter(obj_lab(ic, :, 2), obj_lab(ic, :, 1), [], squeeze(obj_rgb(ic, :, :).^(1/2.2)), 'filled');
            p_obj = plot(mean_obj_lab(ic, 2), mean_obj_lab(ic, 1), 's', 'MarkerFaceColor', mean_obj_rgb(ic, :).^(1/2.2), 'MarkerSize', 15);

            h_leg = legend([p_obj_pts, p_obj], {'Object', 'Mean Object'}, 'Location', 'NorthEastOutside');

            subplot(2,2,3);
            hold on;
            xlabel('a*');
            ylabel('b*');
            axis([-150 150 -150 150]);
            axis square;
            grid;
            p_obj_pts = plot(obj_lab(ic, bnd_idxs_a_b, 2), obj_lab(ic, bnd_idxs_a_b, 3), 'k-');
            p_obj = plot(mean_obj_lab(ic, 2), mean_obj_lab(ic, 3), 's', 'MarkerFaceColor', mean_obj_rgb(ic, :).^(1/2.2), 'MarkerSize', 15);
            p_high = plot(mean_high_lab(ic, 2), mean_high_lab(ic, 3), 'o', 'MarkerFaceColor', mean_high_rgb(ic, :).^(1/2.2), 'MarkerSize', 15);

            for sc = 1:size(avgdata_lab, 1)
                p_matches = scatter(avgdata_lab(sc, ic, 2), avgdata_lab(sc, ic, 3), [], squeeze(avgdata_rgb(sc, ic, :).^(1/2.2))', 'filled');
            end

            p_match = plot(popavgdata_lab(ic, 2), popavgdata_lab(ic, 3), 'rd', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', mean_match_rgb, 'MarkerSize', 15);

            h_leg = legend([p_obj_pts, p_obj, p_high, p_matches, p_match], {'Object', 'Mean Object', 'Mean Highlight', 'Obs. matches', 'Grand mean'}, 'Location', 'NorthEastOutside');

            subplot(2,2,4);
            hold on;
            xlabel('a*');
            ylabel('L*');
            axis([-150 150 0 100]);
            axis square;
            grid;
            p_obj_pts = plot(obj_lab(ic, bnd_idxs_l_a, 2), obj_lab(ic, bnd_idxs_l_a, 1), 'k-');
            p_obj = plot(mean_obj_lab(ic, 2), mean_obj_lab(ic, 1), 's', 'MarkerFaceColor', mean_obj_rgb(ic, :).^(1/2.2), 'MarkerSize', 15);
            p_high = plot(mean_high_lab(ic, 2), mean_high_lab(ic, 1), 'o', 'MarkerFaceColor', mean_high_rgb(ic, :).^(1/2.2), 'MarkerSize', 15);

            for sc = 1:size(avgdata_lab, 1)
                p_matches = scatter(avgdata_lab(sc, ic, 2), avgdata_lab(sc, ic, 1), [], squeeze(avgdata_rgb(sc, ic, :).^(1/2.2))', 'filled');
            end

            p_match = plot(popavgdata_lab(ic, 2), popavgdata_lab(ic, 1), 'rd', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', mean_match_rgb, 'MarkerSize', 15);

            h_leg = legend([p_obj_pts, p_obj, p_high, p_matches, p_match], {'Object', 'Mean Object', 'Mean Highlight', 'Obs. matches', 'Grand mean'}, 'Location', 'NorthEastOutside');

            fnparts = strsplit(ims2{ic}, '_');
            if strcmp(fnparts{end}, 'lightest.png')
               transmission = [fnparts{9} '_less_sat'];
            else
               transmission = [fnparts{9} '_more_sat'];
            end
            export_fig(['~/Dropbox/transparent_color_skin_of_things_nov_2018/' which_exp '/plots/achromatic_transmission_' transmission '_obj_high_obs_match.pdf'], '-pdf');
        end
    end

    %%%% save examples of the masks

    imtom = squeeze(imgs(15,:,:,:));
    imwrite(imtom, ['~/Dropbox/transparent_color_skin_of_things_nov_2018/' which_exp '/images/image_to_mask.png']);

    % show highlight
    immr = imtom(:,:,1); immg = imtom(:,:,2); immb = imtom(:,:,3);
    immr(~bm) = 0.5; immg(~bm) = 0.5; immb(~bm) = 0.5;
    imtos = zeros(size(imtom));
    imtos(:,:,1) = immr; imtos(:,:,2) = immg; imtos(:,:,3) = immb;
    figure(10);
    clf;
    hold on;
    imshow(imtos);
    imwrite(imtos, ['~/Dropbox/transparent_color_skin_of_things_nov_2018/' which_exp '/images/image_masked_high.png']);

    % show object
    immr = imtom(:,:,1); immg = imtom(:,:,2); immb = imtom(:,:,3);
    immr(~pm) = 0.5; immg(~pm) = 0.5; immb(~pm) = 0.5;
    imtos = zeros(size(imtom));
    imtos(:,:,1) = immr; imtos(:,:,2) = immg; imtos(:,:,3) = immb;
    figure(10);
    clf;
    hold on;
    imshow(imtos);
    imwrite(imtos, ['~/Dropbox/transparent_color_skin_of_things_nov_2018/' which_exp '/images/image_masked_obj.png']);

    % show refls
    immr = imtom(:,:,1); immg = imtom(:,:,2); immb = imtom(:,:,3);
    immr(~rm) = 0.5; immg(~rm) = 0.5; immb(~rm) = 0.5;
    imtos = zeros(size(imtom));
    imtos(:,:,1) = immr; imtos(:,:,2) = immg; imtos(:,:,3) = immb;
    figure(10);
    clf;
    hold on;
    imshow(imtos);
    imwrite(imtos, ['~/Dropbox/transparent_color_skin_of_things_nov_2018/' which_exp '/images/image_masked_refl.png']);

    % plot showing observer consistency in a+b and L+a
    for sc = 1:length(fns)
        figure(30);
        clf;
        subplot(1,2,1);
        hold on;
        xlabel("a*");
        ylabel("b*");
        axis([-150 150 -150 150]);
        axis square;
        grid;
        plot(avgdata_lab(sc, :, 2), avgdata_lab(sc, :, 3), 'ko');

        subplot(1,2,2);
        hold on;
        xlabel("b*");
        ylabel("L*");
        axis([-150 150 -150 150]);
        axis square;
        grid;
        plot(avgdata_lab(sc, :, 3), avgdata_lab(sc, :, 1), 'ko');

        export_fig(['~/Dropbox/transparent_color_skin_of_things_nov_2018/' which_exp '/plots/per_obs/' num2str(sc) '_mean_match_obj_correlation.pdf'], '-pdf');

        figure(40);
        clf;
        subplot(1,2,1);
        hold on;
        xlabel("L* obj");
        ylabel("L* obs");
        axis([0 100 0 100]);
        axis square;
        grid;
        plot(mean_obj_lab(:, 1), avgdata_lab(sc, :, 1), 'ko');

        subplot(1,2,2);
        hold on;
        xlabel("a* obj");
        ylabel("a* obs");
        % axis([-150 150 -150 150]);
        axis square;
        grid;
        plot(mean_obj_lab(:, 2), avgdata_lab(sc, :, 2), 'ko');

        fparts = strsplit(fns(sc).name, '_');
        export_fig(['~/Dropbox/transparent_color_skin_of_things_nov_2018/' which_exp '/plots/per_obs/' fparts{4} '_mean_match_obj_correlation.pdf'], '-pdf');
    end

    % show some portrayal of the variability of observer settings
    figure(31);
    clf;
    subplot(1,3,1);
    hold on;
    lsds = reshape(squeeze(stddata_lab(:, :, 1)), 16*length(fns), 1);
    xc = 1;
    ic = 1;
    title('L*');
    xlabel('obs per img.');
    ylabel('stand. dev.');
    for x = 1:16
        bar(xc:xc+length(fns)-1, lsds(ic:ic+length(fns)-1));
        xc = xc + length(fns)+2;
        ic = ic + length(fns);
    end

    subplot(1,3,2);
    hold on;
    asds = reshape(squeeze(stddata_lab(:, :, 2)), 16*length(fns), 1);
    xc = 1;
    ic = 1;
    title('a*');
    xlabel('obs per img.');
    ylabel('stand. dev.');
    for x = 1:16
        bar(xc:xc+length(fns)-1, asds(ic:ic+length(fns)-1));
        xc = xc + length(fns)+2;
        ic = ic + length(fns);
    end

    subplot(1,3,3);
    hold on;
    bsds = reshape(squeeze(stddata_lab(:, :, 3)), 16*length(fns), 1);
    xc = 1;
    ic = 1;
    title('b*');
    xlabel('obs per img.');
    ylabel('stand. dev.');
    for x = 1:16
        bar(xc:xc+length(fns)-1, bsds(ic:ic+length(fns)-1));
        xc = xc + length(fns)+2;
        ic = ic + length(fns);
    end

    export_fig(['~/Dropbox/transparent_color_skin_of_things_nov_2018/' which_exp '/plots/overall_variability.pdf'], '-pdf')
end
