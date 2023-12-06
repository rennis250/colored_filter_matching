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

% close enough to non-colored, completely transparent filter
rg_mid = 0.611;
by_mid = 0.588;
ld_mid = 1.319;

ld_max = 3;

% it goes:
% dark, max red, max green, max blue, max yellow, max white
extreme_lds = [0,      ld_mid,  ld_mid, ld_mid,  ld_mid, ld_max];
extreme_rgs = [rg_mid, 1,       -1,     rg_mid,  rg_mid, rg_mid];
extreme_bys = [by_mid, by_mid,  by_mid, 1,       -1,     by_mid];

% use this and the ld loop below to find the max ld possible,
% which would be a filter of pure max white at all points
% rg_mix = rg_mid;
% by_mix = by_mid;
% for ld = 0:0.2:10
% ld = ld_max;

rmc_extremes = [];
rsd_extremes = [];

for extc = 1:length(extreme_rgs)
    % for extc = 4
    ld     = extreme_lds(extc);
    rg_mix = extreme_rgs(extc);
    by_mix = extreme_bys(extc);
    
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
        extc
        disp('this one goes out of bounds');
    end
    
    % simulate the monitor's action of clamping RGB values to [0, 1]
    rgbs_sent_to_monitor(rgbs_sent_to_monitor > 1) = 1;
    rgbs_sent_to_monitor(rgbs_sent_to_monitor < 0) = 0;
    
    clear rgbs_gc_lin;
    rgbs_gc_lin(:, 1) = rgbs_sent_to_monitor(:, 1).^gR;
    rgbs_gc_lin(:, 2) = rgbs_sent_to_monitor(:, 2).^gG;
    rgbs_gc_lin(:, 3) = rgbs_sent_to_monitor(:, 3).^gB;
    
    %     only useful when searching for max ld value. already done.
    %     rgbs_in_filter = rgbs_gc_lin(filter_idxs, :);
    %     if all(rgbs_in_filter(:) >= 1) || all(rgbs_in_filter(:) < 0)
    %         disp('found it!');
    %         ld
    %         break;
    %     end
    
    lms = real(rgb2lms(RGB2LMS, rgbs_gc_lin'))';
    [rmc, rsd] = rmc_rsd(lms, ~filter_idxs, filter_idxs);
    
    rmc_extremes = [rmc_extremes; rmc];
    rsd_extremes = [rsd_extremes; rsd];
end

% figure;
% imshow(reshape(vor_rgb, 256, 256, 3))
% figure;
% imshow(reshape(lms, 256, 256, 3))

% it goes:
% dark, max red, max green, max blue, max yellow, max white
% colors = {'k', 'r', 'g', 'b', 'y', 'k'};
% symbols = {'o', 'o', 'o', 'o', 'o', 'v'};

% figure(1);
% clf;
% hold on;
% for extc = 1:length(extreme_lds)
%     plot3(rsd_extremes(extc, 1), rsd_extremes(extc, 2), rsd_extremes(extc, 3), [colors{extc} symbols{extc}]);
% end
% % line from red to green
% line([rsd_extremes(2, 1), rsd_extremes(3, 1)], [rsd_extremes(2, 2), rsd_extremes(3, 2)], [rsd_extremes(2, 3), rsd_extremes(3, 3)]);
% % line from blue to yellow
% line([rsd_extremes(4, 1), rsd_extremes(5, 1)], [rsd_extremes(4, 2), rsd_extremes(5, 2)], [rsd_extremes(4, 3), rsd_extremes(5, 3)]);
% % line from black to white
% line([rsd_extremes(1, 1), rsd_extremes(6, 1)], [rsd_extremes(1, 2), rsd_extremes(6, 2)], [rsd_extremes(1, 3), rsd_extremes(6, 3)]);

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
            [rmc, rsd] = rmc_rsd(lms, ~filter_idxs, filter_idxs);
            
            ld_rg_bys = [ld_rg_bys; ld, rg_mix, by_mix];
            rmc_sphere = [rmc_sphere; rmc];
            rsd_sphere = [rsd_sphere; rsd];
        end
    end
end

%% okay, now let's see what happens, when we find the filters that are most similar to 
% the rmc_rsd_comp glavens in the paper

rmc_rsd_comp_ins = dir('../images/exp5/*png');
rmc_rsd_comp_images = zeros(9, 400, 400, 3);
for x = 1:10    
    rmc_rsd_comp_images = im2double(imread(['../images/exp5/' rmc_rsd_comp_ins(x).name]));
end

good_idxs = [1, 2, 3, 4, 6, 7, 8, 9, 10];

glaven_rmcs = [0.479, 0.587, 0.299; ... % 1
               0.563, 0.692, 0.314; ... % 2
               0.184, 0.37, 0.395; ... % 3
               0.217, 0.442, 0.447; ... % 4
               0.107, 0.133, 0.076; ... % 6
               0.628, 0.553, 0.214; ... % 7
               0.463, 0.392, 0.195; ... % 8
               0.15, 0.176, 0.233; ... % 9
               0.147, 0.174, 0.223]; % 10

glaven_rsds = [0.802, 0.981, 0.527; ... % 1
               0.772, 0.965, 0.454; ... % 2
               0.412, 0.909, 0.657; ... % 3
               0.444, 0.943, 0.652; ... % 4
               0.351, 0.403, 0.161; ... % 6
               2.944, 2.685, 0.858; ... % 7
               0.743, 0.629, 0.235; ... % 8
               0.484, 0.563, 0.908; ... % 9
               0.387, 0.451, 0.72]; % 10

figure;
for gc = 1:size(glaven_rsds, 1)
    rmc_eqdiff = sqrt(sum((glaven_rmcs(gc, :) - rmc_sphere).^2, 2));
    [~, rmc_idx] = min(rmc_eqdiff);
    
    rsd_eqdiff = sqrt(sum((glaven_rsds(gc, :) - rsd_sphere).^2, 2));
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
    
    export_fig(['../figures/best_rmc_rsd_filter_' good_idxs(gc) '.png']);
    
    pause;
end