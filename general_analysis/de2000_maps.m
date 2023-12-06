clear; close all; clc;

% de2000 maps
DEDist = 20;

expns = {'1', '2', '3', '5', '11a', '11b'};
patch_dye_exp = 2;
eizo_exps = patch_dye_exp;
gamma_fiasco_exps = patch_dye_exp;

avg_bkgd_color = [0.8, 0.8, 0.8];

wout_high_mask = imread('../images/masks/obj_wout_high_mask_correct.png');
wout_high_mask_r = wout_high_mask(:, :, 1);
wout_high_mask_g = wout_high_mask(:, :, 2);
wout_high_mask_b = wout_high_mask(:, :, 3);
wout_high_mask = wout_high_mask_r;
reject_pxs = wout_high_mask_r == 127 & wout_high_mask_g == 127 & wout_high_mask_b == 127;
reject_pxs = logical(reject_pxs(:));
bkgd_pxs = wout_high_mask_r == 0 & wout_high_mask_g == 0 & wout_high_mask_b == 0;
bkgd_pxs = logical(bkgd_pxs(:));
obj_pxs = ~(reject_pxs | bkgd_pxs);
obj_pxs = obj_pxs(:);
wout_high_mask(obj_pxs) = logical(1);
wout_high_mask = wout_high_mask(:);

expn = expns{patch_dye_exp};
imgns = dir(['../images/exp' expn '/*.png']);

gcRGB = csvread('../calibration/eizo_gamma.csv');
monxyY = csvread('../calibration/eizo_chroma.csv');
monxyz = xyY2XYZ(monxyY);

popstats = readtable('../data/popstats_for_R.csv');
stats_patch_dye = popstats(popstats.exp_name == categorical({'patch_dye'}), :);

tic;
for imgc = 1:length(imgns)
    imggc = gammaCorr(im2double(imread(['../images/exp' expn '/' imgns(imgc).name])), gcRGB);
    
    fparts = strsplit(imgns(imgc).name, '_');
    if strcmp(fparts{end}, 'lightest.png')
        hilo = 1;
        illum = fparts{7};
        if strcmp(illum, 'green') || strcmp(illum, 'red')
            continue;
        end
        body = fparts{10};
    else
        hilo = 0;
        illum = fparts{6};
        if strcmp(illum, 'green') || strcmp(illum, 'red')
            continue;
        end
        body = fparts{9};
    end
    
    pop_stats_for_image = stats_patch_dye(stats_patch_dye.body_glaven == categorical({body}) & ...
        stats_patch_dye.illum_glaven == categorical({illum}) & ...
        stats_patch_dye.hilo_glaven == hilo & ...
        stats_patch_dye.mask_name == categorical({'obj_mask'}), :);
    
    anchorLAB = [pop_stats_for_image.mean_mean_gwL_match, ...
        pop_stats_for_image.mean_mean_gwa_match, ...
        pop_stats_for_image.mean_mean_gwb_match];
    
    deMap = zeros(size(imggc, 1)*size(imggc, 2), 3);
    
    r = imggc(:, :, 1); rs = r(:);
    g = imggc(:, :, 2); gs = g(:);
    b = imggc(:, :, 3); bs = b(:);
    
    rs(bkgd_pxs) = 0;
    gs(bkgd_pxs) = 0;
    bs(bkgd_pxs) = 0;
    
    rs(reject_pxs) = 0;
    gs(reject_pxs) = 0;
    bs(reject_pxs) = 0;
    
    lab = real(rgb2labRob(monxyz, [rs, gs, bs]));
    
    goodidxs = deltaE2000(lab, repmat(anchorLAB, size(lab, 1), 1)) < DEDist;
    deMap(goodidxs, :) = [rs(goodidxs), gs(goodidxs), bs(goodidxs)];
    deMap(~goodidxs, :) = rgb2gray([rs(~goodidxs), gs(~goodidxs), bs(~goodidxs)])./2.5;
    deMap(deMap > 1) = 1;
    deMap(deMap < 0) = 0;
    
    deMap(:, 1) = deMap(:, 1).^(1.0/gcRGB(1));
    deMap(:, 2) = deMap(:, 2).^(1.0/gcRGB(2));
    deMap(:, 3) = deMap(:, 3).^(1.0/gcRGB(3));
    
    deMap = reshape(deMap, size(imggc));
    if strcmp(illum, 'yellow')
        illum = 'white';
    end
    imwrite(deMap, ['../figures/de2000Map_body_', body, '_illum_', illum, '_hilo_', num2str(hilo), '.png']);
end
toc;