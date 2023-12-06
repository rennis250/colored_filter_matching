clear all; close all; clc;

gcRGB = csvread('../calibration/eizo_gamma.csv');

monxyY = csvread('../calibration/eizo_chroma.csv');
monxyz = xyY2XYZ(monxyY);

DKL2RGB = dkl2rgbFromCalib('../calibration/eizo_chroma.csv');
RGB2DKL_T = inv(DKL2RGB);

RGB2LMS = rgb2lmsFromCalib_For_Convergence('../calibration/eizo_mon_spectra.csv');

blue_img_wout_filter = gammaCorr(im2double(imread('../images/empty_scene_blue_illum.png')), gcRGB);
yellow_img_wout_filter = gammaCorr(im2double(imread('../images/empty_scene_yellow_illum.png')), gcRGB);

masks = {'../images/masks/obj_mask_correct.png', '../images/masks/refl_mask_correct.png', '../images/masks/obj_wout_high_mask_correct.png'};

% the dark_thresh is always 5% of the max luminance in the Glaven
% (a more strict threshold; it excludes A LOT of pixels in the dark Glavens),
% but we don't want it to be biased by the highlight, so load that mask into
% a separate array for later
wout_high_mask = imread('../images/masks/obj_wout_high_mask_correct.png');
wout_high_mask_r = wout_high_mask(:, :, 1);
wout_high_mask_g = wout_high_mask(:, :, 2);
wout_high_mask_b = wout_high_mask(:, :, 3);
wout_high_mask = wout_high_mask_r;
reject_pxs = wout_high_mask_r == 127 & wout_high_mask_g == 127 & wout_high_mask_b == 127;
reject_pxs = reject_pxs(:);
bkgd_pxs = wout_high_mask_r == 0 & wout_high_mask_g == 0 & wout_high_mask_b == 0;
bkgd_pxs = bkgd_pxs(:);
obj_pxs = ~(reject_pxs | bkgd_pxs);
obj_pxs = obj_pxs(:);
wout_high_mask(obj_pxs) = 255;
wout_high_mask = wout_high_mask(:);

rmse = {};
first_time_through = true;

ref_lev = {'normal', 'little', 'no'};

% figure(1);

tic;
for refc = 1:length(ref_lev)
    if refc == 1 % normal level of refraction that was tested
        imgs = dir('../images/patterned/*.png');
    else
        imgs = dir(['./' ref_lev{refc} '_refraction/*.png']);
    end

    for exclude_darkpxs = [false, true]
        for ic = 1:length(imgs)
            disp(['image num: ', num2str(ic)]);

            if refc == 1 % normal level of refraction that was tested
                img_with_filter = gammaCorr(im2double(imread(['../images/patterned/' imgs(ic).name])), gcRGB);
                fparts = strsplit(imgs(ic).name, '_');
                ld_parts = strsplit(fparts{end}, '.');
                light_or_dark = ld_parts{1};

                if strcmp(light_or_dark, 'lightest')
                    % we got the lightest version
                    body_color = fparts{10};
                    illum = fparts{7};
                    dark_thresh = 2;
                else
                    % we got the darker version
                    body_color = fparts{9};
                    illum = fparts{6};
                    dark_thresh = 1;
                end
            else
                img_with_filter = gammaCorr(im2double(imread(['./' ref_lev{refc} '_refraction/' imgs(ic).name])), gcRGB);
                fparts = strsplit(imgs(ic).name, '_');
                ld_parts = strsplit(fparts{end}, '.');
                light_or_dark = ld_parts{1};

                % did not make darker versions of these to save time
                % we got the lightest version
                body_color = fparts{12};
                illum = fparts{9};
                dark_thresh = 2;
            end

            if ~strcmp(illum, 'blue') && ~strcmp(illum, 'yellow')
                continue;
            end

            if strcmp(illum, 'blue')
                img_wout_filter = blue_img_wout_filter;
            else
                img_wout_filter = yellow_img_wout_filter;
            end

            for mc = 1:length(masks)
                mask = imread(masks{mc});
                mask_r = mask(:, :, 1);
                mask_g = mask(:, :, 2);
                mask_b = mask(:, :, 3);
                mask = mask_r;
                reject_pxs = mask_r == 127 & mask_g == 127 & mask_b == 127;
                reject_pxs = reject_pxs(:);
                bkgd_pxs = mask_r == 0 & mask_g == 0 & mask_b == 0;
                bkgd_pxs = bkgd_pxs(:);
                obj_pxs = ~(reject_pxs | bkgd_pxs);
                obj_pxs = obj_pxs(:);
                mask(obj_pxs) = 255;
                mask = mask(:);

                fparts = strsplit(masks{mc}, '/');
                fparts = strsplit(fparts{end}, '.');
                mask_name = fparts{1};

                %%%%%%%%%%%%
                % image with filter (convert to LMS, DKL, and LAB)

                r = img_with_filter(:,:,1); r = r(:);
                g = img_with_filter(:,:,2); g = g(:);
                b = img_with_filter(:,:,3); b = b(:);

                xyz_filter = real(rgb2xyzRob(monxyz, [r(:), g(:), b(:)]'))';

                darkpxs = [];
                if exclude_darkpxs
                    % dark excluded stats
					dark_thresh = 0.05*max(squeeze(xyz_filter(wout_high_mask, 2)));
                    % dark_thresh = 3;
                    dark_pxs = squeeze(xyz_filter(:, 2)) < dark_thresh;
                    obj_pxs(dark_pxs) = false;
                    bkgd_pxs(dark_pxs) = false;
                end

                r_filter = r(obj_pxs(:));
                g_filter = g(obj_pxs(:));
                b_filter = b(obj_pxs(:));

                lab_filter = real(rgb2labRob(monxyz, [r_filter(:), g_filter(:), b_filter(:)]));
                dkl_filter = real(rgb2dkl(RGB2DKL_T, [r_filter(:), g_filter(:), b_filter(:)]'))';
                lms_filter = real(rgb2lms(RGB2LMS, [r_filter(:), g_filter(:), b_filter(:)]'))';
                dkl_from_lms_filter = real(lms2dkl(lms_filter));

                %%%%%%%%%%%%
                % image without glass filter (convert to LMS, DKL, and LAB)

                rwof = img_wout_filter(:,:,1); rwof = rwof(:);
                gwof = img_wout_filter(:,:,2); gwof = gwof(:);
                bwof = img_wout_filter(:,:,3); bwof = bwof(:);

                r_wout_filter = rwof(obj_pxs(:));
                g_wout_filter = gwof(obj_pxs(:));
                b_wout_filter = bwof(obj_pxs(:));

                lab_wout_filter = real(rgb2labRob(monxyz, [r_wout_filter(:), g_wout_filter(:), b_wout_filter(:)]));
                dkl_wout_filter = real(rgb2dkl(RGB2DKL_T, [r_wout_filter(:), g_wout_filter(:), b_wout_filter(:)]'))';
                lms_wout_filter = real(rgb2lms(RGB2LMS, [r_wout_filter(:), g_wout_filter(:), b_wout_filter(:)]'))';
                dkl_from_lms_wout_filter = real(lms2dkl(lms_wout_filter));

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%       LAB           %%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%
                % make chroma only version of LAB data and compute forms of convergence for it

                lab_chroma_only_orig_data = [zeros(size(lab_wout_filter, 1), 1), lab_wout_filter(:,2), lab_wout_filter(:,3)];
                lab_chroma_only_filter_data = [zeros(size(lab_filter, 1), 1), lab_filter(:,2), lab_filter(:,3)];
                [rmse_fmin, rmse_dzmura, rmse_rigid, rmse_affinemap, rmse_identity] = compute_forms_of_convergence(lab_chroma_only_orig_data, lab_chroma_only_filter_data, true);
                if first_time_through
                    rmse = {rmse_fmin, rmse_dzmura, rmse_rigid, rmse_affinemap, rmse_identity, 'lab', 'chroma', illum, body_color, light_or_dark, mask_name, ref_lev{refc}, exclude_darkpxs};
                    first_time_through = false;
                else
                    rmse = [rmse; {rmse_fmin, rmse_dzmura, rmse_rigid, rmse_affinemap, rmse_identity, 'lab', 'chroma', illum, body_color, light_or_dark, mask_name, ref_lev{refc}, exclude_darkpxs}];
                end

                %%%%%%%%%%%%
                % compute forms of convergence for full LAB data

                [rmse_fmin, rmse_dzmura, rmse_rigid, rmse_affinemap, rmse_identity] = compute_forms_of_convergence(lab_wout_filter, lab_filter, false);
                rmse = [rmse; {rmse_fmin, rmse_dzmura, rmse_rigid, rmse_affinemap, rmse_identity, 'lab', 'full', illum, body_color, light_or_dark, mask_name, ref_lev{refc}, exclude_darkpxs}];

                %%%%%%%%%%%%
                % now let's try divergence for LAB data (uses chroma only in actual divergence_map routine)

                % div = divergence_map(lab_wout_filter, lab_filter);
                % save(['./output/divergence/LAB_divergence_vals_illum_' illum '_body_' body_color '_' mask_name '_' light_or_dark '_ref_' ref_lev{refc}  '.mat'], 'div');

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%       DKL           %%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%
                % make chroma only version of DKL data and compute forms of convergence for it

                dkl_chroma_only_orig_data = [zeros(size(dkl_wout_filter, 1), 1), dkl_wout_filter(:,2), dkl_wout_filter(:,3)];
                dkl_chroma_only_filter_data = [zeros(size(dkl_filter, 1), 1), dkl_filter(:,2), dkl_filter(:,3)];
                [rmse_fmin, rmse_dzmura, rmse_rigid, rmse_affinemap, rmse_identity] = compute_forms_of_convergence(dkl_chroma_only_orig_data, dkl_chroma_only_filter_data, true);
                rmse = [rmse; {rmse_fmin, rmse_dzmura, rmse_rigid, rmse_affinemap, rmse_identity, 'dkl', 'chroma', illum, body_color, light_or_dark, mask_name, ref_lev{refc}, exclude_darkpxs}];

                %%%%%%%%%%%%
                % compute forms of convergence for full DKL data

                [rmse_fmin, rmse_dzmura, rmse_rigid, rmse_affinemap, rmse_identity] = compute_forms_of_convergence(dkl_wout_filter, dkl_filter, false);
                rmse = [rmse; {rmse_fmin, rmse_dzmura, rmse_rigid, rmse_affinemap, rmse_identity, 'dkl', 'full', illum, body_color, light_or_dark, mask_name, ref_lev{refc}, exclude_darkpxs}];

                %%%%%%%%%%%%
                % now let's try divergence for DKL data (uses chroma only in actual divergence_map routine)

                % div = divergence_map(dkl_wout_filter, dkl_filter);
                % save(['./output/divergence/DKL_divergence_vals_illum_' illum '_body_' body_color '_' mask_name '_' light_or_dark '_ref_' ref_lev{refc}  '.mat'], 'div');

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%       DKL from LMS          %%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%
                % make chroma only version of DKL data and compute forms of convergence for it

                dkl_from_lms_chroma_only_orig_data = [zeros(size(dkl_from_lms_wout_filter, 1), 1), dkl_from_lms_wout_filter(:,2), dkl_from_lms_wout_filter(:,3)];
                dkl_from_lms_chroma_only_filter_data = [zeros(size(dkl_from_lms_filter, 1), 1), dkl_from_lms_filter(:,2), dkl_from_lms_filter(:,3)];
                [rmse_fmin, rmse_dzmura, rmse_rigid, rmse_affinemap, rmse_identity] = compute_forms_of_convergence(dkl_from_lms_chroma_only_orig_data, dkl_from_lms_chroma_only_filter_data, true);
                rmse = [rmse; {rmse_fmin, rmse_dzmura, rmse_rigid, rmse_affinemap, rmse_identity, 'dkl_from_lms', 'chroma', illum, body_color, light_or_dark, mask_name, ref_lev{refc}, exclude_darkpxs}];

                %%%%%%%%%%%%
                % compute forms of convergence for full DKL data

                [rmse_fmin, rmse_dzmura, rmse_rigid, rmse_affinemap, rmse_identity] = compute_forms_of_convergence(dkl_from_lms_wout_filter, dkl_from_lms_filter, false);
                rmse = [rmse; {rmse_fmin, rmse_dzmura, rmse_rigid, rmse_affinemap, rmse_identity, 'dkl_from_lms', 'full', illum, body_color, light_or_dark, mask_name, ref_lev{refc}, exclude_darkpxs}];

                %%%%%%%%%%%%
                % now let's try divergence for DKL data (uses chroma only in actual divergence_map routine)

                % div = divergence_map(dkl_from_lms_wout_filter, dkl_from_lms_filter);
                % save(['./output/divergence/DKL_from_LMS_divergence_vals_illum_' illum '_body_' body_color '_' mask_name '_' light_or_dark '_ref_' ref_lev{refc}  '.mat'], 'div');

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%       LMS           %%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % does not make sense to do a "chroma" only version in LMS, without making a pseudo-MB-DKL space...

                %%%%%%%%%%%%
                % compute forms of convergence for full LMS data

                [rmse_fmin, rmse_dzmura, rmse_rigid, rmse_affinemap, rmse_identity] = compute_forms_of_convergence(lms_wout_filter, lms_filter, false);
                rmse = [rmse; {rmse_fmin, rmse_dzmura, rmse_rigid, rmse_affinemap, rmse_identity, 'lms', 'full', illum, body_color, light_or_dark, mask_name, ref_lev{refc}, exclude_darkpxs}];

                % print out an image of the best fit for LMS most general affine map, for comparison with original
                % D'Zmura report

                % global Bfit_affinemap;

                % clf;
                % hold on;
                % axis square;
                % plot3(lms_filter(:, 1), lms_filter(:, 2), lms_filter(:, 3), 'bo');
                % plot3(Bfit_affinemap(:, 1), Bfit_affinemap(:, 2), Bfit_affinemap(:, 3), 'ro');
                % axis([0 1 0 1 0 0.2]);

                % fn = ['./output/convergence/affine_map_illum_' illum '_body_' body_color '_' mask_name '_' light_or_dark '_ref_' ref_lev{refc} '.png'];
                % print('-dpng', fn);
            end
        end
    end
end
toc;

rmse = [{'rmse_fmin', 'rmse_dzmura', 'rmse_rigid', 'rmse_affinemap', 'rmse_identity', 'cspace', 'chroma_or_full', 'illum', 'body_color', 'light_or_dark', 'mask_name', 'ref_lev', 'exclude_darkpxs'}; rmse];
writecell(rmse);
