clear all; close all; clc;

gcRGB = csvread('../calibration/eizo_gamma.csv');

monxyY = csvread('../calibration/eizo_chroma.csv');
monxyz = xyY2XYZ(monxyY);

DKL2RGB = dkl2rgbFromCalib('../calibration/eizo_chroma.csv');
RGB2DKL_T = inv(DKL2RGB);

blue_img_wout_filter = gammaCorr(im2double(imread('../images/empty_scene_blue_illum.png')), gcRGB);
yellow_img_wout_filter = gammaCorr(im2double(imread('../images/empty_scene_yellow_illum.png')), gcRGB);

imgs = dir('../images/patterned/*.png');

masks = {'../images/masks/obj_mask.png'};

for exclude_darkpxs = [false, true]
    for ic = 1:length(imgs)
        tic;
        disp(['image num: ', num2str(ic)]);
        
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
        
        if ~strcmp(illum, 'blue') && ~strcmp(illum, 'yellow')
            continue;
        end
        
        if strcmp(illum, 'blue')
            img_wout_filter = blue_img_wout_filter;
        else
            img_wout_filter = yellow_img_wout_filter;
        end
        
        for mc = 1:length(masks)
            mask = imread(masks{mc}); mask = logical(mask);
            fparts = strsplit(masks{mc}, '/');
            fparts = strsplit(fparts{end}, '.');
            mask_name = fparts{1};
            
            %%%%%%%%%%%%
            % image with filter (convert to LMS, DKL, and LAB)
            
            r = img_with_filter(:,:,1); r = r(:); r_filter = r(mask(:));
            g = img_with_filter(:,:,2); g = g(:); g_filter = g(mask(:));
            b = img_with_filter(:,:,3); b = b(:); b_filter = b(mask(:));
            
            darkpxs = [];
            if exclude_darkpxs
                xyz_filter = real(rgb2xyz([r_filter(:), g_filter(:), b_filter(:)]', monxyz))';
                darkpxs = xyz_filter(:, 2) < dark_thresh;
                r_filter = r_filter(~darkpxs);
                g_filter = g_filter(~darkpxs);
                b_filter = b_filter(~darkpxs);
            end
            
            lab_filter = real(rgb2lab([r_filter(:), g_filter(:), b_filter(:)], monxyz));
            dkl_filter = real(rgb2dkl(RGB2DKL_T, [r_filter(:), g_filter(:), b_filter(:)]'))';
            
            %%%%%%%%%%%%
            % image without glass filter (convert to LMS, DKL, and LAB)
            
            rwof = img_wout_filter(:,:,1); rwof = rwof(:); r_wout_filter = rwof(mask(:));
            gwof = img_wout_filter(:,:,2); gwof = gwof(:); g_wout_filter = gwof(mask(:));
            bwof = img_wout_filter(:,:,3); bwof = bwof(:); b_wout_filter = bwof(mask(:));
            
            if exclude_darkpxs
                darkpxs = xyz_filter(:, 2) < dark_thresh;
                r_wout_filter = r_wout_filter(~darkpxs);
                g_wout_filter = g_wout_filter(~darkpxs);
                b_wout_filter = b_wout_filter(~darkpxs);
            end
            
            lab_wout_filter = real(rgb2lab([r_wout_filter(:), g_wout_filter(:), b_wout_filter(:)], monxyz));
            dkl_wout_filter = real(rgb2dkl(RGB2DKL_T, [r_wout_filter(:), g_wout_filter(:), b_wout_filter(:)]'))';
            
            %%%%%%%%%%%%
            % now let's try divergence for LAB data

            div = divergence_map(lab_wout_filter, lab_filter);
            save(['./output/divergence/LAB_divergence_vals_illum_' illum '_body_' body_color '_' mask_name '_' light_or_dark '.mat'], 'div');
            
            %%%%%%%%%%%%
            % now let's try divergence for DKL data
            
            div = divergence_map(dkl_wout_filter, dkl_filter);
            save(['./output/divergence/DKL_divergence_vals_illum_' illum '_body_' body_color '_' mask_name '_' light_or_dark '.mat'], 'div');
            
            toc;
        end
    end
end
