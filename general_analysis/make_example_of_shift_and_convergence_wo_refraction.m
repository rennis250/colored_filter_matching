% blc = 2; bc = 1; mc = 1; darker = false; do_lab = false; chroma_only = false;

clear all; close all; clc;

monxyY = csvread('../calibration/eizo_chroma.csv');
monxyz = xyY2XYZ(monxyY);

DKL2RGB = dkl2rgbFromCalib('../calibration/eizo_chroma.csv');
RGB2DKL_T = inv(DKL2RGB);

blue_img_wout_filter = im2double(imread('../images/empty_scene_blue_illum.png')).^(2.2);
yellow_img_wout_filter = im2double(imread('../images/empty_scene_yellow_illum.png')).^(2.2);

blurs = [0.0, 2.0, 4.0];
illum_colors = {'blue', 'yellow'};
body_colors = {'red', 'green', 'blue', 'yellow'};
masks = {'../images/masks/obj_mask.png', '../images/masks/refl_mask.png', '../images/masks/obj_wout_high_mask.png'};

ref_lev = {'no', 'little'};

rmse = struct();

for refc = 1:length(ref_lev)
    for blc = 1:length(blurs)
        blur_lev = blurs(blc);
        
        if blc == 1
            blstr = 'no_blur';
        elseif blc == 2
            blstr = 'mod_blur';
        else
            blstr = 'high_blur';
        end
        
        for ic = 1:length(illum_colors)
            illumc = illum_colors{ic};
            
            if strcmp(illumc, 'blue')
                img_wout_filter = blue_img_wout_filter;
            else
                img_wout_filter = yellow_img_wout_filter;
            end
            
            illumc
            
            for bc = 1:length(body_colors)
                bodyc = body_colors{bc};
                fn = ['./' ref_lev{refc} '_refraction/mitsuba_caustics_' ref_lev{refc} '_refraction_rgb_trans_lightest_illum_' illumc '_refl_munsell_' bodyc '_lightest.png'];
                img_with_filter = im2double(imread(fn)).^(2.2);
                
                %                 dark_fn = ['./' ref_lev{refc} '_refraction/mitsuba_caustics_' ref_lev{refc} '_refraction_rgb_trans_illum_' illumc '_refl_munsell_' bodyc '_lighter.png'];
                %                 img_with_dark_filter = im2double(imread(dark_fn)).^(2.2);
                
                if blc > 1
                    img_with_filter = imgaussfilt(img_with_filter, blur_lev);
                    %                     img_with_dark_filter = imgaussfilt(img_with_dark_filter, blur_lev);
                end
                
                bodyc
                
                for mc = 1:length(masks)
                    mask = imread(masks{mc}); mask = logical(mask);
                    fparts = strsplit(masks{mc}, '/');
                    disp(fparts{end});
                    
                    fparts = strsplit(fparts{end}, '.');
                    mstr = fparts{1};
                    
                    %                     for darker = [true, false]
                    %                         if darker
                    %                             disp('darker');
                    %                             dstr = 'darker';
                    %                             rd = img_with_dark_filter(:,:,1); rd = rd(:); rd_filter = rd(mask(:));
                    %                             gd = img_with_dark_filter(:,:,2); gd = gd(:); gd_filter = gd(mask(:));
                    %                             bd = img_with_dark_filter(:,:,3); bd = bd(:); bd_filter = bd(mask(:));
                    %                             global lab_filter
                    %                             lab_filter = real(rgb2labRob([rd_filter(:) gd_filter(:) bd_filter(:)], monxyz));
                    %                             global dkl_filter;
                    %                             dkl_filter = real(rgb2dkl(RGB2DKL_T, [rd_filter(:) gd_filter(:) bd_filter(:)]'))';
                    % if ~darker && mc == 1 && bc == 3
                    % lab_w_dark_filter = real(rgb2labRob([rd_filter(:) gd_filter(:) bd_filter(:)], monxyz));
                    % end
                    
                    %                             darkerc = 1;
                    %                         else
                    disp('lighter');
                    dstr = 'lighter';
                    r = img_with_filter(:,:,1); r = r(:); r_filter = r(mask(:));
                    g = img_with_filter(:,:,2); g = g(:); g_filter = g(mask(:));
                    b = img_with_filter(:,:,3); b = b(:); b_filter = b(mask(:));
                    global lab_filter
                    lab_filter = real(rgb2labRob([r_filter(:) g_filter(:) b_filter(:)], monxyz));
                    global dkl_filter;
                    dkl_filter = real(rgb2dkl(RGB2DKL_T, [r_filter(:) g_filter(:) b_filter(:)]'))';
                    
                    darkerc = 1;
                    %                         end
                    
                    rwf = img_wout_filter(:,:,1); rwf = rwf(:); r_wout_filter = rwf(mask(:));
                    gwf = img_wout_filter(:,:,2); gwf = gwf(:); g_wout_filter = gwf(mask(:));
                    bwf = img_wout_filter(:,:,3); bwf = bwf(:); b_wout_filter = bwf(mask(:));
                    global lab_wout_filter
                    lab_wout_filter = real(rgb2labRob([r_wout_filter(:) g_wout_filter(:) b_wout_filter(:)], monxyz));
                    global dkl_wout_filter;
                    dkl_wout_filter = real(rgb2dkl(RGB2DKL_T, [r_wout_filter(:) g_wout_filter(:) b_wout_filter(:)]'))';
                    
                    % if ~darker && mc == 1 && bc == 3
                    % csvwrite('../images/patterned/mitsuba_caustics_rgb_trans_lightest_illum_blue_refl_munsell_blue_lightest_lab_coords.csv', [lab_filter lab_wout_filter lab_w_dark_filter]);
                    % end
                    
                    x0 = [1.0 1.0 1.0; 1.0 1.0 1.0; 1.0 1.0 1.0; 0.0 0.0 0.0];
                    for do_lab = [true, false]
                        for chroma_only = [true, false]
                            if do_lab
                                disp('LAB transforms');
                                if chroma_only
                                    disp('Chroma only');
                                    noise = randn(size(lab_wout_filter, 1), 1);
                                    orig_data = [50*ones(size(lab_wout_filter, 1), 1)+noise lab_wout_filter(:,2) lab_wout_filter(:,3)];
                                    filter_data = [50*ones(size(lab_filter, 1), 1)+noise lab_filter(:,2) lab_filter(:,3)];
                                    
                                    cstr = 'chroma_only';
                                    
                                    chromc = 1;
                                else
                                    disp('All colors');
                                    orig_data = lab_wout_filter;
                                    filter_data = lab_filter;
                                    
                                    cstr = 'all_colors';
                                    
                                    chromc = 2;
                                end
                                
                                lstr = 'LAB';
                                
                                n = size(orig_data, 1);
                                
                                x = fminsearch(@(x) eq1_lab(x), x0, optimset('MaxFunEvals', 5000));
                                Bfit = (x(1:3,:)*lab_wout_filter')' + repmat(x(4,:), n, 1);
                                
                                labc = 1;
                            else
                                disp('DKL transforms');
                                if chroma_only
                                    disp('Chroma only');
                                    noise = 0.1*randn(size(dkl_wout_filter, 1), 1);
                                    orig_data = [zeros(size(dkl_wout_filter, 1), 1)+noise dkl_wout_filter(:,2) dkl_wout_filter(:,3)];
                                    filter_data = [zeros(size(dkl_filter, 1), 1)+noise dkl_filter(:,2) dkl_filter(:,3)];
                                    
                                    cstr = 'chroma_only';
                                    
                                    chromc = 1;
                                else
                                    disp('All colors');
                                    orig_data = dkl_wout_filter;
                                    filter_data = dkl_filter;
                                    
                                    cstr = 'all_colors';
                                    
                                    chromc = 2;
                                end
                                
                                lstr = 'DKL';
                                
                                n = size(orig_data, 1);
                                
                                x = fminsearch(@(x) eq1_dkl(x), x0, optimset('MaxFunEvals', 5000));
                                Bfit = (x(1:3,:)*dkl_wout_filter')' + repmat(x(4,:), n, 1);
                                
                                labc = 2;
                            end
                            
                            err = Bfit - filter_data;
                            err = err .* err;
                            err = sum(err(:));
                            rmse.fmin(blc, ic, bc, mc, darkerc, labc, chromc, refc) = sqrt(err/n);
                            
                            [ret_R, ret_t] = rigid_transform_3D(orig_data, filter_data);
                            Bfit = (ret_R*orig_data') + repmat(ret_t, 1, n);
                            Bfit = Bfit';
                            
                            err = Bfit - filter_data;
                            err = err .* err;
                            err = sum(err(:));
                            rmse.rigid(blc, ic, bc, mc, darkerc, labc, chromc, refc) = sqrt(err/n);
                            
                            [A, b] = affinemap(orig_data, filter_data);
                            Bfit = A*orig_data' + b;
                            Bfit = Bfit';
                            
                            err = Bfit - filter_data;
                            err = err .* err;
                            err = sum(err(:));
                            rmse.affinemap(blc, ic, bc, mc, darkerc, labc, chromc, refc) = sqrt(err/n);
                            
                            figure(1);
                            clf;
                            hold on;
                            axis square;
                            plot(filter_data(:,2), filter_data(:,3), 'bo');
                            plot(Bfit(:,2), Bfit(:,3), 'ro');
                            if do_lab
                                axis([-100 100 -100 100]);
                            else
                                axis([-2 2 -2 2]);
                            end
                            grid;
                            
                            fn = ['./output/divergence/linear_map_' bodyc '_' mstr '_' dstr '_' lstr '_' cstr '_' blstr '_.png'];
                            print('-dpng', fn);
                            
                            % now look at divergence
                            
                            vecs = filter_data - orig_data;
                            
                            if do_lab
                                xmin = min(orig_data(:,1)) - 10;
                                ymin = min(orig_data(:,2)) - 10;
                                zmin = min(orig_data(:,3)) - 10;
                                
                                xmax = max(orig_data(:,1)) + 10;
                                ymax = max(orig_data(:,2)) + 10;
                                zmax = max(orig_data(:,3)) + 10;
                            else
                                xmin = min(orig_data(:,1)) - 0.1;
                                ymin = min(orig_data(:,2)) - 0.1;
                                zmin = min(orig_data(:,3)) - 0.1;
                                
                                xmax = max(orig_data(:,1)) + 0.1;
                                ymax = max(orig_data(:,2)) + 0.1;
                                zmax = max(orig_data(:,3)) + 0.1;
                            end
                            
                            if chroma_only
                                [ygr, zgr] = meshgrid(linspace(ymin, ymax, 100), linspace(zmin, zmax, 100));
                                % else
                                % [xgr, ygr, zgr] = meshgrid(linspace(xmin, xmax, 100), linspace(ymin, ymax, 100), linspace(zmin, zmax, 100));
                                % end
                                
                                % try my own "interpolation" approach
                                
                                ys = linspace(ymin, ymax, 100);
                                zs = linspace(zmin, zmax, 100);
                                [ygr, zgr] = meshgrid(ys, zs);
                                ygr = fliplr(ygr);
                                zgr = flipud(zgr);
                                a_interp = zeros(100, 100);
                                b_interp = zeros(100, 100);
                                interp_orig_colors = zeros(1,3);
                                interp_filtered_colors = zeros(1,3);
                                as = squeeze(orig_data(:,2));
                                bs = squeeze(orig_data(:,3));
                                vec_as = squeeze(vecs(:,2));
                                vec_bs = squeeze(vecs(:,3));
                                step_size = 0.0077;
                                cc = 1;
                                for yc = 1:length(ys)
                                    y = ys(yc);
                                    for zc = 1:length(zs)
                                        z = zs(zc);
                                        idxs = (as > y - step_size & as < y + step_size) & (bs > z - step_size & bs < z + step_size);
                                        vec_a = vec_as(idxs);
                                        vec_b = vec_bs(idxs);
                                        if isempty(find(idxs))
                                            a_interp(end - (yc - 1), end - (zc - 1)) = 0;
                                            b_interp(end - (yc - 1), end - (zc - 1)) = 0;
                                        else
                                            interp_orig_colors(cc,2) = y;
                                            interp_filtered_colors(cc,2) = y + mean(vec_a);
                                            a_interp(end - (yc - 1), end - (zc - 1)) = mean(vec_a);
                                            
                                            interp_orig_colors(cc,3) = z;
                                            interp_filtered_colors(cc,3) = z + mean(vec_b);
                                            b_interp(end - (yc - 1), end - (zc - 1)) = mean(vec_b);
                                            
                                            cc = cc + 1;
                                        end
                                    end
                                end
                                a_interp = a_interp.';
                                b_interp = b_interp.';
                                
                                div = divergence(ygr, zgr, a_interp, b_interp);
                                
                                save(['./output/divergence/divergence_vals_' bodyc '_' mstr '_' dstr '_' lstr '_' cstr '_' blstr '_.mat'], 'ygr', 'zgr', 'a_interp', 'b_interp', 'div', 'interp_orig_colors', 'interp_filtered_colors', 'as', 'bs', 'vec_as', 'vec_bs');
                            end
                        end
                    end
                end
            end
        end
    end
end
% end

save('output/divergence/rmse_vals.mat', 'rmse');
