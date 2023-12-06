clear all; close all; clc;

fns = dir('../data/exp5/*.mat');
ims = dir('../images/exp5/*.png');
stats = importdata(['../data/exp5_img_stats_dark_pxs_excl.csv'], ',', 1);
img_names_re = {stats.textdata{2:2:end}};
img_names_others = {stats.textdata{2:2:end}};

good_idxs = ~strcmp(img_names_others, 'mitsuba_comp_rmc_rsd_illum_4_rg_1_by_2_ld_1_bkgd_voronoi_diff_dist.png');

img_names_others = img_names_others(good_idxs);

img_names = img_names_others; % original 're' data is now completely excluded

stats = stats.data(good_idxs, :);

fh = fopen('../data/exp5_dark_pxs_excl_obs_data.csv', 'w');
fprintf(fh, 'sub_id,trial,img,img_name,illum,rg,by,lighter_darker,bkgd,mean_ld,mean_rg,mean_yv,mean_l,mean_a,mean_b,mean_cone_l,mean_cone_m,mean_cone_s,max_l,max_a,max_b,mean_ratio_l,mean_ratio_m,mean_ratio_s,std_ratio_l,std_ratio_m,std_ratio_s,lrc,redness,wpld,wprg,wpyv,gwld,gwrg,gwyv,gwl,gwa,gwb,hv,mean_cone_l,mean_cone_m,mean_cone_s,sd_cone_l,sd_cone_m,sd_cone_s,rmc_l,rmc_m,rmc_s,rsd_l,rsd_m,rsd_s,wpL,wpA,wpB,ld,rg_mix,by_mix,qualityc\n');

total_data = 1;
skipped_data = 1;

for fc = 1:length(fns)
    fn = fns(fc).name;
    load(['../data/exp5/' fn], 'data', 'subID', 'trialOrder', 'ntrials', 'st');

    if strcmp(subID, 're')
        continue
    end

    if strcmp(subID, '02') || strcmp(subID, '03')
        redtrans = '../base_stimuli/spectra/munsell_red_better.spd';
        greentrans = '../base_stimuli/spectra/munsell_green_better.spd';
        bluetrans = '../base_stimuli/spectra/munsell_blue_better.spd';
        yellowtrans = '../base_stimuli/spectra/munsell_yellow_better.spd';
    else
        redtrans = '../base_stimuli/spectra/munsell_red_EXTREEEMMMEE.spd';
        greentrans = '../base_stimuli/spectra/munsell_green_EXTREEEMMMEE.spd';
        bluetrans = '../base_stimuli/spectra/munsell_blue_EXTREEEMMMEE.spd';
        yellowtrans = '../base_stimuli/spectra/munsell_yellow_EXTREEEMMMEE.spd';
    end

    for tc = 1:length(img_names)
        fparts = strsplit(ims(tc).name, '_');

        % illum 3 = red
        % illum 4 = green
        % ld 1 = darker
        % ld 2 = lighter

        if strcmp(fparts{6}, '3')
            illum = 'red';
        else
            illum = 'green';
        end
        rg = fparts{8};
        by = fparts{10};
        if strcmp(fparts{12}, '1')
            lighter_darker = 'darker';
        else
            lighter_darker = 'lighter';
        end
        bkgd = ims(tc).name(50:end-4);

        img_idx = 1;
        fnparts = strsplit(fn, '_');
	    for img_idx = 1:length(ims)
            if strcmp(ims(tc).name, img_names{img_idx})
                break
          	end
    	end

        idxs = find(trialOrder == tc);
        c = 1;
        disp([length(fns), fc, tc])
        for idx = idxs
            rg_mix = data(idx, 1);
            by_mix = data(idx, 2);
            ld = data(idx, 3);
            qualityc = data(idx, 4);

            total_data = total_data + 1;

            [~, output] = system(['voronoi_filters' ...
            	' -gammaFile=../calibration/exp5.gamma.csv' ...
            	' -spectraFile=../calibration/exp5.spectra.csv' ...
            	' -chromaFile=../calibration/exp5.chroma.csv' ...
            	' -ld=' num2str(ld) ...
            	' -rg=' num2str(rg_mix) ...
            	' -by=' num2str(by_mix) ...
            	' -red=' redtrans ...
            	' -green=' greentrans ...
            	' -blue=' bluetrans ...
            	' -yellow=' yellowtrans]);

            output_parts = strsplit(output, ' ');

            if length(output_parts) > 1 && strcmp(output_parts{3}, 'OOG:')
                subID
                disp([ld rg_mix by_mix]);
                disp('voronoi filter went OOG');
                skipped_data = skipped_data + 1;
                continue
            end

            output_parts = strsplit(output, ',');
            
            mean_ld = output_parts{6};
            mean_rg = output_parts{7};
            mean_yv = output_parts{8};

            mean_l = output_parts{9};
            mean_a = output_parts{10};
            mean_b = output_parts{11};
            
            mean_cone_l = output_parts{13};
            mean_cone_m = output_parts{14};
            mean_cone_s = output_parts{15};

            mean_ratio_l = output_parts{19};
            mean_ratio_m = output_parts{20};
            mean_ratio_s = output_parts{21};

            std_ratio_l = output_parts{22};
            std_ratio_m = output_parts{23};
            std_ratio_s = output_parts{24};

            max_l = output_parts{25};
            max_a = output_parts{26};
            max_b = output_parts{27};

            fprintf(fh, [subID ',' ...
            num2str(c) ',' ...
            num2str(tc) ',' ...
            img_names{img_idx} ',' ...
            illum ',' ...
            num2str(rg) ',' ...
            num2str(by) ',' ...
            lighter_darker ',' ...
            bkgd ',' ...
            mean_ld ',' ...
            mean_rg ',' ...
            mean_yv ',' ...
            mean_l ',' ...
            mean_a ',' ...
            mean_b ',' ...
            mean_cone_l ',' ...
            mean_cone_m ',' ...
            mean_cone_s ',' ...
            max_l ',' ...
            max_a ',' ...
            max_b ',' ...
            mean_ratio_l ',' ...
            mean_ratio_m ',' ...
            mean_ratio_s ',' ...
            std_ratio_l ',' ...
            std_ratio_m ',' ...
            std_ratio_s ',' ...
            num2str(stats(img_idx,3-2)) ',' ...   % [1,3] = lrc
            num2str(stats(img_idx,4-2)) ',' ...   % [1,4] = redness
            num2str(stats(img_idx,5-2)) ',' ...   % [1,5] = wpld
            num2str(stats(img_idx,6-2)) ',' ...   % [1,6] = wprg
            num2str(stats(img_idx,7-2)) ',' ...   % [1,7] = wpyv
            num2str(stats(img_idx,8-2)) ',' ...   % [1,8] = gwld
            num2str(stats(img_idx,9-2)) ',' ...   % [1,9] = gwrg
            num2str(stats(img_idx,10-2)) ',' ...  % [1,10] = gwyv
            num2str(stats(img_idx,11-2)) ',' ...  % [1,11] = gwl
            num2str(stats(img_idx,12-2)) ',' ...  % [1,12] = gwa
            num2str(stats(img_idx,13-2)) ',' ...  % [1,13] = gwb
            num2str(stats(img_idx,14-2)) ',' ...  % [1,14] = hv
            num2str(stats(img_idx,15-2)) ',' ...  % [1,15] = mean_cone_l
            num2str(stats(img_idx,16-2)) ',' ...  % [1,16] = mean_cone_m
            num2str(stats(img_idx,17-2)) ',' ...  % [1,17] = mean_cone_s
            num2str(stats(img_idx,18-2)) ',' ...  % [1,18] = sd_cone_l
            num2str(stats(img_idx,19-2)) ',' ...  % [1,19] = sd_cone_m
            num2str(stats(img_idx,20-2)) ',' ...  % [1,20] = sd_cone_s
            num2str(stats(img_idx,21-2)) ',' ...  % [1,21] = rmc_l
            num2str(stats(img_idx,22-2)) ',' ...  % [1,22] = rmc_m
            num2str(stats(img_idx,23-2)) ',' ...  % [1,23] = rmc_s
            num2str(stats(img_idx,24-2)) ',' ...  % [1,24] = rsd_l
            num2str(stats(img_idx,25-2)) ',' ...  % [1,25] = rsd_m
            num2str(stats(img_idx,26-2)) ',' ...  % [1,26] = rsd_s
            num2str(stats(img_idx,27-2)) ',' ...  % [1,27] = wpL
            num2str(stats(img_idx,28-2)) ',' ...  % [1,28] = wpA
            num2str(stats(img_idx,29-2)) ',' ...  % [1,29] = wpB
            num2str(ld) ',' ... % obs ld setting for this trial
            num2str(rg_mix) ',' ... % obs rg_mix setting for filters for this trial
            num2str(by_mix) ',' ... % obs by_mix setting for filters for this trial
            num2str(qualityc) '\n' ]); % obs rating of how good their match was

            c = c + 1;
        end
    end
end

fclose(fh);

disp('exp 5');
skipped_data
total_data
skipped_data/total_data
