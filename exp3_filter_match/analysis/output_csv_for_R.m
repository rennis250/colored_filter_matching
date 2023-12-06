clear all; close all; clc;

ims = dir('../../images/patterned/mitsuba_caus*.png');
fns = dir('../data/*.mat');
stats = importdata(['../../images/patterned/eizo_chrom_stats.csv'], ',', 1);
img_names = {stats.textdata{2:2:end}};
stats = stats.data;

fh = fopen('../data/all_obs_data_together.csv', 'w');
fprintf(fh, 'sub_id,trial,img,illum,body,highlow,mean_l,mean_a,mean_b,max_l,max_a,max_b,mean_ratio_l,mean_ratio_m,mean_ratio_s,std_ratio_l,std_ratio_m,std_ratio_s,lrc,redness,wpld,wprg,wpyv,gwld,gwrg,gwyv,gwl,gwa,gwb,hv,mean_cone_l,mean_cone_m,mean_cone_s,sd_cone_l,sd_cone_m,sd_cone_s,rmc_l,rmc_m,rmc_s,rsd_l,rsd_m,rsd_s,wpL,wpA,wpB\n');

for fc = 1:length(fns)
    fn = fns(fc).name;
    load(['../data/' fn], 'data', 'subID', 'trialOrder', 'ntrials', 'st');

    if strcmp(subID, 'll')
        continue;
        % by_mult = 3;
        % redtrans = '../../base_stimuli/spectra/refl_redBall.spd';
        % greentrans = '../../base_stimuli/spectra/refl_greenBall.spd';
        % bluetrans = '../../base_stimuli/spectra/refl_blueBall.spd';
        % yellowtrans = '../../base_stimuli/spectra/refl_yellowBall.spd';
    else
        % by_mult = 1;
        redtrans = '../../base_stimuli/spectra/munsell_red_better.spd';
        greentrans = '../../base_stimuli/spectra/munsell_green_better.spd';
        bluetrans = '../../base_stimuli/spectra/munsell_blue_better.spd';
        yellowtrans = '../../base_stimuli/spectra/munsell_yellow_better.spd';
    end

    for tc = 1:size(ims, 1)
        fparts = strsplit(ims(tc).name, '_');
        if strcmp(fparts{end}, 'lighter.png')
            if strcmp(fparts{6}, 'red') || strcmp(fparts{6}, 'green')
                continue;
            end
            illum = fparts{6};
            body = fparts{9};
            hilo = 1;
        else
            if strcmp(fparts{7}, 'red') || strcmp(fparts{7}, 'green')
                continue;
            end
            illum = fparts{7};
            body = fparts{10};
            hilo = 0;
        end

        img_idx = 1;
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

            [~, output] = system(['voronoi_filters' ...
            	' -gR=2.1442' ...
            	' -gG=2.1514' ...
            	' -gB=1.9483' ...
            	' -spectraFile=../../calibration/eizo_mon_spectra.csv' ...
            	' -chromaFile=../../calibration/eizo_chroma.csv' ...
            	' -ld=' num2str(ld) ...
            	' -rg=' num2str(rg_mix) ...
            	' -by=' num2str(by_mix) ...
            	' -red=' redtrans ...
            	' -green=' greentrans ...
            	' -blue=' bluetrans ...
            	' -yellow=' yellowtrans]);

            output_parts = strsplit(output, ',');

            mean_l = output_parts{9};
            mean_a = output_parts{10};
            mean_b = output_parts{11};

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
            illum ',' ...
            body ',' ...
            num2str(hilo) ',' ...
            mean_l ',' ...
            mean_a ',' ...
            mean_b ',' ...
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
            num2str(stats(img_idx,29-2)) '\n']);  % [1,29] = wpB

            c = c + 1;
        end
    end
end

fclose(fh);
