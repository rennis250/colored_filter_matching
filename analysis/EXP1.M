clear all; close all; clc;

gR = 2.1102;
gG = 2.1243;
gB = 2.0170;

gR_act = 2.1442;
gG_act = 2.1514;
gB_act = 1.9483;

monxyY = csvread('../calibration/exp1.chroma.csv');
monxyz = xyY2XYZ(monxyY);

DKL2RGB = dkl2rgbFromCalib('../calibration/eizo_chroma.csv');
RGB2DKL_T = inv(DKL2RGB);

RGB2LMS = rgb2lmsFromCalib('../calibration/eizo_mon_spectra.csv');

fns = dir('../data/exp1/*.mat');
ims = dir('../images/exp1/*.png');
stats = importdata(['../images/exp1/exp1_chrom_stats.csv'], ',', 1);
img_names = {stats.textdata{2:2:end}};
stats = stats.data;

fh = fopen('../data/exp1_obs_data.csv', 'w');
fprintf(fh, 'sub_id,trial,img,illum,body,highlow,cone_l,cone_m,cone_s,ld,rg,yv,l,a,b,lrc,redness,wpld,wprg,wpyv,gwld,gwrg,gwyv,gwl,gwa,gwb,hv,mean_cone_l,mean_cone_m,mean_cone_s,sd_cone_l,sd_cone_m,sd_cone_s,wpL,wpA,wpB\n');

total_data = 1;
skipped_data = 1;

for fc = 1:length(fns)
    fn = fns(fc).name;
    load(['../data/exp1/' fn], 'data', 'subID', 'trialOrder');

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
        for x = idxs
            rgb = squeeze(data(x, 1, :))';
            rgb_sent_to_monitor = real(rgb.^[gR gG gB]);
            rgb_gc_lin = real(rgb_sent_to_monitor.^[gR_act gG_act gB_act]);
            
            dkl = real(rgb2dkl(RGB2DKL_T, rgb_gc_lin'));
            lms = real(rgb2lms(rgb_gc_lin, RGB2LMS));
            lab = real(rgb2lab(rgb_gc_lin, monxyz));

            total_data = total_data + 1;

            if any(rgb_sent_to_monitor > 1) || any(rgb_sent_to_monitor < 0)
                subID
                squeeze(data(x, 1, :))
                skipped_data = skipped_data + 1;
                continue
            end

            fprintf(fh, [subID ',' ...
            num2str(c) ',' ...
            num2str(tc) ',' ...
            illum ',' ...
            body ',' ...
            num2str(hilo) ',' ...
            num2str(lms(1)) ',' ...
            num2str(lms(2)) ',' ...
            num2str(lms(3)) ',' ...
            num2str(dkl(1)) ',' ...
            num2str(dkl(2)) ',' ...
            num2str(dkl(3)) ',' ...
            num2str(lab(1)) ',' ...
            num2str(lab(2)) ',' ...
            num2str(lab(3)) ',' ...
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
            num2str(stats(img_idx,27-2)) ',' ...  % [1,27] = wpL
            num2str(stats(img_idx,28-2)) ',' ...  % [1,28] = wpA
            num2str(stats(img_idx,29-2)) '\n']);  % [1,29] = wpB

            c = c + 1;
        end
    end
end

fclose(fh);

disp('exp 1');
skipped_data
total_data
skipped_data/total_data


