clear all; close all; clc;

wlns = 380:800;

bluetrans = dlmread("munsell_blue_extrem.spd");
bluetranstmp = [];
bluetranstmp(:,1) = wlns(1:14:end);
bluetranstmp(:,2) = interp1(bluetrans(:,1), bluetrans(:,2), wlns(1:14:end));
bluetrans = bluetranstmp;

greentrans = dlmread("munsell_green_extrem.spd");
greentranstmp = [];
greentranstmp(:,1) = wlns(1:14:end);
greentranstmp(:,2) = interp1(greentrans(:,1), greentrans(:,2), wlns(1:14:end));
greentrans = greentranstmp;

yellowtrans = dlmread("munsell_yellow_extrem.spd");
yellowtranstmp = [];
yellowtranstmp(:,1) = wlns(1:14:end);
yellowtranstmp(:,2) = interp1(yellowtrans(:,1), yellowtrans(:,2), wlns(1:14:end));
yellowtrans = yellowtranstmp;

redtrans = dlmread("munsell_red_extrem.spd");
redtranstmp = [];
redtranstmp(:,1) = wlns(1:14:end);
redtranstmp(:,2) = interp1(redtrans(:,1), redtrans(:,2), wlns(1:14:end));
redtrans = redtranstmp;

for ld = 0:1/5:1
    for rg_mix = 0:1/10:1
        for by_mix = 0:1/10:1
            mixed_spectrum = ld.*(((rg_mix.*redtrans(:,2) + (1 - rg_mix).*greentrans(:,2)) + (by_mix.*bluetrans(:,2) + (1 - by_mix).*yellowtrans(:,2))));
            dlmwrite(['ld_' num2str(ld) '_rg_' num2str(rg_mix) '_by_' num2str(by_mix) '.spd'], [wlns(1:14:end)' mixed_spectrum(:)], 'delimiter', ' ');
        end
    end
end
