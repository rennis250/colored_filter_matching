clear all; close all; clc;

imgs = dir('*.png');
for ic = 1:length(imgs)
    img = im2double(imread(imgs(ic).name));
    imggc = img.^(2.2);
    if strcmp(imgs(ic).name, 'voronoi_diff_dist.png')
        imggc14 = imggc.*9.0; % since it is such a dark image
    else
        imggc14 = imggc.*1.5;
    end
    img14 = imggc14.^(1.0/2.2);
    imwrite(img14, [imgs(ic).name(1:end-4) '_brighter.png']);
end
