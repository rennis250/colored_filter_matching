clear all; close all; clc;

load('PRIMARIES_04_07_2019_OLED.mat');

monxyY = [data(1).x data(1).y data(1).Y
        data(2).x data(2).y data(2).Y
        data(3).x data(3).y data(3).Y];
csvwrite('oled_chroma.csv', monxyY);

mon_spectra = [data(1).spectralData(:) data(2).spectralData(:) data(3).spectralData(:)];
csvwrite('oled_mon_spectra.csv', mon_spectra);

load('PRIMARIES_04_07_2019_EIZO.mat');

monxyY = [data(1).x data(1).y data(1).Y
        data(2).x data(2).y data(2).Y
        data(3).x data(3).y data(3).Y];
csvwrite('eizo_chroma.csv', monxyY);

mon_spectra = [data(1).spectralData(:) data(2).spectralData(:) data(3).spectralData(:)];
csvwrite('eizo_mon_spectra.csv', mon_spectra);
