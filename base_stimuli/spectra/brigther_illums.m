clear all; close all; clc;

spds = dir('spectra/illum_*.spd');

for x = 1:length(spds)
	illum = dlmread(spds(x).name);
	illum(:,2) = 0.5.*illum(:,2); % 180 - value used on re
	dlmwrite(['spectra/brighter_' spds(x).name], illum, ' ');
end