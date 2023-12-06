clear all; close all; clc;

yellowspectra = dlmread('/Users/rje/my_docs/simpler_pathmarcher/spectra/de_brighter_illum_yellow.spd');

[~,yellowsx] = colormatchPlanck(yellowspectra(:,2), yellowspectra(:,1));
yellowsx = yellowsx./3000;

monxyY = [.6804 .3073 30.94
    .2029 .6968 74.22
    .1527 .0508 6.74];
monxyz = xyY2XYZ(monxyY);

X = monxyz(:, 1);
Y = monxyz(:, 2);
Z = monxyz(:, 3);

Xs = yellowsx(1);
Ys = yellowsx(2);
Zs = yellowsx(3);

illumXYZ = [sum(X) sum(Y) sum(Z)];
illumX = illumXYZ(1); illumY = illumXYZ(2); illumZ = illumXYZ(3);

labnys = Ys./illumY; labnxs = Xs./illumX; labnzs = Zs./illumZ;

labxcubis = find(labnxs > (6/29)^3); labxlinis = find(labnxs <= (6/29)^3);
labycubis = find(labnys > (6/29)^3); labylinis = find(labnys <= (6/29)^3);
labzcubis = find(labnzs > (6/29)^3); labzlinis = find(labnzs <= (6/29)^3);

% first we do LAB, following the formulae on wikipedia
linmt = (1/3)*(29/6)^2; linat = (4/29);
labnxs(labxlinis) = linmt.*(labnxs(labxlinis)) + linat;
labnxs(labxcubis) = labnxs(labxcubis).^(1/3);
labnys(labylinis) = linmt.*(labnys(labylinis)) + linat;
labnys(labycubis) = labnys(labycubis).^(1/3);
labnzs(labzlinis) = linmt.*(labnzs(labzlinis)) + linat;
labnzs(labzcubis) = labnzs(labzcubis).^(1/3);

labls = 116.*labnys - 16; labas = 500.*(labnxs - labnys); labbs = 200.*(labnys - labnzs);

lab = [labls(:) labas(:) labbs(:)];
