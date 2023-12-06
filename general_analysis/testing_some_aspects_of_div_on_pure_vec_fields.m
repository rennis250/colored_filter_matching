clear all; close all; clc;

[xs, ys] = meshgrid(-2.5:.25:2.5, -2.5:.25:2.5);
ls = sqrt(xs.*xs + ys.*ys);
vec_xs = -ls.*xs;
vec_ys = -ls.*ys;

div = divergence(xs, ys, vec_xs, vec_ys);
figure;
pcolor(xs, ys, div);
shading interp
axis square;
colorposneg
colorbar
hold on
quiver(xs, ys, vec_xs, vec_ys, 'Color', 'black');
export_fig pure_converging_vec_field.tiff

vec_xs(10:15, 10:15) = vec_xs(10:15, 10:15) - 0.8;
vec_ys(10:15, 10:15) = vec_ys(10:15, 10:15) + 0.8;

div = divergence(xs, ys, vec_xs, vec_ys);
figure;
pcolor(xs, ys, div);
shading interp
axis square;
colorposneg
colorbar
hold on
quiver(xs, ys, vec_xs, vec_ys, 'Color', 'black');
export_fig slightly_perturbed_converging_vec_field.tiff