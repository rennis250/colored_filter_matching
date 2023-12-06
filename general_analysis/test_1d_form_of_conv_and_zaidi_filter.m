clear all; close all; clc;

x = 0.75;
y = rand(1, 10);
z1 = (x^2.2).*y;
other_col = 0.8;
z2 = x*y + (1 - x)*other_col;
figure;
hold on;
plot(z1, 'k');
plot(z2, 'r');

clear the_orig_data;
global the_orig_data;
the_orig_data = y;

n = size(the_orig_data, 1);

clear the_filter_data;
global the_filter_data;
the_filter_data = z1;
x0 = [1.0, 0.0];
x1 = fminsearch(@(x) oned_trans_to_min(x), x0, optimset('MaxFunEvals', 5000));
Bfit1 = x1(1)*y + x1(2);
err = Bfit1 - the_filter_data;
err = err .* err;
err = sum(err(:));
rmse_x1 = sqrt(err/n);
plot(Bfit1, 'kx');

err = the_orig_data - the_filter_data;
err = err .* err;
err = sum(err(:));
rmse_identity_x1 = sqrt(err/n);

clear the_filter_data;
global the_filter_data;
the_filter_data = z2;
x0 = [1.0, 0.0];
x2 = fminsearch(@(x) oned_trans_to_min(x), x0, optimset('MaxFunEvals', 5000));
Bfit2 = x2(1)*y + x2(2);
err = Bfit2 - the_filter_data;
err = err .* err;
err = sum(err(:));
rmse_x2 = sqrt(err/n);
plot(Bfit2, 'rx');

err = the_orig_data - the_filter_data;
err = err .* err;
err = sum(err(:));
rmse_identity_x2 = sqrt(err/n);

(rmse_identity_x1 - rmse_x1)/rmse_identity_x1
(rmse_identity_x2 - rmse_x2)/rmse_identity_x2