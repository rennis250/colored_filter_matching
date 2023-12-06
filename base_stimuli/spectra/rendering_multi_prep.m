clear all; close all; clc;

global B;
global Uhats;
global csm;
global m;

m = 12; % how many reduced dimensions do we want to work with

% CIE1931 XYZ CMFs
cmfData = csvread('ciexyz31_1.csv');
wavelength_cmf = cmfData(:,1);
x_bar = cmfData(:,2);
y_bar = cmfData(:,3);
z_bar = cmfData(:,4);

% redS = dlmread('refl_redBall.spd');
% redS(:,2) = smoothdata(redS(:,2), 'gaussian', 50);
% greenS = dlmread('refl_greenBall.spd');
% greenS(:,2) = smoothdata(greenS(:,2), 'gaussian', 50);
% blueS = dlmread('refl_blueBall.spd');
% blueS(:,2) = smoothdata(blueS(:,2), 'gaussian', 50);
% yellowS = dlmread('refl_yellowBall.spd');
% yellowS(:,2) = smoothdata(yellowS(:,2), 'gaussian', 50);

% S = [redS(20:end-36,2) greenS(20:end-36,2) blueS(20:end-36,2) yellowS(20:end-36,2)]; % columns contain reflectance spectra

redE = dlmread('illum_redOrtho.spd');
greenE = dlmread('illum_greenOrtho.spd');
blueE = dlmread('illum_blueOrtho.spd');
yellowE = dlmread('illum_yellowOrtho.spd');

E = [redE(20:end-36,2) greenE(20:end-36,2) blueE(20:end-36,2) yellowE(20:end-36,2)]; % columns contain illuminant spectra

wlns = redE(20:end-36,1);

% let's use munsell reflectances instead, which will allow more light through as transparency
redS = dlmread('munsell_red_lighter.spd');
redS(:,2) = smoothdata(redS(:,2), 'gaussian', 20);
redS = [wlns interp1(redS(:,1), redS(:,2), wlns, 'spline')];
greenS = dlmread('munsell_green_lighter.spd');
greenS(:,2) = smoothdata(greenS(:,2), 'gaussian', 20);
greenS = [wlns interp1(greenS(:,1), greenS(:,2), wlns, 'spline')];
blueS = dlmread('munsell_blue_lighter.spd');
blueS(:,2) = smoothdata(blueS(:,2), 'gaussian', 20);
blueS = [wlns interp1(blueS(:,1), blueS(:,2), wlns, 'spline')];
yellowS = dlmread('munsell_yellow_lighter.spd');
yellowS(:,2) = smoothdata(yellowS(:,2), 'gaussian', 20);
yellowS = [wlns interp1(yellowS(:,1), yellowS(:,2), wlns, 'spline')];

S = [redS(:,2) greenS(:,2) blueS(:,2) yellowS(:,2)]; % columns contain reflectance spectra

% wlns = redS(:,1);
cmf = interp1(wavelength_cmf(20:422), [x_bar(20:422) y_bar(20:422) z_bar(20:422)], wlns, 'spline');
% cmf = [x_bar(20:422) y_bar(20:422) z_bar(20:422)];
cmf(cmf < 0) = 0;

XYZ_CMFs = [cmf(:,1) cmf(:,2) cmf(:,3)];

% load('fredSmoothed.mat');

% where do all the fred spectra plot to

% cmf_fred = interp1(wavelength_cmf(20:422), [x_bar(20:422) y_bar(20:422) z_bar(20:422)], wlns, 'spline');
% cmf_fred(cmf_fred < 0) = 0;

% XYZ_CMFs = [cmf_fred(:,1) cmf_fred(:,2) cmf_fred(:,3)];
% xyz = XYZ_CMFs' * fredSmooth;
% xyY = [xyz(1,:)./sum(xyz,1); xyz(2,:)./sum(xyz,1); xyz(2,:)];

% figure(1);
% clf;
% hold on;
% axis([0 1 0 1]);
% plot(xyY(1,:), xyY(2,:), 'ko');

% lasso = csvread('xyYcoords.csv');
% plot(lasso(:,1), lasso(:,2), 'k-');

% plank = csvread('planck.csv');
% plot(plank(:,1), plank(:,2), 'r-', 'LineWidth', 2);

% % find green point
% greenidx = find(xyY(1,:) > 0.3794 & xyY(1,:) < 0.3796 & xyY(2,:) > 0.4639 & xyY(2,:) < 0.4641);

% % find red point
% redidx = find(xyY(1,:) > 0.4927 & xyY(1,:) < 0.4929 & xyY(2,:) > 0.3108 & xyY(2,:) < 0.3110);

% % find blue point
% blueidx = find(xyY(1,:) > 0.2465 & xyY(1,:) < 0.2467 & xyY(2,:) > 0.2420 & xyY(2,:) < 0.2422);

% % find yellow point
% yellowidx = find(xyY(1,:) > 0.4342 & xyY(1,:) < 0.4344 & xyY(2,:) > 0.4257 & xyY(2,:) < 0.4259);

% % find gray point
% grayidx = find(xyY(1,:) > 0.3601 & xyY(1,:) < 0.3603 & xyY(2,:) > 0.362 & xyY(2,:) < 0.364);

% % S = [fredSmooth(:, [greenidx redidx(1) blueidx yellowidx]) 0.8*ones(61,1)]; % columns contain reflectance spectra
% S = fredSmooth(:, [greenidx redidx blueidx yellowidx grayidx]); % columns contain reflectance spectra
% E = fredSmooth(:, [greenidx redidx blueidx yellowidx]); % columns contain illuminant spectra

SE = []; % columns will contain all possible lights that result after 1 bounce
SecondBounce = []; % columns will contain all possible lights that result after 2 bounces
Ck = []; % columns will contain all possible lights that result after
         % 0 (i.e., light direct into sensor for illuminant), 1, or 2 bounces
c = 1;
for x = 1:size(S, 2)
    for y = 1:size(E, 2)
        SE(:, c) = S(:, x) .* E(:, y);
        c = c + 1;
    end
end

c = 1;
for x = 1:size(S, 2)
    for y = 1:size(SE, 2)
        SecondBounce(:, c) = S(:, x) .* SE(:, y);
        c = c + 1;
    end
end
Ck = [E SecondBounce SE]; % union of E, SecondBounce, and SE

[Uck, Sck, ~] = svd(Ck); % getting a basis for the spectra
B = Uck(:, 1:m); % reduced basis
c = B' * Ck;
csm = sqrtm(c * c');

[Us, Ss, ~] = svd(S);
Uhats = Us * Ss;

% T0 = [-0.4161    0.0511   -3.7990   -0.0929    1.3437    0.4094    2.9287   -7.8231
%     6.7000   11.3002   -2.2089   -2.1951   -3.1328    3.3331    4.7723  -10.4459
%    -8.3373    0.4055    1.6709   -0.4671  -11.6566    0.1298    1.3663    8.2074
%   -70.9398   20.9355   -0.8393    2.2923    0.7246   36.0565    1.5087  -24.3881
%    44.4044    4.1085   -8.8419    6.8235    9.5184  -17.3218  -10.5716   -1.2480
%    25.6442    8.9401    6.7877   -5.6207   -1.0011   -1.4562   -0.0403   -9.4320
%    24.6296   -2.7981    0.6280   -0.6663   -7.8591  -23.4720   13.6885   -2.9262
%     0.4819   -0.3216   -0.4761    0.3390    5.5735  -13.6074   -7.8340    2.9183];

T0 = randn(m, m);

options = optimset('MaxFunEvals', 10000000, 'MaxIter', 10000000);
T = fminsearch(@(x) f(x), T0, options); % transformation matrix for giving us the sharpened basis

Ti = inv(T);
s = B' * S;
e = B' * E;
esquiggle = Ti * e; % basis to use for lights in rendering
ssquiggle = Ti * s; % basis to use for surfaces in rendering

% OLED XYZ2RGB
monxyy = [0.6751 0.3224 33.32
        0.1941 0.7253 76.97
        0.1415 0.0511 7.71];

rxyz = xyY2XYZ(monxyy(1,:));
gxyz = xyY2XYZ(monxyy(2,:));
bxyz = xyY2XYZ(monxyy(3,:));

M = [rxyz(1) gxyz(1) bxyz(1)
        rxyz(2) gxyz(2) bxyz(2)
        rxyz(3) gxyz(3) bxyz(3)];

XYZ2RGB = inv(M);

M = XYZ2RGB * XYZ_CMFs';

coeff2rgb = M * B * T; % matrix to take vectors from sharpened basis space
                       % to linear RGB where M = Y * X = XYZ2RGB * XYZ_CMFs

coeff2spect = B * T; % matrix to take vectors from sharpened basis space to spectra

coeff2xyz = XYZ_CMFs' * B * T; % matrix to take vectors from sharpened basis space to CIEXYZ

disp('rough approx. difference between illuminant->rgb and (8-dim basis coeffs)->rgb')
abs((M * E) - (coeff2rgb * esquiggle))

disp('rough approx. difference between surfaces->rgb and (8-dim basis coeffs)->rgb')
abs((M * S) - (coeff2rgb * ssquiggle))

% see how interpolating between two illuminants looks in terms of spectra

figure(2);
clf;
hold on;
plot(E(:,1),'r-');
plot(E(:,2),'b-');

for x = 0:0.1:1
  plot(coeff2spect * (x*esquiggle(:,1) + (1 - x)*esquiggle(:,2)), 'k--');
  plot(x*E(:,1) + (1 - x)*E(:,2), 'gx');
  pause;
end

figure(3);
clf;
hold on;
plot(E(:,3),'r-');
plot(E(:,4),'b-');

for x = 0:0.1:1
  plot(coeff2spect * (x*esquiggle(:,3) + (1 - x)*esquiggle(:,4)), 'k--');
  plot(x*E(:,3) + (1 - x)*E(:,4), 'gx');
  pause;
end

% all spectra on same plot

figure(4);
clf;
hold on;
plot(E(:,1),'r-');
plot(E(:,2),'b-');
plot(E(:,3),'g-');
plot(E(:,4),'k-');

plot(coeff2spect * esquiggle(:,1), 'kx');
plot(coeff2spect * esquiggle(:,2), 'kx');
plot(coeff2spect * esquiggle(:,3), 'kx');
plot(coeff2spect * esquiggle(:,4), 'rx');

% see how interpolating between two illuminants looks in terms of CIExyY

Exyz = XYZ_CMFs' * E;
ExyY = zeros(3,4);
ExyY(:, 1) = [Exyz(1,1)/sum(Exyz(:,1)) Exyz(2,1)/sum(Exyz(:,1)) Exyz(2,1)];
ExyY(:, 2) = [Exyz(1,2)/sum(Exyz(:,2)) Exyz(2,2)/sum(Exyz(:,2)) Exyz(2,2)];
ExyY(:, 3) = [Exyz(1,3)/sum(Exyz(:,3)) Exyz(2,3)/sum(Exyz(:,3)) Exyz(2,3)];
ExyY(:, 4) = [Exyz(1,4)/sum(Exyz(:,4)) Exyz(2,4)/sum(Exyz(:,4)) Exyz(2,4)];

figure(5);
clf;
hold on;
plot(ExyY(1,1), ExyY(2,1), 'rx');
plot(ExyY(1,2), ExyY(2,2), 'bx');
axis([0 1 0 1]);

for x = 0:0.1:1
  ispec = coeff2spect * (x*esquiggle(:,1) + (1 - x)*esquiggle(:,2));
  ixyz = XYZ_CMFs' * ispec;
  ixyY = [ixyz(1)/sum(ixyz) ixyz(2)/sum(ixyz()) ixyz(2)];
  plot(ixyY(1), ixyY(2), 'ko');

  ixyz = x*Exyz(:,1) + (1 - x)*Exyz(:,2);
  ixyY = [ixyz(1)/sum(ixyz) ixyz(2)/sum(ixyz()) ixyz(2)];
  plot(ixyY(1), ixyY(2), 'gx');
  pause;
end

plot(ExyY(1,3), ExyY(2,3), 'rx');
plot(ExyY(1,4), ExyY(2,4), 'bx');

for x = 0:0.1:1
  ispec = coeff2spect * (x*esquiggle(:,3) + (1 - x)*esquiggle(:,4));
  ixyz = XYZ_CMFs' * ispec;
  ixyY = [ixyz(1)/sum(ixyz) ixyz(2)/sum(ixyz()) ixyz(2)];
  plot(ixyY(1), ixyY(2), 'ko');

  ixyz = x*Exyz(:,3) + (1 - x)*Exyz(:,4);
  ixyY = [ixyz(1)/sum(ixyz) ixyz(2)/sum(ixyz()) ixyz(2)];
  plot(ixyY(1), ixyY(2), 'gx');
  pause;
end

% does render output make sense

se = []; % columns will contain all possible coeffs that result after 1 bounce
rgb = []; % columns will contain all linear rgb output for se
rgb_gc = []; % columns will contain all gamma-corrected rgb output for se
c = 1;
for x = 1:size(ssquiggle, 2)
    for y = 1:size(esquiggle, 2)
        se(:, c) = ssquiggle(:, x) .* esquiggle(:, y);
        rgb(:, c) = coeff2rgb * se(:, c);
        rgb_gc(:, c) = real(rgb(:, c).^(1.0/2.2));
        c = c + 1;
    end
end

% produce image containing colors of all first bounces

render_out_img = zeros(4*20, 4*20, 3);
c = 1;
for y = 1:20:4*20
  for x = 1:20:4*20
    render_out_img(y:y+20, x:x+20, 1) = 0.5*rgb_gc(1, c);
    render_out_img(y:y+20, x:x+20, 2) = 0.5*rgb_gc(2, c);
    render_out_img(y:y+20, x:x+20, 3) = 0.5*rgb_gc(3, c);
    c = c + 1;
  end
end

% does direct spectral rendering output make sense

se_spec = []; % columns will contain all possible coeffs that result after 1 bounce
rgb_spec = []; % columns will contain all linear rgb_spec output for se_spec
rgb_gc_spec = []; % columns will contain all gamma-corrected rgb_spec output for se_spec
c = 1;
for x = 1:size(S, 2)
    for y = 1:size(E, 2)
        se_spec(:, c) = S(:, x) .* E(:, y);
        rgb_spec(:, c) = M * se_spec(:, c);
        rgb_gc_spec(:, c) = real(rgb_spec(:, c).^(1.0/2.2));
        c = c + 1;
    end
end

% produce image containing colors of all first bounces

render_out_img_spec = zeros(4*20, 4*20, 3);
c = 1;
for y = 1:20:4*20
  for x = 1:20:4*20
    render_out_img_spec(y:y+20, x:x+20, 1) = 0.5*rgb_gc_spec(1, c);
    render_out_img_spec(y:y+20, x:x+20, 2) = 0.5*rgb_gc_spec(2, c);
    render_out_img_spec(y:y+20, x:x+20, 3) = 0.5*rgb_gc_spec(3, c);
    c = c + 1;
  end
end

figure;
clf;
subplot(2,1,1)
imshow(render_out_img)
subplot(2,1,2)
imshow(render_out_img_spec)

% try getting a reduced basis via the technique in section 5.4 of waddle's thesis
% if it is even necessary

% export data for mitsuba comparison

% % refls
% green = [wlns(:) fredSmooth(:,1064)];
% dlmwrite("green.spd", green, 'delimiter', ' ')
% red = [wlns(:) fredSmooth(:,1025)];
% dlmwrite("red.spd", red, 'delimiter', ' ')
% blue = [wlns(:)  fredSmooth(:,903)];
% dlmwrite("blue.spd", blue, 'delimiter', ' ')
% yellow = [wlns(:) fredSmooth(:,1946)];
% dlmwrite("yellow.spd", yellow, 'delimiter', ' ')
% gray = [wlns(:) fredSmooth(:,1946)];
% dlmwrite("gray.spd", gray, 'delimiter', ' ')

% % illums
% green = [wlns(:) 200.*fredSmooth(:,1064)];
% dlmwrite("illum_green.spd", green, 'delimiter', ' ')
% red = [wlns(:) 200.*fredSmooth(:,1025)];
% dlmwrite("illum_red.spd", red, 'delimiter', ' ')
% blue = [wlns(:)  200.*fredSmooth(:,903)];
% dlmwrite("illum_blue.spd", blue, 'delimiter', ' ')
% yellow = [wlns(:) 200.*fredSmooth(:,1946)];
% dlmwrite("illum_yellow.spd", yellow, 'delimiter', ' ')
% gray = [wlns(:) 200.*fredSmooth(:,1946)];
% dlmwrite("illum_gray.spd", gray, 'delimiter', ' ')

% refls
green = [wlns(:) S(:,1)];
dlmwrite("green.spd", green, 'delimiter', ' ')
red = [wlns(:) S(:,2)];
dlmwrite("red.spd", red, 'delimiter', ' ')
blue = [wlns(:)  S(:,3)];
dlmwrite("blue.spd", blue, 'delimiter', ' ')
yellow = [wlns(:) S(:,4)];
dlmwrite("yellow.spd", yellow, 'delimiter', ' ')

% illums
green = [wlns(:) E(:,1)];
dlmwrite("illum_green.spd", green, 'delimiter', ' ')
red = [wlns(:) E(:,2)];
dlmwrite("illum_red.spd", red, 'delimiter', ' ')
blue = [wlns(:)  E(:,3)];
dlmwrite("illum_blue.spd", blue, 'delimiter', ' ')
yellow = [wlns(:) E(:,4)];
dlmwrite("illum_yellow.spd", yellow, 'delimiter', ' ')

% write out coeff2rgb and coeffs for easy copying to rendering system

csvwrite('coeff.csv', coeff2rgb);
csvwrite('esq.csv', esquiggle');
csvwrite('ssq.csv', ssquiggle');
