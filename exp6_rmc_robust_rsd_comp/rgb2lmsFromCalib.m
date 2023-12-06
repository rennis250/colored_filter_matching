function [RGB2LMS] = rgb2lmsFromCalib(fn)
monspectra = csvread(fn);
rs = monspectra(:, 1);
gs = monspectra(:, 2);
bs = monspectra(:, 3);

lms_absorp = csvread('linss2_10e_1.csv');
wlns = lms_absorp(:, 1);
l = lms_absorp(:, 2);
m = lms_absorp(:, 3);
s = lms_absorp(:, 4);

offs1 = 11:length(rs);
offs2 = 1:(length(rs) - 10);

rL = rs(offs1)'*l(offs2);
rM = rs(offs1)'*m(offs2);
rS = rs(offs1)'*s(offs2);

gL = gs(offs1)'*l(offs2);
gM = gs(offs1)'*m(offs2);
gS = gs(offs1)'*s(offs2);

bL = bs(offs1)'*l(offs2);
bM = bs(offs1)'*m(offs2);
bS = bs(offs1)'*s(offs2);

RGB2LMS = [rL, gL, bL;
           rM, gM, bM;
           rS, gS, bS];
return
end
