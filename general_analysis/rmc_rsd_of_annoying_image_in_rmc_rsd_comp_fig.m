clear all; close all; clc;

expns = {'5_rmc_rsd_comp'};
rmc_rsd_comp_exp = 1;
eizo_exps = [rmc_rsd_comp_exp]; % which exp in the iteration loop is it?
filter_exps = [rmc_rsd_comp_exp];

masks = {'../images/masks/obj_mask_correct.png', '../images/masks/refl_mask_correct.png', '../images/masks/obj_wout_high_mask_correct.png'};
mask_names = {'obj_mask', 'refl_mask', 'obj_wout_high_mask'};
stats = [struct()]; % exps x imgs x masks

tic;
expc = 1;
expn = expns{expc};
imgns = dir(['../images/exp' expn '/*.png']);

if any(expc == eizo_exps)
  gcRGB = csvread('../calibration/eizo_gamma.csv');
  monxyY = csvread('../calibration/eizo_chroma.csv');
  DKL2RGB = dkl2rgbFromCalib('../calibration/eizo_chroma.csv');
  RGB2LMS = rgb2lmsFromCalib('../calibration/eizo_mon_spectra.csv');
end
monxyz = xyY2XYZ(monxyY);
RGB2DKL_T = inv(DKL2RGB);

imgn = 'mitsuba_comp_rmc_rsd_illum_4_rg_1_by_2_ld_1_bkgd_voronoi_diff_dist.png';
img = gammaCorr(im2double(imread(['../images/exp' expn '/' imgn])), gcRGB);

for maskc = 1:length(masks)
  mask = imread(masks{maskc});
  mask_r = mask(:, :, 1);
  mask_g = mask(:, :, 2);
  mask_b = mask(:, :, 3);
  mask = mask_r;
  reject_pxs = mask_r == 127 & mask_g == 127 & mask_b == 127;
  reject_pxs = logical(reject_pxs(:));
  bkgd_pxs = mask_r == 0 & mask_g == 0 & mask_b == 0;
  bkgd_pxs = logical(bkgd_pxs(:));
  obj_pxs = ~(reject_pxs | bkgd_pxs);
  obj_pxs = logical(obj_pxs(:));
  mask(obj_pxs) = 1;
  mask = logical(mask(:));

  r = img(:, :, 1); rs = r(:);
  g = img(:, :, 2); gs = g(:);
  b = img(:, :, 3); bs = b(:);

  rgb_gc_lin = [rs, gs, bs];

  rgb_gc_masked = [rs(obj_pxs(:)), gs(obj_pxs(:)), bs(obj_pxs(:))];

  lms = real(rgb2lms(RGB2LMS, rgb_gc_lin'))';

  [rmc, rsd] = rmc_rsd(lms, bkgd_pxs, obj_pxs);

  disp(maskc);
  disp([rmc; rsd]);
end
