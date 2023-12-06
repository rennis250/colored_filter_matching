clear all; close all; clc;

gcRGB = csvread('../calibration/eizo_gamma.csv');

monxyY = csvread('../calibration/eizo_chroma.csv');
monxyz = xyY2XYZ(monxyY);

DKL2RGB = dkl2rgbFromCalib('../calibration/eizo_chroma.csv');
RGB2DKL_T = inv(DKL2RGB);

RGB2LMS = rgb2lmsFromCalib_For_Convergence('../calibration/eizo_mon_spectra.csv');

blue_img_wout_filter = gammaCorr(im2double(imread('../images/empty_scene_blue_illum.png')), gcRGB);
yellow_img_wout_filter = gammaCorr(im2double(imread('../images/empty_scene_yellow_illum.png')), gcRGB);

mask = imread('../images/masks/obj_mask.png');
mask = logical(mask);

imgs = dir('../images/patterned/*.png');
ic = 2
% for ic = 1:length(imgs)
    disp(['image num: ', num2str(ic)]);
    
    img_with_filter = gammaCorr(im2double(imread(['../images/patterned/' imgs(ic).name])), gcRGB);
    fparts = strsplit(imgs(ic).name, '_');
    ld_parts = strsplit(fparts{end}, '.');
    light_or_dark = ld_parts{1};
    
    if strcmp(light_or_dark, 'lightest')
        % we got the lightest version
        body_color = fparts{10};
        illum = fparts{7};
        dark_thresh = 2;
    else
        % we got the darker version
        body_color = fparts{9};
        illum = fparts{6};
        dark_thresh = 1;
    end
    
    if ~strcmp(illum, 'blue') && ~strcmp(illum, 'yellow')
%         continue;
        error('invalid illuminant!');
    end
    
    if strcmp(illum, 'blue')
        img_wout_filter = blue_img_wout_filter;
    else
        img_wout_filter = yellow_img_wout_filter;
    end
    
    %%%%%%%%%%%%
    % image with filter (convert to LMS, DKL, and LAB)
    
    r = img_with_filter(:,:,1); r = r(:);
    g = img_with_filter(:,:,2); g = g(:);
    b = img_with_filter(:,:,3); b = b(:);
    
    r_filter = r(mask(:));
    g_filter = g(mask(:));
    b_filter = b(mask(:));
    
    dkl_filter = real(rgb2dkl(RGB2DKL_T, [r_filter(:), g_filter(:), b_filter(:)]'))';
    
    %%%%%%%%%%%%
    % image without glass filter (convert to LMS, DKL, and LAB)
    
    rwof = img_wout_filter(:,:,1); rwof = rwof(:);
    gwof = img_wout_filter(:,:,2); gwof = gwof(:);
    bwof = img_wout_filter(:,:,3); bwof = bwof(:);
    
    r_wout_filter = rwof(mask(:));
    g_wout_filter = gwof(mask(:));
    b_wout_filter = bwof(mask(:));
    
    dkl_wout_filter = real(rgb2dkl(RGB2DKL_T, [r_wout_filter(:), g_wout_filter(:), b_wout_filter(:)]'))';
    
    [ld_rg_counts_orig, ld_yv_counts_orig, yv_rg_counts_orig, ld_rg_counts_filter, ld_yv_counts_filter, yv_rg_counts_filter] = maps_to_find_sinks_and_sources(dkl_wout_filter, dkl_filter);
% end

handles.ImgSeq = zeros(size(ld_rg_counts_orig, 1), size(ld_rg_counts_orig, 2), 2);
handles.ImgSeq(:, :, 1) = yv_rg_counts_orig;
handles.ImgSeq(:, :, 2) = yv_rg_counts_filter;
[handles.dim1,handles.dim2,handles.nFrames] = size(handles.ImgSeq);
handles.ImgSeqLoaded = 1;
handles.uvCLGcalculated = 0;
handles.uvHScalculated = 0;
handles.uvTScalculated = 0;

%% load the craniotomy window mask
handles.Mask = ones(size(ld_rg_counts_orig));
[handles.rMask,handles.cMask] = find(handles.Mask);
handles.idxMask = sub2ind(size(handles.Mask),handles.rMask,handles.cMask);

%% optical flow analysis
% run for HS and CLG then save
handles.SaveOF = 1;
% CLG params
handles.runCLG = 1;
handles.saveCLG = 1;
handles.CLGparams.alpha = 0.03;
handles.CLGparams.ratio = 0.5;
handles.CLGparams.minWidth = round(min(handles.dim1,handles.dim2)*handles.CLGparams.ratio/2);
handles.CLGparams.nOuterFPIterations = 7;
handles.CLGparams.nInnerFPIterations = 1;
handles.CLGparams.nSORIterations = 30;
% HS params
handles.runHS = 1;
handles.saveHS = 1;
handles.HSparams.alpha = 0.35;
handles.HSparams.iterations = 2000;
% specify frames that you want to run optical flow on
handles.FstartOF = 1;
handles.FendOF = handles.nFrames;
% run optical flow and save
handles = RunOpticalFlowAnalysisButton(handles);

%% Source-Sink analysis
handles.SaveSS = 1;
method = 'CLG';
handles = calc_save_SourceSink(handles,method);
method = 'HS';
handles = calc_save_SourceSink(handles,method);

