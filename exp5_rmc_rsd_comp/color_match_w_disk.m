clear all; close all; clc

KbName('UnifyKeyNames');

Beeper;

monxyY = [.6804 .3073 30.94
    .2029 .6968 74.22
    .1527 .0508 6.74];
initmon(monxyY);

rgb2lms = rgb2lmsFromCalib('./calibration/mon_spectra.csv');
lms2rgb = inv(rgb2lms);

monitor_spectra = csvread('./calibration/mon_spectra.csv');

lms = csvread('../aux/linss2_10e_1.csv');
wlns = lms(:,1);
lms = lms(1:14:end,2:end);

vor_texture = csvread('./vor.csv');

circ_mask = zeros(256,256);
for r = 1:60
    for theta = 0.001:0.001:2*pi;
        y = r*sin(theta);
        x = r*cos(theta);
        if sqrt(x^2 + y^2) < 256
            y = ceil(y) + 256/2;
            x = ceil(x) + 256/2;
            circ_mask(y, x) = 1.0;
        end
    end
end
idxs = find(circ_mask == 1.0);
[y, x] = ind2sub(size(circ_mask), find(circ_mask == 1.0));

bluetrans = dlmread('../base_stimuli/spectra/munsell_blue_better.spd');
bluetranstmp = [];
bluetranstmp(:,1) = wlns(1:14:end);
bluetranstmp(:,2) = interp1(bluetrans(:,1), bluetrans(:,2), wlns(1:14:end));
bluetrans = bluetranstmp;

greentrans = dlmread('../base_stimuli/spectra/munsell_green_better.spd');
greentranstmp = [];
greentranstmp(:,1) = wlns(1:14:end);
greentranstmp(:,2) = interp1(greentrans(:,1), greentrans(:,2), wlns(1:14:end));
greentrans = greentranstmp;

redtrans = dlmread('../base_stimuli/spectra/munsell_red_better.spd');
redtranstmp = [];
redtranstmp(:,1) = wlns(1:14:end);
redtranstmp(:,2) = interp1(redtrans(:,1), redtrans(:,2), wlns(1:14:end));
redtrans = redtranstmp;

yellowtrans = dlmread('../base_stimuli/spectra/munsell_yellow_better.spd');
yellowtranstmp = [];
yellowtranstmp(:,1) = wlns(1:14:end);
yellowtranstmp(:,2) = interp1(yellowtrans(:,1), yellowtrans(:,2), wlns(1:14:end));
yellowtrans = yellowtranstmp;

whiteillum = max(max(monitor_spectra))*ones(length(wlns(1:14:end)),1);

vor_col = vor_texture.*3;
vor_col = bsxfun(@times, vor_col, whiteillum');
vor_col2 = vor_col;

gR = 2.1102;
gG = 2.1243;
gB = 2.0170;

subID = input('What are the observers initials? ','s');
condition = input('Which condition? (u/p) ','s');

if strcmp(condition, 'u')
    ims = dir('../stimuli/images/*.png');
elseif strcmp(condition, 'p') || strcmp(condition, 'movie')
    ims = dir('../stimuli/images/*.png');
end
imgs = zeros(length(ims),400,400,3);
for x = 1:length(ims)
    imgs(x,:,:,:) = im2double(imread(['../stimuli/images/' ims(x).name]));
end

wd = 400;
he = 400;

st = 1;
tc = 1;
ntrials = 5;

data = zeros(length(ims)*ntrials, 3, 3);

trialOrder = [];
for x = 1:ntrials
    if x == 1
        trialOrder = randperm(length(ims));
    else
        trialOrder = [trialOrder randperm(length(ims))];
    end
end
trialOrder = trialOrder(randperm(length(trialOrder)));

avgBck = [0.8, 0.8, 0.8];

if exist(['./data/color_match_transparency_rmc_rsd_' subID '_' condition '.mat'], 'file') == 2
    load(['./data/color_match_transparency_rmc_rsd_' subID '_' condition '.mat']);
end

st = tc;

try
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask','General','FloatingPoint32BitIfPossible');
    PsychImaging('AddTask','General','EnableNative10BitFramebuffer');
    [w, wRect]=PsychImaging('OpenWindow', max(Screen('Screens')), [0.8 0.8 0.8], [0 0 800 800]);
    Priority(MaxPriority(w));

%     if strcmp(condition, 'movie')
        HideCursor;
%     end

    [width, height]=Screen('WindowSize', w);

    [xCenter, yCenter] = RectCenter(wRect);

    Screen('FillRect', w, [avgBck(1), avgBck(2), avgBck(3)]);
    Screen('FillOval', w, [0, 0, 0], [xCenter-3 yCenter-3 xCenter+3 yCenter+3]);
    Screen('Flip', w);

    if ~strcmp(subID, 'test')
        WaitSecs(60);
    end

    Beeper;

    Screen('FillRect', w, [avgBck(1), avgBck(2), avgBck(3)]);
    Screen('FillOval', w, [1, 0, 0], [xCenter-3 yCenter-3 xCenter+3 yCenter+3]);
    Screen('Flip', w);

    WaitSecs(0.5);

    Beeper;

    Screen('FillRect', w, [avgBck(1), avgBck(2), avgBck(3)]);
    Screen('FillOval', w, [0, 0, 0], [xCenter-3 yCenter-3 xCenter+3 yCenter+3]);
    Screen('Flip', w);

    WaitSecs(0.5);

    Beeper;

    Screen('FillRect', w, [avgBck(1), avgBck(2), avgBck(3)]);
    Screen('FillOval', w, [1, 0, 0], [xCenter-3 yCenter-3 xCenter+3 yCenter+3]);
    Screen('Flip', w);

    WaitSecs(0.5);

    Screen('FillRect', w, [avgBck(1), avgBck(2), avgBck(3)]);
    Screen('FillOval', w, [0, 0, 0], [xCenter-3 yCenter-3 xCenter+3 yCenter+3]);
    Screen('Flip', w);

    WaitSecs(1.5);

    done = false;
    for tc = st:length(trialOrder)
        currt = trialOrder(tc);
        rgb1 = squeeze(imgs(currt,:,:,:));

        tex1 = Screen('MakeTexture', w, rgb1, [], [], 1);

        Beeper;

        oog = false;
        ld = 2*rand();
        rg_mix = 0.5;
        by_mix = 0.5;
        voronoi_img = zeros(256, 256, 3);

        trialDone = false;
        frameC = 1;
        while ~trialDone
            FlushEvents;
            [xm,ym,buttons,focus,valuators,valinfo] = GetMouse(w);
            [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();

            if buttons(3)
                if ~oog
                    ld_old = ld;
                end
                ld = ld + 0.02;
                if ld > 3.0
                    ld = 1.0;
                end
            elseif buttons(1)
                if ~oog
                    ld_old = ld;
                end
                ld = ld - 0.02;
                if ld < -1.0
                    ld = 0.0;
                end
            elseif strcmp(KbName(keyCode), 'q')
                trialDone = true;
                done = true;
            end

            Screen('FillRect', w, [avgBck(1), avgBck(2), avgBck(3)]);
            Screen('DrawTexture', w, tex1, [], [xCenter-wd/2, yCenter-he/2, xCenter+wd/2, yCenter+he/2]);

            if ~oog
                rg_mix_old = rg_mix;
                by_mix_old = by_mix;
                voronoi_img_old = voronoi_img;
            end

            rg_mix = xm/width;
            by_mix = ym/height;

            if ~oog
                filter = ld.*(((rg_mix.*redtrans(:,2) + (1 - rg_mix).*greentrans(:,2)) + (by_mix.*bluetrans(:,2) + (1 - by_mix).*yellowtrans(:,2))));
                vor_col2(idxs, :) = bsxfun(@times, vor_col(idxs, :), filter'.^2.2);
                vor_col2(isnan(vor_col2) | isinf(vor_col2)) = 0.0;
                voronoi_img = reshape(vor_col2*lms*lms2rgb', 256, 256, 3);
            end

            if any(voronoi_img < 0 | voronoi_img > 1)
                oog = true;
                voronoi_tex = Screen('MakeTexture', w, real(voronoi_img_old.^(1/2.2)), [], [], 1);
                Screen('DrawTexture', w, voronoi_tex, [], [xCenter-128 + 340, yCenter-128, xCenter+128 + 340, yCenter+128]);
            else
                oog = false;
                voronoi_tex = Screen('MakeTexture', w, real(voronoi_img.^(1/2.2)), [], [], 1);
                Screen('DrawTexture', w, voronoi_tex, [], [xCenter-128 + 340, yCenter-128, xCenter+128 + 340, yCenter+128]);
            end

            Screen('Flip', w);

            if buttons(2)
                data(tc,1) = rg_mix_old;
                data(tc,2) = by_mix_old;
                data(tc,3) = ld_old;
                trialDone = true;
            end

            if strcmp(condition, 'movie')
                imageArray = Screen('GetImage', w, [xCenter - 210, yCenter - 200, xCenter + 480, yCenter + 200]);
                imwrite(imageArray, ['../images/movie_frames/' num2str(frameC) '.png']);
                frameC = frameC + 1;
            end
        end

        WaitSecs(1);

        Screen('FillRect', w, [avgBck(1), avgBck(2), avgBck(3)]);
        Screen('FillOval', w, [0, 0, 0], [xCenter-3 yCenter-3 xCenter+3 yCenter+3]);
        Screen('Flip', w);

        Screen('Close', tex1);
        Screen('Close', voronoi_tex);

        Beeper;
        WaitSecs(0.2);
        Beeper;

        Screen('FillRect', w, [avgBck(1), avgBck(2), avgBck(3)]);
        DrawFormattedText(w, 'Please rate the quality of your match (1-5).', 'center', 'center');
        Screen('Flip', w);

        answered = false;
        while ~answered
            [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
            if strcmp(KbName(keyCode), '1')
                answered = true;
                data(tc,4) = 1;
            elseif strcmp(KbName(keyCode), '2')
                answered = true;
                data(tc,4) = 2;
            elseif strcmp(KbName(keyCode), '3')
                answered = true;
                data(tc,4) = 3;
            elseif strcmp(KbName(keyCode), '4')
                answered = true;
                data(tc,4) = 4;
            elseif strcmp(KbName(keyCode), '5')
                answered = true;
                data(tc,4) = 5;
            end
        end

        Screen('FillRect', w, [avgBck(1), avgBck(2), avgBck(3)]);
        Screen('Flip', w);

        if ~strcmp(subID, 'test')
            save(['../data/color_match_transparency_rmc_rsd_' subID '_' condition '.mat']);
        end

        if done
            break;
        end
    end

    sca;
    RestoreCluts;
    ShowCursor;

    if ~strcmp(subID, 'test')
        save(['../data/color_match_transparency_rmc_rsd_' subID '_' condition '.mat']);
    end
catch e
    ShowCursor;
    sca;
    RestoreCluts;
    e
end

