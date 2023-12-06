clear all; close all; clc

KbName('UnifyKeyNames');

Beeper;

monxyY = [.6804 .3073 30.94
    .2029 .6968 74.22
    .1527 .0508 6.74];
initmon(monxyY);

gR = 2.1102;
gG = 2.1243;
gB = 2.0170;

subID = input('What are the observers initials? ','s');
condition = input('Which condition? (u/p) ','s');

if strcmp(condition, 'u')
    ims = dir('images/uniform/mitsuba_caus*.png');
elseif strcmp(condition, 'p')
    ims = dir('images/patterned/mitsuba_caus*.png');
end
imgs = zeros(length(ims),400,400,3);
for x = 1:length(ims)
    imgs(x,:,:,:) = im2double(imread(['images/patterned/' ims(x).name]));
end

wd = 400;
he = 400;

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

if exist(['data/color_match_transparency_' subID '_' condition '.mat'], 'file') == 2
    load(['data/color_match_transparency_' subID '_' condition '.mat']);
end

st = 1;
tc = 1;

try
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask','General','FloatingPoint32BitIfPossible');
    PsychImaging('AddTask','General','EnableNative10BitFramebuffer');
    [w, wRect]=PsychImaging('OpenWindow', max(Screen('Screens')), [0.8 0.8 0.8]);
    Priority(MaxPriority(w));

    [width, height]=Screen('WindowSize', w);

    [xCenter, yCenter] = RectCenter(wRect);

    done = false;
    for tc = st:length(trialOrder)
        currt = trialOrder(tc);
        rgb1 = squeeze(imgs(currt,:,:,:));

        tex1 = Screen('MakeTexture', w, rgb1, [], [], 1);

        Beeper;

        trialDone = false;
        while ~trialDone
            FlushEvents;
            [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();

            Screen('FillRect', w, [avgBck(1), avgBck(2), avgBck(3)]);
            Screen('DrawTexture', w, tex1, [], [xCenter-wd/2, yCenter-he/2, xCenter+wd/2, yCenter+he/2]);

            rgb = squeeze(data(tc,1,:));
            rgb_gc = real(rgb'.^[gR gG gB]);

            Screen('FillRect', w, rgb_gc, [xCenter-50 + 300, yCenter-50, xCenter+50 + 300, yCenter+50]);

            Screen('Flip', w);

            if strcmp(KbName(keyCode), 'space')
                trialDone = true;
            elseif strcmp(KbName(keyCode), 'q')
                trialDone = true;
                done = true;
            end
        end

        Screen('FillRect', w, [avgBck(1), avgBck(2), avgBck(3)]);
        Screen('FillOval', w, [0, 0, 0], [xCenter-3 yCenter-3 xCenter+3 yCenter+3]);
        Screen('Flip', w);

        Screen('Close', tex1);

        Screen('FillRect', w, [avgBck(1), avgBck(2), avgBck(3)]);

        if done
            break;
        end
    end

    sca;
    RestoreCluts;
    ShowCursor;
catch e
    ShowCursor;
    sca;
    RestoreCluts;
    e
end

