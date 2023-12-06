function [imgp, imgb, gisp, gisb] = betterMaskingFunction(img)
    mask = img(:,:,1) ~= 0 & img(:,:,2) ~= 0 & img(:,:,3) ~= 0;
    pm = logical(mask); gisp = find(pm);

    bm = ~pm; gisb = find(bm);

    sm = logical(ones(size(pm,1),size(pm,2)));

    xc = 121; yc = 129;
    backx = 212; backy = 61;
    shadx = 600; shady = 281;

    gR = 2.2; gG = 2.2; gB = 2.2;

    [imgp, imgb] = maskAndLinearTransparent(img, pm, bm, sm, xc, yc, backx, backy, shadx, shady, gR, gG, gB);

    return;
end
