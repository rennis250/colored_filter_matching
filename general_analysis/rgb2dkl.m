function dkl = rgb2dkl(T, rgb)
    rgbScaled = 2.0.*(rgb - 0.5);
    dkl = T*rgbScaled;
    return;
end
