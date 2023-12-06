function [v] = hash22(p)
    p3 = fract([p(1) p(2) p(1)] .* [.1031, .1030, .0973]);
    p3 = p3 + dot(p3, [p3(2) p3(3) p3(1)] + 19.19);
    v = fract([(p3(1) + p3(2))*p3(3), (p3(1) + p3(3))*p3(2)]);
    return;
end
