function [m] = voronoi(x)
    n = floor(x);
    f = fract(x);

    m = [8.0, 8.0, 8.0];
    for j = -1:1
        for i = -1:1
            g = [i, j];
            o = hash22(n + g);
            r = g - f + (0.5 + 0.5*sin(6.2831*o));
            d = dot(r, r);
            if d < m(1)
                m = [d, o];
            end
        end
    end

    m = [sqrt(m(1)), m(2)+m(3)];
    return;
end
