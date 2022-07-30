function [radial_avg] = ravg(f, r)
    %ravg Average of field f at radius r
    [ysize xsize] = size(f);
    ycenter = fix(ysize / 2);
    xcenter = fix(xsize / 2);
    r = min([r ycenter xcenter ysize-ycenter xsize-xcenter]);
    circle = false(size(f));
    for theta=0.0:0.1:2*pi
        circle( ...
            fix(ycenter+r*cos(theta)), ...
            fix(xcenter+r*sin(theta))) = true;
    end
    radial_avg = mean(f(circle), 'all');
end