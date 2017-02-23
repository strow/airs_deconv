
% 2-point linear fit for channel widths
function y = airs_fwhm(x)
  x1 =  700;   y1 = 0.5;
  x2 = 2183;   y2 = 1.7;
  y = ((x - x1) ./ (x2 - x1)) .* (y2 - y1) + y1;
end

