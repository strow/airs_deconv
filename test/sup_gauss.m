
% higher order Gaussian
function y = sup_gauss(x, b, c)
  c = c / 2.35482;
% y = exp(-(x - b).^2 / (2*c^2));
  y = exp(-((x - b).^2 / (2*c^2)).^1.5);
end

