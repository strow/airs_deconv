
% generalized bell function
%
% y = bell_fun(x, c, w, p)
% x = input
% c = center
% w = fWMM
% p = boxyness (optional)

function y = bell_fun(x, c, w, p)

% default to rough AIRS SRF fit
if nargin == 3
  p = 1.8;
end

a = w/2;
y = 1 ./ (1 + abs((x-c)./a).^(2*p));

