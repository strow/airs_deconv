%
% pen_lift - add NaNs for pen lift at frequency gaps
%

function[x2, y2] = pen_lift(x1, y1)

x1 = x1(:);
y1 = y1(:);
n = length(x1);
m = length(y1);
if m ~= n, error('x1 and y1 lengths differ'), end

x2 = zeros(n+20, 1);
y2 = zeros(n+20, 1);

j = 1;
x2(j) = x1(j);
y2(j) = y1(j);

for i = 2 : length(x1)

  j = j + 1;

  if x1(i) > x1(i-1) + 10
     x2(j) = NaN;
     y2(j) = NaN;
     j = j + 1;
  end

  x2(j) = x1(i);
  y2(j) = y1(i);
end

x2 = x2(1:j);
y2 = y2(1:j);

% ********* test only *************
% x2 = x1;
% y2 = y1;

