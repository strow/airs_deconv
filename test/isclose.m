
function c = isclose(a, b, k)

if nargin == 2
  k = 1;
end

d = rms(a(:) - b(:)) / rms(b(:));

if d < eps * k
  c = 1;
else
  c = 0;
end

