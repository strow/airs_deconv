%
% spectral space inverse hamming for column-order data
%

function r2 = hamm_app(r1);

[m, n] = size(r1);

r2 = inv(full(mkhamm(m))) * r1;

