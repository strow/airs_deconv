%
% spectral space inverse hamming for column-order data
%

function r2 = hamm_inv(r1);

[m, n] = size(r1);

r2 = inv(full(mkhamm(m))) * r1;

