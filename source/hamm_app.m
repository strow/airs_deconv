%
% spectral space hamming apodization for column-order data
%

function r2 = hamm_app(r1);

[m, n] = size(r1);

r2 = mkhamm(m) * r1;

