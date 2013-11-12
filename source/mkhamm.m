%
% mkhamm - build a sparse matrix for fast hamming convolutions
%

function h = mkhamm(n)

w1 = 0.54;
w2 = 0.23;

s1 = 1 : n;
s2 = 2 : n;
s3 = 1 : n-1;

irow = [s1, s2, s3];
icol = [s1, s3, s2];

v1 = ones(1, n) * w1;
v2 = ones(1, 2*(n-1)) * w2;

h = sparse(irow, icol, [v1, v2]);

h(1, 1) = h(1, 1) + w2;
h(n, n) = h(n, n) + w2;

