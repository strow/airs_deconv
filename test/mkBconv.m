%
% NAME
%   mkBconv - build a sparse gaussian convolution matrix
%
% SYNOPSIS
%   [Bconv, vcol, vrow] = mkBconv(dvb, res, p)
%
% INPUTS
%   dvb    - output grid spacing
%   res    - output resolving power, FWHM = v / res
%   p      - higher-order gaussian exponent
% 
% OUTPUTS
%   Bconv  - m x n sparse convolution matrix
%   vcol   - n-vector, Bconv cols, SRF tabulation grid
%   vrow   - m-vector, Bconv rows, output channel grid
% 
% DESCRIPTION
%   mkBconv returns a tabulation of "generalized Gaussian" SRFs 
%   with channel spacing dvb and FWHM = v / res, in an m x n sparse
%   convolution matrix.  vcol is chosen to span the SRFs of the m
%   output channels.
%
%   typical parameters for convolution of kcarta radiances to the
%   AIRS deconvolution grid are dvb = 0.1 cm-1 and res = 2000
%
% AUTHOR
%   H. Motteler, 3 Mar 2017
%

function [Bconv, vcol, vrow] = mkBconv(dvb, res, p)

% kcarta grid spacing
dvk = 0.0025;

% build row output grid
vb1 = 645;
vb2 = 2670;
nrow = 1 + round((vb2 - vb1) / dvb);
vrow = vb1 + (0 : nrow-1)' * dvb;

% build col kcarta res grid
vk1 = 640;
vk2 = 2675;
ncol = 1 + round((vk2 - vk1) / dvk);
vcol = vk1 + (0 : ncol-1)' * dvk;

% initialize sparse matrix lists
si = []; sj = []; sd = [];

% loop on channels, calculate SRFs at vcol grid
for i = 1 : nrow

  % eval span for SRF i
  vc = vrow(i);         % current channel center
% fwhm = res * dvb;     % mod for fixed FWHM
  fwhm = vc / res;      % current channel FWHM 
  v1 = vc - 4 * fwhm;   % low end of tabulation span
  v2 = vc + 4 * fwhm;   % high end of tabulation span
  jx = find(v1 <= vcol & vcol <= v2);  % tabulation index
  k = length(jx);

  % evaulate and normalize the SRF
  stmp = sup_gauss(vcol(jx), vc, fwhm, p);
  stmp = stmp ./ sum(stmp);

  % save sparse indices and data
  si = [si; ones(k, 1) * i];
  sj = [sj; jx];
  sd = [sd; stmp];
end

% build a sparse matrix from the lists
Bconv = sparse(si, sj, sd, nrow, ncol, length(sd));

end % mkBconv definition

