%
% NAME
%   L1d_conv -- build a sparse AIRS L1d convolution matrix
%
% SYNOPSIS
%   [dconv, v_tab, v_L1d] = L1d_conv(res, dv_tab, v_base)
%
% INPUTS
%   res    - resolving power, v/FWHM
%   dv_tab - step size for the SRF tabulation grid
%   v_base - start of output grid, i.e., v_L1d(1)
% 
% OUTPUTS
%   dconv  - m x n sparse convolution matrix
%   v_tab  - n-vector, sconv cols, SRF tabulation grid
%   v_L1d  - m-vector, sconv rows, output channel grid
% 
% DESCRIPTION
%   L1d_conv returns a tabulation of "generalized Gaussian" SRFs
%   with FWHM = v/res and channel spacing dv = FWHM/2 in an m x n
%   sparse convolution matrix.  v_tab is chosen to span the SRFs 
%   of the requested channels.
%
% AUTHOR
%   H. Motteler, 21 Feb 2017
%

function [dconv, v_tab, v_L1d] = L1d_conv(res, dv_tab, v_base)

% get the L1c channel set
v_L1c = load('freq2645.txt');
n_L1c = length(v_L1c);

% get the two L1c bands
v1b1 = v_L1c(1);
v2b1 = v_L1c(2162);
v1b2 = v_L1c(2163);
v2b2 = v_L1c(end);

% build L1d channel grid
k = 6000; 
v_L1d = zeros(k, 1);
v_L1d(1) = v_base;
for i = 2 : k
  v_L1d(i) = v_L1d(i-1) + v_L1d(i-1) / (2*res);
  if v_L1d(i) > v2b2, break, end
end

% trim v_L1d to the two L1c bands
ix = (v1b1 <= v_L1d & v_L1d <= v2b1) | (v1b2 <= v_L1d & v_L1d <= v2b2);
v_L1d = v_L1d(ix);
n_L1d = length(v_L1d);

% get spanning frequencies for the L1d SRF set
span = 4;
v1 = v_L1d(1) - span * v_L1d(1) / res;
v2 = v_L1d(end) + span * v_L1d(end) / res;

% get the dconv column grid
n = round((v2 - v1) / dv_tab);
v_tab = v1 + (0 : n)' * dv_tab;
n_tab = n + 1;

% initialize sparse matrix lists
si = []; sj = []; sd = [];

% loop on channels, calculate SRFs at v_tab grid
for i = 1 : n_L1d

  % eval span for SRF i
  vc = v_L1d(i);                  % current channel center
  vs = vc / res;                  % current channel FWHM 
  v1 = vc - span * vs;            % low end of tabulation span
  v2 = vc + span * vs;            % high end of tabulation span
  jx = find(v1 <= v_tab & v_tab <= v2);  % tabulation index
  k = length(jx);

  % evaulate and normalize the SRF
  stmp = sup_gauss(v_tab(jx), vc, vs);
  stmp = stmp ./ sum(stmp);

  % save sparse indices and data
  si = [si; ones(k, 1) * i];
  sj = [sj; jx];
  sd = [sd; stmp];
end

% build a sparse matrix from the lists
dconv = sparse(si, sj, sd, n_L1d, n_tab, length(sd));

end % L1d_conv definition

% % higher order Gaussian
% function y = sup_gauss(x, b, c)
%   c = c / 2.35482;
% % y = exp(-(x - b).^2 / (2*c^2));
%   y = exp(-((x - b).^2 / (2*c^2)).^1.5);
% end

