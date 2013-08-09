%
% NAME
%   damp_pinv - produce a set of smoothed pseudo-SRFs
% 
% SYNOPSIS
%   d = damp_pinv(vc, vg) 
% 
% INPUTS
%   vc  n-vector, channel center frequencies
%   vg  m-vector, output frequency grid, spans vc
%
% OUTPUT
%   d  m x n matrix of pseudo-SRFs
%

function d = damp_pinv(vc, vg)

m = length(vc);
n = length(vg);

d = zeros(m,n);

% loop on rows
for i = 1 : m

  v1 = vc(i);
  v2 = vc(i);
% vr = chanstep(vc(i)) * 6;  % for smoothing
  vr = chanstep(vc(i)) * 2.2;  % for fake channels

  % get passband indices in vg
  j1 = max(find(vg <= v1));
% j2 = min(find(v2 <= vg));
  j2 = j1;

  % get indices for filter wings 
  k1 = max([find(vg <= v1 - vr); 1]);
  k2 = min([find(v2 + vr <= vg); n]);

  % get sizes of each segment
  n1 = j1 - k1;      % points in LHS rolloff
  n2 = k2 - j2;      % points in RHS rolloff
  n3 = j2 - j1 + 1;  % points in passband

  % scale cosine for the rolloffs
  f1 = (1+cos(pi+(0:n1-1)*pi/n1))/2;
  f2 = (1+cos((1:n2)*pi/n2))/2;

  % build the filter
  stmp = [f1, ones(1,n3), f2];
  d(i, k1:k2) =  stmp ./ sum(stmp);

end

