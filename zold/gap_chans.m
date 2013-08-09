% 
% NAME
%   gap_chans - make a list of fake channel centers for gap edges
% 
% SYNOPSIS
%   vg = gap_chans(vc)
% 
% INPUT
%   vc  n-vector, channel center frequencies
%
% OUTPUT
%   vg  m-vector, new channel centers near gaps
%

function vg = gap_chans(vc)

k = 4;     % extra channels for each edge of the gap
w = 1.1;   % stretch factor applied to chanstep values

n = length(vc);
dvc = diff(vc);
dmin = chanstep(vc);
vg = [];

% loop on channels
for i = 1 : n - 1
  
  % test for a gap
  if dvc(i) > 4 * dmin(i)

    v1 = vc(i);
    v2 = vc(i+1);

    for j = 1 : k
      v1 = v1 + w * dmin(i);
      v2 = v2 - w * dmin(i+1);
      vg = [vg, v1, v2];
    end
  end
end

vg = sort(vg);

