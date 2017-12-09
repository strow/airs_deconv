%
% dt = dr2dt(v, r, dr), 
%
% dt = B(v, r + dr) - B(v, r);
%

function dt = dr2dt(v, r, dr)

dt = real(rad2bt(v, r + dr)) - real(rad2bt(v, r));

