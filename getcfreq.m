
function cfreq = getcfreq(hfile, cind)

% addpath /asl/matlab/h4tools

% read the srf data
[alist, fattr] = h4sdread(hfile);

for i = 1 : length(alist)
  switch alist{i}{1}
    case 'chanid', chanid = double(alist{i}{2})';
    case 'freq',   freq   = double(alist{i}{2})';
    case 'fwgrid', fwgrid = double(alist{i}{2})';
    case 'srfval', srfval = double(alist{i}{2})';
    case 'width',  width  = double(alist{i}{2})';
  end
end

if nargin > 1 
 % use cind as a channel index
  cfreq = freq(cind);
else
  % use the standard 1b channels
  cfreq = freq(1:2378,1);
end

