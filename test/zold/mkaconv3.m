%
% NAME
%   mkaconv -- build an nchan x nchan AIRS convolution matrix
%
% SYNOPSIS
%   function mkaconv(sfile, afile);
%
% INPUTS
%   sfile  - HDF-format SRF filename
%   afile  - matlab output filename (optional)
% 
% OUTPUT is a Matlab data file, with the variables
%   Amat    - AIRS convolution matrix
%   Afin    - Amat input frequency list (a dv grid spanning all SRFs)
%   Afout   - Amat output frequency list (chan centers from the SRF file)
%   Aiout   - Amat output channel number (chan IDs from the SRF file)
% 
% DESCRIPTION
%
% AUTHOR
%  H. Motteler, 28 Aug 2012
%

% function mkaconv(sfile, afile)

% addpath /home/motteler/mot2008/hdf/h4tools

nargin = 2;
sfile = '/asl/matlab/srftest/srftables_m140f_withfake_mar08.hdf';
afile = 'asrf.mat';

% default SRF source file
if nargin < 1
  sfile = '/asl/data/airs/srf/srftablesV10.hdf';
end

% default matlab output file
if nargin < 2
  [path,name,ext] = fileparts(sfile);
  afile = [name,'.mat'];
end

% read the srf data
[alist, fattr] = h4sdread(sfile);

for i = 1 : length(alist)
  switch alist{i}{1}
    case 'chanid', chanid = double(alist{i}{2})';
    case 'freq',   freq   = double(alist{i}{2})';
    case 'fwgrid', fwgrid = double(alist{i}{2})';
    case 'srfval', srfval = double(alist{i}{2})';
    case 'width',  width  = double(alist{i}{2})';
  end
end
clear alist fattr

[nchan, nspts] = size(srfval);

% for now, use the L1B channel set
nchan = 2378;

% use the sorted channel set internally
[fsrt, isrt] = sort(freq(1:nchan));
chanid = chanid(isrt);
srfval = srfval(isrt, :);
width = width(isrt);

% srfreq is the nchan x nspts frequency grids for srfval
srfreq = (width * fwgrid) + fsrt * ones(1, nspts);

% get spanning frequencies for this SRF set
vmin = min(min(srfreq));
vmax = max(max(srfreq));

% build transform input frequency list Cfin
dv = 0.01;
v1 = ceil(vmin/dv)*dv;
v2 = floor(vmax/dv)*dv;
Cfin = (v1 : dv : v2)';
ncol = length(Cfin);

Amat = zeros(nchan, ncol);

% for ci = 1 : nchan
for ci = 2000 : nchan

   v1 = srfreq(ci, 1);
   v2 = srfreq(ci, nspts);

   i1 = min(find(v1 <= Cfin));
   i2 = max(find(Cfin <= v2));

%  stmp = interp1(srfreq(ci,:), srfval(ci,:), Cfin(i1:i2), 'linear');
   stmp = interp1(srfreq(ci,:), srfval(ci,:), Cfin(i1:i2), 'spline');

   Amat(ci, i1:i2) = stmp ./ sum(stmp);

%  plot(fsrt(i1:i2), stmp)
%  pause

   if i1 < i2
     continue
   end

   keyboard

end

rank(Amat)

% save the results
% eval(sprintf('save %s Amat Cfin Cfout Ciout', afile));

