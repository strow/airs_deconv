
function mksconv(sfile, mfile, dv);

% NAME
%
%   mksconv -- build a sparse SRF convolution matrix
%
% SYNOPSIS
%
%   function mksconv(sfile, mfile, dv);
%
% INPUTS
%
%   sfile  - HDF-format SRF filename
%   mfile  - matlab output filename (optional)
%   dv     - convolution input frequency step (optional)
% 
% OUTPUT is a Matlab data file, with the variables
%
%   Cmat    - sparse convolution matrix
%   Cfin    - Cmat input frequency list (a dv grid spanning all SRFs)
%   Cfout   - Cmat output frequency list (chan centers from the SRF file)
%   Ciout   - Cmat output channel number (chan IDs from the SRF file)
% 
% DESCRIPTION
%
%   mkconv creates an nchan by length(Cfin) sparse matrix Cmat such 
%   that r2 = Cmat * r1 is a convolution of r1 as specified by the HDF 
%   SRF data file sfile.  Cfin are the frequencies associated with the
%   columns of Cmat, and Cfout with the rows; Cfout is taken from the 
%   channel centers of the SRF file.
%   
%   If dv is not specified, a value of 0.0025 1/cm is used.  If mfile 
%   is not specified, <infile>.mat is used, where <infile> is the name 
%   (minus any extension) of the HDF input file.
%
%   mkconv prints a "." to stdout for every 100-th channel processed.
%
% BUGS
%
%   This function is very slow, and uses a lot of memory.  About 75%
%   of the time is spent in building the index lists for the sparse 
%   matrix, which would not be too hard to speed up.  But since mkconv 
%   only needs to be run when the SRF tables are updated, this may not 
%   be worth fixing.
%
% H. Motteler, 15 Jan 02
%


% default SRF source file
if nargin < 1
  sfile = '/asl/data/airs/srf/srftablesV10.hdf';
end

% default matlab output file
if nargin < 2
  [path,name,ext] = fileparts(sfile);
  mfile = [name,'.mat'];
end

% default to a 0.0025 1/cm input grid
if nargin < 3
  dv = 0.0025;
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

% srfreq is an nchan x nspts array of frequency points for srfval
srfreq = (width * fwgrid) + freq * ones(1, nspts);

% get spanning frequencies for this SRF set
vmin = min(min(srfreq));
vmax = max(max(srfreq));

% build transform input frequency list Cfin
v1 = ceil(vmin/dv)*dv;
v2 = floor(vmax/dv)*dv;
Cfin = (v1 : dv : v2)';

% Cfout is channel centers from the SRF file; Ciout is channel IDs
Cfout = freq; 
Ciout = chanid;

% initialize a non-sparse test matrix
% C = zeros(nchan, length(Cfin));

% intitalize the sparse matrix index lists
si = [];
sj = [];
sd = [];

f1 = Cfin(1);
f2 = Cfin(length(Cfin));

% loop on SRFs
for i = 1 : nchan

  % get the frequency span of the current SRF
  v1 = srfreq(i,1);
  v2 = srfreq(i,nspts);

  % if the SRF is outside the specified Cfin, skip this channel
  if v2 <= f1 | f2 <= v1
    % fprintf(1, 'mkconv(): WARNING -- SRF %d outside of input range\n', i)
    continue
  end

  % if the SRF overlaps the specified Cfin, lop it off to fit
  if v1 < f1 
    if v1 < f1 - dv
      fprintf(1, 'mkconv(): WARNING -- truncating LHS of SRF %d\n', j)
    end
    v1 = f1;
  end
  if f2 < v2
    if f2 + dv < v2
      fprintf(1, 'mkconv(): WARNING -- truncating RHS of SRF %d\n', j)
    end
    v2 = f2;
  end

  % find the indices of the current SRF in Cfin 
  v1ind = ceil((v1-f1)/dv) + 1;
  v2ind = floor((v2-f1)/dv) + 1;

  % interpolate the SRF to a subinterval of the Cfin grid
% s1 = interp1(srfreq(i,:), srfval(i,:), Cfin(v1ind:v2ind), 'spline');
  s1 = interp1(srfreq(i,:), srfval(i,:), Cfin(v1ind:v2ind), 'linear');

  % normalize the SRF
  s1 = s1 ./ sum(s1);

  % build a non-sparse test matrix
  % C(i, v1ind:v2ind) = s1';

  % build the index list for the sparse matrix
  si = [si; ones(v2ind-v1ind+1, 1)*i ];
  sj = [sj; (v1ind:v2ind)'];
  sd = [sd; s1];

  if mod(i,100) == 0
    fprintf(1, '.');
  end
end
fprintf(1, '\n');

% create the sparse matrix
Cmat = sparse(si, sj, sd, nchan, length(Cfin), length(sd));

% save the results
eval(sprintf('save %s Cmat Cfin Cfout Ciout', mfile));

