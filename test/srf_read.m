%
% NAME
%   srf_read  - read an hdf4 format SRF file
%
% SYNOPSIS
%   [tfreq, tgrid, srfval, chanid] = srf_read(sfile)
%
% INPUTS
%   sfile  - HEF 4 SRF tabulation file
%   
% OUTPUTS
%   tfreq   - nchan vector of channel center frequencies
%   tgrid   - nchan x nspts frequency grid for SRF tabulations
%   srfval  - nchan x nspts SRF tabulations
%   chanid  - nchan vector of tabulation channel IDs
%
% DISCUSSION
%   uses old matlab HDF 4 libs; to get rid of warnings, try
%   warning('off', 'MATLAB:imagesci:hdf:removalWarningHDFSD')
%
% AUTHOR
%  H. Motteler, 17 Nov 2016
%

function [tfreq, tgrid, srfval, chanid] = srf_read(sfile)

[alist, fattr] = h4sdread(sfile);

for i = 1 : length(alist)
  switch alist{i}{1}
    case 'chanid', chanid = double(alist{i}{2})';
    case 'freq',   tfreq  = double(alist{i}{2})';
    case 'fwgrid', fwgrid = double(alist{i}{2})';
    case 'srfval', srfval = double(alist{i}{2})';
    case 'width',  width  = double(alist{i}{2})';
  end
end

[nchan, nspts] = size(srfval);
tgrid = (width * fwgrid) + tfreq * ones(1, nspts);

