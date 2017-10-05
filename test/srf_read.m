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
% AUTHOR
%  H. Motteler, 17 Nov 2016
%

function [tfreq, tgrid, srfval, chanid] = srf_read(sfile)

chanid = double(hdfread(sfile, 'chanid'));
tfreq  = double(hdfread(sfile, 'freq'));
fwgrid = double(hdfread(sfile, 'fwgrid'));
srfval = double(hdfread(sfile, 'srfval'));
width  = double(hdfread(sfile, 'width'));

[nchan, nspts] = size(srfval);
tgrid = (width * fwgrid) + tfreq * ones(1, nspts);

