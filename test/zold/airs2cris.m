%
% NAME
%   airs2cris -- translate AIRS to CrIS radiances
%
% SYNOPSIS
%   [crad, cfreq, opt2] = airs2cris(arad, afrq, sfile, opt1)
%
% INPUTS
%   arad   - AIRS channel radiances, m x k array
%   afrq   - AIRS channel frequencies, m-vector
%   sfile  - AIRS HDF4 SRF tabulation file
%   opt1   - optional input parameters
%
% OPT1 FIELDS
%   dvb    - deconv grid step, default 0.1
%   bfile  - pinv temp file, default 'bconv.mat'
%   hapod  - 0 = default, 1 = Hamming apodization
%   inst_res - 'lowres' (default), 'hires1-3', 'hi3to2'
%   user_res - 'lowres' (default), 'hires'
%
% OUTPUTS
%   crad   - simulated CrIS radiances, n x k array
%   cfrq   - CrIS channel frequencies, n-vector
%   opt2   - optional, various internal values
%
% DISCUSSION
%   airs2cris translates AIRS to CrIS radiances by deconvolving
%   AIRS to an intermediate grid and then reconvolving to the CrIS
%   user grid.
%
%   output frequencies are the intersection of AIRS input frequencies
%   with the three CrIS bands, and can change if the input frequencies
%   change.  Results with anything less than the full AIRS 1c channel
%   set may be unpredictable.
%
%   The code is vectorized, so k obs can be processed with 1 call.
%   For most applications both input and output opts fields can be
%   omitted.  Dropping output opts may improve performance slightly.
%
%   Calculating the pseudo-inverse matrix for the deconvolution is
%   relatively slow so this data is saved in a file.  The matrix is
%   a function of the SRF tabulation file, channel frequencies, and
%   deconvolution grid step, and is updated if any of these change.
%
% COPYRIGHT
%   Copyright 2012-2013, Atmospheric Spectroscopy Laboratory.  
%   This code is distributed under the terms of the GNU GPL v3.
%
% AUTHOR
%   H. Motteler, 4 Oct 2013
%

function [crad, cfrq, opt2] = airs2cris(arad, afrq, sfile, opt1)

% defaults
dvb = 0.1;            % deconv grid step size
bfile = 'bconv.mat';  % temp file for pinv matrix
hapod = 0;            % no Hamming apodization

% process input options
if nargin == 4
  if isfield(opt1, 'dvb'), dvb = opt1.dvb; end
  if isfield(opt1, 'bfile'), bfile = opt1.bfile; end
  if isfield(opt1, 'hapod'), hapod = opt1.hapod; end
end

% CrIS params
wlaser = 773.1301;    % nominal value
bstr{1} = 'LW';       % band by numbber
bstr{2} = 'MW';
bstr{3} = 'SW';

% initialize outputs
crad = []; cfrq = []; 
bpts = []; bv1 = []; bv2 = [];

% deconvolve the AIRS radiances
% [afrq, ifrq] = trim_chans(afrq);
  [afrq, ifrq] = sort(afrq);
arad = arad(ifrq, :);
[brad, bfrq] = airs_decon(arad, afrq, sfile, bfile, dvb);

% loop on CrIS bands
for bi = 1 : 3

  % get the CrIS user grid
  band = bstr{bi};
  [inst, user] = inst_params(band, wlaser, opt1);

  % find AIRS and CrIS band intersection
  ftmp = user.v1 : user.dv : user.v2;
  [i, j] = seq_match(afrq, ftmp, 2 * user.dv);
  tv1 = ftmp(j(1));  tv2 = ftmp(j(end));

  % filter deconvolved data to band intersection 
  rtmp = bandpass(bfrq, brad, tv1, tv2, user.vr);

  % reconvolve to the CrIS user grid dv
  [rtmp, ftmp] = finterp(rtmp, bfrq, user.dv);
  ftmp = ftmp(:);

  % trim reconvolved data to band intersection 
  ix = find(tv1 <= ftmp & ftmp <= tv2);
  rtmp = rtmp(ix, :);
  ftmp = ftmp(ix);

  % option for hamming apodization
  if hapod
    rtmp = hamm_app(rtmp);
  end

  % concatenate output columns
  crad = [crad; rtmp];
  cfrq = [cfrq; ftmp];

  % save band size and info
  bpts = [bpts, length(ix)];
  bv1 = [bv1, tv1];
  bv2 = [bv2, tv2];
end

% option to return more data
if nargout == 3
  opt2.brad = brad;   % deconvolved radiances
  opt2.bfrq = bfrq;   % deconvolved frequencies
  opt2.afrq = afrq;   % trimmed input frequencies
  opt2.bpts = bpts;   % band size list
  opt2.bv1 = bv1;     % band v1 list
  opt2.bv2 = bv2;     % band v2 list
end

