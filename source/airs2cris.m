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
% opt1 fields (with defaults) include
%   hapod  - 0 = unapodized, 1 = Hamming apodization
%   scorr  - 0 = no correction, 1 = statistical correction
%   bfile  - 'bconv.mat', inverse temp file cache
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
%   Output frequencies are the concatenation of the intersection of
%   the two major AIRS sub-bands with the CrIS bands.  So we lose
%   some of the high end of the MW and low end of the CrIS SW bands.
%   Results with anything less than the full AIRS 1c channel set may
%   be unpredictable.
%
%   The band specs were calculated on the fly in earlier versions,
%   but it seems safer to make this explicit; see the default values
%   below.  To get the old values, set opt1.pH = [1095 1613.8 2550];
%   the other values should be the same.
%
%   The code is vectorized, so k obs can be processed with 1 call.
%   For most applications both input and output opts fields can be
%   omitted.
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

% general defaults
dvb = 0.1;  % deconv grid step size
hapod = 0;  % no Hamming apodization
scorr = 0;  % no statistical correction
bfile = 'bconv.mat';     % inverse matrix cache file
cfile = 'lin_corr.mat';  % linear correction weights

% band and filter parameters
%         LW     MW      SW
  pL = [ 650   1210    2182.5];   % passband low
  pH = [1095   1605    2550  ];   % passband high
  rL = [   5     20      10  ];   % rolloff low
  rH = [  20     10      20  ];   % rolloff high

% process input options
if nargin == 4
  if isfield(opt1, 'dvb'), dvb = opt1.dvb; end
  if isfield(opt1, 'hapod'), hapod = opt1.hapod; end
  if isfield(opt1, 'scorr'), scorr = opt1.scorr; end
  if isfield(opt1, 'bfile'), bfile = opt1.bfile; end
  if isfield(opt1, 'cfile'), cfile = opt1.cfile; end
  if isfield(opt1, 'pL'), pL = opt1.pL; end
  if isfield(opt1, 'pH'), pH = opt1.pH; end
  if isfield(opt1, 'rL'), rL = opt1.rL; end
  if isfield(opt1, 'rH'), rH = opt1.rH; end
else
  opt1 = struct;
end

% setup for statistical correction
if scorr
  hapod = 1; 
  load lin_corr
end

% CrIS params
wlaser = 773.1301;
bstr = {'LW' 'MW' 'SW'};

% initialize outputs
crad = [];  cfrq = []; 

% sort the AIRS channels
[afrq, ifrq] = sort(afrq);
arad = arad(ifrq, :);

% deconvolve the AIRS radiances
[brad, bfrq] = airs_decon(arad, afrq, sfile, bfile, dvb);

% loop on CrIS bands
for j = 1 : 3

  % get the CrIS user grid
  [inst, user] = inst_params(bstr{j}, wlaser, opt1);

  % AIRS spanning band for this channel
  ix = find(pL(j)-rL(j) <= bfrq & bfrq <= pH(j)+rH(j));

  % filter deconvolved data to band intersection 
  rtmp = bandpass(bfrq(ix), brad(ix,:), pL(j), pH(j), rL(j), rH(j));

  % reconvolve to the CrIS user grid dv
  [rtmp, ftmp] = finterp(rtmp, bfrq(ix), user.dv);
  ftmp = ftmp(:);

  % trim reconvolved data to the intersection grid
  ix = find(pL(j) <= ftmp & ftmp <= pH(j));
  rtmp = rtmp(ix, :);
  ftmp = ftmp(ix);

  % option for hamming apodization
  if hapod
    rtmp = hamm_app(rtmp);
  end

  % option for statistical correction
  if scorr
    nchan = length(ftmp);
    switch j
      case 1, Ptmp = Pcor2LW; vtmp = tcfrqLW;
      case 2, Ptmp = Pcor2MW; vtmp = tcfrqMW;
      case 3, Ptmp = Pcor2SW; vtmp = tcfrqSW;
    end
    if length(vtmp) ~= length(ftmp)
      error('statistical correction file array mismatch')
    end
    btmp = real(rad2bt(ftmp, rtmp));
    for i = 1 : nchan
      btmp(i, :) = polyval(Ptmp(i,:), btmp(i, :));
    end
    rtmp = bt2rad(ftmp, btmp);
  end

  % concatenate output columns
  crad = [crad; rtmp];
  cfrq = [cfrq; ftmp];
end

% option to return internal values
if nargout == 3
  opt2.brad = brad;   % deconvolved radiances
  opt2.bfrq = bfrq;   % deconvolved frequencies
  opt2.afrq = afrq;   % sorted AIRS frequencies
  opt2.dvb = dvb;
  opt2.bfile = bfile;
  opt2.hapod = hapod;
  opt2.pL = pL;
  opt2.pH = pH;
  opt2.rL = rL;
  opt2.rH = rH;
end

