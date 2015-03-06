%
% NAME
%   kc2cris - convolve kcarta to CrIS channel radiances
%
% SYNOPSIS
%   [rad2, frq2] = kc2cris(user, rad1, frq1)
%
% INPUTS
%   user  - CrIS user grid params
%   rad1  - kcarta radiances, m x n array
%   frq1  - kcarta frequencies, m-vector
%
% OUTPUTS
%   rad2  - CrIS radiances, k x n array
%   frq2  - CrIS frequency grid, k-vector
%
% DISCUSSION
%   see doc/finterp.pdf for the derivations used here.
%
%   kc2cris uses inst and user structs from the ccast function
%   inst_params.  The relevant fields are
%
%   user fields
%     v1  - user grid start
%     v2  - user grid end
%     vr  - out-of-band rolloff
%     dv  - user grid dv
%
%   note that since we start with kcarta radiances, large n 
%   (number of observations) can quickly reach memory limits.
%
%   kc2cris calls isclose.m from airs_decon/test, bandpass.m from
%   ccast/source, and sets user grid parameters from inst_params.m, 
%   also in ccast/source
%
%   kc2cris and kc2iasi are identical except for the way the
%   parameters v1, v2, vr, and dv2 are set, and for the IASI
%   apodization
%
% HM, 22 Oct 2014
%

function [rad2, frq2] = kc2cris(user, rad1, frq1)

% check that array sizes match
frq1 = frq1(:);
[m, nobs] = size(rad1);
if m ~= length(frq1)
  error('rad1 and frq1 sizes do not match')
end

%-----------------------------------
% set up interferometric parameters
%-----------------------------------

% kcarta param
dv1  = 0.0025;       % kcarta dv

% check input frequency spacing
if abs(dv1 - (frq1(2) - frq1(1))) > 1e-10
  error('input frequency spacing not 0.0025 1/cm')
end

% CrIS params
v1 = user.v1;         % user grid start
v2 = user.v2;         % user grid end
vr = user.vr;         % out-of-band rolloff
vb = v2 + vr;         % transform max
dv2 = user.dv;        % user grid dv

% get rational approx to dv1/dv2
[m1, m2] = rat(dv1/dv2);
if ~isclose(m1/m2, dv1/dv2, 4)
  error('no rational approximation for dv1 / dv2')
end

% get the tranform sizes
for k = 4 : 24
  if m2 * 2^k * dv1 >= vb, break, end
end
N1 = m2 * 2^k;
N2 = m1 * 2^k;

% get (and check) dx
dx1 = 1 / (2*dv1*N1);
dx2 = 1 / (2*dv2*N2);
if ~isclose(dx1, dx2, 4)
  error('dx1 and dx2 are different')
end
dx = dx1;

% fprintf(1, 'kc2cris: N1 = %7d, N2 = %5d, dx = %6.3e\n', N1, N2, dx);

%-------------------------------
% take kcarta to CrIS radiances
%-------------------------------

% set kcarta radiance passband to the user grid
rad1 = bandpass(frq1, rad1, v1, v2, vr);

% embed kcarta radiance in a 0 to Vmax grid
ftmp = (0:N1)' * dv1;
rtmp = zeros(N1+1, nobs);
[ix, jx] = seq_match(ftmp, frq1);
rtmp(ix, :) = rad1(jx, :);

% radiance to interferogram
igm1 = real(ifft([rtmp; flipud(rtmp(2:N1, :))]));
igm1 = igm1(1:N1+1, :);

% apply a time-domain apodization
% dtmp = (0:N2)' * dx;
% apod = gaussapod(dtmp, 2) * ones(1, nobs);
% igm1(1:N2+1, :) = igm1(1:N2+1, :) .* apod;

% interferogram to radiance
rad2 = real(fft([igm1(1:N2+1,:); flipud(igm1(2:N2,:))]));
frq2 = (0:N2)' * dv2;

% return just the CrIS user grid
ix = find(v1 <= frq2 & frq2 <= v2);
rad2 = rad2(ix, :);
frq2 = frq2(ix);

