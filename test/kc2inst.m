%
% NAME
%   kc2inst - convolve kcarta to CrIS sensor grid
%
% SYNOPSIS
%   [rad3, frq3, opt2] = kc2inst(inst, user, rkc, vkc, opt1)
%
% INPUTS
%   inst    - instrument params struct
%   user    - user grid params struct
%   rkc     - kcarta grid radiances
%   vkc     - frequency grid for rkc
%   opt1    - optional input parameters
%
% OUTPUTS
%   rad3    - cris radiances
%   frq3    - radiance frequency grid
%   opt2    - selected internal values
%
% DISCUSSION
%   derived from the setup and reference path parts of cris_igm2.m,
%   from the cris_sim package
% 
% key internal variables
%   user - user grid params
%   inst - instrument grid params 
%   N    - kcarta to interferogram transform size
%   N3   - undecimated half-path size, inst.df * inst.npts / 2
%   igm1 - N2+1 pt single-sided high res igm from kcarta radiances
%   rad3 - instrument grid radiances from igm1
%   frq3 - instrument grid frequencies for rad3
% 
% HM, 4 Mar 2013
%

function [rad3, frq3, opt2] = kc2inst(inst, user, rkc, vkc, opt1)

%------------------------------
% check inputs and set defauts
%------------------------------

% devault values
doplot = 0;        % option for lots of plots
figtype = 'fig';

% option to override defaults with opt1 fields
if nargin == 5
  optvar = fieldnames(opt1);
  for i = 1 : length(optvar)
    vname = optvar{i};
    if exist(vname, 'var')
      % fprintf(1, 'cris_igm: setting %s\n', vname)
      eval(sprintf('%s = opt1.%s;', vname, vname));
    else
      fprintf(1, 'kc2inst: unknown option variable %s\n', vname);
    end
  end
end

%-----------------------------------
% set up interferometric parameters
%-----------------------------------

dvk  = 0.0025;            % kcarta dv
dx = 1 / inst.vlaser;     % undecimated IGM distance step

% get dv and transform size
for k = 20 : 24
  N = 2^k;
  dv = 1 / (2*N*dx);
  if dv <= dvk, break, end
end            

% get Lmax and Vmax
Lmax = N * dx;
Vmax = N * dv;

% single-sided undecimated interferogram size
N3 = inst.df * inst.npts / 2;   

% sanity check for dx and inst.opd
if N3 ~= round(inst.opd / dx)
  error('df * npts / 2 != inst.opd / dx')
end

% fprintf(1, 'kc2inst: N = %7d, N3 = %5d, dx = %6.3e, dv = %6.3e\n', ...
%             N, N3, dx, dv);

%----------------------------------------------
% create a single-sided high res interferogram
%----------------------------------------------

% embed kcarta-grid radiance in 0 to Vmax N+1 point grid
frq2 = (0:N)' * dv;
rad2 = zeros(N+1, 1);

% set kcarta radiance passband to the user grid
rkc = bandpass(vkc, rkc, user.v1, user.v2, user.vr);

% interpolate rkc to the dv grid
ix = find(user.v1 - user.vr < frq2 & frq2 < user.v2 + user.vr);
rad2(ix) = interp1(vkc, rkc, frq2(ix), 'linear');
% rad3(ix) = interp1(vkc, rkc, frq2(ix), 'spline');

% compare interpolated vs actual kcarta radiances
% plot(vkc, rkc, frq2(ix), rad2(ix), frq2(ix), rad3(ix))
% ax=axis; ax(1)=user.v1-user.vr; ax(2)=user.v2+user.vr; axis(ax);
% legend('kcarta', 'linear', 'spline'); grid on; zoom on

% do the N+1 point cosine transform (as a 2*N point FFT)
% igm1 = real(ifft([rad2; rad2(N:-1:2,1)]));
igm1 = real(ifft([rad2; flipud(rad2(2:N,1))]));
% igm_full = igm1;
igm1 = igm1(1:N+1,1);

%---------------------
% reference transform
%---------------------
% truncate the single-sided high res igm1 to instrument OPD and
% transform back to radiance.  This is the usual interferometric
% interpolation as done in finterp or fconvkc
%
rad3 = real(fft([igm1(1:N3+1,1); flipud(igm1(2:N3,1))]));
frq3 = (0:N3)' * inst.dv;
ix = interp1(frq3, (1:N3+1)', inst.freq, 'nearest');
rad3 = rad3(ix);
frq3 = frq3(ix);

% option to return more data
if nargout == 3
  opt2 = struct;
  opt2.N = N;
  opt2.N3 = N3;
  opt2.igm = igm1;
end

%------------------
% option for plots
%------------------
if doplot == 0
  return
end

% plot reference radiance
figure(1); clf
plot(frq3, rad3)
title('CrIS reference radiance')
xlabel('wavenumber, 1/cm')
ylabel('radiance')
grid on; zoom on
saveas(gcf, 'figure_1', figtype)

