% 
% a2cris_test - basic AIRS-to-CrIS tests with 49 fitting profiles
% 
% key variables
%   rtc, vtc - true CrIS rad and freq
%   rta, vta - true AIRS rad and freq
%   rac, vac - AIRS-to-CrIS rad and freq
%
% note that setting opt1.scorr = 1 (statistical correction) forces
% opt1.hapod = 1 (Hamming apodization) in airs2cris.
%

addpath ../source
addpath /asl/packages/ccast/source

%-------------
% test params
%-------------

% CrIS params
band = 'SW';
wlaser = 773.1307;
opt1 = struct;
opt1.user_res = 'hires';
opt1.inst_res = 'hires3';
[inst, user] = inst_params(band, wlaser, opt1);

% test params
opt1.hapod = 1;   % 0 = off, 1 = hamming apodization
opt1.scorr = 1;   % 0 = off, 1 = do linear correction
switch opt1.user_res
  case 'lowres', opt1.cfile = 'corr_lowres.mat';
  case 'midres', opt1.cfile = 'corr_midres.mat';
  case 'hires',  opt1.cfile = 'corr_hires.mat';
end

% AIRS params
sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';
cfreq = load('freq2645.txt');
dvk = 0.0025; 
[sconv, sfreq, ofreq] = mksconv2(sfile, cfreq, dvk);

% kcarta test data
kcdir = '/home/motteler/cris/sergio/JUNK2012/';
flist =  dir(fullfile(kcdir, 'convolved_kcart*.mat'));

%-----------------------------
% get true CrIS and true AIRS
%-----------------------------

% loop on kcarta files
rtc = []; rta = [];
for i = 1 : length(flist)
  d1 = load(fullfile(kcdir, flist(i).name));
  vkc = d1.w(:); rkc = d1.r(:);

  % convolve kcarta radiances to CrIS channels
  [rtmp, ftmp] = kc2cris(user, rkc, vkc);
  rtc = [rtc, rtmp];

  % apply the AIRS convolution
  ix = interp1(vkc, 1:length(rkc), sfreq, 'nearest');
  rtmp = sconv * rkc(ix);
  rta = [rta, rtmp];

  fprintf(1, '.');
end
fprintf(1, '\n')
vtc = ftmp(:);     % from kc2cris
vta = ofreq(:);    % from mksconv
clear d1 vkc rkc

%------------------------
% transform AIRS to CrIS
%------------------------

[rac, vac, opt2] = airs2cris(rta, vta, sfile, opt1);
rdc = opt2.brad;
vdc = opt2.bfrq;

% option to apodize true CrIS
if opt1.hapod
  rtc = hamm_app(rtc);
end

%-----------------
% stats and plots
%-----------------

plot(diff(vac))
axis([0, 1200, 0, 3])

% take radiances to brightness temps
btc = real(rad2bt(vtc, rtc));   % true CrIS
bta = real(rad2bt(vta, rta));   % true AIRS
bdc = real(rad2bt(vdc, rdc));   % deconvolved AIRS
bac = real(rad2bt(vac, rac));   % AIRS-to-CrIS

% AIRS and CrIS spectra
figure(1); clf; j = 1; 
plot(vtc, btc(:,j), vta, bta(:,j), vdc, bdc(:,j), vac, bac(:,j))
switch band
 case 'LW', axis([650, 1100, 200, 300]);
 case 'MW', axis([1200, 1750, 200, 300]);
 case 'SW', axis([2175, 2550, 200, 300]);
end
legend('true CrIS', 'true AIRS', 'AIRS dec', 'AIRS CrIS', ...
       'location', 'southeast')
xlabel('wavenumber (cm-1)'); ylabel('BT (K)')
title(sprintf('AIRS 1C and CrIS %s profile %d', band, j));
grid on; zoom on

% AIRS CrIS minus true CrIS mean
figure(2); clf
subplot(2,1,1)
[itc, iac] = seq_match(vtc, vac);
plot(vtc(itc), mean(bac(iac,:) - btc(itc,:), 2))
switch band
 case 'LW', axis([650, 1100, -0.2, 0.2]);
 case 'MW', axis([1200, 1620, -0.2, 0.2]);
 case 'SW', axis([2175, 2550, -0.2, 0.2]);
end
title(sprintf('AIRS CrIS minus true CrIS %s mean', band));
ylabel('dBT (K)')
grid on; zoom on

% AIRS CrIS minus true CrIS std
subplot(2,1,2)
plot(vtc(itc), std(bac(iac,:) - btc(itc,:), 0, 2))
switch band
 case 'LW', axis([ 650, 1100, 0, 0.1]);
 case 'MW', axis([1200, 1620, 0, 0.1]);
 case 'SW', axis([2175, 2550, 0, 0.1]);
end
title(sprintf('AIRS CrIS minus true CrIS %s std', band));
xlabel('wavenumber (cm-1)'); 
ylabel('dBT (K)')
grid on; zoom on

% AIRS and AIRS-to-CrIS together
return
figure(3)
[x1, y1] = pen_lift(vta, bta(:,1));
[x2, y2] = pen_lift(vac, bac(:,1));
plot(x1, y1, x2, y2)
legend('AIRS L1c', 'AIR-to-CrIS', 'location', 'south')
axis([650, 2700, 200, 300])
title(sprintf('AIRS L1c and %s AIRS-to-CrIS', opt1.user_res))
xlabel('wavenumber (cm-1)'); 
ylabel('BT (K)')
grid on; zoom on

