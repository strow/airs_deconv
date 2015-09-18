% 
% airs_test1 -- compare CrIS AIRS with true AIRS
% 
% true AIRS: kcarta radiances convolved with AIRS SRFs
% 
% true CrIS: kcarta radiances convolved to CrIS channels
%
% CrIS AIRS: interpolate CrIS to an intermediate grid, typically at
% a 0.1 cm-1 spacing, and reconvolve to AIRS.  Compare CrIS AIRS vs
% true AIRS.

%-----------------
% test parameters
%-----------------

addpath /asl/matlib/h4tools
addpath /asl/packages/iasi_decon
addpath /asl/packages/airs_decon/source
addpath /asl/packages/ccast/source

dvb = 0.1;       % deconvolution frequency step
dvk = 0.0025;    % kcarta frequency spacing

% kcarta test data
kcdir = '/home/motteler/cris/sergio/JUNK2012/';
flist =  dir(fullfile(kcdir, 'convolved_kcart*.mat'));

% AIRS 1C channel frequencies
cfrq = load('freq2645.txt');  

% specify an AIRS SRF tabulation
sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';

% build the AIRS convolution matrix
[sconv, sfreq, tfreq] = mksconv1(sfile, cfrq, dvk);

% opts for kc2cris
opt1 = struct;
opt1.ng = 2;

% opts for inst_params
opt2 = struct;
opt2.resmode = 'hires2';
wlaser = 773.1301;

% cris user grid structs
[instLW, userLW] = inst_params('LW', wlaser, opt2);
[instMW, userMW] = inst_params('MW', wlaser, opt2);
[instSW, userSW] = inst_params('SW', wlaser, opt2);

%-----------------------------
% get true CrIS and true AIRS
%-----------------------------

% loop on kcarta files
rad1LW = []; rad1MW = []; rad1SW = []; rad2 = [];
for i = 1 : length(flist)
  d1 = load(fullfile(kcdir, flist(i).name));
  vkc = d1.w(:); rkc = d1.r(:);

  % kcarta to CrIS channel radiances
  [rtmp, frq1LW] = kc2cris(userLW, rkc, vkc, opt1);
  rad1LW = [rad1LW, rtmp];

  [rtmp, frq1MW] = kc2cris(userMW, rkc, vkc, opt1);
  rad1MW = [rad1MW, rtmp];

  [rtmp, frq1SW] = kc2cris(userSW, rkc, vkc, opt1);
  rad1SW = [rad1SW, rtmp];

  % kcarta to AIRS channel radiances
  [ix, jx] = seq_match(sfreq, vkc);
  rtmp = zeros(length(sfreq), 1);
  rtmp(ix) = rkc(jx);
  rtmp = sconv * rtmp;
  rad2 = [rad2, rtmp];

  fprintf(1, '.');
end
fprintf(1, '\n')
rad1 = [rad1LW; rad1MW; rad1SW];
frq1 = [frq1LW; frq1MW; frq1SW];
frq2 = tfreq(:);    % from mksconv1
clear d1 vkc rkc

%------------------------
% transform CrIS to AIRS
%------------------------

[rad4, frq4, rad3, frq3] = ...
  cris2airs(rad1LW, rad1MW, rad1SW, frq1LW, frq1MW, frq1SW, ...
            sfile, cfrq, opt1);
%-----------------
% stats and plots
%-----------------

% take radiances to brightness temps
bt1 = real(rad2bt(frq1, rad1));   % true CrIS
bt2 = real(rad2bt(frq2, rad2));   % true AIRS
bt3 = real(rad2bt(frq3, rad3));   % deconvolved CrIS
bt4 = real(rad2bt(frq4, rad4));   % CrIS to AIRS

% CrIS and AIRS spectra
figure(1); clf; j = 1; 
set(gcf, 'Units','centimeters', 'Position', [4, 10, 24, 16])
plot(frq1, bt1(:,j), frq2, bt2(:,j), frq3, bt3(:,j), frq4, bt4(:,j))
axis([600, 2800, 140, 320])
legend('true CrIS', 'true AIRS', 'CrIS decon', 'CrIS AIRS', ...
       'location', 'south')
xlabel('wavenumber'); ylabel('brighness temp')
title(sprintf('CrIS and AIRS profile %d', j));
grid on; zoom on
saveas(gcf, 'cris2airs_1', 'png')

% match true AIRS and CrIS AIRS channels
[ix, jx] = seq_match(frq2, frq4);
frq2 = frq2(ix);   
frq4 = frq4(jx);
bt2 = bt2(ix, :);  
bt4 = bt4(jx, :);

% CrIS AIRS minus true AIRS mean
figure(2); clf
set(gcf, 'Units','centimeters', 'Position', [4, 10, 24, 16])
subplot(2,1,1)
plot(frq2, mean(bt4 - bt2, 2))
axis([600, 2800, -5, 5])
xlabel('wavenumber'); ylabel('dBT')
title('CrIS AIRS minus true AIRS mean');
grid on; zoom on

% CrIS AIRS minus true AIRS std
subplot(2,1,2)
plot(frq2, std(bt4 - bt2, 0, 2))
axis([600, 2800, 0, 2])
xlabel('wavenumber'); ylabel('dBT')
title('CrIS AIRS minus true AIRS std');
grid on; zoom on
saveas(gcf, 'cris2airs_2', 'png')

