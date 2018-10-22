% 
% cris_test7 -- compare deconvolved AIRS with 0.2 cm-1 sinc 
% 
% includes the overview plot from cris_test5 but no stats, 
% uses radiances from one (typically the first) fitting profile
%
%   bkc = real(rad2bt(vkc, rkc));     % kcarta
%   bt1 = real(rad2bt(frq1, rad1));   % true CrIS
%   bt2 = real(rad2bt(frq2, rad2));   % true AIRS
%   bt3 = real(rad2bt(frq3, rad3));   % deconvolved AIRS
%   bt4 = real(rad2bt(frq4, rad4));   % AIRS CrIS
%   bt5 = real(rad2bt(frq5, rad5));   % decon grid sinc

%-----------------
% test parameters
%-----------------

% set paths to asl libs
addpath ../source
addpath ../h4tools
addpath /asl/packages/ccast/source

% test params
band = 'LW';            % cris band
hapod = 0;              % flag for Hamming apodization
dvb = 0.1;              % deconvolution frequency step
dvs = 0.2;              % karta to sinc ILS resolution
bfile = 'bconv7.mat';   % deconvolution temp file

% kcarta test data
kcdir = '/home/motteler/cris/sergio/JUNK2012/';
flist =  dir(fullfile(kcdir, 'convolved_kcart*.mat'));

% get the kcarta to AIRS convolution matrix
sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';
cfreq = load('freq2645.txt');
% cfreq = trim_chans(cfreq);
dvk = 0.0025; 
[sconv, sfreq, ofreq] = mksconv2(sfile, cfreq, dvk);

% get CrIS inst and user params
opts.resmode = 'lowres';
wlaser = 773.1301;  % nominal value
[inst, user] = inst_params(band, wlaser, opts);

%-----------------------------
% get true CrIS and true AIRS
%-----------------------------

% kcarta
i = 1;
d1 = load(fullfile(kcdir, flist(i).name));
vkc = d1.w(:); rkc = d1.r(:);
clear d1

% true CrIS
[rad1, ftmp] = kc2cris(user, rkc, vkc);
frq1 = ftmp(:);

% true AIRS
ix = interp1(vkc, 1:length(rkc), sfreq, 'nearest');
rad2 = sconv * rkc(ix);
frq2 = ofreq(:);

% decon grid sinc 
ix = find(min(cfreq) - 10 <= vkc & vkc <= max(cfreq) + 10);
vtmp = vkc(ix); rtmp = rkc(ix); 
rtmp = bandpass(vtmp, rtmp, min(cfreq), max(cfreq), 10);
% [rad5, frq5] = finterp(rtmp, vtmp, dvb);
[rad5, frq5] = finterp(rtmp, vtmp, dvs);

%------------------------
% transform AIRS to CrIS
%------------------------

opt1.dvb = dvb;
opt1.bfile = bfile;
opt1.hapod = hapod;

[rad4, frq4, opt2] = airs2cris(rad2, frq2, sfile, opt1);
rad3 = opt2.brad;
frq3 = opt2.bfrq;

% option to apodize true CrIS
if hapod
  rad1 = hamm_app(rad1);
end

%-----------------
% stats and plots
%-----------------

% take radiances to brightness temps
bkc = real(rad2bt(vkc, rkc));     % kcarta
bt1 = real(rad2bt(frq1, rad1));   % true CrIS
bt2 = real(rad2bt(frq2, rad2));   % true AIRS
bt3 = real(rad2bt(frq3, rad3));   % deconvolved AIRS
bt4 = real(rad2bt(frq4, rad4));   % AIRS CrIS
bt5 = real(rad2bt(frq5, rad5));   % decon grid sinc

% plot parameters
[i1, i4] = seq_match(frq1, frq4); 

pv1 = min(frq1(i1)) - 10; 
pv2 = max(frq1(i1)) + 10;
if hapod,  psf = 0.2; app = 'hamm';
else psf = 2.0; app = 'noap'; end

% AIRS and CrIS spectra
figure(1); clf;
% set(gcf, 'Units','centimeters', 'Position', [4, 10, 24, 16])
plot(frq1, bt1, frq2, bt2, frq3, bt3, frq4, bt4)
ax(1)=pv1; ax(2)=pv2; ax(3)=180; ax(4)=320; axis(ax)
legend('true CrIS', 'true AIRS', 'AIRS dec', 'AIRS CrIS', ...
       'location', 'southeast')
xlabel('wavenumber'); ylabel('brighness temp')
title(sprintf('AIRS 1C and CrIS %s profile %d', band, j));
grid on; zoom on
pname = sprintf('airs_cris_spec_%s_%s', band, app);
% saveas(gcf, pname, 'png')
% export_fig([pname, '.pdf'], '-m2', '-transparent')

% deconvolved AIRS and sinc ILS 
figure(2); clf
% set(gcf, 'Units','centimeters', 'Position', [4, 10, 24, 16])
% plot(frq3, bt3, vkc, bkc)
plot(frq3, bt3, frq5, bt5, 'linewidth', 2)
axis([660, 680, 200, 260])
legend('deconvolved AIRS', sprintf('%.2f cm-1 sinc ILS', dvs))
title(sprintf('deconvolved AIRS and sinc ILS at %.2f cm-1', dvs))
xlabel('wavenumber'); ylabel('brighness temp')
grid on; zoom on
% export_fig('airs_decon_res.pdf', '-m2', '-transparent')