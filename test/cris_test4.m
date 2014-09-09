% 
% cris_test4
% 
% reference truth: start with kcarta radiances, convolve these to 
% CrIS radiances at the user grid, and call the result “true CrIS”.
% 
% deconvolution: start with kcarta radiances, convolve to AIRS
% channel radiances (“true AIRS”), deconvolve to an intermediate
% grid, e.g. 0.05 1/cm spacing, and reconvolve to the CrIS user grid
% (“AIRS Cris”).  Then compare AIRS CrIS vs true CrIS for various
% profiles and any of the three CrIS bands
%
% derived frim cris_test3, processes all the fitting profiles and
% plots the means of residuals.  Calls the airs2cris wrapper instead
% of airs_decon, and includes an option for hamming apodization of
% all the test radiances.
%
%-----------------
% test parameters
%-----------------

% use my bcast utils and HDF libs
addpath /home/motteler/cris/ccast/source
addpath /home/motteler/cris/airs_decon/source
addpath /home/motteler/mot2008/hdf/h4tools

% test params
band = 'SW';            % cris band
bfile = 'bconv4.mat';   % deconvolution temp file
dvb = 0.1;              % deconvolution frequency step
fig = 'fig';            % plot type
dohamm = 1;

% kcarta test data
kcdir = '/home/motteler/cris/sergio/JUNK2012/';
flist =  dir(fullfile(kcdir, 'convolved_kcart*.mat'));

% get the kcarta to AIRS convolution matrix
sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';
cfreq = load('freq2645.txt');
cfreq = trim_chans(cfreq);
dvs = 0.0025; 
[sconv, sfreq, ofreq] = mksconv2(sfile, cfreq, dvs);

% get CrIS inst and user params
opts.resmode = 'lowres';
wlaser = 773.1301;  % nominal value
[inst, user] = inst_params(band, wlaser, opts);

% set wlaser so inst grid == user grid
wlaser = 1e7/(inst.df/(2*user.opd/inst.npts));
[inst, user] = inst_params(band, wlaser, opts);

% intersection of AIRS and CrIS bands
switch upper(band)
  case 'LW', tv1 =  650; tv2 = 1095; dt1 = 2; dt2 = 6;
  case 'MW', tv1 = 1210; tv2 = 1614; dt1 = 2; dt2 = 6;
  case 'SW', tv1 = 2185; tv2 = 2550; dt1 = 2; dt2 = 2;
end

%-----------------------------------
% true CrIS and true AIRS radiances
%-----------------------------------

% loop on kcarta files
rad1 = []; rad2 = [];
for i = 1 : length(flist)
  d1 = load(fullfile(kcdir, flist(i).name));
  vkc = d1.w(:); rkc = d1.r(:);

  % convolve kcarta radiances to CrIS channels
  [rtmp, ftmp] = kc2cris(inst, user, rkc, vkc);
  rad1 = [rad1, rtmp];

  % apply the AIRS convolution
  ix = interp1(vkc, 1:length(rkc), sfreq, 'nearest');
  rtmp = sconv * rkc(ix);
  rad2 = [rad2, rtmp];

  fprintf(1, '.');
end
fprintf(1, '\n')
frq1 = ftmp(:);     % from kc2cris
frq2 = ofreq(:);    % from mksconv
clear d1 vkc rkc

%-------------------------------
% convolution and interpolation
%-------------------------------

opt1.dvb = dvb;
opt1.bfile = bfile;

[rad4, frq4, opt2] = airs2cris(rad2, frq2, sfile, opt1);
rad3 = opt2.brad;
bfrq = opt2.bfrq;

% deconvolve the AIRS radiances
% [rad3, bfrq] = airs_decon(rad2, cfreq, sfile, bfile, dvb);

% try an extra smoothing step
% rad3 = mkhamm(length(rad3)) * rad3;

% apply the bandpass filter
% rad3 = bandpass(bfrq, rad3, tv1, tv2, user.vr);

% reconvolve to CrIS
% [rad4, frq4] = finterp(rad3, bfrq, user.dv);
% frq4 = frq4(:);

% AIRS direct interpolation to CrIS
rad5 = interp1(frq2, rad2, frq1, 'spline', 'extrap');

% AIRS interpolation and convolution to CrIS
vx = ceil(frq2(1)/dvb) * dvb;
nx = floor((frq2(end) - frq2(1)) / dvb);
ftmp = vx + (0 : nx - 1) * dvb;
rtmp = interp1(frq2, rad2, ftmp', 'spline', 'extrap');
rtmp = bandpass(ftmp, rtmp, tv1, tv2, user.vr);
[rad6, frq6] = finterp(rtmp, ftmp, user.dv);
frq6 = frq6(:);

%-----------------
% stats and plots
%-----------------

% band-specific plot values
switch upper(band)
  case 'LW', dt1 = 2; dt2 = 6;
  case 'MW', dt1 = 2; dt2 = 6;
  case 'SW', dt1 = 2; dt2 = 2;
end

% option for apodization
if dohamm
  rad1 = hamm_app(rad1);
  rad2 = hamm_app(rad2);
  rad3 = hamm_app(rad3);
  rad4 = hamm_app(rad4);
  rad5 = hamm_app(rad5);
  rad6 = hamm_app(rad6);
end

% take radiances to brightness temps
bt1 = real(rad2bt(frq1, rad1));   % true CrIS
bt2 = real(rad2bt(frq2, rad2));   % true AIRS
bt3 = real(rad2bt(bfrq, rad3));   % deconvolved AIRS
bt4 = real(rad2bt(frq4, rad4));   % AIRS CrIS
bt5 = real(rad2bt(frq1, rad5));   % interp AIRS
bt6 = real(rad2bt(frq6, rad6));   % interp conv AIRS

figure(1); clf; j = 1; 
plot(frq1, bt1(:,j), frq2, bt2(:,j), bfrq, bt3(:,j), frq4, bt4(:,j))
ax(1)=tv1-20; ax(2)=tv2+20; ax(3)=180; ax(4)=320; axis(ax)
legend('true CrIS', 'true AIRS', 'AIRS dec', 'AIRS CrIS', ...
       'location', 'southeast')
xlabel('wavenumber'); ylabel('brighness temp')
title(sprintf('AIRS 1C and CrIS %s profile %d', band, j));
grid on; zoom on
saveas(gcf, sprintf('test4_fig_1_%s', band), fig)

% residuals for real CrIS and AIRS CrIS
figure(2); clf
[i1, i4] = seq_match(frq1, frq4);
plot(frq1(i1), mean(bt4(i4,:) - bt1(i1,:), 2))
ax(1)=tv1; ax(2)=tv2; ax(3)=-dt1; ax(4)=dt1; axis(ax)
xlabel('wavenumber'); ylabel('dBT')
title(sprintf('AIRS CrIS minus true CrIS %s mean', band));
grid on; zoom on
saveas(gcf, sprintf('test4_fig_2_%s', band), fig)

% residuals for real CrIS and interpolated AIRS
figure(3); clf
plot(frq1, mean(bt5 - bt1, 2))
ax(1)=tv1; ax(2)=tv2; ax(3)=-dt2; ax(4)=dt2; axis(ax)
xlabel('wavenumber'); ylabel('dBT')
title(sprintf('interpolated CrIS minus true CrIS %s mean', band));
grid on; zoom on
saveas(gcf, sprintf('test4_fig_3_%s', band), fig)

% residuals for real CrIS and interpolated convolved AIRS
[j1, j6] = seq_match(frq1, frq6);
figure(4); clf
plot(frq1(j1), mean(bt6(j6,:) - bt1(j1,:), 2))
% plot(frq1(j1), mean(bt6(j6,:) - bt1(j1,:), 2), frq1, 0, '+')
ax(1)=tv1; ax(2)=tv2; ax(3)=-dt2; ax(4)=dt2; axis(ax)
xlabel('wavenumber'); ylabel('dBT')
title(sprintf('interpolated convolved CrIS minus true CrIS %s mean', band));
grid on; zoom on
saveas(gcf, sprintf('test4_fig_4_%s', band), fig)

