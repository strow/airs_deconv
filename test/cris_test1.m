% 
% cris_test1
% 
% reference truth: start with kcarta radiances, convolve these to 
% CrIS radiances at the user grid, and call the result “true CrIS”.
% 
% deconvolution: start with kcarta radiances, convolve to AIRS
% channel radiances (“true AIRS”), deconvolve to an intermediate
% grid, e.g. 0.05 1/cm spacing, and reconvolve to the CrIS user grid
% (“AIRS Cris”).  Then compare AIRS CrIS vs true CrIS for various
% profiles and any of the three CrIS bands

% use  my bcast utils and HDF libs
addpath /home/motteler/cris/bcast/source
addpath /home/motteler/cris/bcast/motmsc/utils
addpath /home/motteler/mot2008/hdf/h4tools

% bconv temp file, saves binv, bconv, bfreq, and fx, and depends
% on sfile, cfreq, dvb, so delete it delete if any of these change
bfile = 'bconv1.mat';

% plot type
fig = 'fig';   

% load kcarta test radiance
kcfile = '/home/motteler/cris/sergio/JUNK2012/convolved_kcarta1.mat';
d1 = load(kcfile);
vkc = d1.w;  % frequency grid
rkc = d1.r;  % radiances
clear d1

% get CrIS inst and user params
band = 'LW';   % cris band
wlaser = 773.1301;  % real value
[inst, user] = inst_params(band, wlaser);

% get CrIS channel radiances for rkc
[rad1, frq1] = kc2cris(user, rkc, vkc);

% get the kcarta to AIRS convolution matrix
sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';
cfreq = load('freq2645.txt');
cfreq = trim_chans(cfreq);
dvs = 0.0025; 
[sconv, sfreq, ofreq] = mksconv2(sfile, cfreq, dvs);

% apply the AIRS convolution
ix = interp1(vkc, 1:length(rkc), sfreq, 'nearest');
rad2 = sconv * rkc(ix);
frq2 = ofreq;

% deconvolve the AIRS radiances
% dvb = 0.05;
dvb = 0.1;
% [bconv, bfreq, fx] = mksconv2(sfile, cfreq, dvb);
% tic, rad3 = pinv(full(bconv)) * rad2; toc
if exist(bfile) == 2
  fprintf(1, 'loading %s...\n', bfile)
  load(bfile)
else
  [bconv, bfreq, fx] = mksconv2(sfile, cfreq, dvb);
  binv = pinv(full(bconv));
  save(bfile, 'binv', 'bconv', 'bfreq', 'fx')
end
rad3 = binv * rad2;
clear binv;

% try an extra smoothing step
% rad3 = mkhamm(length(rad3)) * rad3;

% reconvolve to CrIS
[rad4, frq4] = finterp(rad3, bfreq, user.dv);
frq4 = frq4(:);

% AIRS direct interpolation to CrIS
rad5 = interp1(frq2, rad2, frq1, 'spline', 'extrap');

% AIRS interpolation and convolution to CrIS
vx = ceil(frq2(1)/dvb) * dvb;
nx = floor((frq2(end) - frq2(1)) / dvb);
ftmp = vx + (0 : nx - 1) * dvb;
rtmp = interp1(frq2, rad2, ftmp', 'spline', 'extrap');
[rad6, frq6] = finterp(rtmp, ftmp, user.dv);
frq6 = frq6(:);

% compare simulated AIRS 1C and CrIS data
bt1 = real(rad2bt(frq1, rad1));   % true cris
bt2 = real(rad2bt(frq2, rad2));   % true airs
bt3 = real(rad2bt(bfreq, rad3));  % deconvolved airs
bt4 = real(rad2bt(frq4, rad4));   % airs cris
bt5 = real(rad2bt(frq1, rad5));   % interp cris 1
bt6 = real(rad2bt(frq6, rad6));   % interp cris 2

figure(1); clf
plot(frq1, bt1, frq2, bt2, bfreq, bt3, frq4, bt4)
ax(1)=user.v1-50; ax(2)=user.v2+50; ax(3)=180; ax(4)=320; axis(ax)
legend('true CrIS', 'true AIRS', 'AIRS dec', 'AIRS CrIS', ...
       'location', 'southeast')
xlabel('wavenumber'); ylabel('brighness temp')
title(sprintf('AIRS 1C and CrIS %s comparison', band));
grid on; zoom on
% saveas(gcf, sprintf('decon_fig_1_%s', band), fig)

% residuals for real CrIS and AIRS CrIS
[i1, i4] = seq_match(frq1, frq4);
figure(2); clf
plot(frq1(i1), bt4(i4) - bt1(i1))
ax(1)=user.v1; ax(2)=user.v2; ax(3)=-2; ax(4)=2; axis(ax)
xlabel('wavenumber'); ylabel('dBT')
title(sprintf('AIRS CrIS minus true CrIS %s comparison', band));
grid on; zoom on
% saveas(gcf, sprintf('decon_fig_2_%s', band), fig)

% residuals for real CrIS and interpolated CrIS
figure(3); clf
plot(frq1, bt5 - bt1)
ax(1)=user.v1; ax(2)=user.v2; ax(3)=-8; ax(4)=8; axis(ax)
xlabel('wavenumber'); ylabel('dBT')
title(sprintf('interpolated CrIS minus true CrIS %s comparison', band));
grid on; zoom on
% saveas(gcf, sprintf('decon_fig_3_%s', band), fig)

% residuals for real CrIS and AIRS CrIS
[j1, j6] = seq_match(frq1, frq6);
figure(4); clf
plot(frq1(j1), bt6(j6) - bt1(j1))
ax(1)=user.v1; ax(2)=user.v2; ax(3)=-8; ax(4)=8; axis(ax)
xlabel('wavenumber'); ylabel('dBT')
title(sprintf('interp-conv CrIS minus true CrIS %s comparison', band));
grid on; zoom on
% saveas(gcf, sprintf('decon_fig_4_%s', band), fig)

