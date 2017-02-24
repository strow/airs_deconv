%
% jpl_test -- test JPL L1b channel shift
%

% set paths to libs
addpath ./data
addpath ../source
addpath ../h4tools
addpath /asl/packages/ccast/source

% turn off HDF 4 update warnings
warning('off', 'MATLAB:imagesci:hdf:removalWarningHDFSD')

% chose shift coefficients
% coeff = 'jpl_1b';
  coeff = 'umbc_1b';

% kcarta test data
kcdir = '/home/motteler/cris/sergio/JUNK2012/';
flist =  dir(fullfile(kcdir, 'convolved_kcart*.mat'));

% specify SRF tabulations
sdir = '/asl/matlab2012/srftest/';
  srf1 = fullfile(sdir, 'srftables_m130f_withfake_mar08.hdf');
  srf2 = fullfile(sdir, 'srftables_m140f_withfake_mar08.hdf');
% srf2 = fullfile(sdir, 'srftables_m150f_withfake_mar08.hdf');

% get channel frequencies
tf1 = srf_read(srf1);
tf2 = srf_read(srf2);

% use the 1b channel set
  [~, i1b] = sort(tf1(1:2378));
% [~, i1b] = trim_chans(tf1(1:2378));
tf1 = tf1(i1b);
tf2 = tf2(i1b);

% convolution matrices for reference truth
dvk = 0.0025; 
[S1, sfS1, tfS1] = mksconv(srf1, i1b, dvk);
[S2, sfS2, tfS2] = mksconv(srf2, i1b, dvk);

% loop on kcarta files
rad1 = []; rad2 = []; rad3 = []; rad4 = [];
for i = 1 : length(flist)
  d1 = load(fullfile(kcdir, flist(i).name));
  vkc = d1.w(:); rkc = d1.r(:);

  % apply the S1 and S2 convolutions
  ix = interp1(vkc, 1:length(rkc), sfS1, 'nearest');
  r1 = S1 * rkc(ix);  rad1 = [rad1, r1];

  ix = interp1(vkc, 1:length(rkc), sfS2, 'nearest');
  r2 = S2 * rkc(ix);  rad2 = [rad2, r2];

  % apply the JPL shift 
  tb1 = real(rad2bt(tf1, r1));
  tb3 = jpl_shift(tb1, tf1, tf2, coeff, i1b);
  r3 = bt2rad(tf2, tb3);
  rad3 = [rad3, r3];

  % try a simple spline shift
  r4 = interp1(tf1, r1, tf2, 'spline');  
  rad4 = [rad4, r4];

  fprintf(1, '.');
end
fprintf(1, '\n')
frq1 = tf1(:);
frq2 = tf2(:);
clear d1

% take radiances to brightness temps
bt1 = real(rad2bt(frq1, rad1));   % SRF set 1
bt2 = real(rad2bt(frq2, rad2));   % SRF set 2
bt3 = real(rad2bt(frq2, rad3));   % 1 shifted to 2
bt4 = real(rad2bt(frq2, rad4));   % 1 interpolated to 2

% jpl shift residuals
figure(1); clf; 
% set(gcf, 'Units','centimeters', 'Position', [4, 10, 24, 16])
subplot(2,1,1)
plot(frq2, mean(bt3 - bt2, 2))
axis([600, 2700, -0.1, 0.1])
title('JPL shift minus ref, 49 profile mean');
ylabel('dTb')
grid on; zoom on

subplot(2,1,2)
plot(frq2, std(bt3 - bt2, 0, 2))
axis([600, 2700, 0, 0.12])
title('JPL shift minus ref, 49 profile std');
xlabel('wavenumber'); 
ylabel('dTb')
grid on; zoom on

return

% spline fit residuals
figure(2); clf; 
% set(gcf, 'Units','centimeters', 'Position', [4, 10, 24, 16])
subplot(2,1,1)
plot(frq2, mean(bt4 - bt2, 2))
axis([600, 2700, -0.1, 0.1])
title('spline minus ref, 49 profile mean');
ylabel('dTb')
grid on; zoom on

subplot(2,1,2)
plot(frq2, std(bt4 - bt2, 0, 2))
axis([600, 2700, 0, 0.12])
title('spline minus ref, 49 profile std');
xlabel('wavenumber'); 
ylabel('dTb')
grid on; zoom on

