%
% jpl_test3 -- test JPL shift with JPL 1c data
%

% set paths to libs
addpath ../source
addpath ../h4tools
addpath /asl/packages/ccast/source

% turn off HDF 4 update warnings
warning('off', 'MATLAB:imagesci:hdf:removalWarningHDFSD')

% load the JPL data
d1 = load('jpl_test_data');
[nchan, nobs] = size(d1.r130);

% specify SRF tabulations
sdir = '/asl/matlab2012/srftest/';
srf1 = fullfile(sdir, 'srftables_m130f_withfake_mar08.hdf');
srf2 = fullfile(sdir, 'srftables_m140f_withfake_mar08.hdf');

% read the SRF files
[tf1, tg1, sv1, id1] = srf_read(srf1);
[tf2, tg2, sv2, id2] = srf_read(srf2);
[tf1, if1] = sort(tf1); 
[tf2, if2] = sort(tf2);
% isequal(if1, if2)

% use the JPL L1C channel set
cfreq = load('freq2645.txt');

% use the L1B channel set
% cfreq = srf_read(srf1);
% cfreq = trim_chans(cfreq(1:2378));

% match L1C and SRF channel sets
[ix1, jx1] = seq_match(tf1, cfreq, 0.04);
tf1 = tf1(ix1);  if1 = if1(ix1);  
tf2 = tf2(ix1);  if2 = if2(ix1);
cfreq = cfreq(jx1);

% convolution matrices for reference truth
% dvk = 0.0025; 
% [S1, sfS1, tfS1] = mksconv(srf1, if1, dvk);
% [S2, sfS2, tfS2] = mksconv(srf2, if2, dvk);

% loop on JPL radiances
rad1 = []; rad2 = []; rad3 = []; rad4 = [];
for i = 1 : nobs

  r1 = d1.r130(:, i);
  rad1 = [rad1, r1];

  r2 = d1.r140(:, i);
  rad2 = [rad2, r2];

  % apply the JPL shift 
  tb1 = real(rad2bt(tf1, r1));
  tb3 = jpl_shift2(tb1, tf1, tf2);
  r3 = bt2rad(tf2, tb3);
  rad3 = [rad3, r3];

  % try a simple spline shift
  r4 = interp1(tf1, r1, tf2, 'spline');  
  rad4 = [rad4, r4];

end

frq1 = tf1(:);
frq2 = tf2(:);
clear d1

% take radiances to brightness temps
bt1 = real(rad2bt(frq1, rad1));   % SRF set 1
bt2 = real(rad2bt(frq2, rad2));   % SRF set 2
bt3 = real(rad2bt(frq2, rad3));   % 1 shifted to 2
bt4 = real(rad2bt(frq2, rad4));   % 1 interpolated to 2

% mean shift and spline differences
figure(1); clf; 
% set(gcf, 'Units','centimeters', 'Position', [4, 10, 24, 16])
subplot(2,1,1)
plot(frq2, mean(bt3 - bt2, 2))
% axis([600, 2700, -0.1, 0.1])
ylabel('dTb')
title('JPL shift minus ref, JPL test data');
grid on; zoom on

subplot(2,1,2)
plot(frq2, mean(bt4 - bt2, 2))
% axis([600, 2700, -0.1, 0.1])
xlabel('wavenumber'); 
ylabel('dTb')
title('spline minus ref, JPL test data');
grid on; zoom on

return

figure(2); clf; 
% set(gcf, 'Units','centimeters', 'Position', [4, 10, 24, 16])
% subplot(2,1,1)
plot(frq2, mean(bt2 - bt1, 2))
% axis([600, 2700, -0.1, 0.1])
ylabel('dTb')
title('JPL shift difference, JPL test data');
grid on; zoom on

