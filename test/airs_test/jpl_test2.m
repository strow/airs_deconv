%
% jpl_test2c -- test UMBC 1c version of the jpl shift
%

% set paths to libs
addpath ./data
addpath ../source
addpath ../h4tools
addpath /asl/packages/ccast/source

% turn off HDF 4 update warnings
warning('off', 'MATLAB:imagesci:hdf:removalWarningHDFSD')

% kcarta test data
kcdir = '/home/motteler/cris/sergio/JUNK2012/';
flist =  dir(fullfile(kcdir, 'convolved_kcart*.mat'));

% specify SRF tabulations
sdir = '/asl/matlab2012/srftest/';
srf1 = fullfile(sdir, 'srftables_m130f_withfake_mar08.hdf');
srf2 = fullfile(sdir, 'srftables_m140f_withfake_mar08.hdf');

% get SRF channel frequencies
tf1 = srf_read(srf1);
tf2 = srf_read(srf2);

% use the L1b subset
% tf1 = tf1(1:2378);
% tf2 = tf2(1:2378);

% sort the full SRF channel set
[~, i1c] = sort(tf1); 
tf1 = tf1(i1c);
tf2 = tf2(i1c);

% get the JPL L1C channel set
c2645 = load('freq2645.txt');

% match L1C and SRF channel sets
[ix, jind] = seq_match(tf1, c2645, 0.04);
i1c = i1c(ix);  
tf1 = tf1(ix);  
tf2 = tf2(ix);

% convolution matrices for reference truth
dvk = 0.0025; 
[S1, sfS1, tfS1] = mksconv(srf1, i1c, dvk);
[S2, sfS2, tfS2] = mksconv(srf2, i1c, dvk);

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
  tb3 = jpl_shift2(tb1, tf1, tf2, jind);
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

% 49 profile mean difference
figure(1); clf; 
% set(gcf, 'Units','centimeters', 'Position', [4, 10, 24, 16])
subplot(2,1,1)
plot(frq2, mean(bt3 - bt2, 2))
axis([600, 2700, -0.1, 0.1])
ylabel('dTb')
title('JPL shift minus ref, 49 profile mean');
grid on; zoom on

subplot(2,1,2)
plot(frq2, mean(bt4 - bt2, 2))
axis([600, 2700, -0.1, 0.1])
xlabel('wavenumber'); 
ylabel('dTb')
title('spline minus ref, 49 profile mean');
grid on; zoom on

