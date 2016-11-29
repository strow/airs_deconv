%
% jpl_coeffs -- get new "a" and "b" coeffs for jpl_shift2
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
% srf3 = fullfile(sdir, 'srftables_m150f_withfake_mar08.hdf');

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
dvk = 0.0025; 
[S1, sfS1, tfS1] = mksconv(srf1, if1, dvk);
[S2, sfS2, tfS2] = mksconv(srf2, if2, dvk);

% loop on kcarta files
rad1 = []; rad2 = []; rad3 = []; rad4 = []; rad5 = [];
for i = 1 : length(flist)
  d1 = load(fullfile(kcdir, flist(i).name));
  vkc = d1.w(:); rkc = d1.r(:);

  % apply the S1 and S2 convolutions
  ix = interp1(vkc, 1:length(rkc), sfS1, 'nearest');
  r1 = S1 * rkc(ix);  rad1 = [rad1, r1];

  ix = interp1(vkc, 1:length(rkc), sfS2, 'nearest');
  r2 = S2 * rkc(ix);  rad2 = [rad2, r2];

  % spline shift from tf1 to tf2
  r3 = interp1(tf1, r1, tf2, 'spline');  
  rad3 = [rad3, r3];

  fprintf(1, '.');
end
fprintf(1, '\n')
frq1 = tf1(:);
frq2 = tf2(:);
clear d1

% take radiances to brightness temps
bt1 = real(rad2bt(frq1, rad1));   % SRF set 1
bt2 = real(rad2bt(frq2, rad2));   % SRF set 2
bt3 = real(rad2bt(frq2, rad3));   % spline shift of 1 to 2

% regress for "a" and "b" weights
dv = frq2 - frq1;
[m, n] = size(bt1);
a = zeros(m, 1);
b = zeros(m, 1);

% loop on channels
for i = 1 : m
 X = [(bt3(i, :)' - bt1(i, :)') / dv(i), ones(n,1)];
 Z =  (bt2(i, :)' - bt1(i, :)') / dv(i);
 % regress for X * Y = Z
 Y = X \ Z;
 a(i) = Y(1);
 b(i) = Y(2);
end

save jpl_shift2 frq1 frq2 a b

