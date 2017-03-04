%
% jpl_coeffs -- get new "a" and "b" coeffs for jpl_shift
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

% get SRF channel frequencies
tf1 = srf_read(srf1);
tf2 = srf_read(srf2);

% use the L1b subset
  tf1 = tf1(1:2378);
  tf2 = tf2(1:2378);

% sort the full SRF channel set
[~, i1c] = sort(tf1); 
tf1 = tf1(i1c);
tf2 = tf2(i1c);

% get the JPL L1C channel set
% c2645 = load('freq2645.txt');

% match L1C and SRF channel sets
% [ix, jind] = seq_match(tf1, c2645, 0.04);
% i1c = i1c(ix);  
% tf1 = tf1(ix);  
% tf2 = tf2(ix);

% convolution matrices for reference truth
dvk = 0.0025; 
[S1, sfS1, tfS1] = mksconv(srf1, i1c, dvk);
[S2, sfS2, tfS2] = mksconv(srf2, i1c, dvk);

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

% save shift_L1c frq1 frq2 a b
  save shift_L1b frq1 frq2 a b

