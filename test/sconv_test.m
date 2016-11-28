%
% sconv_test - compare mksconv with mksconv1 and mksconv2
%

% set paths to libs
addpath ../source
addpath ../h4tools
addpath /asl/packages/ccast/source

% turn off HDF 4 update warnings
warning('off', 'MATLAB:imagesci:hdf:removalWarningHDFSD')

% specify SRF tabulations
sdir = '/asl/matlab2012/srftest/';
  srf1 = fullfile(sdir, 'srftables_m130f_withfake_mar08.hdf');
  srf2 = fullfile(sdir, 'srftables_m140f_withfake_mar08.hdf');
  srf3 = fullfile(sdir, 'srftables_m150f_withfake_mar08.hdf');

% read the SRF files
[tf1, tg1, sv1, id1] = srf_read(srf1);
[tf2, tg2, sv2, id2] = srf_read(srf2);
[tf3, tg3, sv3, id3] = srf_read(srf3);
[tf1, if1] = sort(tf1);
[tf2, if2] = sort(tf2);
[tf3, if3] = sort(tf3);
[isequal(if1, if2), isequal(if2, if3)]

% read the JPL L1C channel set
cfreq = load('freq2645.txt');
cfreq = trim_chans(cfreq);

% match L1C and SRF channel sets
[ix1, jx1] = seq_match(tf1, cfreq, 0.04);
tf1 = tf1(ix1);  if1 = if1(ix1);  

[ix2, jx2] = seq_match(tf2, cfreq, 0.04);
tf2 = tf2(ix2);  if2 = if2(ix2);

[ix3, jx3] = seq_match(tf3, cfreq, 0.04);
tf3 = tf3(ix3);  if3 = if3(ix3);
[isequal(jx1, jx2), isequal(jx2, jx3)]

cfreq = cfreq(jx3);

dvk = 0.0025; 
[S1, sfrq1, tfrq1] = mksconv2(srf1, cfreq, dvk);
[S2, sfrq2, tfrq2] = mksconv2(srf2, cfreq, dvk);
[S3, sfrq3, tfrq3] = mksconv2(srf3, cfreq, dvk);
[isclose(tf1, tfrq1), isclose(tf2, tfrq2), isclose(tf3, tfrq3)]

[S1x, sfrq1x, tfrq1x] = mksconv(srf1, if1, dvk);
[S2x, sfrq2x, tfrq2x] = mksconv(srf2, if2, dvk);
[S3x, sfrq3x, tfrq3x] = mksconv(srf3, if3, dvk);
[isclose(tf1, tfrq1x), isclose(tf2, tfrq2x), isclose(tf3, tfrq3x)]

% compare the two convolvers
[isclose(sfrq1, sfrq1x), isclose(tfrq1, tfrq1x), isclose(S1, S1x)]
[isclose(sfrq2, sfrq2x), isclose(tfrq2, tfrq2x), isclose(S2, S2x)]
[isclose(sfrq3, sfrq3x), isclose(tfrq3, tfrq3x), isclose(S3, S3x)]

