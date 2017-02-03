%
% airs_diffs -- first and second differences of SRF shifts
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
  srf3 = fullfile(sdir, 'srftables_m150f_withfake_mar08.hdf');

% read the SRF files
[tf1, tg1, sv1, id1] = srf_read(srf1);
[tf2, tg2, sv2, id2] = srf_read(srf2);
[tf3, tg3, sv3, id3] = srf_read(srf2);
[tf1, if1] = sort(tf1); 
[tf2, if2] = sort(tf2);
[tf3, if3] = sort(tf3);
[isequal(if1, if2), isequal(if2, if3)]

% use the JPL L1C channel set
cfreq = load('freq2645.txt');

% use the L1B channel set
% cfreq = srf_read(srf1);
% cfreq = trim_chans(cfreq(1:2378));

% match L1C and SRF channel sets
[ix1, jx1] = seq_match(tf1, cfreq, 0.04);
tf1 = tf1(ix1);  if1 = if1(ix1);  
tf2 = tf2(ix1);  if2 = if2(ix1);
tf3 = tf3(ix1);  if3 = if3(ix1);
cfreq = cfreq(jx1);

% convolution matrices for reference truth
dvk = 0.0025; 
[S1, sfS1, tfS1] = mksconv(srf1, if1, dvk);
[S2, sfS2, tfS2] = mksconv(srf2, if2, dvk);
[S3, sfS3, tfS3] = mksconv(srf3, if3, dvk);

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

  ix = interp1(vkc, 1:length(rkc), sfS3, 'nearest');
  r3 = S3 * rkc(ix);  rad3 = [rad3, r3];

  % spline shift from tf1 to tf2
  r4 = interp1(tf1, r1, tf2, 'spline');  
  rad4 = [rad4, r4];

  % spline shift from tf1 to tf3
  r5 = interp1(tf1, r1, tf3, 'spline');  
  rad5 = [rad5, r5];

  fprintf(1, '.');
end
fprintf(1, '\n')
frq1 = tf1(:);
frq2 = tf2(:);
frq3 = tf3(:);
clear d1

% take radiances to brightness temps
bt1 = real(rad2bt(frq1, rad1));   % SRF set 1
bt2 = real(rad2bt(frq2, rad2));   % SRF set 2
bt3 = real(rad2bt(frq3, rad3));   % SRF set 3
bt4 = real(rad2bt(frq2, rad4));   % spline shift of 1 to 2
bt5 = real(rad2bt(frq3, rad5));   % spline shift of 1 to 2

% dbt1 = (bt2 - bt1) ./ ((frq2 - frq1) * ones(1, 49));
% dbt2 = (bt3 - bt2) ./ ((frq3 - frq2) * ones(1, 49));
dbt1 = bt2 - bt1;
dbt2 = bt3 - bt2;
ddbt = dbt2 - dbt1;

dsp1 = bt4 - bt1;
dsp2 = bt5 - bt4;
ddsp = dsp2 - dsp1;

figure(1); clf; 
subplot(2,1,1)
plot(frq1, mean(dbt1, 2))
axis([600, 2700, -0.4, 0.4])
title('first difference mean')
ylabel('dTb')
grid on

subplot(2,1,2)
plot(frq1, mean(ddbt, 2))
axis([600, 2700, -0.01, 0.01])
title('second difference mean')
xlabel('wavenumber')
ylabel('ddTb')
grid on

figure(2); clf; 
subplot(2,1,1)
plot(frq1, std(dbt1, 0, 2))
axis([600, 2700, 0, 0.12])
title('first difference std')
ylabel('dTb')
grid on

subplot(2,1,2)
plot(frq1, mean(dsp1, 2))
axis([600, 2700, -0.4, 0.4])
title('spline difference mean')
xlabel('wavenumber')
ylabel('dTb')
grid on

figure(3); clf; 
subplot(2,1,1)
plot(frq1, mean(dsp1 - dbt1, 2))
axis([600, 2700, -0.1, 0.1])
title('spline diff minus first diff mean')
ylabel('ddTb')
grid on

subplot(2,1,2)
plot(frq1, std(dsp1 - dbt1, 0, 2))
axis([600, 2700, 0, 0.03])
title('spline diff minus first diff std')
xlabel('wavenumber')
ylabel('ddTb')
grid on

return

% spline diff minus first diff mean
serr =  mean(dsp1 - dbt1, 2);
dfix = dsp1 + serr * ones(1, 49);

figure(4); clf
subplot(2,1,1)
plot(frq1, mean(dbt1 - dfix, 2))

subplot(2,1,2)
plot(frq2, mean(bt1 + dfix - bt2, 2))


