%
% regr_test1 -- chan shift regression with a large dependent set
%

% set paths to libs
addpath ../source
addpath ../h4tools
addpath /asl/packages/ccast/source
addpath /home/motteler/matlab/export_fig

% turn off HDF 4 update warnings
warning('off', 'MATLAB:imagesci:hdf:removalWarningHDFSD')

% kcarta dependent set
kcd1 = '/asl/s1/motteler/kc7377/cloudy';
nkc1 = 7377;

% kcarta independent set
kcd2 = '/home/motteler/cris/sergio/JUNK2012';
nkc2 = 49;

% specify SRF tabulations
sdir = '/asl/matlab2012/srftest/';
  f1srf = fullfile(sdir, 'srftables_m130f_withfake_mar08.hdf');
  f2srf = fullfile(sdir, 'srftables_m140f_withfake_mar08.hdf');
% f3srf = fullfile(sdir, 'srftables_m150f_withfake_mar08.hdf');

% get SRF channel frequencies
v1srf = srf_read(f1srf);

% use the L1b subset
% v1srf = v1srf(1:2378);

% sort the SRF channel set
[~, i1srf] = sort(v1srf); 
v1srf = v1srf(i1srf);

% get the JPL L1C channel set
c2645 = load('freq2645.txt');

% match L1C and SRF channel sets
[ix, ~] = seq_match(v1srf, c2645, 0.04);
i1srf = i1srf(ix);  

% convolution matrices for reference truth
dvk = 0.0025; 
[C1, v1col, v1row] = mksconv(f1srf, i1srf, dvk);
[C2, v2col, v2row] = mksconv(f2srf, i1srf, dvk);

nchan = length(v1row);
rdep1 = zeros(nchan, nkc1);
rdep2 = zeros(nchan, nkc1);
rind1 = zeros(nchan, nkc2);
rind2 = zeros(nchan, nkc2);

% loop on dependent set kcarta files
tic
for i = 1 : nkc1

  kcmat = fullfile(kcd1, sprintf('kc%04d.mat', i));
  d1 = load(kcmat);
  rkc = d1.rad; vkc = d1.frq; clear d1

  ix = interp1(vkc, 1:length(rkc), v1col, 'nearest');
  rdep1(:, i) = C1 * rkc(ix);

  jx = interp1(vkc, 1:length(rkc), v2col, 'nearest');
  rdep2(:, i) = C2 * rkc(jx);

  if mod(i, 100) == 0, fprintf(1, '.'), end
end
fprintf(1, '\n'); toc

% loop on independent set kcarta files
tic
for i = 1 : nkc2

  kcmat = fullfile(kcd2, sprintf('convolved_kcarta%d.mat', i));
  d1 = load(kcmat);
  rkc = d1.r; vkc = d1.w; clear d1

  ix = interp1(vkc, 1:length(rkc), v1col, 'nearest');
  rind1(:, i) = C1 * rkc(ix);

  jx = interp1(vkc, 1:length(rkc), v2col, 'nearest');
  rind2(:, i) = C2 * rkc(jx);

  fprintf(1, '.');
end
fprintf(1, '\n'); toc

% save regr_test1 rdep1 rdep2 rind1 rind2

R12 = (rdep1' \ rdep2')';

rind3 = R12 * rind1;

% get brightness temp
bt1 = real(rad2bt(v1row, rind1));   
bt2 = real(rad2bt(v2row, rind2));
bt3 = real(rad2bt(v2row, rind3));

% get residual stats
mdif = mean(bt3 - bt2, 2);
sdif = std(bt3 - bt2, 0, 2);

% plot residuals
figure(1); clf; 
set(gcf, 'Units','centimeters', 'Position', [4, 10, 24, 16])
subplot(2,1,1)
plot(v2row, mdif);
axis([600, 2700, -0.01, 0.01])
title('shift minus ref, 49 profile mean');
ylabel('dTb')
grid on; zoom on

subplot(2,1,2)
plot(v2row, sdif)
axis([600, 2700, 0, 4e-3])
title('shift minus ref, 49 profile std');
xlabel('wavenumber'); 
ylabel('dTb')
grid on; zoom on

% saveas(gcf, 'regr_shift_L1c', 'png')
