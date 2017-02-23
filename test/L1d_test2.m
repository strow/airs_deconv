%
% L1d_test2 - try deconvolution to the L1d SRF set
%
%  L1d is used here as a new reference truth.  We calculate true
%  L1d AIRS, convolved directly from kcarta, and L1c to L1d AIRS
%  (via de- and reconvolution) and compare them
%

% set paths to libs
addpath ../source
addpath ../h4tools
addpath /asl/packages/ccast/source

% L1d resolution scaling factor
s_fact = 1;      

% turn off HDF 4 update warnings
warning('off', 'MATLAB:imagesci:hdf:removalWarningHDFSD')

% kcarta test data
kcdir = '/home/motteler/cris/sergio/JUNK2012/';
flist =  dir(fullfile(kcdir, 'convolved_kcart*.mat'));

% get the kcarta to L1c AIRS convolution matrix
sdir = '/asl/matlab2012/srftest/';
srf1 = fullfile(sdir, 'srftables_m140f_withfake_mar08.hdf');

v_tmp = load('freq2645.txt');
dvk = 0.0025;    % kcarta dv
dvb = 0.1;       % decon dv

% L1c convolution matrices
[S1, vScol, v_L1c] = mksconv2(srf1, v_tmp, dvk);
[A1, vAcol, vArow] = mksconv2(srf1, v_tmp, dvb);

% L1d convolution matrices
[D1, vDcol, v_L1d] = L1d_conv(s_fact, dvk);
[B1, vBcol, VBrow] = L1d_conv(s_fact, dvb);

% match intermediate grids
[ix, jx] = seq_match(vAcol, vBcol);
vAcol = vAcol(ix);  A1 = A1(:, ix);
vBcol = vBcol(jx);  B1 = B1(:, jx);

% build the decon/recon matrix
fprintf(1, 'inverting A ...\n')
tic
CtoD = B1 * pinv(full(A1));
toc

return

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

  % apply the SRF shift transform
  r3 = Bshift * r1;  rad3 = [rad3, r3];

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
axis([600, 2700, -0.12, 0.12])
ylabel('dTb')
title('decon minus ref, 49 profile mean');
grid on; zoom on

subplot(2,1,2)
plot(frq2, mean(bt4 - bt2, 2))
axis([600, 2700, -0.12, 0.12])
xlabel('wavenumber'); 
ylabel('dTb')
title('spline minus ref, 49 profile mean');
grid on; zoom on

return

% single profile mean difference
figure(1); clf; 
j = 1; 
% set(gcf, 'Units','centimeters', 'Position', [4, 10, 24, 16])
subplot(2,1,1)
plot(frq2, bt3(:, j) - bt2(:,j))
axis([600, 2700, -0.12, 0.12])
ylabel('dTb')
title(sprintf('decon minus ref, profile %d', j));
grid on; zoom on

subplot(2,1,2)
plot(frq2, bt4(:, j) - bt2(:,j))
axis([600, 2700, -0.12, 0.12])
xlabel('wavenumber'); 
ylabel('dTb')
title(sprintf('spline minus ref, profile %d', j));
grid on; zoom on

