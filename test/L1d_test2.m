%
% L1d_test2 - try deconvolution to the L1d SRF set
%
%  calculate true L1c and L1d AIRS, convolved directly from kcarta,
%  and L1c to L1d AIRS "C to D" AIRS (via de- and reconvolution) and
%  compare them
%

% set paths to libs
addpath ../source
addpath ../h4tools
addpath /asl/packages/ccast/source

% L1d resolution, dv = v / res
% res = 1200;  % L1c nominal
  res =  700;  % best for L1d

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
[C1, vCcol, v_L1c] = mksconv2(srf1, v_tmp, dvk);
[A1, vAcol, vArow] = mksconv2(srf1, v_tmp, dvb);

% L1d convolution matrices
[D1, vDcol, v_L1d] = L1d_conv(res, dvk);
[B1, vBcol, VBrow] = L1d_conv(res, dvb);

% match intermediate grids
[ix, jx] = seq_match(vAcol, vBcol);
vAcol = vAcol(ix);  A1 = A1(:, ix);
vBcol = vBcol(jx);  B1 = B1(:, jx);

% build the decon/recon matrix
fprintf(1, 'inverting A ...\n')
tic
CtoD = B1 * pinv(full(A1));
toc

% loop on kcarta files
trueCrad = []; trueDrad = [];
for i = 1 : length(flist)
  d1 = load(fullfile(kcdir, flist(i).name));
  vkc = d1.w(:); rkc = d1.r(:);

  % apply the S1 and D1 convolutions
  ix = interp1(vkc, 1:length(rkc), vCcol, 'nearest');
  r1 = C1 * rkc(ix);  trueCrad = [trueCrad, r1];

  ix = interp1(vkc, 1:length(rkc), vDcol, 'nearest');
  r2 = D1 * rkc(ix);  trueDrad = [trueDrad, r2];

  fprintf(1, '.');
end
fprintf(1, '\n')

% de/reconvolve L1C to L1d
CtoDrad = CtoD * trueCrad;

% take radiances to brightness temps
trueCbt = real(rad2bt(v_L1c, trueCrad));
trueDbt = real(rad2bt(v_L1d, trueDrad));
CtoDbt  = real(rad2bt(v_L1d, CtoDrad));

% plot residuals
figure(1); clf; 
% set(gcf, 'Units','centimeters', 'Position', [4, 10, 24, 16])
subplot(2,1,1)
plot(v_L1d, mean(CtoDbt - trueDbt, 2))
  axis([650, 2700, -0.3, 0.3])
% axis([650, 2700, -1, 1])
% axis([650, 2700, -2, 2])
title('C to D minus true L1d , 49 profile mean');
ylabel('dTb')
grid on; zoom on

subplot(2,1,2)
plot(v_L1d, std(CtoDbt - trueDbt, 0, 2))
  axis([650, 2700, 0, 0.1])
% axis([650, 2700, 0, 0.3])
% axis([650, 2700, 0, 0.6])
title('C to D minus true L1d , 49 profile std');
ylabel('dTb')
xlabel('wavenumber')
grid on; zoom on

