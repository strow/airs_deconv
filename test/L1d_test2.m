%
% L1d_test2 - try deconvolution to the L1d SRF set
%
%   C1 takes kcarta to the L1c channel grid v_L1c
%   A1 takes the 0.1 cm-1 grid to L1c grid vArow
%   D1 takes kcarta to the L1d channel grid v_L1d
%   B1 takes the 0.1 cm-1 grid to L1d grid vBrow
% 
%   true AIRS L1c, rc = C1 * rkc
%   true AIRS L1d, rd = D1 * rkc
%   CtoD AIRS L1d, rd = B1 * inv(A1) * rc
%
% then we compare true AIRS L1d with CtoD AIRS L1d
%

% set paths to libs
addpath ../source
addpath /asl/packages/ccast/source
addpath /home/motteler/matlab/export_fig

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

% set start of output grid
% v_base = v_L1c(1);  % 649.6620
% v_base = 649.6620;  % best for res 1200
  v_base = 649.8220;  % best for res 700

% % temporary loop
% rmin = 1000;
% for v_base = v_L1c(1) : .02 :  v_L1c(2)

% L1d convolution matrices
[D1, vDcol, v_L1d] = L1d_conv(res, dvk, v_base);
[B1, vBcol, VBrow] = L1d_conv(res, dvb, v_base);

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

% direct spline interpolation from L1c to L1d
Ci1Drad = interp1(v_L1c, trueCrad, v_L1d, 'spline');

% spline to intermediate grid followed by conv to L1d
rtmp = interp1(v_L1c, trueCrad, vBcol, 'spline');
Ci2Drad = B1 * rtmp;

% take radiances to brightness temps
trueCbt = real(rad2bt(v_L1c, trueCrad));
trueDbt = real(rad2bt(v_L1d, trueDrad));
CtoDbt = real(rad2bt(v_L1d, CtoDrad));
Ci1Dbt = real(rad2bt(v_L1d, Ci1Drad));
Ci2Dbt = real(rad2bt(v_L1d, Ci2Drad));

rms(mean(CtoDbt - trueDbt, 2))

% rtmp = mean(CtoDbt - trueDbt, 2);
% ix = 650 < v_L1d & v_L1d < 1600;
%   rtmp = rms(rtmp(ix, :));
% % rtmp = max(abs(rtmp(ix, :)));
% if rtmp < rmin
%   rmin = rtmp;
%   vmin = v_base
% end
% end % temporary test loop
% 
% rmin, vmin

% plot decon residuals
figure(1); clf; 
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
% saveas(gcf, 'CtoD_decon_diff', 'fig')

% summary decon and interp residuals
figure(2); clf; 
subplot(2,1,1)
x1 = mean(Ci1Dbt - trueDbt, 2);
x2 = mean(Ci2Dbt - trueDbt, 2);
x3 = mean(CtoDbt - trueDbt, 2);
plot(v_L1d, x1, v_L1d, x2, v_L1d, x3) 
axis([650, 2700, -6, 6])
title('L1c to L1d residual mean');
legend('spline interpolation', 'L1c interp/L1d conv', ...
       'L1c decon/L1d conv', 'location', 'northeast');
ylabel('dTb')
grid on; zoom on

subplot(2,1,2)
x1 = std(Ci1Dbt - trueDbt, 0, 2);
x2 = std(Ci2Dbt - trueDbt, 0, 2);
x3 = std(CtoDbt - trueDbt, 0, 2);
plot(v_L1d, x1, v_L1d, x2, v_L1d, x3) 
axis([650, 2700, 0, 3])
title('L1c to L1d residual std dev');
legend('spline interpolation', 'L1c interp/L1d conv', ...
       'L1c decon/L1d conv', 'location', 'northeast');
ylabel('dTb')
xlabel('wavenumber')
grid on; zoom on
  saveas(gcf, 'CtoD_interp_diff', 'fig')

