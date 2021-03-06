%
% L1d_test3 - compare L1c to L1d de/recon vs regression
%
% test direct regression on both radiance and brightness temps 
% vs the decon/recon CtoD linear transform
%
% transform matrices
%   C1 takes kcarta to the L1c channel grid v_L1c
%   A1 takes the 0.1 cm-1 grid to L1c grid vArow
%   D1 takes kcarta to the L1d channel grid v_L1d
%   B1 takes the 0.1 cm-1 grid to L1d grid vBrow
% 
% radiance sets
%   true AIRS L1c, rc = C1 * rkc
%   true AIRS L1d, rd = D1 * rkc
%   CtoD AIRS L1d, rd = B1 * inv(A1) * rc
%
% then we compare true AIRS L1d with CtoD AIRS L1d
%

% set paths to libs
addpath ../source
addpath ../h4tools
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
  v_base = v_L1c(1);  % 649.6620
% v_base = 649.6620;  % best for res 1200
% v_base = 649.8220;  % best for res 700

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

% solve X * C = D as  C' * X' = D',  X' = C' \ D';
% ix = 1:49; jx = 1:49; 
  ix = 1:2:49; jx = 2:2:49; 
% ix = 1:24; jx = 25:49; 
X = (trueCrad(:,ix)' \ trueDrad(:,ix)')';
Cr1Drad = X * trueCrad;

% % test of "\" independent columns
% [m, n] = size(X);
% Z = zeros(m, n);
% for i = 1 : m
%   Z(i, :) = trueCrad(:,ix)' \ trueDrad(i,ix)';
% end
% rms(Z(:) - X(:))

% take radiances to brightness temps
trueCbt = real(rad2bt(v_L1c, trueCrad));
trueDbt = real(rad2bt(v_L1d, trueDrad));
CtoDbt = real(rad2bt(v_L1d, CtoDrad));
Cr1Dbt = real(rad2bt(v_L1d, Cr1Drad));

% solve Y * trueCbt = trueDbt
% ix = 2:2:49; jx = 1:2:49; 
Y = (trueCbt(:,ix)' \ trueDbt(:,ix)')';
Cr2Dbt = Y * trueCbt;

% summary decon and interp residuals
figure(2); clf; 
% set(gcf, 'Units','centimeters', 'Position', [4, 10, 24, 16])
subplot(2,1,1)
x0 = mean(CtoDbt - trueDbt, 2);
x1 = mean(Cr1Dbt(:,jx) - trueDbt(:,jx), 2);
x2 = mean(Cr2Dbt(:,jx) - trueDbt(:,jx), 2);
plot(v_L1d, x0, v_L1d, x1, v_L1d, x2)
% axis([650, 2700, -0.3, 0.3])
  axis([650, 2700, -1, 1])
title('L1c to L1d residual mean');
legend('decon/recon', 'rad regression', 'bt regression', ...
       'location', 'northeast');
ylabel('dTb')
grid on; zoom on

subplot(2,1,2)
x0 = std(CtoDbt - trueDbt, 0, 2);
x1 = std(Cr1Dbt(:,jx) - trueDbt(:,jx), 0, 2);
x2 = std(Cr2Dbt(:,jx) - trueDbt(:,jx), 0, 2);
plot(v_L1d, x0, v_L1d, x1, v_L1d, x2)
% axis([650, 2700, 0, 0.1])
  axis([650, 2700, 0, 1])
title('L1c to L1d residual std dev');
legend('decon/recon', 'rad regression', 'bt regression', ...
       'location', 'northeast');
ylabel('dTb')
xlabel('wavenumber')
grid on; zoom on

return

% plot decon residuals
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
% export_fig('CtoD_decon_diff.pdf', '-m2', '-transparent')

