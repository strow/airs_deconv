%
% plot_Binv - plot sample deconvolution and C to D rows
%

% set paths to libs
addpath ../source
addpath ../h4tools
addpath /asl/packages/ccast/source
addpath /home/motteler/matlab/export_fig

% L1d resolution, dv = v / res
% res = 2000;  % for direct regression test
% res = 1200;  % L1c nominal
  res =  700;  % best for L1d

% set start of output grid
% v_base = v_L1c(1);  % default, 649.6620
% v_base = 649.6620;  % best for res 1200
  v_base = 649.8220;  % best for res 700

% turn off HDF 4 update warnings
warning('off', 'MATLAB:imagesci:hdf:removalWarningHDFSD')

% get the kcarta to L1c AIRS convolution matrix
sdir = '/asl/matlab2012/srftest/';
srf1 = fullfile(sdir, 'srftables_m140f_withfake_mar08.hdf');

v_tmp = load('freq2645.txt');
dvk = 0.0025;    % kcarta dv
dvb = 0.1;       % decon dv

% L1c convolution matrix
[A1, vAcol, vArow] = mksconv2(srf1, v_tmp, dvb);

% L1d convolution matrix
[B1, vBcol, vBrow] = L1d_conv(res, dvb, v_base);

% match intermediate grids
[ix, jx] = seq_match(vAcol, vBcol);
vAcol = vAcol(ix);  A1 = A1(:, ix);
vBcol = vBcol(jx);  B1 = B1(:, jx);

fprintf(1, 'inverting A ...\n')
tic
Ainv =  pinv(full(A1));
toc

% build the decon/recon matrix
CtoD = B1 * Ainv;

figure(1); clf
subplot(2,1,1)
plot(vArow, Ainv(1001:1002, :), 'linewidth', 2)
axis([740, 752, -2.4, 2.9])
title('sample deconvolution matrix rows')
legend('745.5 cm-1', '745.6 cm-1', 'location', 'northeast')
ylabel('weight')
grid on; zoom on

subplot(2,1,2)
plot(vArow, CtoD(193:196, :), 'linewidth', 2)
axis([740, 752, -0.1, 0.4])
title('sample C to D matrix rows')
legend('745.30 cm-1', '745.84 cm-1', '746.37 cm-1', ...
       '746.90 cm-1', 'location', 'northeast')
xlabel('AIRS L1c wavenumber')
ylabel('weight')
grid on; zoom on
saveas(gcf, 'airs_decon_basis', 'fig')

return

% L1c synthetic channel list
% load L1c_synth

figure(1); clf
  plot(vAcol, Ainv(:, 201:202), vArow, 0, 'ok', 'linewidth', 2)
  axis([697, 706, -2, 3])
% plot(bfrq, binv(:, 2401:2402), afrq, 0, 'ok', 'linewidth', 2)
% plot(bfrq, binv(:, 2401:2402), L1c_syn, 0, '*g', 'linewidth', 2)
% axis([2370, 2460, -4, 4])
title('sample deconvolution basis functions')
  legend('701.06 cm-1', '701.34 cm-1', 'AIRS channels')
% legend('column 201', 'column 202', 'AIRS channels')
% legend('column 2401', 'column 2402', 'AIRS channels')
xlabel('wavenumber')
ylabel('weight')
grid on; zoom on


