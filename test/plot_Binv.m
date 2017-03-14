%
% plot sample SRFs and columns of inverse SRF matrix
%

% set paths to libs
addpath ../source
addpath ../h4tools
addpath /asl/packages/ccast/source
addpath /home/motteler/matlab/export_fig

% turn off HDF 4 update warnings
warning('off', 'MATLAB:imagesci:hdf:removalWarningHDFSD')

% get the kcarta to L1c AIRS convolution matrix
sdir = '/asl/matlab2012/srftest/';
srf1 = fullfile(sdir, 'srftables_m140f_withfake_mar08.hdf');

v_tmp = load('freq2645.txt');
dvk = 0.0025;    % kcarta dv
dvb = 0.1;       % decon dv

% L1c convolution matrices
[A1, vAcol, vArow] = mksconv2(srf1, v_tmp, dvb);
fprintf(1, 'inverting A ...\n')
tic
Ainv =  pinv(full(A1));
toc

% L1c synthetic channel list
% load L1c_synth

figure(1); clf
  set(gcf, 'Units','centimeters', 'Position', [4, 10, 24, 16])
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
% export_fig('airs_decon_basis.pdf', '-m2', '-transparent')

