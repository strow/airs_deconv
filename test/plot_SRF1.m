%
% plot sample SRFs and columns of inverse SRF matrix
%

addpath ../h4tools
addpath ../source

% turn off HDF 4 update warnings
warning('off', 'MATLAB:imagesci:hdf:removalWarningHDFSD')

% sample decon basis functions
load bconv.mat
load L1c_synth
figure(1); clf
% set(gcf, 'Units','centimeters', 'Position', [4, 10, 24, 16])
% plot(bfrq, binv(:, 201:202), afrq, 0, 'ok', 'linewidth', 2)
% axis([697, 706, -2, 3])
% plot(bfrq, binv(:, 2401:2402), afrq, 0, 'ok', 'linewidth', 2)
  plot(bfrq, binv(:, 2401:2402), L1c_syn, 0, '*g', 'linewidth', 2)
  axis([2370, 2460, -4, 4])
title('sample deconvolution basis functions')
% legend('column 201', 'column 202', 'AIRS channels')
  legend('column 2401', 'column 2402', 'AIRS channels')
xlabel('wavenumber')
ylabel('weight')
grid on; zoom on
% export_fig('airs_decon_basis.pdf', '-m2', '-transparent')

% sample AIRS srfs
sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';
cfreq = load('freq2645.txt');
% cfreq = trim_chans(cfreq);
dvk = 0.0025; 
[sconv, sfreq, ofreq] = mksconv2(sfile, cfreq, dvk);
figure(2)
% set(gcf, 'Units','centimeters', 'Position', [4, 10, 24, 16])
% plot(sfreq, sconv(201:207, :), 'linewidth', 2)
% axis([700, 704, 0, 5e-3])
  plot(sfreq, sconv(2414:2417, :), 'linewidth', 2)
  axis([2418, 2428, 0, 1.4e-3])
title('sample AIRS SRFs')
xlabel('wavenumber')
ylabel('weight')
grid on; zoom on
%  export_fig('airs_sample_SRFs.pdf', '-m2', '-transparent')

