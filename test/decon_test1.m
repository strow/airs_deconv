%
% decon_test1 - direct vs decon to intermediate grid
%
%   C1 takes kcarta to the L1c channel grid,
%   B1 takes kcarta to the 0.1 cm-1 decon grid
%   A1 takes the 0.1 cm-1 grid to L1c grid
% 
%   true AIRS = C1 * rkc
%   true decon = B1 * rkc
%   decon AIRS = inv(A1) * C1 * rkc
%
% then we compare true decon with decon AIRS
%

% set paths to libs
addpath ../source
addpath ../h4tools
addpath /asl/packages/ccast/source

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
dvb = 0.1;       % direct convolution and decon spacing
span = 4;        % direct convolution FWHM = span * dvb 

% L1c convolution matrices
[C1, vCcol, v_L1c] = mksconv2(srf1, v_tmp, dvk);
[A1, vAcol, vArow] = mksconv2(srf1, v_tmp, dvb);

fprintf(1, 'inverting A ...\n')
tic
A1inv = pinv(full(A1));
toc

% intermediate grid convolution matrix
[B1, vBcol, vBrow] = mkBconv(dvb, span);

% loop on kcarta files
trueCrad = []; trueBrad = [];
for i = 1 : length(flist)
  d1 = load(fullfile(kcdir, flist(i).name));
  vkc = d1.w(:); rkc = d1.r(:);

  % convolve kcarta to AIRS L1c
  ix = interp1(vkc, 1:length(rkc), vCcol, 'nearest');
  r1 = C1 * rkc(ix);  trueCrad = [trueCrad, r1];

  % convolve kcarta to intermediate "B" grid
  ix = interp1(vkc, 1:length(rkc), vBcol, 'nearest');
  r2 = B1 * rkc(ix);  trueBrad = [trueBrad, r2];

  fprintf(1, '.');
end
fprintf(1, '\n')

% deconvolve L1C to intermediate "B" grid
CtoBrad = A1inv * trueCrad;

% take radiances to brightness temps
trueCbt = real(rad2bt(v_L1c, trueCrad));
trueBbt = real(rad2bt(vBrow, trueBrad));
CtoBbt  = real(rad2bt(vAcol, CtoBrad));

% match frequency grids
[ix, jx] = seq_match(vBrow, vAcol);
vBrow = vBrow(ix); trueBbt = trueBbt(ix, :);
vAcol = vAcol(jx); CtoBbt = CtoBbt(jx, :);

% profile 1 spectra and zoom
figure(1); clf
subplot(2,1,1)
plot(vBrow, trueBbt(:, 1), vAcol, CtoBbt(:, 1))
axis([650, 2650, 200, 310])
title('direct and deconvolution comparison')
legend('gauss', 'decon', 'location', 'north')
% xlabel('wavenumber')
ylabel('Tb, K')
grid on; zoom on

subplot(2,1,2)
plot(vBrow, trueBbt(:, 1), vAcol, CtoBbt(:, 1))
axis([660, 680, 200, 260])
title('direct and deconvolution detail')
legend('gauss', 'decon', 'location', 'northeast')
xlabel('wavenumber')
ylabel('Tb, K')
grid on; zoom on

% mean and std residuals
figure(2); clf
subplot(2,1,1)
plot(vBrow, mean(CtoBbt - trueBbt, 2))
axis([650, 2650, -15, 15])
title('mean decon minus gauss')
% xlabel('wavenumber')
ylabel('dTb, K')
grid on; zoom on

subplot(2,1,2)
plot(vBrow, std(CtoBbt - trueBbt, 0, 2))
axis([650, 2650, 0, 5])
title('std decon minus gauss')
xlabel('wavenumber')
ylabel('dTb, K')
grid on; zoom on

