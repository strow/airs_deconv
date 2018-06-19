%
% decon_test1 - compare AIRS decon with true decon
%
% transform matrices
%   C1 takes kcarta to the L1c channel grid, 
%   B1 takes kcarta to the 0.1 cm-1 decon grid
%   A1 takes the 0.1 cm-1 grid to L1c grid
% 
% radiance sets
%   true AIRS = C1 * rkc
%   true decon = B1 * rkc
%   decon AIRS = inv(A1) * C1 * rkc
%

% set paths to libs
addpath ../source
addpath /asl/packages/ccast/source
addpath /home/motteler/matlab/export_fig

% kcarta test data
kcdir = '/home/motteler/cris/sergio/JUNK2012/';
flist =  dir(fullfile(kcdir, 'convolved_kcart*.mat'));

% get the kcarta to L1c AIRS convolution matrix
sdir = '/asl/matlab2012/srftest/';
srf1 = fullfile(sdir, 'srftables_m140f_withfake_mar08.hdf');

v_tmp = load('freq2645.txt');
dvk = 0.0025;    % kcarta dv
dvb = 0.1;       % deconvolution spacing
res = 2000;      % set resolving power
p = 1.4;         % sup gauss exponent

% L1c convolution matrices
[C1, vCcol, v_L1c] = mksconv2(srf1, v_tmp, dvk);
[A1, vAcol, vArow] = mksconv2(srf1, v_tmp, dvb);

fprintf(1, 'inverting A ...\n')
tic
A1inv = pinv(full(A1));
toc

% intermediate grid convolution matrix
[B1, vBcol, vBrow] = mkBconv(dvb, res, p);

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
subplot(3,1,1)
plot(vBrow, trueBbt(:, 1), vAcol, CtoBbt(:, 1))
axis([650, 2650, 200, 310])
title('deconvolution and reference spectra')
legend('decon ref', 'AIRS decon', 'location', 'north')
% xlabel('wavenumber (cm^{-1})')
ylabel('BT (K)')
grid on; zoom on

subplot(3,1,2)
plot(vBrow, trueBbt(:, 1), vAcol, CtoBbt(:, 1))
axis([660, 680, 200, 260])
% title('deconvolution and reference detail')
legend('decon ref', 'AIRS decon', 'location', 'northeast')
% xlabel('wavenumber (cm^{-1})')
ylabel('BT (K)')
grid on; zoom on

subplot(3,1,3)
plot(vBrow, trueBbt(:, 1), vAcol, CtoBbt(:, 1))
axis([2280, 2340, 220, 260])
% title('deconvolution and reference detail')
legend('decon ref', 'AIRS decon', 'location', 'northwest')
xlabel('wavenumber (cm^{-1})')
ylabel('BT (K)')
grid on; zoom on
% saveas(gcf, 'airs_decon_spec', 'fig')

% profile 1 decon comparison with AIRS
d1 = load(fullfile(kcdir, flist(1).name));
vkc = d1.w(:); rkc = d1.r(:);
bkc = real(rad2bt(vkc, rkc));

figure(2); clf
subplot(2,1,1)
plot(vkc, bkc, vBrow, trueBbt(:, 1), vAcol, CtoBbt(:, 1), ...
     v_L1c, trueCbt(:,1), 'linewidth', 2)
% axis([650, 2650, 180, 310])   % used for overview
% axis([1012, 1016, 230, 300])
% axis([669, 672, 180, 280])
  axis([1312, 1316, 220, 280])  % used for sample zoom 1
% axis([812, 817, 260, 300])    % used for sample zoom 2

title('kcarta, AIRS, and deconvolution')
legend('kcarta','decon ref','AIRS decon','true AIRS','location','southeast')
% xlabel('wavenumber')
ylabel('BT (K)')
grid on; zoom on

subplot(2,1,2)
plot(vkc, bkc, vBrow, trueBbt(:, 1), vAcol, CtoBbt(:, 1), ...
     v_L1c, trueCbt(:,1), 'linewidth', 2)
% axis([650, 2650, 180, 310])   % used for overview
% axis([1012, 1016, 230, 300])
% axis([669, 672, 180, 280])
% axis([1312, 1316, 220, 280])  % used for sample zoom 1
  axis([812, 816, 260, 300])    % used for sample zoom 2

% title('kcarta, AIRS, and deconvolutions')
legend('kcarta','decon ref','AIRS decon','true AIRS','location','southeast')
xlabel('wavenumber (cm^{-1})')
ylabel('BT (K)')
grid on; zoom on
% saveas(gcf, 'airs_decon_zoom', 'fig')

% 49 profile mean and std residuals
figure(3); clf
subplot(2,1,1)
plot(vBrow, mean(CtoBbt - trueBbt, 2))
axis([650, 2650, -15, 15])
title('mean decon minus ref')
% xlabel('wavenumber (cm^{-1})')
ylabel('\Delta BT (K)')
grid on; zoom on

subplot(2,1,2)
plot(vBrow, std(CtoBbt - trueBbt, 0, 2))
axis([650, 2650, 0, 5])
title('std decon minus ref')
xlabel('wavenumber (cm^{-1})')
ylabel('\Delta BT (K)')
grid on; zoom on
% saveas(gcf, 'airs_decon_diff', 'fig')

