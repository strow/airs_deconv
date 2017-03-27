%
% conv_loop3 - generate a few black body spectra
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
% derived from conv_loop1 and L1d_test2
%

% set paths to libs
addpath ../source
addpath ../h4tools
addpath /asl/packages/ccast/source
addpath /home/motteler/matlab/export_fig

% L1d resolution, dv = v / res
% res = 2000;  % for direct regression test
  res = 1200;  % L1c nominal
% res =  700;  % best for L1d

% turn off HDF 4 update warnings
warning('off', 'MATLAB:imagesci:hdf:removalWarningHDFSD')

% get the kcarta to L1c AIRS convolution matrix
sdir = '/asl/matlab2012/srftest/';
srf1 = fullfile(sdir, 'srftables_m140f_withfake_mar08.hdf');

v_tmp = load('freq2645.txt');
dvk = 0.0025;    % kcarta dv
dvb = 0.1;       % decon dv

% black body temperature list and freq grid
nkc = round((3000 - 600) / dvk);
vkc = 600 + (0 : nkc-1)' * dvk;
bbt = [220 240 260 280, 300];
nbbt = length(bbt);

% L1c convolution matrices
[C1, vCcol, v_L1c] = mksconv2(srf1, v_tmp, dvk);
[A1, vAcol, vArow] = mksconv2(srf1, v_tmp, dvb);

% set start of output grid
% v_base = v_L1c(1);  % default, 649.6620
  v_base = 649.6620;  % best for res 1200
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

% initialize output arrays
nc = length(v_L1c);
nd = length(v_L1d);
CtoDrad  = zeros(nd, nbbt);
trueDrad = zeros(nd, nbbt);
trueCrad = zeros(nc, nbbt);

% loop on kcarta files
for i = 1 : nbbt

  % get kcarta radiances
  rkc = real(bt2rad(vkc, bbt(i)));
 
  % kcarta to L1c
  ix = interp1(vkc, 1:length(rkc), vCcol, 'nearest');
  r1 = C1 * rkc(ix);
  trueCrad(:, i) = r1;

  % L1c to L1d
  r2 = CtoD * r1; 
  CtoDrad(:, i) = r2;

  % kcarta to L1d
  ix = interp1(vkc, 1:length(rkc), vDcol, 'nearest');
  r3 = D1 * rkc(ix);  
  trueDrad(:, i) = r3;

end

% save the convolutions
save L1d1200_bbt v_L1c v_L1d CtoDrad trueDrad trueCrad bbt res

