%
% AIRS channel deconvolution test
%
% Let S1 and S2 be two kcarta to channel grid convolutions, and A1
% and A2 the corresponding channel convolution matrices, so that S1
% and A1 have channel grid C1, and S2 and A2 channel grid C2.
%
% We compare the following values:
%
%   1. kcarta rad -> S2 convolution
%
%   2. kcarta rad -> S1 convolution -> interpolate to S2 grid
%
%   3. kcarta rad -> S1 convolution -> A1 deconvolution ->
%                    interpolate to A2 grid -> A2 convolution
%
% this runs OK but the results are similar to or slightly worse
% than with spline interpolation.

addpath /asl/matlab/h4tools

% specify HDF SRF tabulation files
hfile1 = '/asl/matlab/srftest/srftables_m140f_withfake_mar08.hdf';
hfile2 = 'fake_srf1.hdf';

% specify the convolution mat files
sfile1 = 'test2_S1';
sfile2 = 'test2_S2';
afile1 = 'test2_A1';
afile2 = 'test2_A2';

% get the 1b channel frequencies
cfreq1 = getcfreq(hfile1);
cfreq2 = getcfreq(hfile2);

% load some sample kcarta data
%   d1.r  - radiance
%   d1.w  - frequency
kcfile = '/home/motteler/cris/sergio/JUNK2012/convolved_kcarta13.mat';
d1 = load(kcfile);

% call mkaconv2 to build Sa, convolution matrix at channel grid
%   Amat   - AIRS convolution matrix
%   afrq   - sorted trimmed AIRS channel set
%   aind   - index of afrq in cfreq
%   sind   - index of afrq in SRF tabulation
mkaconv2(cfreq1, hfile1, afile1)
mkaconv2(cfreq2, hfile2, afile2)
A1 = load(afile1);
A2 = load(afile2);

% call mksconv to build Sk, convolution matrix at kcarta grid
%   Cmat    - sparse convolution matrix
%   Cfin    - Cmat input frequency list (a dv grid spanning all SRFs)
%   Cfout   - Cmat output frequency list (chan centers from the SRF file)
%   Ciout   - Cmat output channel number (chan IDs from the SRF file)
mksconv(hfile1, sfile1);
mksconv(hfile2, sfile2);
S1 = load(sfile1);
S2 = load(sfile2);

% match kcarta and S1 convolution matrix column indices
Cfin = S1.Cfin;
dvk = 0.0025;
vk1 = max(d1.w(1), Cfin(1));
vk2 = min(d1.w(end), Cfin(end));
ng1 = round((vk2-vk1)/dvk);
vg1 = vk1 + (0:ng1-1)*dvk;
ixk = interp1(d1.w, 1:length(d1.w), vg1, 'nearest');
ixs = interp1(Cfin, 1:length(Cfin), vg1, 'nearest');

vkc = d1.w(ixk);   % kcarta common frequencies
rkc = d1.r(ixk);   % kcarta common radiances

% get channel radiances from S1
rc1 = S1.Cmat(A1.sind, ixs) * rkc;
vc1 = S1.Cfout(A1.sind);
clear S1

% indexing sanity check
max(abs(vkc - Cfin(ixs))) / mean(vkc)
isequal(A1.afrq, vc1)

% match kcarta and S2 convolution matrix column indices
Cfin = S2.Cfin;
dvk = 0.0025;
vk1 = max(d1.w(1), Cfin(1));
vk2 = min(d1.w(end), Cfin(end));
ng1 = round((vk2-vk1)/dvk);
vg1 = vk1 + (0:ng1-1)*dvk;
ixk = interp1(d1.w, 1:length(d1.w), vg1, 'nearest');
ixs = interp1(Cfin, 1:length(Cfin), vg1, 'nearest');

vkc = d1.w(ixk);   % kcarta common frequencies
rkc = d1.r(ixk);   % kcarta common radiances

% get channel radiances from S2
rc2 = S2.Cmat(A2.sind, ixs) * rkc;
vc2 = S2.Cfout(A2.sind);
clear S2

% indexing sanity check
max(abs(vkc - Cfin(ixs))) / mean(vkc)
isequal(A2.afrq, vc2)

% direct interpolation of rc1 to the vc2 grid
rc2i = interp1(vc1, rc1, vc2, 'spline', 'extrap');

% deconvolution, interpolation, and reconvolution
rc1t = inv(A1.Amat) * rc1;
rc2t = interp1(vc1, rc1t, vc2, 'spline', 'extrap');
rc2d = A2.Amat * rc2t;

% get brighness temps for plots and stats
bt1 = real(rad2bt(vc1, rc1));     % c1 true
bt2 = real(rad2bt(vc2, rc2));     % c2 true
bt2i = real(rad2bt(vc2, rc2i));   % c2 interpolated
bt2d = real(rad2bt(vc2, rc2d));   % c2 final deconvolved
bt1t = real(rad2bt(vc2, rc1t));   % c1 deconvolved

% plot transform stages
figure(1); clf
plot(vc1, bt1, vc1, bt1t)
title('S1 transform stages')
xlabel('wavenumber')
ylabel('dBt')
legend('channel', 'deconv')
grid; zoom on

% plot residuals 
figure(2); clf
plot(vc2, bt2i - bt2, vc2, bt2d - bt2)
title('S2 transform residuals')
xlabel('wavenumber')
ylabel('dBt')
legend('interp', 'deconv')
grid; zoom on

[rms(bt2d-bt2), rms(bt2i-bt2)]

