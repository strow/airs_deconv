%
% AIRS SRF tabulation grid deconvolution test
%
% Let S1 and S2 be two kcarta to channel grid convolutions, and B1
% and B2 the corresponding SRF-grid to channel convolutions so that
% S1 and B1 have channel grid C1, and S2 and B2 channel grid C2.
%
% this version does deconvolution with the matlab pinv function
%
% We compare the following values:
%
%   1. kcarta rad -> S2 convolution
%
%   2. kcarta rad -> S1 convolution -> shifted to S2 grid
%
%   3. kcarta rad -> S1 convolution -> B1 deconvolution -> B2 convolution
%                    
% the A1 and A2 (channel convolution) files are generated here only
% to get the trimmed channel sets
%
% in a separate final section we compare the column sum of B1 vs the
% row sum of B1 inv
%

% plot setup
figtype = 'fig';  % plot data type
fignote = '1c';   % plot file suffix

% use my versions of HDF tools
addpath /home/motteler/mot2008/hdf/h4tools

% specify HDF SRF tabulation files
hfile1 = '/asl/matlab/srftest/srftables_m140f_withfake_mar08.hdf';
hfile2 = 'fake_srf1.hdf';

% specify the convolution mat files
sfile1 = 'test2_S1';
sfile2 = 'test2_S2';
afile1 = 'test2_A1';
afile2 = 'test2_A2';
bfile1 = 'test2_B1';
bfile2 = 'test2_B2';

% get the SRF channels frequencies
% cind = 1:2378;    % for 1b channels
load index1c
cind = index1c;   % for 1c channels
cfreq1 = getcfreq(hfile1, cind);
cfreq2 = getcfreq(hfile2, cind);

% load some sample kcarta data
%   d1.r  - radiance
%   d1.w  - frequency
kcfile = '/home/motteler/cris/sergio/JUNK2012/convolved_kcarta13.mat';
d1 = load(kcfile);

% call mkaconv2 to build Sa, convolution matrix at channel grid
% (this used here for the cleaned up 1b channel set, not Amat)
%   Amat   - AIRS convolution matrix
%   afrq   - sorted trimmed AIRS channel set
%   aind   - index of afrq in cfreq
%   sind   - index of afrq in SRF tabulation
mkaconv2(cfreq1, hfile1, afile1)
mkaconv2(cfreq2, hfile2, afile2)
A1 = load(afile1);
A2 = load(afile2);

% call mksconv to build convolution matrices at kcarta and SRF grids
%   Cmat    - sparse convolution matrix
%   Cfin    - Cmat input frequency list (a dv grid spanning all SRFs)
%   Cfout   - Cmat output frequency list (chan centers from the SRF file)
%   Ciout   - Cmat output channel number (chan IDs from the SRF file)
dvk = 0.0025;  % kcarta grid step
dvs = 0.05;    % SRF tabulation grid step
mksconv(hfile1, sfile1, dvk);
mksconv(hfile2, sfile2, dvk);
mksconv(hfile1, bfile1, dvs);
mksconv(hfile2, bfile2, dvs);
S1 = load(sfile1);
S2 = load(sfile2);
B1 = load(bfile1);
B2 = load(bfile2);

% match kcarta and S1 convolution matrix column indices
[ixk, ixs] = seq_isect(d1.w, S1.Cfin);
vkc = d1.w(ixk);   % kcarta common frequencies
rkc = d1.r(ixk);   % kcarta common radiances

% get channel radiances from S1
rc1 = S1.Cmat(A1.sind, ixs) * rkc;
vc1 = S1.Cfout(A1.sind);

% indexing sanity check
max(abs(vkc - S1.Cfin(ixs))) / mean(vkc)
isequal(A1.afrq, vc1)
clear S1

% match kcarta and S2 convolution matrix column indices
[ixk, ixs] = seq_isect(d1.w, S2.Cfin);
vkc = d1.w(ixk);   % kcarta common frequencies
rkc = d1.r(ixk);   % kcarta common radiances

% get channel radiances from S2
rc2 = S2.Cmat(A2.sind, ixs) * rkc;
vc2 = S2.Cfout(A2.sind);

% indexing sanity check
max(abs(vkc - S2.Cfin(ixs))) / mean(vkc)
isequal(A2.afrq, vc2)
clear S2

% direct interpolation of rc1 to the vc2 grid
rc2i = interp1(vc1, rc1, vc2, 'spline', 'extrap');

% deconvolution, match grids, and reconvolution
% rc1d = pinv(full(B1.Cmat(A1.sind,:))) * rc1;
B = full(B1.Cmat(A1.sind,:));
Binv = pinv(B);
rc1d = Binv * rc1;
vc1d = B1.Cfin;
vc2d = B2.Cfin;
[ix1,ix2] = seq_isect(vc1d, vc2d);
rc2r = B2.Cmat(A2.sind, ix2) * rc1d(ix1);

% get brighness temps for plots and stats
bt1 = real(rad2bt(vc1, rc1));     % c1 true
bt2 = real(rad2bt(vc2, rc2));     % c2 true
bt2i = real(rad2bt(vc2, rc2i));   % c2 interpolated
bt1d = real(rad2bt(vc1d, rc1d));  % c1 deconvolved
bt2r = real(rad2bt(vc2, rc2r));   % c2 reconvolved

% plot transform stages
figure(1); clf
plot(vc1, bt1, vc1d, bt1d)
title('SRF transform stages')
xlabel('wavenumber')
ylabel('dBt')
legend('channel', 'deconv')
grid; zoom on
saveas(gcf, ['xform_stages_',fignote], figtype)

% plot residuals 
figure(2); clf
plot(vc2, bt2i - bt2, vc2, bt2r - bt2)
title('SRF transform residuals')
xlabel('wavenumber')
ylabel('dBt')
legend('interp', 'decon')
grid; zoom on
saveas(gcf, ['xform_resids_',fignote], figtype)

% summary stats for some sub-bands without gaps
ix = 140:430;
fprintf(1, '%g to %g, deconv %.4g, interp %.4g\n', ...
        vc2(ix(1)), vc2(ix(end)), ...
        rms(bt2r(ix)-bt2(ix)), rms(bt2i(ix)-bt2(ix)));

ix = 1360:1640;
fprintf(1, '%g to %g, deconv %.4g, interp %.4g\n', ...
        vc2(ix(1)), vc2(ix(end)), ...
        rms(bt2r(ix)-bt2(ix)), rms(bt2i(ix)-bt2(ix)));

ix = 1860:2060;
fprintf(1, '%g to %g, deconv %.4g, interp %.4g\n', ...
        vc2(ix(1)), vc2(ix(end)), ...
        rms(bt2r(ix)-bt2(ix)), rms(bt2i(ix)-bt2(ix)));

% plot a row of B (a single SRF) vs a column of Binv
figure(3); clf

vc = 1000;  % select a channel frequency
ic = min(find(vc <= vc1));

plot(vc1d, 50*B(ic, :), vc1d, Binv(:,ic), vc1, 0, 'k+')
ax = axis; ax(1) = round(vc-6); ax(2) = round(vc+6); axis(ax);
title('B row (SRF) and B inv column for one channel')
xlabel('wavenumber')
ylabel('normalized weights')
legend('B row (SRF)', 'B inv column', 'channels', 'location', 'best')
grid on ; zoom on
saveas(gcf, ['xform_SRF_',fignote], figtype)

% plot column sums of B vs row sums of B inv
figure(4); clf

yy = sum(B);       % sum columns of B   
zz = sum(Binv,2);  % sum rows of Binv

plot(vc1d, 10*yy, vc1d, zz, vc1, 0, 'k+')
% plot(vc1d, 10*yy, vc1d, zz)
% ax=axis; ax(1)=1020; ax(2)=1080; ax(3)=0; ax(4) = 1.5; axis(ax);
title('B column and Binv row sums')
xlabel('wavenumber')
ylabel('normalized weights')
legend('B column sum', 'B inv row sum', 'AIRS channels', 'location', 'best')
% legend('B column sum', 'B inv row sum', 'location', 'best')
grid on; zoom on
saveas(gcf, ['xform_sums_',fignote], figtype)



