%
% demo1 -- AIRS to CrIS translation demo
% 

% local libs
addpath ../source
addpath /asl/packages/ccast/source

% AIRS SRF tabulation file
sfile = 'data/airs_demo_srf.hdf';

% AIRS sample radiance data
% arad - nchan x nobs AIRS radiances
% afrq - nchan AIRS channel frequencies
load data/airs_demo_rad

% basic translation 
opt1 = struct;
opt1.user_res = 'midres';  % target resolution
[crad1, cfrq1] = airs2cris(arad, afrq, sfile, opt1);

% apodized translation
opt1.hapod = 1;  % Hamming apodization
opt1.scorr = 1;  % statistical correction
opt1.cfile = 'corr_midres.mat';  % correction weights
[crad2, cfrq2] = airs2cris(arad, afrq, sfile, opt1);

% plot a selected obs
iobs = 101;
abt = real(rad2bt(afrq, arad(:, iobs)));
cbt1 = real(rad2bt(cfrq1, crad1(:, iobs)));
cbt2 = real(rad2bt(cfrq2, crad2(:, iobs)));
figure(1); clf
[x1, y1] = pen_lift(afrq, abt);
[x2, y2] = pen_lift(cfrq1, cbt1);
[x3, y3] = pen_lift(cfrq2, cbt2);
plot(x1, y1, x2, y2, x3, y3)
axis([600, 2500, 200, 300])
title('AIRS to CrIS demo')
legend('AIRS', 'AIRS to CrIS', 'apodized AIRS to CrIS', ...
       'location', 'southeast')
xlabel('wavenumber (cm-1)')
ylabel('BT (K)')
grid on; zoom on

