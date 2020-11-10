% 
% demo2 -- CrIS to AIRS translation demo
% 

% local libs
addpath ../data
addpath ../source
addpath /asl/packages/ccast/source

% sample CrIS SDR granule
cpath = '/asl/cris/ccast/sdr45_npp_HR/2018/091';
cgran = 'CrIS_SDR_npp_s45_d20180401_t2006080_g202_v20d.mat';
d1 = load(fullfile(cpath, cgran));
rLW = squeeze(d1.rLW(:, 5, 15, 21:30));
rMW = squeeze(d1.rMW(:, 5, 15, 21:30));
rSW = squeeze(d1.rSW(:, 5, 15, 21:30));
vLW = d1.vLW; vMW = d1.vMW; vSW = d1.vSW;
clear d1

% AIRS SRF tabulation file
sfile = 'airs_l1c_srf_tables_lls_20181205.hdf';
frq1c = load('data/freq2645.txt');

% translate CrIS to AIRS
opt1 = struct;
opt1.user_res = 'hires'; 
[arad, afrq, brad, bfrq] = ...
       cris2airs(rLW, rMW, rSW, vLW, vMW, vSW, sfile, frq1c, opt1);

% select and plot a sample obs
ix = 4;
crad = [rLW; rMW; rSW];
cfrq = [vLW; vMW; vSW];
cbt = real(rad2bt(cfrq, crad(:, ix)));
bbt = real(rad2bt(bfrq, brad(:, ix)));
abt = real(rad2bt(afrq, arad(:, ix)));
figure(1); clf
[x1, y1] = pen_lift(cfrq, cbt);
[x2, y2] = pen_lift(bfrq, bbt);
[x3, y3] = pen_lift(afrq, abt);
plot(x1, y1, x2, y2, x3, y3)
axis([600, 2500, 200, 300])
legend('CrIS', 'CrIS decon', 'CrIS to AIRS', 'location', 'north')
title('CrIS to AIRS demo')
xlabel('wavenumber (cm-1)'); 
ylabel('BT (K)')
grid on; zoom on

