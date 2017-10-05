% 
% demo2 -- CrIS to AIRS translation demo
% 

% set paths to asl libs
addpath ../source
addpath /asl/packages/ccast/source

% get sample CrIS radiances
f1 = '/asl/data/cris/ccast/sdr60_hr/2015/231/SDR_d20150819_t0821345.mat';
d1 = load(f1);
rLW = squeeze(d1.rLW(:, 5, 15, 21:30));
rMW = squeeze(d1.rMW(:, 5, 15, 21:30));
rSW = squeeze(d1.rSW(:, 5, 15, 21:30));
vLW = d1.vLW; vMW = d1.vMW; vSW = d1.vSW;
clear d1

% specify AIRS SRF and channel set
sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';
cfrq = load('freq2645.txt');
opt1 = struct;
opt1.resmode = 'hires2'; 

% translate CrIS to AIRS
[arad, afrq, brad, bfrq] = ...
       cris2airs(rLW, rMW, rSW, vLW, vMW, vSW, sfile, cfrq, opt1);

% select and plot a sample obs
ix = 4;
abt = real(rad2bt(afrq, arad(:, ix)));
bbt = real(rad2bt(bfrq, brad(:, ix)));
crad = [rLW; rMW; rSW];
cfrq = [vLW; vMW; vSW];
cbt = real(rad2bt(cfrq, crad(:, ix)));
figure(1); clf
plot(cfrq, cbt, bfrq, bbt, afrq, abt)
legend('true cris', 'cris decon', 'cris to airs')
title('CrIS to AIRS demo')
xlabel('wavenumber'); 
ylabel('brightness temp')
grid on; zoom on

