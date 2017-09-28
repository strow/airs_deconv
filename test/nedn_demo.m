%
% nedn_demo - ccast NEdN estimate plots
%

addpath /asl/packages/ccast/source
addpath /asl/packages/ccast/motmsc/utils
addpath ../source

% sample ccast CrIS NEdN estimate
% d1 = load('/asl/data/cris/ccast/sdr60_hr/2016/018/SDR_d20160118_t0801033.mat');
% d1 = load('/asl/data/cris/ccast/sdr60_hr/2017/031/SDR_d20170131_t0749521.mat');
% d1 = load('/asl/data/cris/ccast/sdr60_hr/2017/031/SDR_d20170131_t2149474.mat');
% d1 = load('/asl/data/cris/ccast/sdr60_hr/2016/091/SDR_d20160331_t0807119.mat');
  d1 = load('/asl/data/cris/ccast/sdr60_hr/2016/301/SDR_d20161027_t0650511.mat');
nednLW = mean(d1.nLW, 3);
nednMW = mean(d1.nMW, 3);
nednSW = mean(d1.nSW, 3);
nedn_ccast = [nednLW; nan(1,9); nednMW; nan(1,9); nednSW];
freq_ccast = [d1.vLW; nan; d1.vMW; nan; d1.vSW];

% option for apodized noise 
nedn_ccast = 0.63 * nedn_ccast;

figure(1); clf
% set(gcf, 'Units','centimeters', 'Position', [4, 10, 24, 16])
set(gcf, 'DefaultAxesColorOrder', fovcolors);
semilogy(freq_ccast, nedn_ccast);
axis([600, 2600, 0, 1])
title('CrIS full res NEdN estimates')
legend(fovnames, 'location', 'northeast')
xlabel('wavenumber')
ylabel('NEdN, mw sr-1 m-2')
grid on; zoom on

return

figure(1); clf
% set(gcf, 'Units','centimeters', 'Position', [4, 10, 24, 16])
set(gcf, 'DefaultAxesColorOrder', fovcolors);
plot(d1.vLW, mean(d1.nLW, 3))
legend(fovnames, 'location', 'north')
xlabel('wavenumber')
ylabel('NEdN, mw sr-1 m-2')
grid on

figure(2); clf
% set(gcf, 'Units','centimeters', 'Position', [4, 10, 24, 16])
set(gcf, 'DefaultAxesColorOrder', fovcolors);
plot(d1.vMW, mean(d1.nMW, 3))
legend(fovnames, 'location', 'north')
xlabel('wavenumber')
ylabel('NEdN, mw sr-1 m-2')
grid on

figure(3); clf
% set(gcf, 'Units','centimeters', 'Position', [4, 10, 24, 16])
set(gcf, 'DefaultAxesColorOrder', fovcolors);
plot(d1.vSW, mean(d1.nSW, 3))
legend(fovnames, 'location', 'north')
xlabel('wavenumber')
ylabel('NEdN, mw sr-1 m-2')
grid on

