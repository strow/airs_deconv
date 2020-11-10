%
% nedn_demo - plot sample ccast unsmoothed measured NEdN 
%
% the ccast values are simply the std of the measured ICT,
% without any added smoothing, and so will look noisier than 
% other estimates.
%

addpath /asl/packages/ccast/source
addpath /asl/packages/ccast/motmsc/utils
addpath ../source

% sample ccast CrIS NEdN estimate
p1 = '/asl/cris/ccast/sdr45_npp_HR/2015/121';
g1 = 'CrIS_SDR_npp_s45_d20150501_t2106020_g212_v20a.mat';
d1 = load(fullfile(p1,g1));

%---------------------
% summary 3-band plot
%---------------------
ap_ind = 1;  % 1 = unapodized, 2 = apodized
nednLW = d1.nLW(:,:,ap_ind);
nednMW = d1.nMW(:,:,ap_ind);
nednSW = d1.nSW(:,:,ap_ind);
nedn_ccast = [nednLW; nan(1,9); nednMW; nan(1,9); nednSW];
freq_ccast = [d1.vLW; nan; d1.vMW; nan; d1.vSW];

figure(1); clf
set(gcf, 'DefaultAxesColorOrder', fovcolors);
semilogy(freq_ccast, nedn_ccast);
axis([600, 2600, 0, 1])
title('CrIS ccast unapodized full res NEdN estimates')
legend(fovnames, 'location', 'northeast')
xlabel('wavenumber')
ylabel('NEdN, mw sr-1 m-2')
grid on; zoom on

%---------------------
% single band and FOV
%---------------------
iFOV = 5;

% LW detail
nedn_unap = d1.nLW(:,iFOV,1);
nedn_apod = d1.nLW(:,iFOV,2);
nedn_apxx = d1.nLW(:,iFOV,1) * 0.63;
freq = d1.vLW;

figure(2); clf
semilogy(freq, nedn_unap, freq, nedn_apod, freq,nedn_apxx);
ax = axis; axis([650. 1100, ax(3), ax(4)])
title('CrIS ccast LW unapodized and apodized NEdN')
legend('unapodized', 'apodized', 'apod est', 'location', 'north')
xlabel('wavenumber')
ylabel('NEdN, mw sr-1 m-2')
grid on; zoom on

% MW detail
nedn_unap = d1.nMW(:,iFOV,1);
nedn_apod = d1.nMW(:,iFOV,2);
nedn_apxx = d1.nMW(:,iFOV,1) * 0.63;
freq = d1.vMW;

figure(3); clf
semilogy(freq, nedn_unap, freq, nedn_apod, freq,nedn_apxx);
ax = axis; axis([1200. 1750, ax(3), ax(4)])
title('CrIS ccast MW unapodized and apodized NEdN')
legend('unapodized', 'apodized', 'apod est', 'location', 'north')
xlabel('wavenumber')
ylabel('NEdN, mw sr-1 m-2')
grid on; zoom on

% SW detail
nedn_unap = d1.nSW(:,iFOV,1);
nedn_apod = d1.nSW(:,iFOV,2);
nedn_apxx = d1.nSW(:,iFOV,1) * 0.63;
freq = d1.vSW;

figure(4); clf
semilogy(freq, nedn_unap, freq, nedn_apod, freq,nedn_apxx);
ax = axis; axis([2150. 2550, ax(3), ax(4)])
title('CrIS ccast SW unapodized and apodized NEdN')
legend('unapodized', 'apodized', 'apod est', 'location', 'north')
xlabel('wavenumber')
ylabel('NEdN, mw sr-1 m-2')
grid on; zoom on

