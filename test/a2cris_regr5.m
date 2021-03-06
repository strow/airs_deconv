%
% a2cris_regr5 - AIRS L1c to CrIS PC regression test and plots
%
% uses data from conv_loop4
%   nkcd  - dependent set size
%   nkci  - independent set size
%   va1C  - AIRS channel frequency
%   na1C  - number of AIRS channels
%   a1Crd - AIRS dep set radiance
%   a1Cri - AIRS ind set radiance
%   vcLW,  vcMW,  vcSW  - CrIS channel frequency
%   ncLW,  ncMW,  ncSW  - number of CrIS channels
%   cLWrd, cMWrd, cSWrd - CrIS dep set radiance
%   cLWri, cMWri, cSWri - CrIS ind set radiance
%

addpath /asl/packages/ccast/source
addpath ../source

% colormap for 2D plots
load llsmap5

% correlation color axis
cax = [-0.4, 0.4];

% get radiance data
load('conv_loop4')

% sample ccast CrIS NEdN estimate
% d = load('/asl/data/cris/ccast/sdr60_npp_LR/2016/018/SDR_d20160118_t0801033.mat');
% nedn_crisLW = squeeze(d.nLW(3:end-2, 5, 1));
% nedn_crisMW = squeeze(d.nMW(3:end-2, 5, 1));
% nedn_crisSW = squeeze(d.nSW(3:end-2, 5, 1));
% clear d

% apodized CrIS noise
% nedn_cris = nedn_cris * 0.63;

% add noise to the true CrIS dependent set
% cLWrd = cLWrd + randn(ncLW, nkcd) .* (nedn_crisLW * ones(1, nkcd));
% cMWrd = cMWrd + randn(ncMW, nkcd) .* (nedn_crisMW * ones(1, nkcd));
% cSWrd = cSWrd + randn(ncSW, nkcd) .* (nedn_crisSW * ones(1, nkcd));

% load the AIRS L1c noise spec
d1 = load('int_L1c_NEdN.mat');
nedn_airs = double(d1.nedn);
% freq_airs = d1.v1c;

% add noise to the dependent set
% a1Crd = a1Crd + randn(na1C, nkcd) .* (nedn_airs * ones(1, nkcd));

% trim CrIS bands to the intersection 
jMW = find(vcMW < 1614);
jSW = find(2182 < vcSW);
vcMW = vcMW(jMW); vcSW = vcSW(jSW);
ncMW = length(vcMW); ncSW = length(vcSW);
cMWrd = cMWrd(jMW, :); cMWri = cMWri(jMW, :);
cSWrd = cSWrd(jSW, :); cSWri = cSWri(jSW, :);

% AIRS spans for CrIS bands
iLW = find(vcLW(1)-4 <= va1C & va1C <= vcLW(end)+4);
iMW = find(vcMW(1)-4 <= va1C & va1C <= vcMW(end)+4);
iSW = find(vcSW(1)-4 <= va1C & va1C <= vcSW(end)+4);

% LW SVDs
[aU,~,~] = svd(a1Crd(iLW,:),0);
[cU,~,~] = svd(cLWrd, 0);
aU = aU(:, 1:500); % max 1266
cU = cU(:, 1:500); % max 713
QLW = (cU' * cLWrd) / (aU' * a1Crd(iLW,:));
RLW = cU * QLW * aU';
% RLW =  (cU * cU' * cLWrd) / (aU * aU' * a1Crd(iLW,:));

% MW SVDs
[aU,~,~] = svd(a1Crd(iMW,:),0);
[cU,~,~] = svd(cMWrd, 0);
aU = aU(:, 1:500); % max 692
cU = cU(:, 1:320); % max 324
QMW = (cU' * cMWrd) / (aU' * a1Crd(iMW,:));
RMW = cU * QMW * aU';
% RMW =  (cU * cU' * cMWrd) / (aU * aU' * a1Crd(iMW,:));

% SW SVDs
[aU,~,~] = svd(a1Crd(iSW,:),0);
[cU,~,~] = svd(cSWrd, 0);
aU = aU(:, 1:100); % max 377
cU = cU(:, 1:100); % max 148
QSW = (cU' * cSWrd) / (aU' * a1Crd(iSW,:));
RSW = cU * QSW * aU';
% RSW =  (cU * cU' * cSWrd) / (aU * aU' * a1Crd(iSW,:));

% do the regression
% RLW = band_regr(a1Crd(iLW,:), cLWrd, va1C(iLW), vcLW, 120);
% RLW = cLWrd / a1Crd(iLW,:);
% RMW = cMWrd / a1Crd(iMW,:);
% RSW = cSWrd / a1Crd(iSW,:);

% apply the regression
acLWrd = RLW * a1Crd(iLW,:);  acLWri = RLW * a1Cri(iLW,:);
acMWrd = RMW * a1Crd(iMW,:);  acMWri = RMW * a1Cri(iMW,:);
acSWrd = RSW * a1Crd(iSW,:);  acSWri = RSW * a1Cri(iSW,:);

% unapodized brightness temps
cLWbd = real(rad2bt(vcLW, cLWrd));  cLWbi = real(rad2bt(vcLW, cLWri));
cMWbd = real(rad2bt(vcMW, cMWrd));  cMWbi = real(rad2bt(vcMW, cMWri));
cSWbd = real(rad2bt(vcSW, cSWrd));  cSWbi = real(rad2bt(vcSW, cSWri));

acLWbd = real(rad2bt(vcLW, acLWrd)); acLWbi = real(rad2bt(vcLW, acLWri));
acMWbd = real(rad2bt(vcMW, acMWrd)); acMWbi = real(rad2bt(vcMW, acMWri));
acSWbd = real(rad2bt(vcSW, acSWrd)); acSWbi = real(rad2bt(vcSW, acSWri));

% apodized brightness temps
cLWai =  real(rad2bt(vcLW, hamm_app(cLWri)));
cMWai =  real(rad2bt(vcMW, hamm_app(cMWri)));
cSWai =  real(rad2bt(vcSW, hamm_app(cSWri)));

acLWai = real(rad2bt(vcLW, hamm_app(acLWri)));
acMWai = real(rad2bt(vcMW, hamm_app(acMWri)));
acSWai = real(rad2bt(vcSW, hamm_app(acSWri)));

% unapodized stats
mdifLWbd = mean(acLWbd - cLWbd, 2);    
mdifMWbd = mean(acMWbd - cMWbd, 2);   
mdifSWbd = mean(acSWbd - cSWbd, 2);

sdifLWbd = std(acLWbd - cLWbd, 0, 2); 
sdifMWbd = std(acMWbd - cMWbd, 0, 2); 
sdifSWbd = std(acSWbd - cSWbd, 0, 2);

mdifLWbi = mean(acLWbi - cLWbi, 2);
mdifMWbi = mean(acMWbi - cMWbi, 2);
mdifSWbi = mean(acSWbi - cSWbi, 2);

sdifLWbi = std(acLWbi - cLWbi, 0, 2);
sdifMWbi = std(acMWbi - cMWbi, 0, 2);
sdifSWbi = std(acSWbi - cSWbi, 0, 2);

% apodized stats
mdifLWai = mean(acLWai - cLWai, 2);
mdifMWai = mean(acMWai - cMWai, 2);
mdifSWai = mean(acSWai - cSWai, 2);

sdifLWai = std(acLWai - cLWai, 0, 2);
sdifMWai = std(acMWai - cMWai, 0, 2);
sdifSWai = std(acSWai - cSWai, 0, 2);

%---------------------
% mean residual stats
%---------------------
% dep and ind mean unapodized residuals
figure(1); clf
subplot(3,1,1)
plot(vcLW, mdifLWbi, vcLW, mdifLWbd)
axis([650, 1100, -0.1, 0.1])
title('unapodized direct regression mean residuals')
legend('ind set', 'dep set', 'location', 'north')
ylabel('\Delta BT(K)')
grid on; zoom on

subplot(3,1,2)
plot(vcMW, mdifMWbi, vcMW, mdifMWbd)
axis([1210, 1610, -0.1, 0.1])
legend('ind set', 'dep set', 'location', 'north')
ylabel('\Delta BT (K)')
grid on; zoom on

subplot(3,1,3)
plot(vcSW, mdifSWbi, vcSW, mdifSWbd)
axis([2180, 2550, -0.4, 0.4])
legend('ind set', 'dep set', 'location', 'north')
xlabel('wavenumber (cm^{-1})')
ylabel('\Delta BT(K)')
grid on; zoom on

% dep set max values
% max(abs(mdifLWbd)), max(abs(mdifMWbd)), max(abs(mdifSWbd))

% ind set mean apodized residuals
figure(2); clf
subplot(3,1,1)
plot(vcLW, mdifLWai)
axis([650, 1100, -5e-2, 5e-2])
title('apodized direct regression mean residuals')
ylabel('\Delta BT (K)')
grid on; zoom on

subplot(3,1,2)
plot(vcMW, mdifMWai)
axis([1210, 1610, -5e-2, 5e-2])
ylabel('\Delta BT (K)')
grid on; zoom on

subplot(3,1,3)
plot(vcSW, mdifSWai)
axis([2180, 2550, -2e-1, 2e-1])
xlabel('wavenumber (cm^{-1})')
ylabel('\Delta BT (K)')
grid on; zoom on
% saveas(gcf, 'ap_pc_regr', 'fig')

%--------------------------
% plot regression matrices
%--------------------------
figure(3); clf
pcolor(va1C(iLW), vcLW, RLW)
shading flat
caxis(cax)
colormap(llsmap5)
colorbar
title('AIRS to CrIS LW regression matrix')
xlabel('AIRS wavenumber (cm^{-1})')
ylabel('CrIS wavenumber (cm^{-1})')
grid on
% saveas(gcf, 'LW_pc_regr_mat', 'png')

figure(4); clf
pcolor(va1C(iMW), vcMW, RMW)
shading flat
caxis(cax)
colormap(llsmap5)
colorbar
title('AIRS to CrIS MW regression matrix')
xlabel('AIRS wavenumber (cm^{-1})')
ylabel('CrIS wavenumber (cm^{-1})')
grid on
% saveas(gcf, 'MW_pc_regr_mat', 'png')

figure(5); clf
pcolor(va1C(iSW), vcSW, RSW)
shading flat
caxis(cax)
colormap(llsmap5)
colorbar
title('AIRS to CrIS SW regression matrix')
xlabel('AIRS wavenumber (cm^{-1})')
ylabel('CrIS wavenumber (cm^{-1})')
grid on
% saveas(gcf, 'SW_pc_regr_mat', 'png')


