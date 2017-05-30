%
% a2cris_regr2 - AIRS L1c to CrIS direct regression tests and plots
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

% get radiance data
load('conv_loop4')

% trim CrIS bands to the intersection 
jMW = find(vcMW < 1614);
jSW = find(2182 < vcSW);
vcMW = vcMW(jMW); vcSW = vcSW(jSW);
ncMW = length(vcMW); ncSW = length(vcSW);
cMWrd = cMWrd(jMW, :); cMWri = cMWri(jMW, :);
cSWrd = cSWrd(jSW, :); cSWri = cSWri(jSW, :);

% % split the big dependent set 
% % rstate = rng; save rstate rstate
% load rstate
% rng(rstate)
% % ixd = rand(1, nkcd) > 0.1;
% ixd = logical(randi([0,1], [1, nkcd]));
% ixi = ~ixd;
% nkcd = sum(ixd);
% nkci = sum(ixi); 
% 
% a1Cr = a1Crd; a1Crd = a1Cr(:, ixd);  a1Cri = a1Cr(:, ixi);
% cLWr = cLWrd; cLWrd = cLWr(:, ixd);  cLWri = cLWr(:, ixi);
% cMWr = cMWrd; cMWrd = cMWr(:, ixd);  cMWri = cMWr(:, ixi);
% cSWr = cSWrd; cSWrd = cSWr(:, ixd);  cSWri = cSWr(:, ixi); 
% clear a1Cr cLWr cMWr cSWr

% AIRS spans for CrIS bands
iLW = find(vcLW(1)-4 <= va1C & va1C <= vcLW(end)+4);
iMW = find(vcMW(1)-4 <= va1C & va1C <= vcMW(end)+4);
iSW = find(vcSW(1)-4 <= va1C & va1C <= vcSW(end)+4);

% do the regression
% RLW = band_regr(a1Crd(iLW,:), cLWrd, va1C(iLW), vcLW, 120);
RLW = (a1Crd(iLW,:)' \ cLWrd')';
RMW = (a1Crd(iMW,:)' \ cMWrd')';
RSW = (a1Crd(iSW,:)' \ cSWrd')';

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

% dep/ind residuals for 7377 set split
figure(1); clf
% set(gcf, 'Units','centimeters', 'Position', [4, 10, 24, 16])
subplot(3,1,1)
plot(vcLW, mdifLWbi, vcLW, mdifLWbd)
% axis([650, 1100, -5e-7, 5e-7])
% axis([650, 1100, -1e-6, 1e-6])
title('direct regression unapodized mean residuals')
legend('ind set', 'dep set', 'location', 'north')
ylabel('dTb, K')
grid on; zoom on

subplot(3,1,2)
plot(vcMW, mdifMWbi, vcMW, mdifMWbd)
% axis([1200, 1620, -5e-7, 5e-7])
% axis([1200, 1620, -1e-6, 1e-6])
legend('ind set', 'dep set', 'location', 'north')
ylabel('dTb, K')
grid on; zoom on

subplot(3,1,3)
plot(vcSW, mdifSWbi, vcSW, mdifSWbd)
% axis([2180, 2550, -5e-6, 5e-6])
% axis([2180, 2550, -1e-4, 1e-4])
legend('ind set', 'dep set', 'location', 'north')
xlabel('wavenumber')
ylabel('dTb, K')
grid on; zoom on

% dep set max values
% max(abs(mdifLWbd)), max(abs(mdifMWbd)), max(abs(mdifSWbd))

figure(2); clf
% set(gcf, 'Units','centimeters', 'Position', [4, 10, 24, 16])
pcolor(va1C(iLW), vcLW, RLW)
shading flat
caxis([-0.5, 0.5])
colormap(llsmap5)
colorbar
title('AIRS to CrIS LW regression matrix')
xlabel('AIRS channels')
ylabel('CrIS channels')
grid on

figure(3); clf
% set(gcf, 'Units','centimeters', 'Position', [4, 10, 24, 16])
pcolor(va1C(iMW), vcMW, RMW)
shading flat
caxis([-1, 1])
colormap(llsmap5)
colorbar
title('AIRS to CrIS MW regression matrix')
xlabel('AIRS channels')
ylabel('CrIS channels')
grid on
% saveas(gcf, 'split_7377_MW_regr_mat', 'png')

figure(4); clf
% set(gcf, 'Units','centimeters', 'Position', [4, 10, 24, 16])
pcolor(va1C(iSW), vcSW, RSW)
shading flat
caxis([-0.5, 0.5])
colormap(llsmap5)
colorbar
title('AIRS to CrIS SW regression matrix')
xlabel('AIRS channels')
ylabel('CrIS channels')
grid on

return

% basic plots
figure(2); clf
subplot(2,1,1)
plot(vcLW, mdifLWbi, vcLW, mdifLWai)
% axis([650, 1100, -0.1, 0.1])
title('direct regression mean residuals')
legend('unapodized', 'Hamming', 'location', 'north')
ylabel('dTb, K')
grid on; zoom on

subplot(2,1,2)
plot(vcLW, sdifLWbi, vcLW, sdifLWai)
% axis([650, 1100, 0, 0.03])
title('direct regression std residuals')
legend('unapodized', 'Hamming', 'location', 'north')
xlabel('wavenumber')
ylabel('dTb, K')
grid on; zoom on

