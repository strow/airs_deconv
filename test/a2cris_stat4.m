%
% AIRS to CrIS regression tests
%

addpath ../source
addpath /asl/matlib/h4tools

% load radiance data
d1 = load('cris_cloudy');   % true cris
d2 = load('acris_cloudy');  % airs cris

tcrad = d1.radLW;  tcfrq = d1.frqLW;
acrad = d2.crad;   acfrq = d2.cfrq;
clear d1 d2

% get the cnannel intersection
[tci, aci] = seq_match(tcfrq, acfrq);
tcrad = tcrad(tci, :);  tcfrq = tcfrq(tci);
acrad = acrad(aci, :);  acfrq = acfrq(aci);
[nchan, nobs] = size(tcrad);

% optional apodization
tcrad = hamm_app(tcrad);
acrad = hamm_app(acrad);

% get brightness temps
tcbt = real(rad2bt(tcfrq, tcrad));
acbt = real(rad2bt(acfrq, acrad));

% split into dependent and independent sets
% load rstate
rstate = rng;
rng(rstate)
idep = logical(randi([0,1], [1, nobs]));
iind = ~idep;

acind = acbt(:, iind);
acdep = acbt(:, idep);
tcind = tcbt(:, iind);
tcdep = tcbt(:, idep);

% mean and std of differences
mdifdep = mean(acdep - tcdep, 2);
mdifind = mean(acind - tcind, 2);
sdifdep = std(acdep - tcdep, 0, 2);
sdifind = std(acind - tcind, 0, 2);

% use mdifdep as a correction
[~, nind] = size(acind);
[~, ndep] = size(acdep);

% corrected dependend and independent sets
cordep1 = acdep - (mdifdep * ones(1, ndep));
corind1 = acind - (mdifdep * ones(1, nind));

% residuals of bias correction (dep is zero)
mcordep1 = mean(cordep1 - tcdep, 2);
mcorind1 = mean(corind1 - tcind, 2);
scorind1 =  std(corind1 - tcind, 0, 2);

% per-channel polynomial fits
cordep2 = zeros(nchan, ndep);
cordep3 = zeros(nchan, ndep);
corind2 = zeros(nchan, nind);
corind3 = zeros(nchan, nind);
Pcor2 = zeros(nchan, 2);
Pcor3 = zeros(nchan, 3);

for i = 1 : nchan
  Ptmp = polyfit(acdep(i, :), tcdep(i, :), 1);
  cordep2(i, :) = polyval(Ptmp, acdep(i, :));
  corind2(i, :) = polyval(Ptmp, acind(i, :));
  Pcor2(i, :) = Ptmp;

  Ptmp = polyfit(acdep(i, :), tcdep(i, :), 2);
  cordep3(i, :) = polyval(Ptmp, acdep(i, :));
  corind3(i, :) = polyval(Ptmp, acind(i, :));
  Pcor3(i, :) = Ptmp;
end

% residuals of linear correction (dep is zero)
mcordep2 = mean(cordep2 - tcdep, 2);
mcorind2 = mean(corind2 - tcind, 2);
scorind2 =  std(corind2 - tcind, 0, 2);

% residuals of quadratic correction (dep is zero)
mcordep3 = mean(cordep3 - tcdep, 2);
mcorind3 = mean(corind3 - tcind, 2);
scorind3 =  std(corind3 - tcind, 0, 2);

figure(1); clf
subplot(2,1,1)
plot(tcfrq, mcorind1, tcfrq, mcorind2, tcfrq, mcorind3)
axis([650, 1100, -2e-3, 2e-3])
title('mean residual corrected independent set')
legend('bias correction', 'linear correction', 'quadratic correction', ...
       'location', 'north')
ylabel('dTb')
grid on

subplot(2,1,2)
plot(tcfrq, scorind1, tcfrq, scorind2, tcfrq, scorind3)
axis([650, 1100, 0, 0.1])
legend('bias correction', 'linear correction', 'quadratic correction', ...
       'location', 'north')
title('std residual corrected independent set')
ylabel('dTb')
xlabel('wavenumber')
grid on
% saveas(gcf, 'cor_ind_1', 'png')

figure(2)
subplot(2,1,1)
plot(tcfrq, Pcor2(:, 1))
axis([650, 1100, 0.99, 1.01])
title('"a" (scaling) weights')
grid on

subplot(2,1,2)
plot(tcfrq, Pcor2(:, 2))
axis([650, 1100, -1, 1])
title('"b" (bias) weights')
xlabel('wavenumber')
grid on
% saveas(gcf, 'cor_ind_2', 'png')

