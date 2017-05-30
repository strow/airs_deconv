%
% a2cris_regr2 - AIRS to CrIS Tb regression corrrection
%
% uses data from conv_loop1 and a2cris_loop;
%
% dependent set 7377 profile cloudy, independent set 49 fitting
% profiles
%

addpath ../source
addpath /asl/packages/ccast/source
addpath /home/motteler/matlab/export_fig
% addpath /asl/matlib/h4tools

% load radiance data
d1 = load('cris_cloudy');   % true cris big cloudy set
d2 = load('acris_cloudy');  % airs cris big cloudy set
d3 = load('cris_fit49');    % true cris 49 fitting profiles
d4 = load('acris_fit49');   % airs cris 49 fitting profiles

% set band
band = 'SW';
switch band
  case 'LW', tcfrq = d1.frqLW; tcdep = d1.radLW; tcind = d3.radLW;
  case 'MW', tcfrq = d1.frqMW; tcdep = d1.radMW; tcind = d3.radMW;
  case 'SW', tcfrq = d1.frqSW; tcdep = d1.radSW; tcind = d3.radSW;  
end
acfrq = d2.cfrq;
acdep = d2.crad;   
acind = d4.crad;   
clear d1 d2 d3 d4

% get the channel intersection
[tci, aci] = seq_match(tcfrq, acfrq);
tcfrq = tcfrq(tci); tcdep = tcdep(tci, :); tcind = tcind(tci, :);  
acfrq = acfrq(aci); acdep = acdep(aci, :); acind = acind(aci, :);  
[nchan, ndep] = size(tcdep);
[nchan, nind] = size(tcind);

% optional apodization
  tcdep = hamm_app(tcdep);
  acdep = hamm_app(acdep);
  tcind = hamm_app(tcind);
  acind = hamm_app(acind);

% get brightness temps
tcdep = real(rad2bt(tcfrq, tcdep));
acdep = real(rad2bt(acfrq, acdep));
tcind = real(rad2bt(tcfrq, tcind));
acind = real(rad2bt(acfrq, acind));

% mean and std of differences
mdifdep = mean(acdep - tcdep, 2);
mdifind = mean(acind - tcind, 2);
sdifdep = std(acdep - tcdep, 0, 2);
sdifind = std(acind - tcind, 0, 2);

% basic bias correction
cordep1 = acdep - (mdifdep * ones(1, ndep));
corind1 = acind - (mdifdep * ones(1, nind));

% linear and quadratic corrections
cordep2 = zeros(nchan, ndep);
corind2 = zeros(nchan, nind);
cordep3 = zeros(nchan, ndep);
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

% residuals of bias correction (dep is zero)
mcordep1 = mean(cordep1 - tcdep, 2);
mcorind1 = mean(corind1 - tcind, 2);
scorind1 =  std(corind1 - tcind, 0, 2);

% residuals of linear correction (dep is zero)
mcordep2 = mean(cordep2 - tcdep, 2);
mcorind2 = mean(corind2 - tcind, 2);
scorind2 =  std(corind2 - tcind, 0, 2);

% residuals of quadratic correction (dep is zero)
mcordep3 = mean(cordep3 - tcdep, 2);
mcorind3 = mean(corind3 - tcind, 2);
scorind3 =  std(corind3 - tcind, 0, 2);

figure(1); clf
% set(gcf, 'Units','centimeters', 'Position', [4, 10, 24, 16])
subplot(2,1,1)
plot(tcfrq, mcorind1, tcfrq, mcorind2, tcfrq, mcorind3)
switch band
  case 'LW', axis([650, 1100, -6e-2, 6e-2])
  case 'MW', axis([1200, 1620, -0.02, 0.02])
  case 'SW', axis([2180, 2550, -0.1, 0.1])
end
title('mean residual corrected independent set')
legend('bias correction', 'linear correction', 'quadratic correction', ...
       'location', 'north')
ylabel('dTb')
grid on

subplot(2,1,2)
plot(tcfrq, scorind1, tcfrq, scorind2, tcfrq, scorind3)
switch band
  case 'LW', axis([650, 1100, 0, 0.08])
  case 'MW', axis([1200, 1620, 0, 0.05])
  case 'SW', axis([2180, 2550, 0, 0.2])
end
legend('bias correction', 'linear correction', 'quadratic correction', ...
       'location', 'north')
title('std residual corrected independent set')
ylabel('dTb')
xlabel('wavenumber')
grid on
% saveas(gcf, sprintf('a2cris_regr_%s', band), 'png')
% export_fig(sprintf('a2cris_regr_%s.pdf', band), '-m2', '-transparent')

return

figure(2)
set(gcf, 'Units','centimeters', 'Position', [4, 10, 24, 16])
subplot(2,1,1)
plot(tcfrq, Pcor2(:, 1))
% axis([650, 1100, 0.995, 1.005])
title('"a" (scaling) weights')
grid on

subplot(2,1,2)
plot(tcfrq, Pcor2(:, 2))
% axis([650, 1100, -1, 1])
title('"b" (bias) weights')
xlabel('wavenumber')
grid on
% saveas(gcf, 'cor_ind_2', 'png')
% export_fig(sprintf('a2cris_coef_%s.pdf', band), '-m2', '-transparent')
