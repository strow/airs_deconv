%
% AIRS to CrIS stats on the big cloudy and clear sets
%

addpath ../source

% load the data
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
cordep = acdep - (mdifdep * ones(1, ndep));
corind = acind - (mdifdep * ones(1, nind));

% residuals of corrected sets (dep is zero)
mcordep = mean(cordep - tcdep, 2);
mcorind = mean(corind - tcind, 2);
scorind =  std(corind - tcind, 0, 2);

figure(1);clf
subplot(2,1,1)
plot(tcfrq, mcorind)
axis([650, 1100, -2e-3, 2e-3])
title('mean residual corrected independent set')
ylabel('dTb')
grid on

subplot(2,1,2)
plot(tcfrq, scorind)
axis([650, 1100, 0, 0.1])
title('std residual corrected independent set')
ylabel('dTb')
grid on
% saveas(gcf, 'fig_3_ham_cloudy', 'png')

return

figure(2);clf
subplot(2,1,1)
plot(tcfrq, mdifind)
% axis([650, 1100, -1.5, 1.5])
  axis([650, 1100, -0.2, 0.2])
title('mean residual independent set')
ylabel('dTb')
grid on

subplot(2,1,2)
plot(tcfrq, mdifind - mdifdep)
axis([650, 1100, -0.005, 0.005])
title('independent minus dependent residuals')
xlabel('wavenumber')
ylabel('ddTb')
grid on
% saveas(gcf, 'fig_3_ham_cloudy', 'png')

figure(3); clf
subplot(2,1,1)
plot(tcfrq, mean(tcbt, 2))
axis([650, 1100, 210 270])
title('mean true CrIS all obs')
ylabel('Tb, K')
grid on

subplot(2,1,2)
plot(tcfrq, std(tcbt, 0, 2))
axis([650, 1100, 0, 30])
title('standard deviation all obs')
xlabel('wavenumber')
ylabel('Tb, K')
grid on
% saveas(gcf, 'fig_1_ham_cloudy', 'png')

figure(4); clf
subplot(2,1,1)
plot(tcfrq, mdifdep)
% axis([650, 1100, -1.5, 1.5])
  axis([650, 1100, -0.2, 0.2])
title('mean residual dependent set')
ylabel('dTb')
grid on

subplot(2,1,2)
plot(tcfrq, sdifdep)
% axis([650, 1100, 0, 0.6])
  axis([650, 1100, 0, 0.1])
title('std residual dependent set')
xlabel('wavenumber')
ylabel('dTb')
grid on
% saveas(gcf, 'fig_2_ham_cloudy', 'png')

