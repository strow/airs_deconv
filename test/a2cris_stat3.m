%
% AIRS to CrIS stats ozone vs residuals
%

addpath ../source
addpath /asl/matlib/h4tools

% load profiles
frtp = '/asl/s1/motteler/kc7377/for_kcarta_241_pcrtm_clouds7.op.rtp';
[head, hattr, prof, pattr] = xrtpread(frtp);

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
cordep = acdep - (mdifdep * ones(1, ndep));
corind = acind - (mdifdep * ones(1, nind));

% residuals of corrected sets (dep is zero)
mcordep = mean(cordep - tcdep, 2);
mcorind = mean(corind - tcind, 2);
scorind =  std(corind - tcind, 0, 2);

% residuals for ozone band
o3ch = 561 : 710;
% o3rms = rms(acind(o3ch, :) -  tcind(o3ch, :));
o3rms = rms(corind(o3ch, :) -  tcind(o3ch, :));

% total column values
o3tot = sum(prof.gas_3(:, iind));
h2otot = sum(prof.gas_1(:, iind));

figure(1); clf
subplot(3,1,1)
plot(prof.rlat(iind))
title('latitude')
grid on;

subplot(3,1,2)
% plot(o3tot)
% title('total colunn ozone')
plot(h2otot)
title('total colunn water')
grid on;

subplot(3,1,3)
plot(o3rms)
title('1000 to 1090 cm-1 residual')
xlabel('obs index')
grid on

% saveas(gcf, 'h2o_resids_1', 'png')

return

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

