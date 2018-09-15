%
% AIRS to CrIS stats on Tsurf - Tb900cm-1 subsets
%

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

% index of 900 cm=1
i900 = 400;

% get brightness temps
tcbt = real(rad2bt(tcfrq, tcrad));
acbt = real(rad2bt(acfrq, acrad));

sdiff = prof.stemp(:) - tcbt(i900, :)';

figure(1); clf
histogram(sdiff)
title('Tsurf - Tb900cm-1')
% axis([-2, 18, 0, 800])
xlabel('dTb, K')
ylabel('obs')
grid on;
% saveas(gcf, 'window_diff', 'png')

ix1 = 0 <= sdiff & sdiff <= 6;
ix2 = 40 <= sdiff & sdiff <= 60;
ix3 = rand([1, nobs]) > 0.86;
ix4 = rand([1, nobs]) < 0.14;

mdif1 = mean(acbt(:, ix1) - tcbt(:, ix1), 2);
mdif2 = mean(acbt(:, ix2) - tcbt(:, ix2), 2);
mdif3 = mean(acbt(:, ix3) - tcbt(:, ix3), 2);
mdif4 = mean(acbt(:, ix4) - tcbt(:, ix4), 2);

sdif1 = std(acbt(:, ix1) - tcbt(:, ix1), 0, 2);
sdif2 = std(acbt(:, ix2) - tcbt(:, ix2), 0, 2);

[sum(ix1), sum(ix2), sum(ix3), sum(ix4)]
[rms(mdif1), rms(mdif2), rms(mdif3), rms(mdif4)]

figure(2); clf
subplot(2,1,1)
plot(tcfrq, mean(tcbt(:, ix1), 2), tcfrq, mean(tcbt(:, ix2), 2))
axis([650, 1100, 200, 280])
title('mean of window diff subsets')
legend('0 to 6 K', '40 to 60 K', 'location', 'best')
ylabel('Tb, K')
grid on;

subplot(2,1,2)
plot(tcfrq, std(tcbt(:, ix1), 0, 2), tcfrq, std(tcbt(:, ix2), 0, 2))
axis([650, 1100, 0, 30])
title('std of window diff subsets')
legend('0 to 6 K', '40 to 60 K', 'location', 'best')
ylabel('Tb, K')
grid on;
% saveas(gcf, 'fig2_wind_diff', 'png')

figure(3)
subplot(2,1,1)
plot(tcfrq, mdif1, tcfrq, mdif2)
axis([650, 1100, -2, 2])
title('means of residuals for two sets')
legend('0 to 6 K', '40 to 60 K', 'location', 'best')
ylabel('dTb, K')
grid on;

subplot(2,1,2)
plot(tcfrq, sdif1, tcfrq, sdif2)
axis([650, 1100, 0, 0.6])
title('stds of residuals for two sets')
legend('0 to 6 K', '40 to 60 K', 'location', 'best')
ylabel('dTb, K')
grid on;
% saveas(gcf, 'fig3_wind_diff', 'png')

figure(4); clf
subplot(2,1,1)
plot(tcfrq, mdif1 - mdif2)
axis([650, 1100, -0.25, 0.25])
title('0 to 6 K minus 40 to 60 K')
ylabel('ddTb, K')
grid on; zoom on

subplot(2,1,2)
plot(tcfrq, mdif3 - mdif4)
axis([650, 1100, -0.025, 0.025])
title('random set 1 minus random set 2')
ylabel('ddTb, K')
grid on; zoom on
% saveas(gcf, 'fig4_wind_diff', 'png')

