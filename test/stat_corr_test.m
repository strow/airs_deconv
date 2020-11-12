%
% some more statistical correction tests
%

addpath ../data

% midres apodized, with old stat corr, new SRFs
d1 = load('ac_c1_MR_srf_49');

% midres apodized, with new stat corr, new SRFs
d2 = load('ac_c2_MR_srf_49');

% get apodized BT
cfrq = d1.cfrq;
bt1 = rad2bt(cfrq, d1.crad);
bt2 = rad2bt(cfrq, d2.crad);

figure(1)
plot(cfrq, bt2(:,1) - bt1(:,1))
axis([600, 2600, -0.09, 0.03])
title('new minus old stat corr with new srfs')
xlabel('wavenumber (cm-1)')
ylabel('dBT (K)')
grid on
% saveas(gcf, 'new_minus_old_stat_corr', 'fig')

% return

% old and new corrections
c1 = load('corr_midres_v1');
c2 = load('corr_midres_v2');

% coeff diffs, index 1 is "a", index 2 is "b"
figure(2)
subplot(2,1,1)
plot(c1.tcfrqSW, c2.Pcor2SW(:,1) - c1.Pcor2SW(:,1))
xlim([2180,2550])
title('"a" coeff difference, new minus old')
grid on

subplot(2,1,2)
plot(c1.tcfrqSW, c2.Pcor2SW(:,2) - c1.Pcor2SW(:,2))
xlim([2180,2550])
title('"b" coeff difference, new minus old')
xlabel('wavenumber (cm-1)')
grid on

