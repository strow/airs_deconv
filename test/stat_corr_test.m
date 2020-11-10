%
% some more statistical correction tests
%

% new and old corrections
c1 = load('corr_midres');
c2 = load('corr_midres_old');

% midres apodized, with new stat corr, new SRFs
d1 = load('ac_cc_MR_srf_49');

% midres apodized, with old stat corr, new SRFs
d2 = load('ac_cc_MR_old_49');

% get apodized BT
cfrq = d1.cfrq;
bt1 = rad2bt(cfrq, d1.crad);
bt2 = rad2bt(cfrq, d2.crad);

figure(1)
plot(cfrq, bt1(:,1) - bt2(:,1))
axis([600, 2600, -0.09, 0.03])
title('old and new stat corr with new srfs')
xlabel('wavenumber (cm-1)')
ylabel('dBT (K)')
grid on
saveas(gcf, 'old_and_new_stat_corr', 'fig')

return

% coeff diffs, index 1 is "a", index 2 is "b"
plot(c1.tcfrqSW, c1.Pcor2SW(:,1) - c2.Pcor2SW(:,1))
plot(c1.tcfrqSW, c1.Pcor2SW(:,2) - c2.Pcor2SW(:,2))

