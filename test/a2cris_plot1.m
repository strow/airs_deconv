%
% combined apodized resiual plots from a2cris_tes1
%

d1 = load('diff_ap_LW');
d2 = load('diff_ap_MW');
d3 = load('diff_ap_SW');

figure(1); clf
set(gcf, 'Units','centimeters', 'Position', [4, 10, 24, 16])

% LW residual mean
acbt_ap = d1.acbt_ap;
tcbt_ap = d1.tcbt_ap;
tcfrq = d1.tcfrq;
mean(acbt_ap(:) - tcbt_ap(:))

subplot(3,1,1)
plot(tcfrq, mean(acbt_ap-tcbt_ap,2))
axis([650, 1100, -0.2, 0.2])
title('apodized AIRS CrIS minus true CrIS LW mean')
ylabel('dBT')
grid on; zoom on

% MW residual mean
acbt_ap = d2.acbt_ap;
tcbt_ap = d2.tcbt_ap;
tcfrq = d2.tcfrq;
mean(acbt_ap(:) - tcbt_ap(:))

subplot(3,1,2)
plot(tcfrq, mean(acbt_ap-tcbt_ap,2))
axis([1200, 1620, -0.2, 0.2])
title('apodized AIRS CrIS minus true CrIS MW mean')
ylabel('dBT')
grid on; zoom on

% SW residual mean 
acbt_ap = d3.acbt_ap;
tcbt_ap = d3.tcbt_ap;
tcfrq = d3.tcfrq;
mean(acbt_ap(:) - tcbt_ap(:))

subplot(3,1,3)
plot(tcfrq, mean(acbt_ap-tcbt_ap,2))
axis([2180, 2550, -0.2, 0.2])
title('apodized AIRS CrIS minus true CrIS SW mean')
xlabel('wavenumber')
ylabel('dBT')
grid on; zoom on

% export_fig('combo_ap_dif_mean.pdf', '-m2', '-transparent')

figure(2); clf
set(gcf, 'Units','centimeters', 'Position', [4, 10, 24, 16])

% LW residual std 
acbt_ap = d1.acbt_ap;
tcbt_ap = d1.tcbt_ap;
tcfrq = d1.tcfrq;
subplot(3,1,1)
plot(tcfrq, std(acbt_ap-tcbt_ap, 0, 2))
axis([650, 1100, 0, 0.1])
title('apodized AIRS CrIS minus true CrIS LW mean')
ylabel('dBT')
grid on; zoom on

% MW residual std 
acbt_ap = d2.acbt_ap;
tcbt_ap = d2.tcbt_ap;
tcfrq = d2.tcfrq;
subplot(3,1,2)
plot(tcfrq, std(acbt_ap-tcbt_ap, 0, 2))
axis([1200, 1620, 0, 0.1])
title('apodized AIRS CrIS minus true CrIS MW mean')
ylabel('dBT')
grid on; zoom on

% SW residual std 
acbt_ap = d3.acbt_ap;
tcbt_ap = d3.tcbt_ap;
tcfrq = d3.tcfrq;
subplot(3,1,3)
plot(tcfrq, std(acbt_ap-tcbt_ap, 0, 2))
axis([2180, 2550, 0, 0.1])
title('apodized AIRS CrIS minus true CrIS SW mean')
xlabel('wavenumber')
ylabel('dBT')
grid on; zoom on

% export_fig('combo_ap_dif_std.pdf', '-m2', '-transparent')

