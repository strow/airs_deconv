%
% nedn_test2 - compare apodized and unapodized NEdN 
%
% NEdN for true AIRS, AIRS to CrIS, and true CrIS
% subplot 1: unapodized CrIS
% subplot 2: apoddized CrIS
%

d1 = load('nedn_noap');
d2 = load('nedn_hamm');

% plot apodized an unapodized NEdN
figure(1); clf
subplot(2,1,1)

[x1, y1] = pen_lift(d2.freq_airs, d2.nedn_airs);
[x2, y2] = pen_lift(d2.freq_tran, d2.nedn_tran);
[x3, y3] = pen_lift(d2.freq_cris, d2.nedn_cris);
semilogy(x1, y1, x2, y2, x3, y3)
axis([600, 2600, 0, 1])
title('AIRS to apodized CrIS NEdN')
legend('mean AIRS', 'AIRS to CrIS', 'mean CrIS')
ylabel('NEdN')
grid on; zoom on

subplot(2,1,2)
[x1, y1] = pen_lift(d1.freq_airs, d1.nedn_airs);
[x2, y2] = pen_lift(d1.freq_tran, d1.nedn_tran);
[x3, y3] = pen_lift(d1.freq_cris, d1.nedn_cris);
semilogy(x1, y1, x2, y2, x3, y3)
axis([600, 2600, 0, 1])
title('AIRS to unapodized CrIS NEdN')
legend('mean AIRS', 'AIRS to CrIS', 'mean CrIS')
xlabel('wavenumber')
ylabel('NEdN')
grid on; zoom on

saveas(gcf, 'a2cris_nedn', 'fig')

return

% plot max of true CrIS NEdN and AIRS CrIS NEdn
[ix, jx] = seq_match(d2.freq_tran, d2.freq_cris);
[x1, y1] = pen_lift(d2.freq_tran(ix), ...
                max(d2.nedn_tran(ix), d2.nedn_cris(jx)));
[x2, y2] = pen_lift(d2.freq_cris, d2.nedn_cris);
semilogy(x1, y1, x2, y2)
axis([600, 2600, 0, 1])
title('AIRS to apodized CrIS NEdN')
legend('max (mean cris, airs-to-cris)', 'mean cris')
xlabel('wavenumber')
ylabel('NEdN')
grid on; zoom on
