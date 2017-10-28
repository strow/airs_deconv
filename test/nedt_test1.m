%
% nedt_test1 - AIRS to CrIS NEdT stats
%
% subplot 1: NEdT for true AIRS, AIRS to CrIS, and true Cris
% subplot 2: max of true CrIS and AIRS-to-CrIS NEdN
%

addpath /asl/packages/ccast/source
addpath ../source

% load NEdN specs
c1 = load('nedn_noap');
c2 = load('nedn_hamm');

% load radiances
d1 = load('cris_fit49');   % true cris
d2 = load('airs_fit49');   % true airs
d3 = load('ac_ap_fit49');  % airs cris apodized
d4 = load('ac_cc_fit49');  % airs cris corrected 

% unapodized radiances
tarad = d2.arad;  tafrq = d2.afrq;  % true airs
adrad = d3.brad;  adfrq = d3.bfrq;  % airs decon

% apodized radiances
tcrad = [hamm_app(d1.radLW); hamm_app(d1.radMW); hamm_app(d1.radSW)];
tcfrq = [d1.frqLW; d1.frqMW; d1.frqSW];
acrad = d3.crad;  acfrq = d3.cfrq;  % apodized airs cris
ccrad = d4.crad;  ccfrq = d4.cfrq;  % corrected airs cris
% clear d1 d2 d3 d4

% match true CrIS grids
[ix, jx] = seq_match(c2.freq_cris, tcfrq);
nedn_cris = c2.nedn_cris(ix);
tcfrq = tcfrq(jx);
tcrad = tcrad(jx, :);

[x1, y1] = pen_lift(tafrq, dr2dt(tafrq, tarad(:, 1), c2.nedn_airs));
[x2, y2] = pen_lift(acfrq, dr2dt(acfrq, acrad(:, 1), c2.nedn_tran));
[x3, y3] = pen_lift(tcfrq, dr2dt(tcfrq, tcrad(:, 1), nedn_cris));

figure(1)
subplot(2,1,1)
plot(x1, y1, x2, y2, x3, y3)
axis([650, 2550, 0, 1])
title('AIRS to apodized CrIS NEdT')
legend('mean AIRS', 'AIRS to CrIS', 'mean CrIS')
ylabel('NEdT, K')
% xlabel('wavenumber')
grid on

% saveas(gcf, 'a2cris_nedt', 'fig')
% return

% plot max of true CrIS NEdN and AIRS CrIS NEdn
subplot(2,1,2)
d1 = load('nedn_noap');
d2 = load('nedn_hamm');
[ix, jx] = seq_match(d2.freq_tran, d2.freq_cris);
[x1, y1] = pen_lift(d2.freq_tran(ix), ...
                max(d2.nedn_tran(ix), d2.nedn_cris(jx)));
[x2, y2] = pen_lift(d2.freq_cris, d2.nedn_cris);
semilogy(x2, y2, x1, y1, 'linewidth', 2)
axis([650, 2550, 0, 1])
title('common record NEdN')
legend('mean CrIS', 'max (CrIS, AIRS to CrIS)')
xlabel('wavenumber')
ylabel('NEdN')
grid on; zoom on

saveas(gcf, 'a2cris_nedt', 'fig')

