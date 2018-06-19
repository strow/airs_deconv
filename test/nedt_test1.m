%
% nedt_test1 - AIRS to CrIS NEdT stats
%
% subplot 1: NEdT for true AIRS, AIRS to CrIS, and true Cris
% subplot 2: NEdT for max of true CrIS and AIRS-to-CrIS NEdN
%
% uses data from conv_loop1, a2cris_loop, and apodized data from
% nedn_test1
%

addpath /asl/packages/ccast/source
addpath ../source

% load NEdN specs
c2 = load('nedn_hamm');    % AIRS, apodized CrIS and AIRS-to-CrIS

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

% match CrIS NEdN and radiance grids
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
axis([650, 2550, 0, 0.8])
title('apodized AIRS to CrIS NEdT')
legend('mean AIRS', 'AIRS to CrIS', 'mean CrIS', 'location', 'north')
ylabel('NEdT (K)')
grid on

% plot max of true CrIS NEdN and AIRS CrIS NEdn
subplot(2,1,2)
[ix, jx] = seq_match(c2.freq_tran, c2.freq_cris);
nedn_max = max(c2.nedn_tran(ix), c2.nedn_cris(jx));
freq_max = c2.freq_tran(ix);

[ix, jx] = seq_match(freq_max, tcfrq);
nedn_max = nedn_max(ix);
freq_max = freq_max(ix);
rad_max = tcrad(jx, :);

[x1, y1] = pen_lift(freq_max, dr2dt(freq_max, rad_max(:, 1), nedn_max));
[x2, y2] = pen_lift(tcfrq, dr2dt(tcfrq, tcrad(:, 1), nedn_cris));

plot(x2, y2, x1, y1)
axis([650, 2550, 0, 0.8])
title('apodized common record NEdT')
legend('mean CrIS', 'max CrIS, AIRS-to-CrIS', 'location', 'north')
xlabel('wavenumber (cm^{-1})')
ylabel('NEdT (K)')
grid on; zoom on

saveas(gcf, 'a2cris_nedt', 'fig')

