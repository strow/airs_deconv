%
% a2cris_test2 - compare AIRS CrIS with true CrIS
%
% 3-band summary plots of the airs2cris statistical correction,
% uses convolved data from conv_loop1 and a2cris_loop
%

addpath ../source
addpath /asl/packages/ccast/source
addpath /home/motteler/matlab/export_fig

% load the data
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
clear d1 d2 d3

% true/airs cris channel intersection
[tci, aci] = seq_match(tcfrq, acfrq);
tcrad = tcrad(tci, :);  tcfrq = tcfrq(tci);
acrad = acrad(aci, :);  acfrq = acfrq(aci);
[nchan, nobs] = size(tcrad);

% get brightness temps
tcbt = real(rad2bt(tcfrq, tcrad));
tabt = real(rad2bt(tafrq, tarad));
acbt = real(rad2bt(acfrq, acrad));
adbt = real(rad2bt(adfrq, adrad));
ccbt = real(rad2bt(ccfrq, ccrad));

% 3-band corrected and uncorrected mean residuals 
figure(1); clf
subplot(3,1,1)
plot(tcfrq, mean(acbt-tcbt,2), tcfrq, mean(ccbt-tcbt,2))
axis([650, 1100, -0.2, 0.2])
title('LW AIRS CrIS minus true CrIS mean')
legend('uncorrected', 'corrected', 'location', 'north')
ylabel('dBT')
grid on; zoom on

subplot(3,1,2)
plot(tcfrq, mean(acbt-tcbt,2), tcfrq, mean(ccbt-tcbt,2))
axis([1210, 1610, -0.2, 0.2])
title('MW AIRS CrIS minus true CrIS mean')
legend('uncorrected', 'corrected', 'location', 'north')
ylabel('dBT')
grid on; zoom on

subplot(3,1,3)
plot(tcfrq, mean(acbt-tcbt,2), tcfrq, mean(ccbt-tcbt,2))
axis([2180, 2550, -0.2, 0.2])
title('SW AIRS CrIS minus true CrIS mean')
legend('uncorrected', 'corrected', 'location', 'northeast')
xlabel('wavenumber')
ylabel('dBT')
grid on; zoom on
saveas(gcf, 'a2cris_regr_all', 'fig')

% 3-band corrected mean residual zoom
figure(2); clf
subplot(3,1,1)
plot(tcfrq, mean(ccbt-tcbt,2))
axis([650, 1100, -0.04, 0.04])
title('corrected AIRS CrIS minus true CrIS mean')
ylabel('dBT')
grid on; zoom on

subplot(3,1,2)
plot(tcfrq, mean(ccbt-tcbt,2))
axis([1210, 1610, -0.01, 0.01])
ylabel('dBT')
grid on; zoom on

subplot(3,1,3)
plot(tcfrq, mean(ccbt-tcbt,2))
axis([2180, 2550, -0.06, 0.06])
xlabel('wavenumber')
ylabel('dBT')
grid on; zoom on
saveas(gcf, 'ap_decon_corr', 'fig')

return

% 3-band std of AIRS CrIS minus true CrIS
figure(2); clf
subplot(3,1,1)
plot(tcfrq, std(acbt-tcbt,0,2), tcfrq, std(ccbt-tcbt,0,2))
axis([650, 1095, 0, 0.1])
title('LW AIRS CrIS minus true CrIS std')
legend('uncorrected', 'corrected', 'location', 'north')
ylabel('dBT')
grid on; zoom on

subplot(3,1,2)
plot(tcfrq, std(acbt-tcbt,0,2), tcfrq, std(ccbt-tcbt,0,2))
axis([1210, 1605, 0, 0.1])
title('MW AIRS CrIS minus true CrIS std')
legend('uncorrected', 'corrected', 'location', 'north')
ylabel('dBT')
grid on; zoom on

subplot(3,1,3)
plot(tcfrq, std(acbt-tcbt,0,2), tcfrq, std(ccbt-tcbt,0,2))
axis([2180, 2550, 0, 0.1])
title('SW AIRS CrIS minus true CrIS std')
legend('uncorrected', 'corrected', 'location', 'northeast')
xlabel('wavenumber')
ylabel('dBT')
grid on; zoom on

