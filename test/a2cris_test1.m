%
% a2cris_test1 - compare AIRS CrIS with true CrIS
%
% derived from a2cris_stat1 and cris_test5, uses convolved data from
% conv_loop and a2cris_loop
%

addpath ../source
addpath /asl/packages/ccast/source
addpath /home/motteler/matlab/export_fig

% load the data
d1 = load('cris_fit49');   % true cris
d2 = load('airs_fit49');   % true airs
d3 = load('acris_fit49');  % airs cris

% set band
band = 'LW';
switch band % select true cris by band
  case 'LW', tcrad = d1.radLW;  tcfrq = d1.frqLW;
  case 'MW', tcrad = d1.radMW;  tcfrq = d1.frqMW;
  case 'SW', tcrad = d1.radSW;  tcfrq = d1.frqSW;
end
tarad = d2.arad;   tafrq = d2.afrq;    % true airs
acrad = d3.crad;   acfrq = d3.cfrq;    % airs cris
adrad = d3.brad;   adfrq = d3.bfrq;    % airs decon
clear d1 d2 d3

% get the channel intersection
[tci, aci] = seq_match(tcfrq, acfrq);
tcrad = tcrad(tci, :);  tcfrq = tcfrq(tci);
acrad = acrad(aci, :);  acfrq = acfrq(aci);
[nchan, nobs] = size(tcrad);

% hamming apodization
tcrad_ap = hamm_app(tcrad);
acrad_ap = hamm_app(acrad);

% get brightness temps
tcbt = real(rad2bt(tcfrq, tcrad));
tabt = real(rad2bt(tafrq, tarad));
acbt = real(rad2bt(acfrq, acrad));
adbt = real(rad2bt(adfrq, adrad));
tcbt_ap = real(rad2bt(tcfrq, tcrad_ap));
acbt_ap = real(rad2bt(tcfrq, acrad_ap));

%----------------------------
% plot spectra and residuals
%----------------------------

% plot all data for a selected spectra
figure(1); clf; j = 1;
set(gcf, 'Units','centimeters', 'Position', [4, 10, 24, 16])
subplot(2,1,1)
plot(tafrq, tabt(:,j), adfrq, adbt(:,j), tcfrq, tcbt(:,j), acfrq, acbt(:,j))
  axis([650, 1100, 200, 310])
% axis([1200, 1620, 210, 310])
% axis([2180, 2550, 210, 310])
title(sprintf('AIRS and unapodized CrIS %s profile %d', band, j));
legend('true AIRS', 'AIRS decon', 'true CrIS', 'AIRS CrIS', ...
       'location', 'south')
ylabel('Tb, K')
grid on; zoom on

subplot(2,1,2)
plot(tafrq, tabt(:,j), adfrq, adbt(:,j), tcfrq, tcbt(:,j), acfrq, acbt(:,j))
  axis([660, 680, 200, 260])
% axis([1320, 1350, 210, 290])
% axis([2320, 2360, 210, 260])
title(sprintf('AIRS and CrIS %s profile %d, detail', band, j));
legend('true AIRS', 'AIRS decon', 'true CrIS', 'AIRS CrIS', ...
       'location', 'northeast')
xlabel('wavenumber'); ylabel('Tb, K')
grid on; zoom on
% export_fig(sprintf('a2cris_spec_%s.pdf', band), '-m2', '-transparent')

% plot mean and std of AIRS CrIS minus true CrIS
figure(2); clf
set(gcf, 'Units','centimeters', 'Position', [4, 10, 24, 16])
subplot(2,1,1)
plot(tcfrq, mean(acbt-tcbt,2), tcfrq, mean(acbt_ap-tcbt_ap,2))
  axis([650, 1100, -1.0, 1.0])
% axis([1200, 1620, -0.3, 0.3])
% axis([2180, 2550, -10, 15])
title(sprintf('AIRS CrIS minus true CrIS %s mean', band));
  legend('unapodized', 'Hamming ap.', 'location', 'north') 
% legend('unapodized', 'Hamming ap.', 'location', 'northeast') 
ylabel('dBT')
grid on; zoom on

% AIRS CrIS minus true CrIS std
subplot(2,1,2)
plot(tcfrq, std(acbt-tcbt,0,2), tcfrq, std(acbt_ap-tcbt_ap,0,2))
  axis([650, 1100, 0, 0.5])
% axis([1200, 1620, 0, 0.12])
% axis([2180, 2550, 0, 10])
title(sprintf('AIRS CrIS minus true CrIS %s std dev', band));
  legend('unapodized', 'Hamming ap.', 'location', 'north')
% legend('unapodized', 'Hamming ap.', 'location', 'northeast') 
xlabel('wavenumber')
ylabel('dBT')
grid on; zoom on
% export_fig(sprintf('a2cris_diff_%s.pdf', band), '-m2', '-transparent')

% save basic info for combined plots
save(sprintf('diff_ap_%s', band), 'tcfrq', 'acbt_ap', 'tcbt_ap');

return

%---------------------
% interpolation tests
%---------------------

% basic interpolation
i1rad = interp1(tafrq, tarad, tcfrq, 'spline', 'extrap');

% call modified version of airs2cris that does interpolation rather
% than deconvolution to the intermediate grid.  choice of sfile does
% not matter here
sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';
[i2rad, i2frq] = airs2crisX(tarad, tafrq, sfile);

% get the channel intersection
[tci, ici] = seq_match(tcfrq, i2frq);
i2rad = i2rad(ici, :);  i2frq = i2frq(ici);
% isclose(i2frq, tcfrq)

% get brightness temps
i1bt = real(rad2bt(tcfrq, i1rad));
i2bt = real(rad2bt(tcfrq, i2rad));

figure(3); clf
set(gcf, 'Units','centimeters', 'Position', [4, 10, 24, 16])

subplot(2,1,1)
y1 = mean(i1bt - tcbt, 2); 
y2 = mean(i2bt - tcbt, 2); 
y3 = mean(acbt - tcbt, 2);
plot(tcfrq, y1, 'b', tcfrq, y2, 'g', tcfrq, y3, 'r')
axis([650, 1100, -6, 6])
title(sprintf('AIRS to unapodized CrIS %s translation residual mean', band));
legend('interpolation', 'AIRS interp/CrIS conv', 'AIRS decon/CrIS conv', ...
       'location', 'north');
ylabel('dTb')
grid on; zoom on

subplot(2,1,2)
z1 = std(i1bt - tcbt, 0, 2); 
z2 = std(i2bt - tcbt, 0, 2); 
z3 = std(acbt - tcbt, 0, 2);
plot(tcfrq, z1, 'b', tcfrq, z2, 'g', tcfrq, z3, 'r')
axis([650, 1100, 0, 2])
title(sprintf('AIRS to unapodized CrIS %s translation residual std dev', band));
legend('interpolation', 'AIRS interp/CrIS conv', 'AIRS decon/CrIS conv', ...
       'location', 'north');
xlabel('wavenumber')
ylabel('dTb')
grid on; zoom on
% export_fig(sprintf('a2cris_interp_%s.pdf', band), '-m2', '-transparent')
% saveas(gcf, sprintf('a2cris_interp_%s.pdf', band), 'png')
