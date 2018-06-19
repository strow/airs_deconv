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
d4 = load('ac_ap_fit49');  % apodized airs cris

% unapodized radiance
tcrad = [d1.radLW; d1.radMW; d1.radSW];
tcfrq = [d1.frqLW; d1.frqMW; d1.frqSW];
tarad = d2.arad;  tafrq = d2.afrq;  % true airs
acrad = d3.crad;  acfrq = d3.cfrq;  % airs cris
adrad = d3.brad;  adfrq = d3.bfrq;  % airs decon

% apodized radiance
tcrad_ap = [hamm_app(d1.radLW); hamm_app(d1.radMW); hamm_app(d1.radSW)];
acrad_ap = d4.crad; % apodized airs cris
clear d1 d2 d3 d4

% get the channel intersection
[tci, aci] = seq_match(tcfrq, acfrq);
tcrad = tcrad(tci, :);  tcfrq = tcfrq(tci);
acrad = acrad(aci, :);  acfrq = acfrq(aci);
tcrad_ap = tcrad_ap(tci, :);
acrad_ap = acrad_ap(aci, :);
[nchan, nobs] = size(tcrad);

% get brightness temps
tcbt = real(rad2bt(tcfrq, tcrad));
tabt = real(rad2bt(tafrq, tarad));
acbt = real(rad2bt(acfrq, acrad));
adbt = real(rad2bt(adfrq, adrad));
tcbt_ap = real(rad2bt(tcfrq, tcrad_ap));
acbt_ap = real(rad2bt(tcfrq, acrad_ap));

% get band for plots
band = upper(input('band > ', 's'));

%-----------------------------------
% apodized and unapodized residuals 
%-----------------------------------
figure(1); clf
subplot(2,1,1)
plot(tcfrq, mean(acbt-tcbt,2), tcfrq, mean(acbt_ap-tcbt_ap,2))
switch band
 case 'LW', axis([650, 1095, -1.0, 1.0]);   loc = 'north';
 case 'MW', axis([1210, 1605, -0.3, 0.3]);  loc = 'north';
 case 'SW', axis([2180, 2550, -2, 2]);      loc = 'northeast';
end
title(sprintf('AIRS CrIS minus true CrIS %s mean', band));
legend('unapodized', 'Hamming ap.', 'location', loc)
ylabel('\Delta BT (K)')
grid on; zoom on

% AIRS CrIS minus true CrIS std
subplot(2,1,2)
plot(tcfrq, std(acbt-tcbt,0,2), tcfrq, std(acbt_ap-tcbt_ap,0,2))
switch band
  case 'LW', axis([650, 1095, 0, 0.5]);   loc = 'north';
  case 'MW', axis([1210, 1605, 0, 0.12]); loc = 'north';
  case 'SW', axis([2180, 2550, 0, 1.2]);  loc = 'northeast';
end
title(sprintf('AIRS CrIS minus true CrIS %s std dev', band));
  legend('unapodized', 'Hamming ap.', 'location', loc)
xlabel('wavenumber (cm^{-1})')
ylabel('\Delta BT (K)')
grid on; zoom on
% saveas(gcf, sprintf('a2cris_diff_%s', band), 'fig')

%----------------------------------
% 3-band apodized residual summary
%----------------------------------
figure(2); clf
subplot(3,1,1)
plot(tcfrq, mean(acbt_ap-tcbt_ap,2))
axis([650, 1095, -0.2, 0.2]);
title('apodized AIRS CrIS minus true CrIS mean')
ylabel('\Delta BT (K)')
grid on; zoom on

subplot(3,1,2)
plot(tcfrq, mean(acbt_ap-tcbt_ap,2))
axis([1210, 1605, -0.2, 0.2]); 
ylabel('\Delta BT (K)')
grid on; zoom on

subplot(3,1,3)
plot(tcfrq, mean(acbt_ap-tcbt_ap,2))
axis([2180, 2550, -0.2, 0.2]);
xlabel('wavenumber (cm^{-1})')
ylabel('\Delta BT (K)')
grid on; zoom on

grid on; zoom on
% saveas(gcf, 'a2cris_diff_all', 'fig')

%---------------------------------
% all data for a selected spectra
%---------------------------------
figure(3); clf; j = 1;
subplot(2,1,1)
plot(tafrq, tabt(:,j), adfrq, adbt(:,j), tcfrq, tcbt(:,j), acfrq, acbt(:,j))
switch band
  case 'LW', axis([650, 1095, 200, 310]); loc = 'south';
  case 'MW', axis([1210, 1605, 210, 310]); loc = 'northeast';
  case 'SW', axis([2180, 2550, 210, 310]); loc = 'southeast';
end
title(sprintf('AIRS and unapodized CrIS %s profile %d', band, j));
legend('true AIRS', 'AIRS decon', 'true CrIS', 'AIRS CrIS', 'location', loc)
ylabel('BT (K)')
grid on; zoom on

subplot(2,1,2)
plot(tafrq, tabt(:,j), adfrq, adbt(:,j), tcfrq, tcbt(:,j), acfrq, acbt(:,j))
switch band
  case 'LW', axis([660, 680, 200, 260]); loc = 'northeast';
  case 'MW', axis([1320, 1350, 210, 290]); loc = 'southwest';
  case 'SW', axis([2320, 2360, 210, 260]); loc = 'south';
end
title(sprintf('AIRS and CrIS %s profile %d, detail', band, j));
legend('true AIRS', 'AIRS decon', 'true CrIS', 'AIRS CrIS', 'location', loc)
xlabel('wavenumber (cm^{-1})'); 
ylabel('BT (K)')
grid on; zoom on
% saveas(gcf, sprintf('a2cris_spec_%s', band), 'fig')

%---------------------
% interpolation tests
%---------------------
if ~strcmp(band, 'LW'), return, end

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

% apodize the interpolations
i1rad = hamm_app(i1rad);
i2rad = hamm_app(i2rad);

% get brightness temps
i1bt = real(rad2bt(tcfrq, i1rad));
i2bt = real(rad2bt(tcfrq, i2rad));

figure(4); clf
subplot(2,1,1)
y1 = mean(i1bt - tcbt_ap, 2); 
y2 = mean(i2bt - tcbt_ap, 2); 
y3 = mean(acbt_ap - tcbt_ap, 2);
plot(tcfrq, y1, 'b', tcfrq, y2, 'g', tcfrq, y3, 'r')
axis([650, 1100, -2, 2])
title(sprintf('AIRS to apodized CrIS %s residual mean', band));
legend('spline interpolation', 'AIRS interp/CrIS conv', ...
       'AIRS decon/CrIS conv', 'location', 'north');
ylabel('\Delta BT (K)')
grid on; zoom on

subplot(2,1,2)
z1 = std(i1bt - tcbt_ap, 0, 2); 
z2 = std(i2bt - tcbt_ap, 0, 2); 
z3 = std(acbt_ap - tcbt_ap, 0, 2);
plot(tcfrq, z1, 'b', tcfrq, z2, 'g', tcfrq, z3, 'r')
axis([650, 1100, 0, 1])
title(sprintf('AIRS to apodized CrIS %s residual std dev', band));
legend('spline interpolation', 'AIRS interp/CrIS conv', ...
       'AIRS decon/CrIS conv', 'location', 'north');
xlabel('wavenumber (cm^{-1})')
ylabel('\Delta BT (K)')
grid on; zoom on
  saveas(gcf, sprintf('a2cris_interp_%s', band), 'fig')

