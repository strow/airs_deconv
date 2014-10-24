%
% test crisg4_unapod vs CrIS reference truth
%
% main steps
%   - start with a kcarta radiance, 
%   - convolve to CrIS reference truth
%   - do a spectral-space hamming apodization
%   - apply the function crisg4_unapod, and
%   - compare with the reference truth.

% path for Scott & Larrabee's stuff
addpath /asl/rtp_prod/cris/unapod

% path for inst_params and rad2bt
addpath /home/motteler/cris/bcast/source

% path for mkhamm, seq_match, isclose, etc.
addpath /home/motteler/cris/bcast/motmsc/utils

% sarta CrIS channels
fm1 = 650:0.625:1095;         % LW user grid
fm2 = 1210:1.25:1750;         % MW user grid
fm3 = 2155:2.5:2550;          % SW user grid
fm4 = 647.5:0.625:649.375;    % LW lower guard
fm5 = 1095.625:0.625:1097.5;  % LW upper guard
fm6 = 1205.00:1.25:1208.75;   % MW lower guard
fm7 = 1751.25:1.25:1755;      % MW upper guard
fm8 = 2145.00:2.5:2153.50;    % SW lower guard
fm9 = 2552.5:2.5:2560;        % SW upper guard
fm = [fm1';fm2';fm3';fm4';fm5';fm6';fm7';fm8';fm9'];

% is2f takes sarta to freq order, if2s takes freq to sarta order
[vcris, is2f] = sort(fm);  
[xx, if2s] = sort(is2f);   
% isequal(vcris, fm(is2f))
% isequal(vcris(if2s), fm)

% load kcarta test radiance
kcfile = '/home/motteler/cris/sergio/JUNK2012/convolved_kcarta6.mat';
d1 = load(kcfile);
vkc = d1.w;  % frequency grid
rkc = d1.r;  % radiances
clear d1

% get CrIS reference truth for all 3 bands
band = 'LW';
wlaser = 773.1301;
[inst, user] = inst_params(band, wlaser);
wlaser = 1e7/(inst.df/(2*user.opd/inst.npts));
[inst, user] = inst_params(band, wlaser);
[rLW, vLW] = kc2cris(user, rkc, vkc);

band = 'MW';
[inst, user] = inst_params(band, wlaser);
wlaser = 1e7/(inst.df/(2*user.opd/inst.npts));
[inst, user] = inst_params(band, wlaser);
[rMW, vMW] = kc2cris(user, rkc, vkc);

band = 'SW';
[inst, user] = inst_params(band, wlaser);
wlaser = 1e7/(inst.df/(2*user.opd/inst.npts));
[inst, user] = inst_params(band, wlaser);
[rSW, vSW] = kc2cris(user, rkc, vkc);

% match frequency grids
rALL = [rLW; rMW; rSW];
vALL = [vLW; vMW; vSW];
[ix, jx] = seq_match(vALL, vcris);

rsinc = rALL(ix);
vsinc = vALL(ix);
% isclose(vsinc, vcris)

% apply hamming to the extended CrIS grid
rtmp = mkhamm(length(rALL)) * rALL;
rhamm = rtmp(ix);

figure(1); clf
plot(vsinc, rsinc, vsinc, rhamm)
legend('sinc', 'hamm')
title('true CrIS sinc and hamming')

% translate to sarta indices
rsinc_s = rsinc(if2s);
rhamm_s = rhamm(if2s);
vsinc_s = vsinc(if2s);
% isequal(vsinc_s1, fm)

% call Scott's unapodizer
ichan = (1:length(fm))';
rtest = crisg4_unapod(ichan, rhamm_s);
vtest = [fm1';fm2';fm3'];

% compare
[iy, jy] = seq_match(vsinc, vtest);
vsinc2 = vsinc(iy);
rsinc2 = rsinc(iy);
vtest2 = vtest(jy);
rtest2 = rtest(jy);
% isclose(vsinc2, vtest2)

btsinc2 = rad2bt(vsinc2, rsinc2);
bttest2 = rad2bt(vtest2, rtest2);

figure(2); clf
% plot(vsinc2, rsinc2, vtest2, rtest2)
% legend('sinc', 'test')
% title('true CrIS sinc and test')
plot(vsinc2, rtest2 - rsinc2)
title('crisg4 unapod minus reference truth')
ylabel('brightness temp diff')
xlabel('wavenumber');
grid on; zoom on
saveas(gcf, 'hamming_test', 'fig')

