%
% demo interpolation of CrIS high res to mid or low res
%

addpath /asl/packages/ccast/source

% read a ccast high res granule
cdir = '/asl/cris/ccast/sdr45_j01_HR/2018/091';
gran = 'CrIS_SDR_j01_s45_d20180401_t2148080_g219_v20a.mat';
nscan = 45;

f1 = fullfile(cdir, gran);
d1 = load(f1);
v1MW = d1.vMW;
r1MW = d1.rMW;
v1SW = d1.vSW;
r1SW = d1.rSW;

% specify the target res
opt1 = struct;
opt1.user_res = 'midres';  % 'lowres' or 'midres'
opt1.inst_res = 'hires3';  % nominal inst res
wlaser = 773.1301;         % nominal wlaser

% interpolate the MW user grid
[instMW, userMW] = inst_params('MW', wlaser, opt1);
[r2MW, v2MW] = finterp(r1MW(:,:), v1MW, userMW.dv);
ix = find(userMW.v1 <= v2MW & v2MW <= userMW.v2);
v2MW = v2MW(ix);
r2MW = r2MW(ix, :); 
r2MW = reshape(r2MW, length(v2MW), 9, 30, nscan);
if userMW.v1 ~= v2MW(1) | userMW.v2 ~= v2MW(end)
  error('MW grid mismatch')
end

% interpolate the SW user grid
[instSW, userSW] = inst_params('SW', wlaser, opt1);
[r2SW, v2SW] = finterp(r1SW(:,:), v1SW, userSW.dv);
ix = find(userSW.v1 <= v2SW & v2SW <= userSW.v2);
v2SW = v2SW(ix);
r2SW = r2SW(ix, :); 
r2SW = reshape(r2SW, length(v2SW), 9, 30, nscan);
if userSW.v1 ~= v2SW(1) | userSW.v2 ~= v2SW(end)
  error('SW grid mismatch')
end

% spot check of the interpolation
bMW1 = real(rad2bt(v1MW, r1MW(:,5,15,23)));
bSW1 = real(rad2bt(v1SW, r1SW(:,5,15,23)));
bMW2 = real(rad2bt(v2MW, r2MW(:,5,15,23)));
bSW2 = real(rad2bt(v2SW, r2SW(:,5,15,23)));

figure(1); clf
subplot(2,1,1)
plot(v1MW, bMW1, v2MW, bMW2)
axis([1200, 1750, 200, 300]);
title('MW interpolation test')
legend('hires', opt1.user_res)
grid on; zoom on

subplot(2,1,2)
plot(v1SW, bSW1, v2SW, bSW2)
axis([2150, 2550, 200, 300])
title('SW interpolation test')
legend('hires', opt1.user_res)
grid on; zoom on

