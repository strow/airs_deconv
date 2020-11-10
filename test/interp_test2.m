%
% test interpolation with and without guard channels
%

addpath /asl/packages/ccast/source
addpath /asl/packages/ccast/motmsc/time
addpath /asl/packages/ccast/motmsc/utils
addpath /asl/packages/airs_decon/source

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
opt1.user_res = 'lowres';  % 'lowres' or 'midres'
opt1.inst_res = 'hires3';  % nominal inst res
wlaser = 773.1301;         % nominal wlaser

% interpolate with guard channels
[instMW, userMW] = inst_params('MW', wlaser, opt1);
[r2MW, v2MW] = finterp(r1MW(:,:), v1MW, userMW.dv);
ix = find(userMW.v1 <= v2MW & v2MW <= userMW.v2);
v2MW = v2MW(ix);
r2MW = r2MW(ix, :); 
r2MW = hamm_app(double(r2MW));
r2MW = reshape(r2MW, length(v2MW), 9, 30, nscan);
if userMW.v1 ~= v2MW(1) | userMW.v2 ~= v2MW(end)
  error('MW grid mismatch')
end

% interpolate without guard channels
iy = userMW.v1 <= v1MW & v1MW <= userMW.v2;
[r2MWa, v2MWa] = finterp(r1MW(iy,:), v1MW(iy), userMW.dv);
ixa = find(userMW.v1 <= v2MWa & v2MWa <= userMW.v2);
v2MWa = v2MWa(ixa);
r2MWa = r2MWa(ixa, :); 
r2MWa = hamm_app(double(r2MWa));
r2MWa = reshape(r2MWa, length(v2MWa), 9, 30, nscan);
if userMW.v1 ~= v2MWa(1) | userMW.v2 ~= v2MWa(end)
  error('MWa grid mismatch')
end

bt  = real(rad2bt(v2MW, r2MW));
bta = real(rad2bt(v2MWa, r2MWa));

btm = mean(bt(:,:),2);
btma = mean(bta(:,:),2);
plot(v2MW, btma - btm)
title('interpolation difference with and without guard channels')
% axis([1200, 1750, -2, 2])
  axis([1200, 1750, -0.1, 0.1])
xlabel('wavenumber (cm-1)')
ylabel('dBT (K)')
grid on
saveas(gcf, 'guard_diff_MW', 'png')

