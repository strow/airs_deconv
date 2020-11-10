%
% compare ccast low res SDRs with high res interpolated to low
%

addpath /asl/packages/ccast/source
addpath /asl/packages/ccast/motmsc/time
addpath /asl/packages/ccast/motmsc/utils
addpath /asl/packages/airs_decon/source

p1 = '/asl/cris/ccast/sdr45_npp_LR/2019/196';
p2 = '/asl/cris/ccast/sdr45_npp_HR/2019/196';
gran = 'CrIS_SDR_npp_s45_d20190715_t0218080_g024_v20a.mat';
nscan = 45;

f1 = fullfile(p1, gran);
f2 = fullfile(p2, gran);

d1 = load(f1);
d2 = load(f2);

t1 = d1.geo.FORTime;
[datestr(iet2dnum(t1(1))), '  ', datestr(iet2dnum(t1(end)))]
[d1.geo.Latitude(1), d1.geo.Longitude(1)]

v1MW = d1.vMW; r1MW = d1.rMW;  % low res SDR data
v2MW = d2.vMW; r2MW = d2.rMW;  % high res SDR data

% specify the target res
opt1 = struct;
opt1.user_res = 'lowres';  % 'lowres' or 'midres'
opt1.inst_res = 'hires3';  % nominal inst res
wlaser = 773.1301;         % nominal wlaser

% interpolate high res to the MW low user grid
[instMW, userMW] = inst_params('MW', wlaser, opt1);
[riMW, viMW] = finterp(r2MW(:,:), v2MW, userMW.dv);
ix = find(userMW.v1 <= viMW & viMW <= userMW.v2);
viMW = viMW(ix);
riMW = riMW(ix, :); 
riMW = hamm_app(double(riMW));
riMW = reshape(riMW, length(viMW), 9, 30, nscan);
if userMW.v1 ~= viMW(1) | userMW.v2 ~= viMW(end)
  error('MW grid mismatch')
end

iy = find(userMW.v1 <= v1MW & v1MW <= userMW.v2);
v1MW = v1MW(iy);
r1MW = r1MW(iy, :); 
r1MW = hamm_app(double(r1MW));
r1MW = reshape(r1MW, length(v1MW), 9, 30, nscan);
if userMW.v1 ~= v1MW(1) | userMW.v2 ~= v1MW(end)
  error('MW gr1d mismatch')
end

% compare interpolated and regular low res
b1 = real(rad2bt(v1MW, r1MW));
bi = real(rad2bt(viMW, riMW));
bim = mean(bi(:, :, :), 3);
b1m = mean(b1(:, :, :), 3);

figure(1)
set(gcf, 'DefaultAxesColorOrder', fovcolors);
plot(v1MW, bim - b1m)
title('interpolated minus production mean')
% axis([1200, 1750, -0.5, 0.5])
  axis([1200, 1750, -0.1, 0.1])
legend(fovnames, 'location', 'south')
xlabel('wavenumber (cm-1)')
ylabel('dBT (K)')
grid on
% saveas(gcf, 'interp_resid_MW', 'png')

v1SW = d1.vSW; r1SW = d1.rSW;  % low res SDR data
v2SW = d2.vSW; r2SW = d2.rSW;  % high res SDR data

% specify the target res
opt1 = struct;
opt1.user_res = 'lowres';  % 'lowres' or 'midres'
opt1.inst_res = 'hires3';  % nominal inst res
wlaser = 773.1301;         % nominal wlaser

% interpolate high res to the SW low user grid
[instSW, userSW] = inst_params('SW', wlaser, opt1);
[riSW, viSW] = finterp(r2SW(:,:), v2SW, userSW.dv);
ix = find(userSW.v1 <= viSW & viSW <= userSW.v2);
viSW = viSW(ix);
riSW = riSW(ix, :); 
riSW = hamm_app(double(riSW));
riSW = reshape(riSW, length(viSW), 9, 30, nscan);
if userSW.v1 ~= viSW(1) | userSW.v2 ~= viSW(end)
  error('SW grid mismatch')
end

iy = find(userSW.v1 <= v1SW & v1SW <= userSW.v2);
v1SW = v1SW(iy);
r1SW = r1SW(iy, :); 
r1SW = hamm_app(double(r1SW));
r1SW = reshape(r1SW, length(v1SW), 9, 30, nscan);
if userSW.v1 ~= v1SW(1) | userSW.v2 ~= v1SW(end)
  error('SW gr1d mismatch')
end

% compare interpolated and regular low res
b1 = real(rad2bt(v1SW, r1SW));
bi = real(rad2bt(viSW, riSW));
bim = mean(bi(:, :, :), 3);
b1m = mean(b1(:, :, :), 3);

figure(2)
set(gcf, 'DefaultAxesColorOrder', fovcolors);
plot(v1SW, bim - b1m)
title('interpolated minus production mean')
% axis([2150, 2550, -0.25, 0.25])
  axis([2150, 2550, -0.1, 0.1])
legend(fovnames, 'location', 'southeast')
xlabel('wavenumber (cm-1)')
ylabel('dBT (K)')
grid on
% saveas(gcf, 'interp_resid_SW', 'png')


