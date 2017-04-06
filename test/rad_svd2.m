%
% rad_svd2 -- SVD and associated stats on sample radiance sets
%
% key variables from conv_loop4
%   nkcd  - dependent set size
%   nkci  - independent set size
%   va1C  - AIRS channel frequency
%   na1C  - number of AIRS channels
%   a1Crd - AIRS dep set radiance
%   a1Cri - AIRS ind set radiance
%   vcLW,  vcMW,  vcSW  - CrIS channel frequency
%   ncLW,  ncMW,  ncSW  - number of CrIS channels
%   cLWrd, cMWrd, cSWrd - CrIS dep set radiance
%   cLWri, cMWri, cSWri - CrIS ind set radiance

addpath ../source
addpath /asl/packages/ccast/source

% get radiance data
load('conv_loop4')

v = va1C;
rd = a1Crd;
ri = a1Cri;  
bd = real(rad2bt(v, rd));
bi = real(rad2bt(v, ri));

d1 = load('/home/motteler/shome/sergio/tigr_airs.mat');
rt = squeeze(d1.dKcall(:,:, 1))';
vt = d1.fKc;
clear d1
bt = real(rad2bt(vt, rt));

[ud, sd, ~] = svd(rd, 0);
[ui, si, ~] = svd(ri, 0);
[ut, st, ~] = svd(rt, 0);

sd = diag(sd);
si = diag(si);
st = diag(st);

nd = length(sd);
ni = length(si);
nt = length(st);

figure(1); clf
loglog(1:ni, si(1:ni), 1:nd, sd(1:nd), 1:nt, st(1:nt), 'linewidth', 2)
title('singular values for 3 radiance sets')
legend('49 fitting', '7377 cloudy', '1761 TIGR')
xlabel('singular value index')
ylabel('singular values')
grid on; zoom on

% ind set residual
ni = 48;
ri2 = ui(:, 1:ni) * ui(:, 1:ni)' * ri;
bi2 = real(rad2bt(v, ri2));
rms(bi(:) - bi2(:))

% dep set residual
nd = 260;
rd2 = ud(:, 1:nd) * ud(:, 1:nd)' * rd;
bd2 = real(rad2bt(v, rd2));
rms(bd(:) - bd2(:))

% TIGR residual
nt = 55;
rt2 = ut(:, 1:nt) * ut(:, 1:nt)' * rt;
bt2 = real(rad2bt(vt, rt2));
rms(bt(:) - bt2(:))

% ind with dep basis
nx = 2400;
ri3 = ud(:, 1:nx) * ud(:, 1:nx)' * ri;
bi3 = real(rad2bt(v, ri3));
rms(bi(:) - bi3(:))

% dep with ind basis
ny = 49;
rd3 = ui(:, 1:ny) * ui(:, 1:ny)' * rd;
bd3 = real(rad2bt(v, rd3));
rms(bd(:) - bd3(:))

