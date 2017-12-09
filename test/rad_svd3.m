%
% rad_svd2 -- compare row and column reconstruction scores
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

[ud0, sd0, vd0] = svd(rd, 0);
[ud1, sd1, vd1] = svd(rd', 0);

% rd is nchan x nobs
n0 = 260;
rd2 = ud0(:, 1:n0) * ud0(:, 1:n0)' * rd;
bd2 = real(rad2bt(v, rd2));
rms(bd(:) - bd2(:))

% rd' is nobs x nchan
n1 = 260;
rd2 = ud1(:, 1:n1) * ud1(:, 1:n1)' * rd';
rd2 = rd2';
bd2 = real(rad2bt(v, rd2));
rms(bd(:) - bd2(:))

