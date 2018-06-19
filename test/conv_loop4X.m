%
% conv_loop4 - convolutions for AIRS to CrIS direct regression tests
%
% key variables
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
%

% set paths to libs
addpath ../source
addpath /asl/packages/ccast/source

% kcarta dependent set
kcdd = '/asl/s1/motteler/kc7377/cloudy';
nkcd = 7377;

% kcarta independent set
kcdi = '/home/motteler/cris/sergio/JUNK2012';
nkci = 49;

%-------------
% AIRS params
%-------------

% specify SRF tabulations
sdir = '/asl/matlab2012/srftest/';
  f1srf = fullfile(sdir, 'srftables_m130f_withfake_mar08.hdf');
% f2srf = fullfile(sdir, 'srftables_m140f_withfake_mar08.hdf');
% f3srf = fullfile(sdir, 'srftables_m150f_withfake_mar08.hdf');

% get SRF channel frequencies
v1srf = srf_read(f1srf);

% use the L1b subset
% v1srf = v1srf(1:2378);

% sort the SRF channel set
[~, i1srf] = sort(v1srf); 
v1srf = v1srf(i1srf);

% get the JPL L1C channel set
c2645 = load('freq2645.txt');

% match L1C and SRF channel sets
[ix, ~] = seq_match(v1srf, c2645, 0.04);
i1srf = i1srf(ix);  

% convolution matrices for reference truth
dvk = 0.0025; 
[C1, v1col, v1row] = mksconv(f1srf, i1srf, dvk);

% initialize AIRS radiance arrays
va1C = v1row;
na1C = length(va1C);
a1Cri = zeros(na1C, nkci);
a1Crd = zeros(na1C, nkcd);

%-------------
% CrIS params
%-------------

ng = 0;
opt1 = struct;
opt1.inst_res = 'hires3';
opt1.user_res = 'lowres';
wlaser = 773.1307;
[instLW, userLW] = inst_params('LW', wlaser, opt1);
[instMW, userMW] = inst_params('MW', wlaser, opt1);
[instSW, userSW] = inst_params('SW', wlaser, opt1);

optLW = struct;  optLW.ng = ng;
optMW = struct;  optMW.ng = ng;
optSW = struct;  optSW.ng = ng;

% ccast e5 filters
optLW.pL =  650; optLW.pH = 1100; optLW.rL = 15; optLW.rH = 20;
optMW.pL = 1200; optMW.pH = 1760; optMW.rL = 30; optMW.rH = 30;
optSW.pL = 2145; optSW.pH = 2560; optSW.rL = 30; optSW.rH = 30;

% initialize CrIS radiance arrays
ncLW = length(cris_ugrid(userLW, ng));
ncMW = length(cris_ugrid(userMW, ng));
ncSW = length(cris_ugrid(userSW, ng));
cLWrd = zeros(ncLW, nkcd); cLWri = zeros(ncLW, nkci);
cMWrd = zeros(ncMW, nkcd); cMWri = zeros(ncMW, nkci);
cSWrd = zeros(ncSW, nkcd); cSWri = zeros(ncSW, nkci);

%-------------------------------
% loop on kcarta radiance files
%-------------------------------

% loop on dependent set 
tic
for i = 1 : nkcd

  kcmat = fullfile(kcdd, sprintf('kc%04d.mat', i));
  d1 = load(kcmat);
  rkc = d1.rad; vkc = d1.frq; clear d1

  % CrIS convolutions
  [rtmp, vcLW] = kc2cris(userLW, rkc, vkc, optLW);
  rtmp = hamm_app(rtmp);
  cLWrd(:, i) = rtmp;

  [rtmp, vcMW] = kc2cris(userMW, rkc, vkc, optMW);
  rtmp = hamm_app(rtmp);
  cMWrd(:, i) = rtmp;

  [rtmp, vcSW] = kc2cris(userSW, rkc, vkc, optSW);
  rtmp = hamm_app(rtmp);
  cSWrd(:, i) = rtmp;

 % AIRS convolution
  ix = interp1(vkc, 1:length(rkc), v1col, 'nearest');
  a1Crd(:, i) = C1 * rkc(ix);

  if mod(i, 100) == 0, fprintf(1, '.'), end
end
fprintf(1, '\n'); toc

% loop on independent set
tic
for i = 1 : nkci

  kcmat = fullfile(kcdi, sprintf('convolved_kcarta%d.mat', i));
  d1 = load(kcmat);
  rkc = d1.r; vkc = d1.w; clear d1

  % CrIS convolutions
  [rtmp, vcLW] = kc2cris(userLW, rkc, vkc, optLW);
  rtmp = hamm_app(rtmp);
  cLWri(:, i) = rtmp;

  [rtmp, vcMW] = kc2cris(userMW, rkc, vkc, optMW);
  rtmp = hamm_app(rtmp);
  cMWri(:, i) = rtmp;

  [rtmp, vcSW] = kc2cris(userSW, rkc, vkc, optSW);
  rtmp = hamm_app(rtmp);
  cSWri(:, i) = rtmp;

 % AIRS convolution
  ix = interp1(vkc, 1:length(rkc), v1col, 'nearest');
  a1Cri(:, i) = C1 * rkc(ix);

  fprintf(1, '.');
end
fprintf(1, '\n'); toc

% save key variables
save conv_loop4X ...
  nkcd nkci va1C na1C a1Crd a1Cri ...
  vcLW vcMW vcSW ncLW ncMW ncSW ...
  cLWrd cMWrd cSWrd cLWri cMWri cSWri 

