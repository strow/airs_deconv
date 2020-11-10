%
% conv_loop1 - AIRS and CrIS convolution of kcarta radiances
% 
% edit as needed
%

addpath /asl/packages/ccast/source
addpath ../source

% path to kcarta data
% kcdir = '/asl/s1/motteler/kc7377/cloudy';
% nkcrad = 7377;
  kcdir = '/home/motteler/cris/sergio/JUNK2012';
  nkcrad = 49;

%-------------
% CrIS params
%-------------

ng = 0;
opt1 = struct;
opt1.user_res = 'hires';
wlaser = 773.1307;
[~, userLW] = inst_params('LW', wlaser, opt1);
[~, userMW] = inst_params('MW', wlaser, opt1);
[~, userSW] = inst_params('SW', wlaser, opt1);

optLW = struct;  optLW.ng = ng;
optMW = struct;  optMW.ng = ng;
optSW = struct;  optSW.ng = ng;

% ccast e5 filters
% probably we could just pass on the inst struct values here...
optLW.pL =  650; optLW.pH = 1100; optLW.rL = 15; optLW.rH = 20;
optMW.pL = 1200; optMW.pH = 1760; optMW.rL = 30; optMW.rH = 30;
optSW.pL = 2145; optSW.pH = 2560; optSW.rL = 30; optSW.rH = 30;

% tabulated CrIS radiances
nLW = length(cris_ugrid(userLW, ng));
nMW = length(cris_ugrid(userMW, ng));
nSW = length(cris_ugrid(userSW, ng));
radLW = zeros(nLW, nkcrad);
radMW = zeros(nMW, nkcrad);
radSW = zeros(nSW, nkcrad);

%-------------
% AIRS params
%-------------

% specify an AIRS SRF file
% sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';
  sfile = '/home/sergio/MATLABCODE/airs_l1c_srf_tables_lls_20181205.hdf';

% get nominal AIRS channel centers from a recent L1c file.  These
% are matched against SRF centers from the SRF file, and the latter
% values are used for subsequent processing.
w1 = load('/home/motteler/shome/airs_tiling/airs_l1c_wnum.mat');
cfreq = w1.wnum;

% call mksconv2 to get the SRF convolution matrix
dvk = 0.0025; 
[sconv, sfreq, ofreq] = mksconv2(sfile, cfreq, dvk);

% tabulated AIRS radiances
afrq = ofreq(:);
arad = zeros(length(afrq), nkcrad);

%-------------------------------
% loop on kcarta radiance files
%-------------------------------

for i = 1 : nkcrad

% kcarta radiances for the big cloudy set
% kcmat = fullfile(kcdir, sprintf('kc%04d.mat', i));
% d1 = load(kcmat);
% rkc = d1.rad; vkc = d1.frq; clear d1

% kcarta radiances for the 49 fitting profiles
  kcmat = fullfile(kcdir, sprintf('convolved_kcarta%d.mat', i));
  d1 = load(kcmat);
  rkc = d1.r; vkc = d1.w; clear d1

  % CrIS convolutions
  [rtmp, frqLW] = kc2cris(userLW, rkc, vkc, optLW);
  radLW(:, i) = rtmp;

  [rtmp, frqMW] = kc2cris(userMW, rkc, vkc, optMW);
  radMW(:, i) = rtmp;

  [rtmp, frqSW] = kc2cris(userSW, rkc, vkc, optSW);
  radSW(:, i) = rtmp;

  % AIRS convolution
  ix = interp1(vkc, 1:length(rkc), sfreq, 'nearest');
  rtmp = sconv * rkc(ix);
  arad(:, i) = rtmp;
  
  if mod(i, 100) == 0, fprintf(1, '.'), end
end
fprintf(1, '\n')

%-----------------------
% save the convolutions
%-----------------------

save airs_srf_49 arad afrq sfile

save crisHR_49 radLW radMW radSW frqLW frqMW frqSW ...
     userLW userMW userSW optLW optMW optSW opt1

