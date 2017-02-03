%
% conv_loop - AIRS and CrIS convolution of kcarta radiances
% 

addpath /asl/packages/ccast/source
addpath ../h4tools
addpath ../source

% turn off HDF 4 update warnings
warning('off', 'MATLAB:imagesci:hdf:removalWarningHDFSD')

% path to kcarta data
kcdir = '/asl/s1/motteler/kc7377/cloudy';
nkcrad = 7377;

%-------------
% CrIS params
%-------------

ng = 0;
opt1 = struct;
opt1.inst_res = 'hires3';
opt1.user_res = 'hires';
wlaser = 773.1307;
[instLW, userLW] = inst_params('LW', wlaser, opt1);
[instMW, userMW] = inst_params('MW', wlaser, opt1);
[instSW, userSW] = inst_params('SW', wlaser, opt1);

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

sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';
cfreq = load('freq2645.txt');
% cfreq = trim_chans(cfreq);
dvk = 0.0025; 
[sconv, sfreq, ofreq] = mksconv2(sfile, cfreq, dvk);

% tabulated AIRS radiances
afrq = ofreq(:);
arad = zeros(length(afrq), nkcrad);

%-------------------------------
% loop on kcarta radiance files
%-------------------------------

for i = 1 : nkcrad

  % CrIS convolutions
  kcmat = fullfile(kcdir, sprintf('kc%04d.mat', i));
  d1 = load(kcmat);
  rkc = d1.rad; vkc = d1.frq; clear d1

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

save cris_cloudy radLW radMW radSW frqLW frqMW frqSW ...
            userLW userMW userSW optLW optMW optSW

save airs_cloudy arad afrq sfile

