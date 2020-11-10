%
% simple kc2cris demo
% 

addpath /asl/packages/ccast/source
addpath /asl/packages/airs_decon/source

% sample kcarta radiances
f = '/umbc/xfs2/strow/asl/s1/sergio/home/MATLABCODE/QUICKTASKS_TELECON/Add_2_FFTs/CRIS_CONV_for_Xianglei/testkcrads.mat';

% local radiance names
rkc = rad;
vkc = w;

% number of guard channels
ng = 2;  

% options for inst_params
opt1 = struct;
opt1.user_res = 'lowres';
wlaser = 773.1307;  % purely nominal value
[instLW, userLW] = inst_params('LW', wlaser, opt1);
[instMW, userMW] = inst_params('MW', wlaser, opt1);
[instSW, userSW] = inst_params('SW', wlaser, opt1);

% options for kc2cris
optLW = struct;  optLW.ng = ng;
optMW = struct;  optMW.ng = ng;
optSW = struct;  optSW.ng = ng;

% add ccast e5 filters
optLW.pL =  650; optLW.pH = 1100; optLW.rL = 15; optLW.rH = 20;
optMW.pL = 1200; optMW.pH = 1760; optMW.rL = 30; optMW.rH = 30;
optSW.pL = 2145; optSW.pH = 2560; optSW.rL = 30; optSW.rH = 30;

% CrIS convolutions
[radLW, frqLW] = kc2cris(userLW, rkc, vkc, optLW);
[radMW, frqMW] = kc2cris(userMW, rkc, vkc, optMW);
[radSW, frqSW] = kc2cris(userSW, rkc, vkc, optSW);

