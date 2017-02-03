%
% translate CrIS high res to AIRS L1b
%

addpath /asl/matlib/h4tools
addpath /asl/packages/ccast/source
addpath /asl/packages/airs_decon/source

% specify an AIRS SRF tabulation
sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';

% AIRS 1b channel frequencies
d2 = load('freqL1b');
cfrq = sort(d2.freqL1b);

% CrIS params
wlaser = 773.13;
opt1 = struct;
opt1.resmode = 'hires2';
[instLW, userLW] = inst_params('LW', wlaser, opt1);
[instMW, userMW] = inst_params('MW', wlaser, opt1);
[instSW, userSW] = inst_params('SW', wlaser, opt1);
vLW = cris_ugrid(userLW, 2)';
vMW = cris_ugrid(userMW, 2)';
vSW = cris_ugrid(userSW, 2)';

iLW =    1 :  717;
iMW =  718 : 1586;
iSW = 1587 : 2223;

k = 240;  % max calls in a chunk

flist = dir(fullfile('tsno', 'tsno.*.nc'));

% profile clear
% profile on

for i = 1 : length(flist)

  tfile = fullfile('tsno', flist(i).name);
  rad =  h5read(tfile, '/L1bCrIS/rad');
  [m, nobs] = size(rad);

  for j = 1 : k : nobs
    ix = j : min(j + k - 1, nobs);

    rLW = rad(iLW, ix);
    rMW = rad(iMW, ix);
    rSW = rad(iSW, ix);

    [rtmp, afrq] = ...
         cris2airs(rLW, rMW, rSW, vLW, vMW, vSW, sfile, cfrq, opt1);

    if j == 1
      arad = zeros(length(afrq), nobs);
    end

   arad(:, ix) = rtmp;

   [sp, sn, se] = fileparts(flist(i).name);
   mfile = fullfile('c2airs', [sn, '.mat']);
%  save(mfile, 'arad', 'afrq', '-v7.3')

 end
 fprintf(1, '.')
end
fprintf(1, '\n')

% profile report

