% 
% guard_test -- test kc2cris with added guard channels
% 

addpath /home/motteler/cris/ccast/source
addpath /home/motteler/cris/airs_decon/source

% test params
band = 'SW';
ngc = 4;

% kcarta test data
kcdir = '/home/motteler/cris/sergio/JUNK2012/';
flist =  dir(fullfile(kcdir, 'convolved_kcart*.mat'));

% get CrIS inst and user params
opts.resmode = 'hires2';
wlaser = 773.1301;  % nominal value
[inst, user] = inst_params(band, wlaser, opts);

% loop on kcarta files
rad1 = []; rad2 = [];
for i = 1 : length(flist)
  d1 = load(fullfile(kcdir, flist(i).name));
  vkc = d1.w(:); rkc = d1.r(:);

  % convolve kcarta radiances to CrIS channels
  [rtmp, ftmp] = kc2cris_old(user, rkc, vkc);
  rad1 = [rad1, rtmp];
  frq1 = ftmp;

  % convolve kcarta radiances to CrIS channels
% [rtmp, ftmp] = kc2cris(user, rkc, vkc);
  [rtmp, ftmp] = kc2cris(user, rkc, vkc, ngc);
  rad2 = [rad2, rtmp];
  frq2 = ftmp;

  fprintf(1, '.');
end
fprintf(1, '\n')

isequal(rad1, rad2(ngc+1:end-ngc, :))

plot(frq2, rad2(:, 1), 'r', frq1, rad1(:, 1), 'g')
legend('guard chans', 'user grid')
grid on; zoom on
