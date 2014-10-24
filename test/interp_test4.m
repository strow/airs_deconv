%
% interp_test4 -- check interpolation from sensor to user gird
% 
% main steps
%   1. use kc2inst to take kcarta to CrIS sensor grid
%   2a use kc2inst to take kcarta to CrIS user grid
%   2b use kc2cris to take kcarta to CrIS user grid
%   3. use new finterp to take sensor grid to user grid
%   4. use old finterp to take sensor grid to user grid
%   5. compare 2 and 3
%   6. compare 3 and 4
%

% use my libs for tests
addpath /home/motteler/cris/ccast/source
addpath /home/motteler/cris/airs_decon/source

% kcarta test data
kcdir = '/home/motteler/cris/sergio/JUNK2012/';
flist =  dir(fullfile(kcdir, 'convolved_kcart*.mat'));

% get CrIS inst and user params
band = 'SW';
wlaser = 773.13;
opt1 = struct;
opt1.resmode = 'hires2';
[inst1, user1] = inst_params(band, wlaser, opt1);

% set wlaser so inst grid == user grid
wlaser = 1e7/(inst1.df/(2*user1.opd/inst1.npts));
[inst2, user2] = inst_params(band, wlaser, opt1);

% loop on kcarta files
rad1 = []; rad2 = [];
for i = 1 : length(flist)
  d1 = load(fullfile(kcdir, flist(i).name));
  vkc = d1.w(:); rkc = d1.r(:);

  % 1. convolve kcarta radiances to the CrIS sensor grid
  [rtmp1, ftmp1] = kc2inst(inst1, user1, rkc, vkc);
  rad1 = [rad1, rtmp1];
  frq1 = ftmp1(:);

  % 2. convolve kcarta radiances to the CrIS user grid
% [rtmp2, ftmp2] = kc2inst(inst2, user2, rkc, vkc);
  [rtmp2, ftmp2] = kc2cris(user2, rkc, vkc);
  rad2 = [rad2, rtmp2];
  frq2 = ftmp2(:);

end

% 3. take sensor to user grid with new finterp
tol = 1e-6;
opt2 = struct; opt2.tol = tol; opt2.info = 1
[rad3, frq3] = finterp(rad1, frq1, user2.dv, opt2);

% 4. take sensor to user grid with old finterp
[rad4, frq4] = finterp_old(rad1, frq1, user2.dv, opt2);
frq4 = frq4(:);

% compare new finterp and user grid
[ix, jx] = seq_match(frq2, frq3);
r2 = rad2(ix, :); f2 = frq2(ix, 1);
r3 = rad3(jx, :); f3 = frq3(jx, 1);
bt2 = real(rad2bt(f2, r2));
bt3 = real(rad2bt(f3, r3));
isclose(f2, f3) 

figure(1); clf
plot(f2, bt3 - bt2)
axis([user1.v1, user2.v2, -0.2, 0.2])
title('finterp minus user grid')
xlabel('wavenumber')
ylabel('dBT')
grid on; zoom on

% compare old and new finterp
[ix, jx] = seq_match(frq3, frq4);
r3 = rad3(ix, :); f3 = frq3(ix, 1);
r4 = rad4(jx, :); f4 = frq4(jx, 1);
bt3 = real(rad2bt(f3, r3));
bt4 = real(rad2bt(f4, r4));
isclose(f3, f4) 

figure(2); clf
plot(f3, bt3 - bt4)
% axis([user1.v1, user2.v2, -3e-4, 3e-4])
axis([user1.v1, user2.v2, -0.02, 0.02])
title('new minus old finterp')
xlabel('wavenumber')
ylabel('dBT')
grid on; zoom on

