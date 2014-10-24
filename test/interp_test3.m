%
% interp_test3 -- compare old and new kc2cris with finterp
%
% rad1 - old kc2cris (now kc2inst)
% rad2 - old finterp
% rad3 - new kc2cris
% rad4 - new finterp

% use my libs for tests
addpath /home/motteler/cris/ccast/source
addpath /home/motteler/cris/airs_decon/source

% kcarta test data
kcdir = '/home/motteler/cris/sergio/JUNK2012/';
flist =  dir(fullfile(kcdir, 'convolved_kcart*.mat'));

% get CrIS inst and user params
band = 'LW';              % CrIS band
wlaser = 773.1301;        % nominal wlaser value
opt1 = struct;
opt1.resmode = 'hires2';  % resolution mode
[inst, user] = inst_params(band, wlaser, opt1);

% set wlaser so inst grid == user grid
wlaser = 1e7/(inst.df/(2*user.opd/inst.npts));
[inst, user] = inst_params(band, wlaser, opt1);

% loop on kcarta files
rad1 = []; rad2 = []; rad_kc = [];
for i = 1 : 4
  d1 = load(fullfile(kcdir, flist(i).name));
  vkc = d1.w(:); rkc = d1.r(:);

  % trim kcarta radiances to cris sensor grid
  ix = find(inst.freq(1) <= vkc & vkc <= inst.freq(end));
  vkc =  vkc(ix); rkc = rkc(ix); 

  % save kcarta radiances in column order
  rad_kc = [rad_kc, rkc];

  % convolve kcarta radiances with kc2inst
  [rtmp1, ftmp1] = kc2inst(inst, user, rkc, vkc);
  rad1 = [rad1, rtmp1];
  frq1 = ftmp1(:);

% fprintf(1, '.');
end
% fprintf(1, '\n')

% convolve kcarta radiances with kc2cris
[rad3, frq3] = kc2cris(user, rad_kc, vkc);

% finterp setup
opt2 = struct; opt2.info = 1; opt2.tol = 1e-6;
rad_kc = bandpass(vkc, rad_kc, user.v1, user.v2, user.vr);

% convolve kcarta radiances with old finterp
[rad2, frq2] = finterp_old(rad_kc, vkc, user.dv, opt2);
frq2 = frq2(:);

% convolve kcarta radiances with the new finterp
[rad4, frq4] = finterp(rad_kc, vkc, user.dv, opt2);

% compare old and new versions of finterp
[ix, jx] = seq_match(frq2, frq4);
r2 = rad2(ix, :); f2 = frq2(ix, 1);
r4 = rad4(jx, :); f4 = frq4(jx, 1);
bt2 = real(rad2bt(f2, r2));
bt4 = real(rad2bt(f4, r4));
isclose(f2, f4) 

figure(1); clf
plot(f2, bt4 - bt2)
title('new minus old finterp')
xlabel('wavenumber')
ylabel('dBT')
grid on; zoom on

% compare old and new versions of kc2cris
[ix, jx] = seq_match(frq1, frq3);
r1 = rad1(ix, :); f1 = frq1(ix, 1);
r3 = rad3(jx, :); f3 = frq3(jx, 1);
bt1 = real(rad2bt(f1, r1));
bt3 = real(rad2bt(f3, r3));
isclose(f1, f3) 

figure(2); clf
plot(f1, bt3 - bt1)
ax(1) = user.v1; ax(2) = user.v2; 
ax(3) = -0.15; ax(4) = 0.15;  axis(ax)
title('new minus old kc2cris')
xlabel('wavenumber')
ylabel('dBT')
grid on; zoom on

% compare new kc2cris with finterp
[ix, jx] = seq_match(frq3, frq4);
r3 = rad3(ix, :); f3 = frq3(ix, :);
r4 = rad4(jx, :); f4 = frq4(jx, :);
bt3 = real(rad2bt(f3, r3));
bt4 = real(rad2bt(f4, r4));
isclose(f3, f4)

figure(3); clf
plot(f3, bt3 - bt4)
title('new kc2cris minus finterp')
xlabel('wavenumber')
ylabel('dBT')
grid on; zoom on

% compare old kc2cris with finterp
[ix, jx] = seq_match(frq1, frq4);
r1 = rad1(ix, :); f1 = frq1(ix, :);
r4 = rad4(jx, :); f4 = frq4(jx, :);
bt1 = real(rad2bt(f1, r1));
bt4 = real(rad2bt(f4, r4));
isclose(f1, f4)

figure(4); clf
plot(f1, bt1 - bt4)
ax(1) = user.v1; ax(2) = user.v2; 
ax(3) = -0.15; ax(4) = 0.15;  axis(ax)
title('old kc2cris minus finterp')
xlabel('wavenumber')
ylabel('dBT')
grid on; zoom on

