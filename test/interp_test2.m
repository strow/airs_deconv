
% demo using finterp to take AIRS 1c to CrIS radiances
%
% the key idea is that finterp does a linear interpolation 
% to the dv1 spacing and then a fourier interpolation to dv2.

% path to bcast finterp and inst_params
addpath ../source

% get Larrabee's sample data
d1 = load('airs/motteler_filled.mat');  % 1c demo with grid
db = load('airs/sno_cris_airs_1b.mat'); % 1b data with grids
dc = load('airs/sno_cris_airs_1c.mat'); % 1c data, no grids
fcris = db.fcris; 
% fairs = d1.f;   % for airs 1b grid
fairs = d1.ff;    % for airs 1c grid

% choose AIRS and CrIS samples
bta = dc.bta_mid_nh;
btc = dc.btc_mid_nh;
% bta = db.bta_mid_nh;
% btc = db.btc_mid_nh;

% get CrIS user params, specify band 
[inst, user] = inst_params('LW', 773.13);

% sort the AIRS channel set
[fs, ix] = sort(fairs);

% permute AIRS data accordingly
bta = bta(ix);

% switch to radiances for the convolution
rada = bt2rad(fs, bta);

% the finterp default is to set dv1 = fs(2) - fs(1), but this can be
% changed with the opt parameter
opt = struct;
opt.dv1 = 0.1;  

% do the interpolations
[rad2, frq2] = finterp(rada, fs, user.dv, opt);

% match user grid with the finterp output
ugrid = user.v1 : user.dv : user.v2;
i2 = interp1(frq2, 1:length(frq2), ugrid, 'nearest');

% take AIRS radiances back to brightness temps
bta = rad2bt(ugrid, rad2(i2));

% match user grid with Larrabee's cris freq list
i3 = interp1(fcris, 1:length(fcris), ugrid, 'nearest');
btc = btc(i3);

figure(1); clf
subplot(2,1,1)
plot(ugrid, bta, ugrid, btc)
ax = axis; ax(3) = 200; ax(4) = 300; axis(ax);
legend('airs', 'cris', 'location', 'best')
title('CrIS and interpolated AIRS 1c comparison')
grid on

subplot(2,1,2)
plot(ugrid, bta - btc)
ax = axis; ax(3) = -10; ax(4) = 10; axis(ax);
legend('airs - cris', 'location', 'best')
grid on
zoom on

saveas(gcf, 'CrIS_AIRS_interp_1c', 'fig')

