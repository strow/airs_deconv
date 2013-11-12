
% demo using finterp to take AIRS to CrIS radiances
%
% the key idea is that finterp does a linear interpolation 
% to the dv1 spacing and then a fourier interpolation to dv2.

% path to bcast finterp and inst_params
addpath ../source

% get Larrabee's sample data
load airs/motteler_samples.mat

% choose AIRS and CrIS samples
bta = bta_mid_nh;
btc = btc_mid_nh;

% get CrIS user params, specify band 
[inst, user] = inst_params('LW', 773.13);

% sort the AIRS channel set
[fs, ix] = sort(f);

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

subplot(2,1,1)
plot(ugrid, bta, ugrid, btc)
legend('airs', 'cris')

subplot(2,1,2)
plot(ugrid, bta - btc)
legend('airs - cris')

