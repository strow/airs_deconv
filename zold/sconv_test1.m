
% first test of sconv transform
%

% path to bcast finterp and inst_params
addpath ../source

% get Larrabee's sample data
load airs/motteler_samples.mat

% choose AIRS and CrIS samples
abt = double(bta_mid_nh);
cbt = double(btc_mid_nh);

% sort the AIRS channel set
[fs, ix] = sort(f);

% permute AIRS data accordingly
abt = abt(ix);

% switch to radiances for the convolution
arad = bt2rad(fs, abt);

% load transform matrix
load /home/motteler/sconv/cmat_small

% sort (and take the 2378 chan subset)
S = Cmat(ix, :);  
vchan = Cfout(ix);
vgrid = Cfin;
clear Cmat Cfout Cfin

Sinv = S' * inv(S * S');  % explict right-inverse

% Sinv = pinv(full(S));  % matlab pseudo-inverse, SLOW

asrf = Sinv * arad;

% plot(vgrid, asrf)

plot(vgrid, real(rad2bt(vgrid, asrf)))

return

% get CrIS user params, specify band here
[inst, user] = inst_params('MW', 773.13);

% fourier interpolation of rad w.r.t. SRF basis
[rad1, frq1] = finterp(asrf, vout, user.dv);

% fourier interpolation after simple linear interp
opt = struct;
opt.dv1 = 0.1;  
[rad2, frq2] = finterp(arad, fs, user.dv, opt);

% change AIRS results back to brightness temps
frq1 = frq1(:);
frq2 = frq2(:);
bt1 = real(rad2bt(frq1, rad1));
bt2 = real(rad2bt(frq2, rad2));

% plot(frq1, bt1);
plot(frq1, bt1, frq2, bt2);

return

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

