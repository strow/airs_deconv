
% translate AIRS to CrIS with an inverse convolution matrix
%
% main processing steps
%  - load AIRS and CrIS obs from Larrabee's test data
%  - load an nchan x nchan convolution matrix from mkaconv1
%  - apply the inverse convolution to the AIRS data
%  - fourier interpolate the result to CrIS band specs
%  - compare results

% band spec for convolution and plot labels
band = 'LW';

% path to bcast finterp and inst_params
addpath ../source

% get Larrabee's sample data
load airs/sno_cris_airs_1b.mat

% choose AIRS and CrIS samples
abt1 = double(bta_mid_nh);  test = 'mid_nh';
cbt1 = double(btc_mid_nh);
% abt1 = double(bta_mid_sh);  test = 'mid_sh';
% cbt1 = double(btc_mid_sh);
% abt1 = double(bta_trop);  test = 'trop';
% cbt1 = double(btc_trop);
% abt1 = double(bta_trop_cold);  test = 'trop_cold';
% cbt1 = double(btc_trop_cold);
% abt1 = double(bta_polar_nh);  % lots of NaNs here
% cbt1 = double(btc_polar_nh);
% abt1 = double(bta_polar_sh);  test = 'polar_sh';
% cbt1 = double(btc_polar_sh);

% load transform matrix: defines Amat, afrq, and aind
load airs/SRFtest.mat

% permute AIRS data accordingly
abt1 = abt1(aind);

% switch to radiances for the convolution
arad1 = bt2rad(afrq, abt1);

% apply the inverse AIRS convolution 
Ainv = inv(Amat);
arad2 = Ainv * arad1;
abt2 = real(rad2bt(afrq, arad2)); 

% get CrIS user grid params by band
[inst, user] = inst_params(band, 773.13);
opt = struct;
opt.dv1 = 0.05;  
[arad3, frq3] = finterp(arad2, afrq, user.dv, opt);
abt3 = real(rad2bt(frq3, arad3)); 

% 4-way comparison representing our processing steps: true AIRS,
% AIRS deconvolved, AIRS deconvolved and interpolated to CrIS, and
% true CrIS
figure(1); clf
plot(afrq, abt1, afrq, abt2, frq3, abt3, fcris, cbt1)
axis([600, 2800, 200, 300])
xlabel('frequency');
ylabel('brightness temp');
title('AIRS and CrIS co-location')
legend('true AIRS', 'AIRS deconv', ['AIRS CrIS ', band], ...
       'true CrIS', 'location', 'north')
grid; zoom on
saveas(gcf, 'AIRS_CrIS_1', 'fig')

% match AIRS CrIS with true CrIS user grid
ugrid = user.v1 : user.dv : user.v2;
ix = interp1(frq3, 1:length(frq3), ugrid, 'nearest');
abt4 = abt3(ix);
% frq4 = frq3(ix);  
% max(abs(frq4 - ugrid)) % should be zero

% match Larrabee's CrIS with true CrIS user grid
ix = interp1(fcris, 1:length(fcris), ugrid, 'nearest');
cbt2 = cbt1(ix);
% cfr2 = fcris(ix);
% max(abs(cfr2 - ugrid)) % should be zero

% 2-way comparison of true CrIS and AIRS CrIS, for 1 band
figure(2); clf
subplot(2,1,1)
plot(ugrid, cbt2, ugrid, abt4)
ax = axis; ax(3) = 200; ax(4) = 300;
axis(ax)
xlabel('frequency');
ylabel('brightness temp');
title(['true and AIRS CrIS ', band])
legend('true CrIS', 'AIRS CrIS', 'location', 'best')
grid; zoom on

subplot(2,1,2)
plot(ugrid, cbt2 - abt4)
ax = axis; ax(3) = -10; ax(4) = 10;
axis(ax)
xlabel('frequency');
ylabel('brightness temp');
% title(['true CrIS minus AIRS CrIS ', band])
legend('true CrIS minus AIRS CrIS ', 'location', 'best')
grid; zoom on
saveas(gcf, ['AIRS_CrIS_',band,'_',test,'_1b'], 'fig')

