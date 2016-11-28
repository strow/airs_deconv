%
% translate AIRS to CrIS with an inverse convolution matrix
%
% this version does the deconvolution to the SRF tabulation grid
%
% main processing steps
%  - load Larrabee's AIRS and CrIS colocation data
%  - generate the AIRS "A" and "B" convolution matrices
%  - apply B inverse to the AIRS radiance data
%  - interferometric interpolation to CrIS band specs
%  - compare results
%

% band spec for convolution and plot labels
band = 'SW';

% path to bcast finterp and inst_params
addpath /home/motteler/cris/bcast/source

% load the new matchup data
% d1 = load('data/may1to1000sno_mean.mat');
d1 = load('data/sno_global_mean_night');

% select the airs and cris data
abt1 = double(d1.btal1cf_mean);   fignote = 'global_mean_night';
cbt1 = double(d1.btc_mean);  

% fc is unapodized CrIS frequencies
% d2 = load('data/fc.mat');
% [fcris, ix] = sort(d2.fc);
% cfreq = d1.freq_fill;
[fcris, ix] = sort(d1.fc);
cbt1 = cbt1(ix);
cfreq = d1.ff;

% plot setup
figtype = 'fig';    % plot data type
figtstr = fignote;
figtstr(find(fignote == '_')) = ' ';

% SRF tabulation and convolution parameters
hfile = '/asl/matlab/srftest/srftables_m140f_withfake_mar08.hdf';
afile = 'SRF_1c_A1.mat';   % channel grid deconvolution file
bfile = 'SRF_1c_B1.mat';   % SRF tab grid deconvolution file

dvs = 0.05;                % SRF tab grid step

% do the convolutions
% mkaconv2(cfreq, hfile, afile)   % only call when params change
% mksconv(hfile, bfile, dvs)      % only call when params change
A1 = load(afile);
B1 = load(bfile);

% permute AIRS data accordingly
abt1 = abt1(A1.aind);

% switch to radiances for the convolution
arad1 = bt2rad(A1.afrq, abt1);

% apply the inverse AIRS convolution 
arad2 = pinv(full(B1.Cmat(A1.sind,:))) * arad1;
abt2 = real(rad2bt(B1.Cfin, arad2)); 

% get CrIS user grid params by band
[inst, user] = inst_params(band, 773.13);
[arad3, frq3] = finterp(arad2, B1.Cfin, user.dv);
abt3 = real(rad2bt(frq3, arad3)); 

% 4-way comparison representing our processing steps: true AIRS,
% AIRS deconvolved, AIRS deconvolved and interpolated to CrIS, and
% true CrIS
figure(1); clf
plot(A1.afrq, abt1, B1.Cfin, abt2, frq3, abt3, fcris, cbt1)
axis([600, 2800, 200, 300])
xlabel('frequency');
ylabel('brightness temp');
title(['AIRS and CrIS co-location ',band,' ',figtstr])
legend('true AIRS', 'AIRS deconv', ['AIRS CrIS ', band], ...
       'true CrIS', 'location', 'north')
grid; zoom on
saveas(gcf, ['bconv_all_',band,'_',fignote], figtype)

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
plot(ugrid, abt4, ugrid, cbt2)
ax = axis; ax(3) = 200; ax(4) = 300;
axis(ax)
xlabel('frequency');
ylabel('brightness temp');
title(['AIRS and true CrIS ',band,' ',figtstr])
legend('AIRS CrIS', 'true CrIS', 'location', 'best')
grid; zoom on

subplot(2,1,2)
plot(ugrid, abt4 - cbt2)
ax = axis; ax(3) = -10; ax(4) = 10;
axis(ax)
xlabel('frequency');
ylabel('brightness temp');
% title(['AIRS CrIS minus true CrIS ',band,' ',figtstr])
legend('AIRS CrIS minus true CrIS ', 'location', 'best')
grid; zoom on
saveas(gcf, ['bconv_cmp_',band,'_',fignote], figtype)

