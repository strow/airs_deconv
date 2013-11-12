%
% translate AIRS to CrIS with an inverse convolution matrix
%
% this version is for any AIRS 1b data

% main processing steps
%  - load the AIRS data
%  - load an AIRS nchan x nchan convolution matrix
%  - apply the inverse convolution to the AIRS data
%  - fourier interpolate the result to CrIS band specs
%  - plot results

% band spec for convolution and plot labels
band = 'LW';

% path to bcast finterp and inst_params
addpath ../source

% load the data
% d1 = load('airs/mean_airs.mat');
% abt1 = [d1.mbtobs, d1.mbtcal];
d1 = load('/home/imbiriba/projects/L1c/new_airs_night_tropics_gstats.mat');
arad1 = d1.dat1.robs1_avg;

% load transform matrix: defines Amat, afrq, and aind
load airs/SRFtest1B.mat
% load airs/SRFtest1C.mat

% permute AIRS data accordingly
% abt1 = abt1(aind,:);
arad1 = arad1(aind,:);

% switch to radiances for the convolution
% arad1 = bt2rad(afrq, abt1);
abt1 = real(rad2bt(afrq, arad1));

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

% 3-way comparison representing our processing steps: true AIRS,
% AIRS deconvolved, and AIRS deconvolved and interpolated to CrIS
figure(1); clf
j = 2;
plot(afrq, abt1(:,j), afrq, abt2(:,j), frq3, abt3(:,j))
axis([600, 2800, 200, 300])
xlabel('frequency');
ylabel('brightness temp');
title(['AIRS ',band])
legend('true AIRS', 'AIRS deconv', ['AIRS CrIS ', band], ...
       'location', 'north')
grid; zoom on
% saveas(gcf, ['AIRS_CrIS_all_',band,'_',test], 'fig')

