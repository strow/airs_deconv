
% call mkaconv2 with some recent tests parameters

% load 1b channel frequencies
% d1 = load('airs/freq1b');  
% cfreq = d1.freq1b;

% load 1c channel frequencies
% d1 = load('freq1c');  
% cfreq = d1.freq1c;
d1 = load('data/may1to1000sno_mean.mat');
cfreq = d1.freq_fill;

% specify an SRF tabulation file
sfile = '/asl/matlab/srftest/srftables_m140f_withfake_mar08.hdf';

% set the output mat file
% afile = 'airs/SRFtest1B';
% afile = 'airs/SRFtest1C';
afile= 'SRFtest1Cb';

% call mkaconv2 to build the convolution matrix
mkaconv2(cfreq, sfile, afile)

