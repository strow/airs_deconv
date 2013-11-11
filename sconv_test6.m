

addpath /home/motteler/mot2008/hdf/h4tools
addpath /home/motteler/cris/bcast/motmsc/utils

sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';
cfreq = load('freq2645.txt');
dvs = 0.0025; 

tic,
[sconv, sfreq, tfreq] = mksconv2(sfile, cfreq, dvs);
toc

% now, call mksconv for a check
tic
mksconv(sfile, 'tempx', dvs)
toc
d1 = load('tempx');

ix1 = interp1(d1.Cfin, 1:length(d1.Cfin), sfreq, 'nearest');
jx1 = interp1(d1.Cfout, 1:length(d1.Cfout), tfreq, 'nearest');

sfreq1 = d1.Cfin(ix1);
streq1 = d1.Cfout(jx1);
sconv1 = d1.Cmat(jx1, ix1);

rms(sconv1(:) - sconv(:)) / rms(sconv)
max(sconv1(:) - sconv(:))
