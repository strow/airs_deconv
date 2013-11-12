
% find a useful lower bound for AIRS 1b channel step size
%
% some regular AIRS 1b channels, looked at in sorted order, are too
% close together; this makes the AIRS SRF convolution matrix poorly
% conditioned.
%
% - figure 1 shows a proposed linear function to trim channels
% - figure 2 shows airs channels and diffs by airs index
% - figure 3 shows airs channels and diffs by sorted frequency
% 

addpath /home/motteler/mot2008/hdf/h4tools

sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';

% read the srf data
[alist, fattr] = h4sdread(sfile);

for i = 1 : length(alist)
  switch alist{i}{1}
    case 'chanid', chanid = double(alist{i}{2})';
    case 'freq',   freq   = double(alist{i}{2})';
    case 'fwgrid', fwgrid = double(alist{i}{2})';
    case 'srfval', srfval = double(alist{i}{2})';
    case 'width',  width  = double(alist{i}{2})';
  end
end
clear alist fattr

% start with the L1B channel set
nchan = 2378;
freq = freq(1:nchan);

% look at sorted frequencies
fs = sort(freq);
dfs = diff(fs);

% get a reasonable lower bound on channel step size
% as a function a*x + b of frequency
a = 4e-4;
b = -0.04;

figure(1); clf
plot(fs(1:nchan-1), dfs, fs, a*fs+b);
axis([600, 2700, 0, 1.2])
xlabel('frequency')
ylabel('step size')
title('AIRS channel step lower bound')
legend('channel steps', 'lower bound', 'location', 'south')
grid; zoom on
saveas(gcf, 'AIRS_chans_bound', 'fig')

% return

figure(2); clf
subplot(2,1,1)
plot(1:nchan, freq);
axis([1, 2400, 640, 2700])
xlabel('index')
ylabel('frequency')
title('AIRS channel frequency by index')
grid; zoom on

subplot(2,1,2)
plot(1:nchan-1, diff(freq));
axis([1, 2400, 0, 1.2])
xlabel('index')
ylabel('difference')
title('AIRS channel diffs by index')
grid; zoom on
saveas(gcf, 'AIRS_chans_by_ind', 'fig')

figure(3); clf
subplot(2,1,1)
plot(1:nchan, fs);
axis([1, 2400, 640, 2700])
xlabel('index')
ylabel('freqency')
title('AIRS sorted channel freq by index')
grid; zoom on

subplot(2,1,2)
plot(1:nchan-1, diff(fs));
axis([1, 2400, 0, 1.2])
xlabel('index')
ylabel('difference')
title('AIRS sorted channel diffs by index')
grid; zoom on
saveas(gcf, 'AIRS_sorted_chans', 'fig')

