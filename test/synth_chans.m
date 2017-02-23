%
% plot and save a list of AIRS synthetic channels
%

addpath ../h4tools
addpath ../source

% turn off HDF 4 update warnings
warning('off', 'MATLAB:imagesci:hdf:removalWarningHDFSD')

% get the L1c channel set
L1c_frq = load('freq2645.txt');
nL1c = length(L1c_frq);

% get L1b channel set from an SRF tabulation
sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';

% read the srf data
[alist, fattr] = h4sdread(sfile);

for i = 1 : length(alist)
  switch alist{i}{1}
    case 'chanid', chanid = double(alist{i}{2})';
    case 'freq',   tfreq   = double(alist{i}{2})';
    case 'fwgrid', fwgrid = double(alist{i}{2})';
    case 'srfval', srfval = double(alist{i}{2})';
    case 'width',  width  = double(alist{i}{2})';
  end
end

% 1:2378 are the L1B channels, we want them sorted
L1b_frq = sort(tfreq(1:2378));
nL1b = length(L1b_frq);

% take the intersection of L1b and L1c channels
[ib, ic] = seq_match(L1b_frq, L1c_frq);

% any L1c channel not in the intersection is synthetic
ind_syn = setdiff((1:nL1c)', ic);
L1c_syn = L1c_frq(ind_syn);

save L1c_syn L1b_frq L1c_frq L1c_syn ind_syn

figure(1); clf
plot(1:nL1c, L1c_frq, ind_syn, L1c_syn, 'o')
axis([0, 2700, 600, 2800])
title('AIRS L1c channels')
legend('full L1c set', 'synthetic channels', 'location', 'northwest')
xlabel('index')
ylabel('frequency')
grid on;
saveas(gcf, 'AIRS_L1c_synth', 'png')



