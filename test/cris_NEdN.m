%
% one-day mean ccast CrIS NEdN
%
% notes:
%   for now these are for a single FOV, selected below
%   checking a few days, the means seem to be consistent
%

iFOV = 5;

% cdir = '/asl/data/cris/ccast/sdr60/2016/201'
  cdir = '/asl/data/cris/ccast/sdr60/2016/339'
flist = dir(fullfile(cdir, 'SDR*.mat'));

for i = 1 : length(flist)

  cfile = fullfile(cdir, flist(i).name);
  d1 = load(cfile, 'nLW', 'nMW', 'nSW', 'vLW', 'vMW', 'vSW');

  ntmp = [squeeze(d1.nLW(:, iFOV, 1)); ...
          squeeze(d1.nMW(:, iFOV, 1)); ...
          squeeze(d1.nSW(:, iFOV, 1))];

  if i == 1, 
    ncnt = 0;
    nsum = zeros(length(ntmp), 1);
  end

  ncnt = ncnt + 1;
  nsum = nsum + ntmp;

  if mod(i, 10) == 0, fprintf(1, '.'), end
end
fprintf(1, '\n')

nedn = nsum / ncnt;
freq = [d1.vLW; d1.vMW; d1.vSW];

save cris_NEdN cdir iFOV ncnt nedn freq

figure(1); clf
semilogy(freq, nedn, 'linewidth', 2)
axis([650, 2600, 0, 1])
title('CrIS 1-day mean NEdN')
xlabel('wavenumber')
ylabel('NEdN')
grid on

