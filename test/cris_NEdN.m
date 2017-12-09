%
% one-day mean ccast CrIS NEdN
%

% cdir = '/asl/data/cris/ccast/sdr60/2016/201';
  cdir = '/asl/data/cris/ccast/sdr60/2016/339';  % 4 Dec 2016
% cdir = '/asl/data/cris/ccast/sdr60_hr/2016/339';
flist = dir(fullfile(cdir, 'SDR*.mat'));

for i = 1 : length(flist)

  cfile = fullfile(cdir, flist(i).name);
  d1 = load(cfile, 'nLW', 'nMW', 'nSW', 'vLW', 'vMW', 'vSW');

  ntmp = [mean(d1.nLW, 3); mean(d1.nMW, 3); mean(d1.nSW, 3)];

  if i == 1, 
    ncnt = 0;
    [m, n] = size(ntmp);
    nsum = zeros(m, n);
  end

  ncnt = ncnt + 1;
  nsum = nsum + ntmp;

  if mod(i, 10) == 0, fprintf(1, '.'), end
end
fprintf(1, '\n')

nedn = nsum / ncnt;
freq = [d1.vLW; d1.vMW; d1.vSW];

  save cris_LR_NEdN cdir ncnt nedn freq
% save cris_HR_NEdN cdir ncnt nedn freq

figure(1); clf
[x, y] = pen_lift(freq, nedn);
semilogy(x, y)
axis([650, 2600, 0, 1])
title('CrIS 1-day mean NEdN')
legend(fovnames)
set(gcf, 'DefaultAxesColorOrder', fovcolors);
xlabel('wavenumber, cm-1')
ylabel('NEdN, mw sr-1 m-2')
grid on

  saveas(gcf, 'CrIS_LR_NEdN', 'png')
% saveas(gcf, 'CrIS_HR_NEdN', 'png')

