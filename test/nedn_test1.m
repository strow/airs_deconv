%
% nedn_test1 -- noise estimate for the AIRS-to-CrIS translation
%
% uses 24-hr NEdN means from cris_NEdN.m and int_L1c_NEdN.m
%
% saves apodized and unapodized NEdN estimates for AIRS, CrIS, 
% and AIRS-to-CrIS 

addpath /asl/packages/ccast/source
addpath ../source

% AIRS to CrIS translation options
opt1 = struct;
opt1.hapod = 1;   % 0 = no apod, 1 = hamming apod
opt1.scorr = 0;   % 0 = no corr, 1 = statistical corr
opt1.user_res = 'lowres';   % 'lowres' or 'hires'

% AIRS SRF file (for airs2cris)
sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';

% test params
nrad = 100;    % simulated radiance obs per set
nset = 10;     % number of simulated radiance sets

% CrIS NEdN estimate
d1 = load('cris_NEdN');
nedn_cris = d1.nedn;
freq_cris = d1.freq;

% apodized CrIS NEdN estimate
if opt1.hapod 
  nedn_cris = nedn_cris * 0.63;
end

% AIRS L1c noise estimate
d1 = load('int_L1c_NEdN.mat');
nedn_airs = d1.nedn;
freq_airs = d1.v1c;
nchan = length(freq_airs);

% AIRS radiance at 280K
r_280K = bt2rad(freq_airs, 280) * ones(1, nrad);

% loop on radiance sets
for i = 1 : nset

  % add noise scaled to the AIRS NEdN spec
  r_airs = r_280K + randn(nchan, nrad) .* (nedn_airs * ones(1, nrad));

  % AIRS to CrIS translation
  [r_tran, freq_tran] = airs2cris(r_airs, freq_airs, sfile, opt1);
  r_tran = real(r_tran);

  % initialize tables on the first iteration
  if i == 1
    tab_airs = zeros(length(freq_airs), nset);
    tab_tran = zeros(length(freq_tran), nset);
  end

  % measure and save the simulated noise (as a check)
  tab_airs(:, i) = std(r_airs, 0, 2);

  % measure and save the translated simulated noise
  tab_tran(:, i) = std(r_tran, 0, 2);

  fprintf(1, '.')
end
fprintf(1, '\n')

% take means over the std sets
nedn_asim = mean(tab_airs, 2);
nedn_tran = mean(tab_tran, 2);

% plot the results
figure(1); clf
[x1, y1] = pen_lift(freq_airs, nedn_airs);
[x2, y2] = pen_lift(freq_airs, nedn_asim);
[x3, y3] = pen_lift(freq_tran, nedn_tran);
[x4, y4] = pen_lift(freq_cris, nedn_cris);
semilogy(x1, y1, x2, y2, x3, y3, x4, y4);
axis([600, 2600, 0, 1])
title('AIRS to CrIS NEdN estimates')
legend('AIRS spec', 'AIRS test', 'AIRS to CrIS', 'CrIS ccast')
xlabel('wavenumber')
ylabel('NEdN')
grid on; zoom on

return

if opt1.hapod 
  nfile = 'nedn_hamm';
else
  nfile = 'nedn_noap';
end

save(nfile, 'freq_airs', 'nedn_airs', ...
            'freq_tran', 'nedn_tran', ...
            'freq_cris', 'nedn_cris', 'opt1', 'opt2')

