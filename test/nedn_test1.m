%
% nedn_test1 -- noise estimate for the AIRS to CrIS translation
%

addpath /asl/packages/ccast/source
addpath ../source

% test params
nrad = 100;          % simulated radiance obs per set
nset = 10;           % number of simulated radiance sets

% sample ccast CrIS NEdN estimate
d = load('/asl/data/cris/ccast/sdr60/2016/018/SDR_d20160118_t0801033.mat');
nedn_ccast = [squeeze(d.nLW(:,5,1));squeeze(d.nMW(:,5,1));squeeze(d.nSW(:,5,1))];
freq_ccast = [d.vLW; d.vMW; d.vSW];
clear d

% AIRS to CrIS translation options
opt1 = struct;
opt1.hapod = 1;
opt1.resmode = 'lowres';

% AIRS SRF file (for airs2cris)
sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';

% load the AIRS L1c noise spec
d1 = load('int_L1c_NEdN.mat');
nedn_airs = d1.nedn;
freq_airs = d1.v1c;
nchan = length(freq_airs);

% get AIRS expected ICT radiance
r_280K = bt2rad(freq_airs, 280) * ones(1, nrad);

% loop on radiance sets
for i = 1 : nset

  % add noise scaled by the AIRS NEdN spec
  r_airs = r_280K + randn(nchan, nrad) .* (nedn_airs * ones(1, nrad));

  % AIRS to CrIS translation
  [r_cris, freq_cris] = airs2cris(r_airs, freq_airs, sfile, opt1);
  r_cris = real(r_cris);

  % initialize tables on the first iteration
  if i == 1
    tab_airs = zeros(length(freq_airs), nset);
    tab_cris = zeros(length(freq_cris), nset);
  end

  % measure and save the simulated noise (as a check)
  tab_airs(:, i) = std(r_airs, 0, 2);

  % measure and save the translated simulated noise
  tab_cris(:, i) = std(r_cris, 0, 2);

  fprintf(1, '.')
end
fprintf(1, '\n')

% take means over the std sets
nedn_asim = mean(tab_airs, 2);
nedn_cris = mean(tab_cris, 2);

% plot the results
semilogy(freq_airs, nedn_airs, freq_airs, nedn_asim, ...
         freq_cris, nedn_cris, freq_ccast, nedn_ccast);

axis([600, 2600, 0, 1])
title('AIRS to CrIS NEdN estimates')
legend('AIRS spec', 'AIRS test', 'AIRS to CrIS', 'CrIS ccast')
xlabel('wavenumber')
ylabel('NEdN')
grid on; zoom on

% saveas(gcf, 'a2cris_nedn_noap', 'png')
% save a2cris_nedn_noap freq_airs nedn_airs freq_cris nedn_cris ...
%      freq_ccast nedn_ccast opt1

% save a2cris_nedn_hamm freq_airs nedn_airs freq_cris nedn_cris ...
%      freq_ccast nedn_ccast opt1

