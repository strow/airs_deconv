%
% demo1 -- AIRS to CrIS translation demo
% 

% set paths to asl libs
addpath ../source
addpath /asl/packages/ccast/source

% specify an SRF tabulation file
sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';

% load some sample AIRS data
afile = 'data/demo_L1c';
d1 = load(afile);
afrq = d1.arsNomFr;     % AIRS channel frequencies, an m-vector
arad = d1.catRad';      % AIRS channel radiances, an m x k array
arad = arad(:, 1:20);   % optional subsetting

% do the basic translation
[crad, cfrq, opt2] = airs2cris(arad, afrq, sfile);

% translation with apodization
opt1 = struct;
opt1.hapod = 1;
[crad_hamm, cfrq, opt2] = airs2cris(arad, afrq, sfile, opt1);

% select and plot a sample obs
iobs = 3;
abt = real(rad2bt(afrq, arad(:,iobs)));
cbt = real(rad2bt(cfrq, crad(:,iobs)));
cbt_hamm = real(rad2bt(cfrq, crad_hamm(:,iobs)));
figure(1); clf
plot(afrq, abt, cfrq, cbt, cfrq, cbt_hamm)
title('AIRS to CrIS demo')
legend('true airs', 'airs to cris', 'airs cris hamm', 'location', 'best')
xlabel('wavenumber')
ylabel('brightness temp')
grid on; zoom on

