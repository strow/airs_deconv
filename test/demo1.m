%
% demo1 - AIRS to CrIS translation demo
% 

% set paths to standard libs
addpath /asl/matlib/h4tools
addpath /asl/packages/ccast/source
addpath /asl/packages/airs_decon/source

% specify an SRF tabulation file
sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';

% specify and load the AIRS data
afile = '/home/chepplew/projects/airs/AIRStoCRIS/airsL1c-Sno-Rad_20121001';
d1 = load(afile);
afrq = d1.arsNomFr;   % AIRS channel frequencies, an m-vector
arad = d1.catRad';    % AIRS channel radiances, an m x k array
% arad = arad(:, 1:20);   % optional subsetting

% do the basic translation
[crad, cfrq, opt2] = airs2cris(arad, afrq, sfile);

% translation with apodization
opt1 = struct;
opt1.hapod = 1;
[crad_hamm, cfrq, opt2] = airs2cris(arad, afrq, sfile, opt1);

% plot sample obs
iobs = 1;
abt = real(rad2bt(afrq, arad(:,iobs)));
cbt = real(rad2bt(cfrq, crad(:,iobs)));
cbt_hamm = real(rad2bt(cfrq, crad_hamm(:,iobs)));
figure(1); clf
plot(afrq, abt, cfrq, cbt, cfrq, cbt_hamm)
title('AIRS L1c and AIRS-to-CrIS sample obs')
legend('AIRS', 'CrIS', 'CrIS Hamm', 'location', 'best')
xlabel('wavenumber')
ylabel('brightness temp')
grid on; zoom on

% save the results
save demo1 crad crad_hamm cfrq opt2

