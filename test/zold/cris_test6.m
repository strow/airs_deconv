%
% cris_test6 - AIRS to CrIS translation test
% 

% set paths to asl libs
addpath ../source
addpath ../h4tools
addpath /asl/packages/ccast/source

% specify an SRF tabulation file
sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';

% specify and load the AIRS data
afile = '/asl/s1/chepplew/projects/airs/AIRStoCRIS/airsL1c-Sno-Rad_20121001';
d1 = load(afile);
afrq = d1.arsNomFr;     % AIRS channel frequencies, an m-vector
arad = d1.catRad';      % AIRS channel radiances, an m x k array
arad = arad(:, 1:1000);  % optional subsetting

% do the basic translation
[crad, cfrq, opt2] = airs2cris(arad, afrq, sfile);

% translation with apodization
opt1 = struct;
opt1.hapod = 1;

profile on

[crad_hamm, cfrq, opt2] = airs2cris(arad, afrq, sfile, opt1);

profile report
profile off
