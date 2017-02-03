
addpath /home/strow/Work/Airs/L1c/Sim

cd /home/strow/Data/L1c/Sim/2008/12/02/airicrad
fn130 = 'AIRS.2008.12.02.001.L1C.AIRS_Rad.v6.1.2.0.SIM_rtp0905_cldym130_l1c.T16284141105.hdf';
fn140 = 'AIRS.2008.12.02.001.L1C.AIRS_Rad.v6.1.2.0.SIM_rtp0905_cldym140_l1c.T16284140806.hdf';
fn150 = 'AIRS.2008.12.02.001.L1C.AIRS_Rad.v6.1.2.0.SIM_rtp0905_cldym150_l1c.T16278141012.hdf';

% Read in radiances from each file
d130 = read_airs_l1c(fn130);
d140 = read_airs_l1c(fn140);
d150 = read_airs_l1c(fn150);

cd /home/motteler/cris/airs_decon/test

r130 = d130.radiances(64:73, :, :); 
r140 = d140.radiances(64:73, :, :);
r150 = d150.radiances(64:73, :, :);

whos r130 r140 r150

r130 = reshape(r130, 900, 2645)';
r140 = reshape(r140, 900, 2645)';
r150 = reshape(r150, 900, 2645)';

whos r130 r140 r150

save jpl_test_data r130 r140 r150

