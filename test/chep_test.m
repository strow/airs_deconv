
addpath ../source
addpath /asl/packages/ccast/source

load('/home/chepplew/myLib/data/airs_f.mat');
% Load up a SNO file
a = dir('/home/chepplew/data/sno/airs_cris/ASL/LR/2017/sno_airs_cris_*frmL1c_v20a.mat');
sn = load([a(1).folder '/' a(1).name]);
opts =struct;
opts.hapod=1;
opts.nguard=2;
opts.dvb=0.1000;
opts.inst_res='lowres';
opts.user_res='lowres';
opts.scorr=0;

% sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';
% sfile = '/home/chepplew/projects/airs/L1C/airs_l1c_srf_tables.hdf'
  sfile = '/home/chepplew/projects/airs/srf_study/airs_2834_srf_table_pert_t_gm.hdf';

[crad, cfrq, opt2] = airs2cris(sn.ra, fairs, sfile, opts);

