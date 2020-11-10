%
% a2cris_loop -- call airs2cris in chunks, for large files
%
% edit as neededed.  often called with apodization off for older
% tests, since that can be applied later as needed.  but statistical
% correction (scorr = 1) always includes hamming apodization, even
% if hapod is not set.
%

addpath ../source
addpath /asl/packages/ccast/source

% AIRS SRFs for deconvolution
% sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';
  sfile = '/home/sergio/MATLABCODE/airs_l1c_srf_tables_lls_20181205.hdf';

% load the AIRS data
% d1 = load('airs_srf_7377');
  d1 = load('airs_srf_49');
[~, nobs] = size(d1.arad);

% chunk size
k = 200;

opt1 = struct;
opt1.hapod = 1;
opt1.scorr = 1;
opt1.user_res = 'midres';
opt1.cfile = 'corr_midres.mat';

for j = 1 : k : nobs
  
  ix = j : min(j+k-1, nobs);

  [rtmp, cfrq, opt2] = airs2cris(d1.arad(:, ix), d1.afrq, sfile, opt1);

  if j == 1
    [m, ~] = size(rtmp);
    crad = zeros(m, nobs);

    [m, ~] = size(opt2.brad);
    brad = zeros(m, nobs);
  end

  crad(:, ix) = rtmp;
  brad(:, ix) = opt2.brad;

  fprintf(1, '.')
end
fprintf(1, '\n')

bfrq = opt2.bfrq;

% save ac_HR_srf_7377 crad cfrq brad bfrq opt1
% save ac_HR_srf_49   crad cfrq brad bfrq opt1
% save ac_ap_LR_srf_49  crad cfrq brad bfrq opt1
% save ac_cc_LR_srf_49  crad cfrq brad bfrq opt1
% save ac_ap_MR_srf_49  crad cfrq brad bfrq opt1
  save ac_cc_MR_srf_49  crad cfrq brad bfrq opt1
% save ac_cc_MR_old_49  crad cfrq brad bfrq opt1
