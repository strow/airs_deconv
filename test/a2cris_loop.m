%
% a2cris_loop -- call airs2cris in chunks, for large files
%
% typically called with hamming apodization off, since that can be
% applied later as needed.  But statistical correction (scorr = 1)
% always includes hamming apodization, even if hapod is not set.
%

addpath ../source
addpath /asl/packages/ccast/source

% AIRS SRFs for deconvolution
sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';

% load the AIRS data
  d1 = load('airs_cloudy');
% d1 = load('airs_fit49');
[~, nobs] = size(d1.arad);

% chunk size
k = 200;

opt1 = struct;
opt1.hapod = 0;
opt1.scorr = 0;
opt1.inst_res = 'hires3';
opt1.user_res = 'hires';

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

% save acris_cloudy crad cfrq brad bfrq
% save acris_fit49 crad cfrq brad bfrq
  save ac_HR_cloudy crad cfrq brad bfrq
% save ac_HR_fit49 crad cfrq brad bfrq opt1

