%
% a2cris_loop -- call airs2cris in chunks, for large files
%

addpath ../source
addpath ../h4tools
addpath /asl/packages/ccast/source

% turn off HDF 4 update warnings
warning('off', 'MATLAB:imagesci:hdf:removalWarningHDFSD')

% AIRS SRFs for deconvolution
sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';

% load the AIRS data
d1 = load('airs_cloudy');
[~, nobs] = size(d1.arad);

% chunk size
k = 200;

opt1 = struct;
opt1.hapod = 0;                % flag for Hamming apodization
opt1.dvb = 0.1;                % deconvolution frequency step
opt1.bfile = 'bconv_tmp.mat';  % deconvolution temp file

for j = 1 : k : nobs
  
  ix = j : min(j+k-1, nobs);

  [rtmp, cfrq, opt2] = airs2cris(d1.arad(:, ix), d1.afrq, sfile, opt1);

  if j == 1
    [m, ~] = size(rtmp);
    crad = zeros(m, nobs);
  end

  crad(:, ix) = rtmp;

  fprintf(1, '.')
end
fprintf(1, '\n')

save acris_cloudy crad cfrq

