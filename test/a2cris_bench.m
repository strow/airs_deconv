%
% a2cris_bench -- call airs2cris in chunks, for large files
%
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
% d1 = load('airs_fit49');
[~, nobs] = size(d1.arad);

% chunk size
k = 200;

opt1 = struct;
opt1.hapod = 0;             % flag for Hamming apodization
opt1.scorr = 0;             % flag for statistical correction
opt1.dvb = 0.1;             % deconvolution frequency step
opt1.bfile = 'bconv.mat';   % deconvolution temp file

opt1.inst_res = 'lowres';   % default is lowres
opt1.user_res = 'lowres';   % default is lowres

profile clear
profile on

for j = 1 : k : nobs
  
  ix = j : min(j+k-1, nobs);

  [rtmp, cfrq, opt2] = airs2cris(d1.arad(:, ix), d1.afrq, sfile, opt1);

% [rtmp2, cfrq2, opt3] = airs2cris(d1.arad(:, ix), d1.afrq, sfile, opt1);

% if ~isequal(rtmp, rtmp2) || ~isequal(cfrq, cfrq2)
%   keyboard
% end

% if ~isequal(opt2.brad, opt3.brad) || ...
%    ~isequal(opt2.bfrq, opt3.bfrq) || ...
%    ~isequal(opt2.afrq, opt3.afrq)
%   keyboard
% end

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

profile report

bfrq = opt2.bfrq;

% save acris_cloudy crad cfrq brad bfrq
% save acrHR_fit49 crad cfrq brad bfrq
% save acris_fit49 crad cfrq brad bfrq
% save ac_ap_fit49 crad cfrq brad bfrq
% save ac_cc_fit49 crad cfrq brad bfrq

