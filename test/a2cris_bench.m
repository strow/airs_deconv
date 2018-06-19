%
% a2cris_bench -- call airs2cris in chunks, for large files
%
%

addpath /asl/packages/ccast/source
addpath ../source

% AIRS SRFs for deconvolution
sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';

adir = '/asl/data/airs/L1C/2016/339';
flist = dir(fullfile(adir, 'AIRS*L1C*.hdf'));
ix = 20;
afile = fullfile(adir, flist(ix).name);
arad = hdfread(afile, 'radiances');
arad = permute(arad, [3,2,1]);
arad = reshape(arad, 2645, 90*135);
afrq = load('freq2645.txt');
[~, nobs] = size(arad);

% chunk size
k = 400;

opt1 = struct;
opt1.hapod = 0;             % flag for Hamming apodization
opt1.scorr = 0;             % flag for statistical correction
opt1.dvb = 0.2;             % deconvolution frequency step

tic
profile clear
profile on

for j = 1 : k : nobs
  
  ix = j : min(j+k-1, nobs);

  [rtmp, cfrq, opt2] = airs2cris(arad(:, ix), afrq, sfile, opt1);

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
toc

