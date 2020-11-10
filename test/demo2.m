%
% demo2 -- call airs2cris in chunks, for larger data sets
%

% local libs
addpath ../data
addpath ../source
addpath /asl/packages/ccast/source

% AIRS SRF tabulation file
sfile = 'airs_l1c_srf_tables_lls_20181205.hdf';

% load an AIRS granule
apath = '/asl/airs/l1c_v672/2018/091';
agran = 'AIRS.2018.04.01.222.L1C.AIRS_Rad.v6.7.2.0.G20008080043.hdf';
afile = fullfile(apath, agran);
arad = hdfread(afile, 'radiances');
arad = permute(arad, [3,2,1]);
arad = reshape(arad, 2645, 90*135);
afrq = load('freq2645.txt');
[~, nobs] = size(arad);

% translation options
opt1 = struct;
opt1.user_res = 'lowres';   % target resolution
opt1.hapod = 0;             % no Hamming apodization
opt1.scorr = 0;             % no statistical correction
k = 400;                    % translation chunk size

% profile clear
% profile on
tic

% loop on chunks
for j = 1 : k : nobs
  
  % indices for current chunk
  ix = j : min(j+k-1, nobs);

  % call airs2cris on the chunk
  [rtmp, cfrq, opt2] = airs2cris(arad(:, ix), afrq, sfile, opt1);

  % initialize output after first obs
  if j == 1
    [m, ~] = size(rtmp);
    [n, ~] = size(opt2.brad);
    crad = zeros(m, nobs);
    brad = zeros(n, nobs);
  end

  % save the current chunk
  crad(:, ix) = rtmp;
  brad(:, ix) = opt2.brad;

% fprintf(1, '.')
end
fprintf(1, '\n')

toc
% profile report

