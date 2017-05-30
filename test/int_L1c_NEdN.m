%
% take AIRS L1c NEdN averages and interpolate missing values
%

adir = '/asl/data/airs/L1C/2016/339';
flist = dir(fullfile(adir, 'AIRS*L1C*.hdf'));

% gran = 'AIRS.2016.12.04.229.L1C.AIRS_Rad.v6.1.2.0.G16340122043.hdf';
% afile = fullfile(adir, gran);

% d1 = hdfinfo(afile);
% [alist, fattr] = h4sdread(afile);

% running mean
nOK = zeros(2645, 1);
sOK = zeros(2645, 1);

for i = 1 : length(flist)

  afile = fullfile(adir, flist(i).name);
  NeN = hdfread(afile, 'NeN');

  % reshape to nchan x nobs array
  NeN = reshape(NeN, 90 * 135, 2645)';

  for j = 1 : 90 * 135
    iOK = NeN(:, j) < 2;
    nOK = nOK + iOK;
    sOK = sOK + iOK .* NeN(:, j);
  end

  if mod(i, 10) == 0, fprintf(1, '.'), end
end
fprintf(1, '\n')

mNeN = sOK ./ nOK;
iOK = nOK > 0;
v1c = load('freq2645.txt');
nedn = interp1(v1c(iOK), mNeN(iOK), v1c, 'linear', 'extrap');
save int_L1c_NEdN v1c iOK mNeN nedn 

figure(1); clf
ix = ~iOK;
semilogy(v1c, mNeN, 'o', v1c(ix), nedn(ix), '+')
axis([650, 2750, 0, 1])
title('measured and interpolated L1c NEdN')
legend('measured', 'interpolated')
xlabel('wavenumber')
ylabel('NEdN')
grid on

% saveas(gcf, 'int_L1c_NEdN', 'png')

