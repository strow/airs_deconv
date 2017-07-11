
% sample L1c radiance data
adir = '/asl/data/airs/L1C/2016/339';
flist = dir(fullfile(adir, 'AIRS*L1C*.hdf'));
ix = 20;
afile = fullfile(adir, flist(ix).name);
r = double(hdfread(afile, 'radiances'));

% L1c nominal frequencies
v = load('freq2645.txt');

% plot radiance max by channel
r0 = reshape(r, 135 * 90, 2645);
r0 = max(r0);
plot(r0)
grid on

% max index 440, just over 100
% mid index 1500, just under 40
% min index 2350, below 0.05

% channel index by frequency
[~, ix] = min(abs(v - 902.0436));

% channel subset
% r1 = r(:, :, ix); 
  r1 = r(:, :, 440); 
% r1 = r(:, :, 1500); 
% r1 = r(:, :, 2350); 
r2 = unique(r1(:));

% show smallest and largest values
format long
[r2(1:24), r2(end-23:end)]

return

% Sample dithering for radiances from channel at 902.04 cm-1.
% This appears to be discretized in 125-sized chunks starting in 
% the 0.1 position, that is there are only 3 sighificant figures
n = length(r2);
r3 = r2 + 0.1 * randn(n,1);
plot(1:n, r3 - r2)

% factor r = m * 10^e, for 1 <= m < 10
% e = floor(log10(r));
% m = r ./ 10.^e;

% % sanity check
% r2 = m .* 10.^e;
% isclose(r, r2)

% m2 = round2n(m, 11);
% [m3, jx] = unique(m2);
% e3 = e(jx);
% [m3, e3]

