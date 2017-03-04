%
% L1d_test1 -- basic test of L1d_conv
%

% build the convolution matrix
res = 1400;
dv_tab = 0.1;
[dconv, v_tab, v_L1d] = L1d_conv(res, dv_tab);

% select a channel for the plot
v0 = input('plot chan > ');
n_L1d = length(v_L1d);
ix = interp1(v_L1d, 1:n_L1d, v0, 'nearest');
jx = ix-2 : ix+2;
v1 = v_L1d(ix-5);
v2 = v_L1d(ix+5);
dtmp = dconv(jx, :);
dmax = max(dtmp(:));

% index for directly calculated SRF
vx = v_L1d(ix);
z = sup_gauss(v_tab, vx, 2*vx/res);
z = z * dmax;

figure(1); clf
plot(v_tab, dconv(jx, :), v_tab, z, 'o')
ax = axis; ax(1) = v1; ax(2) = v2; axis(ax)
title('AIRS L1d SRFs')
grid on; zoom on
xlabel('wavenumber')
ylabel('weight')

