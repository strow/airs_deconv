%
% browse SRF tabulations
%

addpath ../source
addpath ../h4tools

warning('off', 'MATLAB:imagesci:hdf:removalWarningHDFSD')

sdir = '/asl/matlab2012/srftest/'
s3 = 'srftables_m130f_withfake_mar08.hdf';
s4 = 'srftables_m140f_withfake_mar08.hdf';
s5 = 'srftables_m150f_withfake_mar08.hdf';

[tf3, tg3, sv3] = srf_read(fullfile(sdir, s3));
[tf4, tg4, sv4] = srf_read(fullfile(sdir, s4));
[tf5, tg5, sv5] = srf_read(fullfile(sdir, s5));

% channel index
ix = 1000;
[tf3(ix), tf4(ix), tf5(ix)]

plot(tg3(ix, :), sv3(ix, :), ...
     tg4(ix, :), sv4(ix, :), ...
     tg5(ix, :), sv5(ix, :));

vx1 = floor(tf3(ix) - 1);
vx2 = ceil(tf3(ix) + 1);
axis([vx1, vx2, 0, 1])
title('sample AIRS SRF tabulations')
legend('m130', 'm140', 'm150')
xlabel('wavenumber')
ylabel('weight')
grid on; zoom on


