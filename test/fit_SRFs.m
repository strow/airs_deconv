%
% try fitting some real SRFs with z_srf
%

addpath ../h4tools
addpath ../source

% turn off HDF 4 update warnings
warning('off', 'MATLAB:imagesci:hdf:removalWarningHDFSD')

% sample AIRS srfs
sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';
cfreq = load('freq2645.txt');
dvk = 0.0025; 
[sconv, sfreq, ofreq] = mksconv2(sfile, cfreq, dvk);

%------------------
% fitted AIRS SRFs
%------------------

% SRF set 1
ix = 201 : 206;
vwid = 0.48;

vtmp = ofreq(ix)
dv = vtmp(2) - vtmp(1);
v1 = vtmp(1) - 2 * dv;
v2 = vtmp(end) + 2 * dv;

stmp = sconv(ix, :);
stmp = stmp / max(stmp(:));

% select an SRF from the set
vcen = vtmp(3);

vg = v1 : 0.02 : v2;
z = sup_gauss(vg, vcen, vwid);

figure(1); clf
subplot(2,1,1)
plot(sfreq, stmp, vg, z, '-.', 'linewidth', 2)
axis([v1, v2, 0, 1.2])
title('sample AIRS SRFs')
ylabel('weight')
grid on; zoom on

% SRF set 2
ix = 2201 : 2206;
vwid = 1.7;

vtmp = ofreq(ix)
dv = vtmp(2) - vtmp(1);
v1 = vtmp(1) - 2 * dv;
v2 = vtmp(end) + 2 * dv;

stmp = sconv(ix, :);
stmp = stmp / max(stmp(:));

% select an SRF from the set
vcen = vtmp(3);

vg = v1 : 0.02 : v2;
z = sup_gauss(vg, vcen, vwid);

subplot(2,1,2)
plot(sfreq, stmp, vg, z, '-.', 'linewidth', 2)
axis([v1, v2, 0, 1.2])
xlabel('wavenumber')
ylabel('weight')
grid on; zoom on
% saveas(gcf, 'AIRS_sample_SRFs', 'png')

%-------------------------------
% AIRS FWHM and channel spacing
%-------------------------------

% linearized AIRS L1c channel spacing
x = cfreq(2:end);
y = diff(cfreq);
x1 =  700;   y1 = 0.24;
x2 = 2183;   y2 = 0.88;
z1 = ((x - x1) ./ (x2 - x1)) .* (y2 - y1) + y1;
z2 = x / 2500;

figure(2); clf
subplot(2,1,1)
plot(x, y, x, z1, x, z2)
axis([650, 2800, 0.2, 1.2])   
title('AIRS L1c channel spacing')
legend('L1c', '2 pt fit', 'v/2500', 'location', 'southeast')
ylabel('spacing')
grid on

% approximately linearized AIRS L1c channel widths
x1 =  700;   y1 = 0.5;
x2 = 2183;   y2 = 1.7;
w1 = ((x - x1) ./ (x2 - x1)) .* (y2 - y1) + y1;
w2 = x / 1250;

subplot(2,1,2)
plot(x, w1, x, w2)
title('AIRS L1c FWHM')
axis([650, 2800, 0.2, 2.5])   
legend('2 pt fit', '2*v/2500', 'location', 'southeast')
xlabel('wavenumber')
ylabel('width')
grid on
% saveas(gcf, 'AIRS_chan_spacing', 'png')

