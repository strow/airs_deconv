%
% plot AIRS L1c SRFs, channel spacing, and resolving poser
%
% sample SRF plots are from fit_SRFs.m and include a higher
% order Gaussian fit for two of the SRFs
%

addpath ../source

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
vwid = 0.492;
p = 1.4;

vtmp = ofreq(ix)
dv = vtmp(2) - vtmp(1);
v1 = vtmp(1) - 2 * dv;
v2 = vtmp(end) + 2 * dv;

stmp = sconv(ix, :);
stmp = stmp / max(stmp(:));

% select an SRF from the set
vcen = vtmp(3);

vg = v1 : 0.02 : v2;
  z = sup_gauss(vg, vcen, vwid, p);
% z = bell_fun(vg, vcen, vwid);

figure(1); clf
subplot(2,1,1)
% semilogy(sfreq, stmp, vg, z, 'k--', 'linewidth', 2)
% axis([v1, v2, 5e-4, 1.2])
  plot(sfreq, stmp, vg, z, 'k--', 'linewidth', 2)
  axis([v1, v2, -0.1, 1.1])
title('sample AIRS SRFs')
ylabel('weight')
grid on; zoom on

% SRF set 2
ix = 2201 : 2206;
vwid = 1.792;

vtmp = ofreq(ix)
dv = vtmp(2) - vtmp(1);
v1 = vtmp(1) - 2 * dv;
v2 = vtmp(end) + 2 * dv;

stmp = sconv(ix, :);
stmp = stmp / max(stmp(:));

% select an SRF from the set
vcen = vtmp(3);

vg = v1 : 0.02 : v2;
  z = sup_gauss(vg, vcen, vwid, p);
% z = bell_fun(vg, vcen, vwid);

subplot(2,1,2)
  plot(sfreq, stmp, vg, z, 'k--', 'linewidth', 2)
  axis([v1, v2, -0.1, 1.1])
% semilogy(sfreq, stmp, vg, z, 'k--', 'linewidth', 2)
% axis([v1, v2, 5e-4, 1.2])
xlabel('wavenumber (cm^{-1})')
ylabel('weight')
grid on; zoom on
% saveas(gcf, 'airs_sample_SRF', 'fig')

%-------------------------------
% AIRS FWHM and channel spacing
%-------------------------------

% AIRS L1b channel spacing
% vtmp = srf_read(sfile);
% vtmp = sort(vtmp(1:2378));
% v_L1b = vtmp(2:end);
% dv_L1b = diff(vtmp);

% AIRS L1c channel spacing
v_L1c = cfreq(2:end);
dv_L1c = diff(cfreq);
fwhm_L1c = 2 * dv_L1c;
res_L1c = v_L1c ./ fwhm_L1c;

fwhm1200 = (v_L1c / 1200);
dv1200 = fwhm1200 / 2;

figure(2); clf
% set(gcf, 'Units','centimeters', 'Position', [4, 10, 24, 16])
subplot(2,1,1)
plot(v_L1c, dv_L1c, 'o', v_L1c, dv1200)
axis([650, 2800, 0.2, 1.2])   
title('AIRS L1c channel spacing')
legend('L1c', 'res 1200', 'location', 'southeast')
ylabel('wavenumber (cm^{-1})')
grid on

subplot(2,1,2)
plot(v_L1c, res_L1c, 'o')
axis([650, 2800, 1000, 1600])
title('AIRS L1c approximate resolving power')
xlabel('wavenumber (cm^{-1})')
ylabel('resolving power')
grid on
% saveas(gcf, 'airs_L1c_res', 'fig')

% show L1b channel spacing
% subplot(2,1,2)
% plot(v_L1b, dv_L1b, 'o')
% axis([650, 2800, 0, 1])
% title('AIRS L1b channel spacing')
% xlabel('wavenumber')
% grid on

