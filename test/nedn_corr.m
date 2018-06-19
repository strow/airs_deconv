%
% nedn_corr - NEdN correlation test, run after nedn_test1
%

% d_airs = r_airs - r_280K;
% xx = corrcoef(d_airs');

b_airs = rad2bt(freq_airs, r_airs);
b_tran = rad2bt(freq_tran, r_tran);

d_airs = b_airs - mean(b_airs, 2);
d_tran = b_tran - mean(b_tran, 2);

acorr = corrcoef(d_airs');
tcorr = corrcoef(d_tran');

load llsmap5;
cax = [-0.4, 0.4];
% v1 = 650; v2 = 2550;
v1 = 1000; v2 = 1040;
vax = [v1, v2, v1, v2];

% plot uncorrelated AIRS data
figure(1)
imagesc(freq_airs, freq_airs, acorr)
set(gca,'YDir','normal')
axis(vax)
colormap(llsmap5)
caxis(cax)
colorbar
title('AIRS noise correlation')
xlabel('wavenumber')
ylabel('wavenumber')

% plot translation correlation
figure(2)
imagesc(freq_tran, freq_tran, tcorr)
set(gca,'YDir','normal')
axis(vax)
colormap(llsmap5)
caxis(cax)
colorbar
title('translation noise correlation')
xlabel('wavenumber')
ylabel('wavenumber')

