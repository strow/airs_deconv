
  sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';
% sfile = '/home/chepplew/projects/airs/srf_study/airs_2834_srf_table_pert_t_gm.hdf';

chanid = double(hdfread(sfile, 'chanid'));
tfreq  = double(hdfread(sfile, 'freq'));
fwgrid = double(hdfread(sfile, 'fwgrid'));
srfval = double(hdfread(sfile, 'srfval'));
width  = double(hdfread(sfile, 'width'));

[nchan, nspts] = size(srfval);

% alternate 5 and 10 pct tweaks
tweak = repmat([1.05; 1.1], nchan/2, 1);

% untweaked SRFs
tgrid1 = (width * fwgrid) + tfreq * ones(1, nspts);

% tweaked SRFs
tgrid2 = ((width .* tweak) * fwgrid) + tfreq * ones(1, nspts);

% plot x-axis
jx = floor(nspts/2);
jx = jx-100: jx+100;

% channel index
ix = 300;
subplot(2,1,1)
plot(tgrid1(ix,jx), srfval(ix,jx), tgrid2(ix,jx), srfval(ix,jx))
legend('original', 'tweak')

ix = 301;
subplot(2,1,2)
plot(tgrid1(ix,jx), srfval(ix,jx), tgrid2(ix,jx), srfval(ix,jx))
legend('original', 'tweak')

