
% load an AIRS granule
apath = '/asl/data/airs/L1C/2018/091';
agran = 'AIRS.2018.04.01.226.L1C.AIRS_Rad.v6.1.2.0.G18092114905.hdf';
afile = fullfile(apath, agran);
arad = hdfread(afile, 'radiances');
arad = permute(arad, [3,2,1]);
arad = squeeze(arad(:, 45, :));
afrq = load('freq2645.txt');
save data/airs_demo_rad arad afrq
