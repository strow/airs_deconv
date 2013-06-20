%
% make a fake AIRS SRF file from a real one
%

% use my hdf libs
addpath /home/motteler/mot2008/hdf/h4tools

% real SRF file
sfile = '/asl/matlab/srftest/srftables_m140f_withfake_mar08.hdf';

% fake SRF file
ffile = 'fake_srf1.hdf';

% read the srf data
[alist, fattr] = h4sdread(sfile);

for i = 1 : length(alist)
  switch alist{i}{1}
    case 'chanid', chanid = alist{i}{2};
    case 'freq',   freq   = alist{i}{2};
    case 'fwgrid', fwgrid = alist{i}{2};
    case 'srfval', srfval = alist{i}{2};
    case 'width',  width  = alist{i}{2};
  end
end
% clear alist fattr

[nspts, nchan] = size(srfval);

% tweak the SRF parameters
shift = 0.2;  % 1/cm shift for the entire channel grid
w1 = 0.00;    % randn scaling for individual channel centers
w2 = 0.05;    % randn scaling for individual channel widths

freq = freq + randn(1,nchan) .* w1 + shift;
width = width .* (1 + randn(1,nchan) .* w2);

% save the modified SRFs
slist = {{'chanid', chanid}, {'freq', freq}, {'fwgrid', fwgrid}, ...
         {'srfval', srfval}, {'width', width}};

h4sdwrite(ffile, slist)

% [slist1, fattr1] = h4sdread(ffile);
% for i = 1:5
%   isequal(slist{i}{2}, slist1{i}{2})
% end

