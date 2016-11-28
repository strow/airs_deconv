
k = 240;  % this is the "chunk" size

for i = 1 : length(flist)   % loop on netCDF files

  % read some CrIS data
  tfile = fullfile('tsno', flist(i).name);
  rad =  h5read(tfile, '/L1bCrIS/rad');
  [m, nobs] = size(rad);

  for j = 1 : k : nobs     % loop on chunks

    ix = j : min(j + k - 1, nobs);   % chunk index

    rLW = rad(iLW, ix);
    rMW = rad(iMW, ix);
    rSW = rad(iSW, ix);

    [rtmp, afrq] = ...
         cris2airs(rLW, rMW, rSW, vLW, vMW, vSW, sfile, cfrq, opt1);

    if j == 1
      arad = zeros(length(afrq), nobs);
    end

    arad(:, ix) = rtmp;

    % save the CrIS to AIRS translation
    [sp, sn, se] = fileparts(flist(i).name);
    mfile = fullfile('c2airs', [sn, '.mat']);
    save(mfile, 'arad', 'afrq', '-v7.3')

  end
end

% The above is a simple example of calling cris2airs in chunks of 240
% obs.  cris2airs is vectorized and will run much faster if you pass
% it multiple obs, but runs into memory limits after a few hunderd
% obs.  Note that line 12 is the key, this sets the current chunk as a
% matlab index.  This sort of chunking is only one or two lines longer
% than calling airs2cris in a simple one-at-a time loop.  Lines 21--23
% initialize the output array the first time through the loop, since
% we don't know the size of arad ahead of time.

