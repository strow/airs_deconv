%
% NAME
%   trim_chans - drop close neighbors from an AIRS channel set
%
% SYNOPSIS
%   [ctrim, itrim] = trim_chans(cfreq)
%
% INPUTS
%   cfreq  - AIRS channel frequencies
% 
% OUTPUT
%   ctrim  - trimmed sorted version of cfreq
%   itrim  - optional indices of ctrim in cfreq
%
% NOTES
%   derived from part of mkaconv2
%
% AUTHOR
%  H. Motteler, 22 June 2013
%

function [ctrim, itrim] = trim_chans(cfreq)

%  cfreq  - input AIRS channel set
%  fsrt   - cfreq in guaranteed sorted order
%  isrt   - indices for fsrt in cfreq

[fsrt, isrt] = sort(cfreq(:));
nchan = length(fsrt);

% set a lower bound (from plot_chans1.m) on channel step size 
% as a linear function a*x + b of frequency x
a = 4e-4;
b = -0.04;

% set up the channel trimming loop
ftmp = zeros(nchan,1);
jsrt = zeros(nchan,1);

ftmp(1) = fsrt(1);
jsrt(1) = isrt(1);
i = 1;
j = 2;

while 1

  while j <= nchan && fsrt(j) <= ftmp(i) + a*ftmp(i) + b
    j = j + 1;
  end

  % after the inner loop either j > nchan and we're done, 
  % or fsrt(j) > ftmp(i) + a*ftmp(i) + b, and we have a good 
  % next channel freq

  if j > nchan
    break
  end

  ftmp(i+1) = fsrt(j);
  jsrt(i+1) = isrt(j);

  i = i + 1;
  j = j + 1;
end

% switch to the trimmed grid
nchan = i;
ctrim = ftmp(1:nchan);   % final frequency grid
itrim = jsrt(1:nchan);   % index of ctrim in cfreq

