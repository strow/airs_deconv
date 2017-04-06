%
% band_regr - banded regression, solve R * B = C for banded R
%
% SYNOPSIS
%   R = band_regr(B, C, vB, vC, s)
% 
% INPUTS
%   B  -  m chan x k obs
%   C  -  n chan x k obs
%   vB -  chan freq for B
%   vC -  chan freq for C
%   s  -  band span is 2*s + 1
%
% OUTPUTS
%   R = n x m, banded
%
% DISCUSSION
%   matlab unbanded soln is R' = B' \ C'
%

function R = band_regr(B, C, vB, vC, s)

% B = a1Crd(iLW,:);
% C = cLWrd;
% vB = va1C(iLW);
% vC = vcLW;
% s = 60;

% try with brightness temp
% B = real(rad2bt(vB, B));
% C = real(rad2bt(vC, C));

[m, k] = size(B);
[n, t] = size(C);
if k ~= t, error('B and C must have the same number of columns'), end

R = zeros(n, m);

% loop on output channels
tic
for i = 1 : n

  % find closest index in vB for vC(i)
  j0 = interp1(vB, 1:length(vB), vC(i), 'nearest');

  % [vC(i), vB(j0)]

  % input channel span
  j1 = max(1, j0 - s);
  j2 = min(m, j0 + s);
  jx = j1 : j2;

% [j1, j2]

  R(i, jx) = B(jx, :)' \ C(i, :)';

end
toc

% return

% R2 = (B' \ C')';
% mean(abs(R(:) - R2(:)))

% figure(1); clf
% pcolor(vB, vC, R); 
% shading flat
% caxis([-0.5, 0.5])
% load llsmap5
% colormap(llsmap5)
% colorbar
% title('AIRS to CrIS LW regression matrix')
% xlabel('AIRS channels')
% ylabel('CrIS channels')
% grid on
% 
% figure(2); clf
% ix = 400;  vC(ix)
% plot(vB, R(ix, :))

