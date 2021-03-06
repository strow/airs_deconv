%
% L1d_regr1 - regression corrrection for L1c to L1d decon/recon
%
% uses data from conv_loop2
%
% radiance data (from load)
%   v_L1c, v_L1d, res
%   trueCrad, true AIRS L1c, rc = C1 * rkc
%   trueDrad, true AIRS L1d, rd = D1 * rkc
%   CtoDrad, CtoD AIRS L1d, rd = B1 * inv(A1) * rc
%
% key local variables:
% tCri, tCrd - true L1c ind and dep set radiance
% tDri, tCrd - true L1d ind and dep set radiance
% CDri, CDrd - L1c to L1d ind and dep set radiance
%
% tCbi, tCbd - true L1c ind and dep set brightness temps
% tDbi, tCbd - true L1d ind and dep set brightness temps
% CDbi, CDbd - L1c to L1d ind and dep set brightness temps
%
% derived from a2cris_regr1 and L1d_test3
%

% set paths to libs
addpath ../source
addpath /asl/packages/ccast/source
addpath /home/motteler/matlab/export_fig

% set resolving power (700 or 1200)
R = 700;

%-------------------
% get radiance data
%-------------------

switch (R)
  case 1200, d1 = load('L1d1200_cldy'); d2 = load('L1d1200_fit49');
  case 700,  d1 = load('L1d700_cldy');  d2 = load('L1d700_fit49');
end
CtoDrad  = [d1.CtoDrad, d2.CtoDrad];
trueCrad = [d1.trueCrad, d2.trueCrad];
trueDrad = [d1.trueDrad, d2.trueDrad];
v_L1c = d1.v_L1c;
v_L1d = d1.v_L1d;
res = d1.res;
clear d1 d2

[nchan, nobs] = size(trueDrad);

%-------------
% basic stats
%-------------

% split into dependent and independent sets
% rstate = rng; save rstate rstate
% load rstate
% rng(rstate)
% ixi = logical(randi([0,1], [1, nobs]));
% ixd = ~ixi;
% nind = sum(ixi); 
% ndep = sum(ixd);

ndep = 7377;
nind = 49;
ixd = 1:ndep; ixi = ndep+1 : ndep+nind;

tCri = trueCrad(:, ixi);  tCrd = trueCrad(:, ixd);
tDri = trueDrad(:, ixi);  tDrd = trueDrad(:, ixd);
CDri =  CtoDrad(:, ixi);  CDrd =  CtoDrad(:, ixd);
clear trueCrad trueDrad CtoCrad

% get brightness temps
tCbi = real(rad2bt(v_L1c, tCri)); tCbd = real(rad2bt(v_L1c, tCrd)); 
tDbi = real(rad2bt(v_L1d, tDri)); tDbd = real(rad2bt(v_L1d, tDrd));
CDbi = real(rad2bt(v_L1d, CDri)); CDbd = real(rad2bt(v_L1d, CDrd)); 

% mean and std of C to D minus true D
mdifCDbd = mean(CDbd - tDbd, 2);
mdifCDbi = mean(CDbi - tDbi, 2);
sdifCDbd = std(CDbd - tDbd, 0, 2);
sdifCDbi = std(CDbi - tDbi, 0, 2);

%-----------------------------------------
% bias, linear, and quadratic corrections
%-----------------------------------------

% corrected dependent and independent sets
cor1CDbd = CDbd - (mdifCDbd * ones(1, ndep));
cor1CDbi = CDbi - (mdifCDbd * ones(1, nind));

% residuals of bias correction (dep is zero)
mcor1CDbd = mean(cor1CDbd - tDbd, 2);
mcor1CDbi = mean(cor1CDbi - tDbi, 2);
scor1CDbi = std(cor1CDbi - tDbi, 0, 2);

% per-channel polynomial fits
cor2CDbd = zeros(nchan, ndep);
cor2CDbi = zeros(nchan, nind);
cor3CDbd = zeros(nchan, ndep);
cor3CDbi = zeros(nchan, nind);
Pcor2 = zeros(nchan, 2);
Pcor3 = zeros(nchan, 3);

for i = 1 : nchan
  Ptmp = polyfit(CDbd(i, :), tDbd(i, :), 1);
  cor2CDbd(i, :) = polyval(Ptmp, CDbd(i, :));
  cor2CDbi(i, :) = polyval(Ptmp, CDbi(i, :));
  Pcor2(i, :) = Ptmp;

  Ptmp = polyfit(CDbd(i, :), tDbd(i, :), 2);
  cor3CDbd(i, :) = polyval(Ptmp, CDbd(i, :));
  cor3CDbi(i, :) = polyval(Ptmp, CDbi(i, :));
  Pcor3(i, :) = Ptmp;
end

% residuals of linear correction (dep is zero)
mcor2CDbd = mean(cor2CDbd - tDbd, 2);
mcor2CDbi = mean(cor2CDbi - tDbi, 2);
scor2CDbi = std(cor2CDbi - tDbi, 0, 2);

% residuals of quadratic correction (dep is zero)
mcor3CDbd = mean(cor3CDbd - tDbd, 2);
mcor3CDbi = mean(cor3CDbi - tDbi, 2);
scor3CDbi = std(cor3CDbi - tDbi, 0, 2);

figure(1); clf
subplot(2,1,1)
[x1, y1] = pen_lift(v_L1d, mdifCDbi);
[x2, y2] = pen_lift(v_L1d, mcor2CDbi);
plot(x1, y1, x2, y2);
switch(R)
 case 700,
   axis([650, 2650, -0.2, 0.2]); % 700 res, fit49 ind
   title('AIRS-to-L1d minus true L1d mean, R = 700');
 case 1200
   axis([650, 2650, -1, 1]);     % 1200 res, fit49 ind
   title('AIRS-to-L1d minus true L1d mean, R = 1200');
end
legend('no correction', 'linear correction', 'location', 'north')
ylabel('\Delta BT (K)')
grid on

subplot(2,1,2)
[x1, y1] = pen_lift(v_L1d, sdifCDbi);
[x2, y2] = pen_lift(v_L1d, scor2CDbi);
plot(x1, y1, x2, y2);
switch (R)
  case 700,
    axis([650, 2650, 0, 0.1]);  % 700 res, fit49 ind
    title('AIRS-to-L1d minus true L1d std dev, R = 700');
  case 1200,
    axis([650, 2650, 0, 0.3]);   % 1200 res, fit49 ind
    title('AIRS-to-L1d minus true L1d std dev, R = 1200');
end
legend('no correction', 'linear correction', 'location', 'north')
ylabel('\Delta BT (K)')
xlabel('wavenumber (cm^{-1})')
grid on
switch (R)
  case 700, saveas(gcf, 'L1d_cor1_700', 'fig');
  case 1200, saveas(gcf, 'L1d_cor1_1200', 'fig');
end

return

figure(2); clf
subplot(2,1,1)
plot(v_L1d, mcor1CDbi, v_L1d, mcor2CDbi, v_L1d, mcor3CDbi);
% axis([650, 2650, -0.1, 0.1])  %  700 res, fit49 ind
  axis([650, 2650, -0.3, 0.3])  % 1200 res, fit49 ind
title('mean residual corrected independent set')
legend('bias correction', 'linear correction', 'quadratic correction', ...
       'location', 'north')
ylabel('\Delta BT (K)')
grid on

subplot(2,1,2)
plot(v_L1d, scor1CDbi, v_L1d, scor2CDbi, v_L1d, scor3CDbi);
% axis([650, 2650, 0, 0.2])   %  700 res, fit49 ind
  axis([650, 2650, 0, 0.6])   % 1200 res, fit49 ind
legend('bias correction', 'linear correction', 'quadratic correction', ...
       'location', 'north')
title('std residual corrected independent set')
ylabel('\Delta BT (K)')
xlabel('wavenumber (cm^{-1})')
grid on
% saveas(gcf, 'L1d_cor2_700', 'fig')
% saveas(gcf, 'L1d_cor2_1200', 'fig')

return

%-------------------------
% direct regression tests
%-------------------------

% direct regression, L1c to L1d radiances
Creg1D = (tCrd' \ tDrd')';
reg1Drd = Creg1D * tCrd;
reg1Dri = Creg1D * tCri;
reg1Dbd = real(rad2bt(v_L1d, reg1Drd));
reg1Dbi = real(rad2bt(v_L1d, reg1Dri));

% direct regression, L1c to L1d brightness temps
Creg2D = (tCbd' \ tDbd')';
reg2Dbd = Creg2D * tCbd;
reg2Dbi = Creg2D * tCbi;

% direct regression residuals
mreg1Dbd = mean(reg1Dbd - tDbd, 2);
mreg1Dbi = mean(reg1Dbi - tDbi, 2);
mreg2Dbd = mean(reg2Dbd - tDbd, 2);
mreg2Dbi = mean(reg2Dbi - tDbi, 2);

sreg1Dbd = std(reg1Dbd - tDbd, 0, 2);
sreg1Dbi = std(reg1Dbi - tDbi, 0, 2);
sreg2Dbd = std(reg2Dbd - tDbd, 0, 2);
sreg2Dbi = std(reg2Dbi - tDbi, 0, 2);

figure(3); clf
subplot(2,1,1)
% plot(v_L1d, mreg2Dbi, v_L1d, mreg1Dbi)
  plot(v_L1d, mreg1Dbi)
% axis([650, 2650, -0.02, 0.02])   % 700 res, fit49 ind
  axis([650, 2650, -0.15, 0.15])   % 1200 res, fit49 ind
title('mean regression residual independent set')
% legend('Tb regression', 'rad regression')
  legend('radiance regression')
ylabel('\Delta BT (K)')
grid on

subplot(2,1,2)
plot(v_L1d, sreg1Dbi)
% axis([650, 2650, 0, 0.01])   % 700 res, fit49 ind
  axis([650, 2650, 0, 0.04])   % 1200 res, fit49 ind
title('std regression residual independent set')
  legend('radiance regression')
ylabel('\Delta BT (K)')
xlabel('wavenumber (cm^{-1})')
grid on
% saveas(gcf, 'L1d_regr_700', 'fig')
% saveas(gcf, 'L1d_regr_1200', 'fig')
  
figure(4)
subplot(2,1,1)
plot(tcfrq, Pcor2(:, 1))
% axis([650, 1100, 0.995, 1.005])
title('"a" (scaling) weights')
grid on

subplot(2,1,2)
plot(tcfrq, Pcor2(:, 2))
% axis([650, 1100, -1, 1])
title('"b" (bias) weights')
xlabel('wavenumber (cm^{-1})')
grid on
% saveas(gcf, 'cor_ind_2', 'png')

