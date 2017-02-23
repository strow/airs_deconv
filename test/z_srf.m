% function [yout]=z_srf(v,v0,fwhm,gfrac,gslope,lexpval);
%
% Experiment with SRF's using emodified version of George Aumann's
% Gauss+Lorentz SRF.  This version does the guassian portion with the
% power chaning linearly with x.
%
% Input:
%    v = frequency array
%    v0 = channel center freq
%    fwhm = channel FWHM
%    gfrac=gauss fraction (lorentz fraction is 1-gfrac)
%    gslope=guass linear slope
%    lexpval=lorentz exponent value
%
% Output:
%    yout = SRF
%
% The equation is:
% x=abs(v - v0)/(0.5*fwhm);
% yout=gfrac*(exp(-a*x.^(2+x)) + (1-gfrac)*(1./(1 + x.^lexpval));
% with nominal values:
%    gfrac=0.95;
%    gexpval=0.5;
%    lexpval=1.8;
%

function [yout]=x_srf(v,v0,fwhm,gfrac,gslope,lexpval);

x=abs(v - v0)/(0.5*fwhm);
% Note that x is twice the distance from chan center in units of FWHM
a=log(2);
yout=gfrac*(exp(-a*x.^(2 + gslope*x))) + (1-gfrac)*(1./(1 + x.^lexpval));
