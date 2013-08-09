
% chanstep - set a lower bound on AIRS channel step size 
%
% the step size is a linear funcion of frequency a*f + b, derived
% from tests in plot_chans1.m.

function s = chanstep(f)

a = 4e-4;
b = -0.04;

s = a.*f + b;

