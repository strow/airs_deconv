%
% basic JPL shift
%

function Tb_resamp = jpl_shift(Tb_in, v_in, v_nom, ix)

d1 = load('L1C.airs_resample.v1.0.0.anc');
a = d1(ix, 1);
b = d1(ix, 2);

dv = v_nom - v_in;

Tb_spline = interp1(v_in, Tb_in, v_nom, 'spline');

Tb_resamp = Tb_in + (a .* (Tb_spline - Tb_in) ./ dv + b) .* dv;

