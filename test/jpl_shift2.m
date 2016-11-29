%
% JPL shift for 1c data with new "a" and "b" weights
%

function Tb_resamp = jpl_shift(Tb_in, v_in, v_nom)

d1 = load('jpl_shift2');

if ~(isclose(v_in, d1.frq1) && isclose(v_nom, d1.frq2))
  error('frequency mismatch')
end

dv = v_nom - v_in;

Tb_spline = interp1(v_in, Tb_in, v_nom, 'spline');

Tb_resamp = Tb_in + (d1.a .* (Tb_spline - Tb_in) ./ dv + d1.b) .* dv;

