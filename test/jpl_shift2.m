%
% JPL shift for 1c data with new "a" and "b" weights
%

function Tb_resamp = jpl_shift(Tb_in, v_in, v_nom, jx)

if nargin == 3
  jx = 1 : length(v_in);
end

d1 = load('jpl_shift2');

a = d1.a(jx);
b = d1.b(jx);

if ~(isclose(v_in, d1.frq1(jx)) && isclose(v_nom, d1.frq2(jx)))
  keyboard
  error('frequency mismatch')
end

dv = v_nom - v_in;

Tb_spline = interp1(v_in, Tb_in, v_nom, 'spline');

Tb_resamp = Tb_in + (a .* (Tb_spline - Tb_in) ./ dv + b) .* dv;

