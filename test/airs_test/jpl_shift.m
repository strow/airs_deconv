%
% basic JPL shift function from the L1C ATBD
%

function Tb_resamp = jpl_shift(Tb_in, v_in, v_nom, coeff, ix)

if nargin < 5
  ix = 1 : length(v_in);
end

switch coeff
  case 'jpl_1b'
    d1 = load('L1C.airs_resample.v1.0.0.anc');
    a = d1(ix, 1); b = d1(ix, 2);
  case 'umbc_1b'
    d1 = load('umbc_shift_1b');
    a = d1.a(ix); b = d1.b(ix);
  case 'umbc_1c'
    d1 = load('umbc_shift_1c');
    a = d1.a(ix); b = d1.b(ix);
  otherwise
    error('bad value for coeff')
end

dv = v_nom - v_in;

Tb_spline = interp1(v_in, Tb_in, v_nom, 'spline');

Tb_resamp = Tb_in + (a .* (Tb_spline - Tb_in) ./ dv + b) .* dv;

% check vs umbc tabulation frequency
% if ~(isclose(v_in, d1.frq1(jx)) && isclose(v_nom, d1.frq2(jx)))
%   fprintf(1, 'frequency mismatch\n')
%   keyboard
%   error('frequency mismatch')
% end
