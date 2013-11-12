
% return the indices of the intersection of 2 regular grids
% grid values should be close but don't have to be identical
%
% example
%
%  s1 = [2 3 4 5 6 7]
%  s2 = [4 5 6 7 8 9]
%
% [i1, i2] = seq_isect(s1, s2)
%
%  i1 == [3 4 5 6]
%  i2 == [1 2 3 4]

function [i1, i2] = seq_isect(s1, s2)

v1 = max(s1(1), s2(1));
v2 = min(s1(end), s2(end));

t1 = interp1(s1, 1:length(s1), [v1, v2], 'nearest');
i1 = t1(1) : t1(2);

t2 = interp1(s2, 1:length(s2), [v1, v2], 'nearest');
i2 = t2(1) : t2(2);

