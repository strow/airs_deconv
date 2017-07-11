% Do the *whole* svd thing from scratch on even TIGR
% subset & use odd subset strictly for testing.
% USE original 40 level profiles with 64 level Tb's

load /carrot/scratch3/motteler/tigr/tigr40
load /carrot/scratch3/motteler/chansets/ind728
load /carrot/scratch3/motteler/tigr/bt64

bt = bt(ind728,:);

%% select training set
temp = [temp; gtemp];

[m,n] = size(bt);
train_set  = (rem(1:n, 2) == 0);
extrap_set = ~train_set;

IN = bt(:, train_set);
OUT = temp(:, train_set);

EXIN = bt(:, extrap_set);
EXOUT = temp(:, extrap_set);

clear bt
clear temp

[u,s,v] = svd(IN,0);
s = s(:,1:35);
B1 = u * s;
B1inv = pinv(B1);
bt35 = B1inv * IN;

[u,s,v] = svd(OUT,0);
s = s(:,1:35);
B2 = u * s;
B2inv = pinv(B2);
p35  = B2inv * OUT;

R = p35 / bt35;
R1 = B2 * R * B1inv;

save svdset B1 B1inv B2 B2inv IN OUT EXIN EXOUT p35 bt35 R R1

