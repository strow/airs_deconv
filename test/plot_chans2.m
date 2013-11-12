%
% specify an AIRS channel set and show channels step sizes before
% and after trimming neighbots that are too close
%

% channel lower bound as a function a*x + b of frequency
a = 4e-4;
b = -0.04;

% old 1B channel set
% f1 = load('freq1b');
% f1 = f1.freq1b;

% old 1C channel set
% f1 = load('freq1c');
% f1 = f1.freq1c;

% new 1C channel set
f1 = load('freq2645.txt');

f2 = sort(f1);
n2 = length(f2);
df2 = diff(f2);

f3 = trim_chans(f2);
n3 = length(f3);
df3 = diff(f3);

figure(1); clf
plot(f2(1:n2-1), df2, f2, a*f2+b);
axis([600, 2700, 0, 1.2])
xlabel('wavenumber')
ylabel('step size')
title('AIRS channel steps before trimming')
legend('channel steps', 'lower bound', 'location', 'southeast')
grid; zoom on

figure(2); clf
plot(f3(1:n3-1), df3, f3, a*f3+b);
axis([600, 2700, 0, 1.2])
xlabel('wavenumber')
ylabel('step size')
title('AIRS channel steps after trimming')
legend('channel steps', 'lower bound', 'location', 'southeast')
grid; zoom on






