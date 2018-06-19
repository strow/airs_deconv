
b = 0;
w = 2;
x1 = -2.5*w:w/40:2.5*w;

y1 = sup_gauss(x1, b, w, 1);
y2 = sup_gauss(x1, b, w, 1.2);
y3 = sup_gauss(x1, b, w, 1.4);
y4 = sup_gauss(x1, b, w, 1.6);

plot(x1, y1, x1, y2, x1, y3, x1, y4)
axis([x1(1), x1(end), -0.1, 1.1])
title('standard and higher order Gaussian functions')
legend('p = 1', 'p = 1.2', 'p = 1.4','p = 1.6')
grid on; zoom on

return

y1 = sup_gauss(x1, b, w, 1.5);
y2 = bell_fun(x1, b, w, 1.8);

plot(x1, y1, x1, y2, 'linewidth', 2)
axis([x1(1), x1(end), -0.1, 1.1])
legend('gauss p = 1.5', 'bell p = xx')
grid on; zoom on

