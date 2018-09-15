% 
% cris_test8 -- plot sample airs decon basis functio
%

load bconv.mat
figure(1); clf
set(gcf, 'Units','centimeters', 'Position', [4, 10, 24, 16])
plot(bfrq, binv(:, 1000)) 
axis([972, 988, -2.5, 3.5])
title('sample AIRS decon basis function')
xlabel('wavenumber')
ylabel('weight')
grid on; zoom on
% saveas(gcf, 'airs_decon_basis', 'png')
  export_fig('airs_decon_basis.pdf', '-m2', '-transparent')

