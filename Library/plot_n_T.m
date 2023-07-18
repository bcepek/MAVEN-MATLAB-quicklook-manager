function plot_n_T (timefrom, timeto, products_O, products_H, products_O2)

% n_p = products_H.concentration;
% T_p = products_H.temp;
n_O = products_O.concentration;
T_O = products_O.temp;
% n_O2 = products_O2.concentration;
% T_O2 = products_O2.temp;

epoch = products_O.epoch;

v = epoch>=timefrom & epoch<=timeto;

figure()
% loglog(n_p(v), T_p(v), n_O(v), T_O(v), n_O2(v), T_O2(v),...
%     'marker', '.', 'markersize', 14, 'linestyle', 'none');
loglog(T_O(v), n_O(v), 'marker', '.', 'markersize', 14, 'linestyle', 'none');
%legend('p', 'O^+', 'O_2^+')
legend('O^+')
grid on
ylabel('Number density, cm^{-3}')
xlabel('Temperature, eV')
title([datestr(timefrom, 'dd mmm yyyy HH:MM:SS') '-' datestr(timeto, 'HH:MM:SS')])
% 
% xlim([1 1e4])
% ylim([1e-2 1e3])

end