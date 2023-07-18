function plot_altit_dependence (X_grid, Y_grid, Z, SZA_min, SZA_max, xmin, xmax, ymin, ymax, text_x, text_y)

test_alt = Y_grid(:,1);
test1 = mean(reshape(Z( X_grid >= -SZA_max &  X_grid <= -SZA_min), 100, []), 2, 'omitnan');
test2 = mean(reshape(Z( X_grid >=  SZA_min &  X_grid <=  SZA_max), 100, []), 2, 'omitnan');
plot(test_alt, test1, 'r', test_alt, test2, 'k')
legend('Southern hemisphere', 'Northern hemisphere','Location','northwest')
xlim([xmin xmax])
ylim([ymin ymax])
% yrange = ylim;
% text_y = yrange(2) - yrange(2)/4;
text(text_x, text_y, ['SZA = ', num2str(SZA_min), ' - ', num2str(SZA_max)], 'FontSize', 14)
%xlabel('Altitude, km')
set(gca, 'FontSize', 16)
 set(gca, 'XTick', [200 600 1000 1400 1800])