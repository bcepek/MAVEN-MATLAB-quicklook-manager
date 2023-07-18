function plot_ion_velocity_MSO(epoch, products_O, products_H, products_O2, windowpos)
h = figure();
windowpos(4) = floor(windowpos(4)/2);
set(h, 'position', windowpos);

verify = products_O.epoch>=epoch(1) & products_O.epoch<=epoch(end);
v_p = products_H.v_mso(verify, :);
v_O = products_O.v_mso(verify, :);
v_O2 = products_O2.v_mso(verify, :);
ymin = min(min([v_p; v_O; v_O2]));
ymax = max(max([v_p; v_O; v_O2]));

h = subplot(3,1,1);
plot(epoch, v_p(:,1), epoch, v_p(:,2), epoch, v_p(:,3),...
    'linewidth', 1.5)
datetick('x')
ylabel('H^+ v_{MSO}, km/s')
legend('v_x', 'v_y', 'v_z')
set(h, 'xticklabel', [], 'fontsize', 12)
grid on
xlim([epoch(1) epoch(end)])
%ylim([ymin ymax])

h = subplot(3,1,2);
plot(epoch, v_O(:,1), epoch, v_O(:,2), epoch, v_O(:,3),...
    'linewidth', 1.5)
datetick('x')
ylabel('O^+ v_{MSO}, km/s')
legend('v_x', 'v_y', 'v_z')
set(h, 'xticklabel', [], 'fontsize', 12)
grid on
xlim([epoch(1) epoch(end)])
ylim([ymin ymax])

h = subplot(3,1,3);
plot(epoch, v_O2(:,1), epoch, v_O2(:,2), epoch, v_O2(:,3),...
    'linewidth', 1.5)
datetick('x')
ylabel('O_2^+ v_{MSO}, km/s')
legend('v_x', 'v_y', 'v_z')
set(h, 'fontsize', 12)
grid on
xlim([epoch(1) epoch(end)])
ylim([ymin ymax])

end