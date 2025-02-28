function plot_STA_orientation(epoch_d1_disp, pos_sc_mso_d1_disp, quat_mso_d1_disp)
cur_sec = hour(epoch_d1_disp)*3600+...
    minute(epoch_d1_disp)*60+...
    second(epoch_d1_disp);
dx = [diff(pos_sc_mso_d1_disp(:, 1)),...
    diff(pos_sc_mso_d1_disp(:, 2)),...
    diff(pos_sc_mso_d1_disp(:, 3))];
dt = diff(cur_sec);
sc_vel = dx./dt;
sc_vel(end+1,:) = sc_vel(end,:);
sc_vel = sc_vel./sqrt(sum(sc_vel.^2,2));

l = size(epoch_d1_disp, 1);
Z_sta = [zeros(l, 1), zeros(l, 1), ones(l, 1)];
Z_sta = quatrotate(quatinv(quat_mso_d1_disp), Z_sta);

sundir = [ones(l, 1), zeros(l, 1), zeros(l, 1)];
normdir = [pos_sc_mso_d1_disp(:, 1),...
    pos_sc_mso_d1_disp(:, 2),...
    pos_sc_mso_d1_disp(:, 3)];
normdir = normdir./sqrt(sum(normdir.^2,2));

angle_vel = acos(sum(Z_sta.*sc_vel,2))*180/pi;
angle_sun = acos(sum(Z_sta.*sundir,2))*180/pi;
angle_norm = acos(sum(Z_sta.*normdir,2))*180/pi;

h = figure();
pos = get(h, 'position');
pos(1) = 0;
pos(2) = 250;
pos(3) = 1100;
pos(4) = pos(4)*1.3;
set(h, 'position', pos)

h = subplot(3,1,1);
plot(epoch_d1_disp, angle_vel, 'linewidth', 1.5)
datetick('x', 'MM:SS')
ylim([0 180])
ylabel('angle with velocity')
pos = get(h, 'position');
set(h, 'position', [0.4*pos(1) pos(2) 1.2*pos(3) 1.2*pos(4)],...
    'xticklabel', [], 'fontsize', 11,...
    'ytick', [0 45 90 135 180], 'ylim', [0 180])
xlim([epoch_d1_disp(1) epoch_d1_disp(end)])
grid on

h = subplot(3,1,2);
plot(epoch_d1_disp, angle_sun, 'linewidth', 1.5)
datetick('x', 'MM:SS')
ylim([0 180])
ylabel('angle with sun direction')
pos = get(h, 'position');
set(h, 'position', [0.4*pos(1) pos(2) 1.2*pos(3) 1.2*pos(4)],...
    'xticklabel', [], 'fontsize', 11,...
    'ytick', [0 45 90 135 180], 'ylim', [0 180])
xlim([epoch_d1_disp(1) epoch_d1_disp(end)])
grid on

h = subplot(3,1,3);
plot(epoch_d1_disp, angle_norm, 'linewidth', 1.5)
datetick('x', 'HH:MM:SS')
ylim([0 180])
ylabel('angle with surface normal')
pos = get(h, 'position');
set(h, 'position', [0.4*pos(1) pos(2) 1.2*pos(3) 1.2*pos(4)], 'fontsize', 11,...
     'ytick', [0 45 90 135 180], 'ylim', [0 180])
xlim([epoch_d1_disp(1) epoch_d1_disp(end)])
grid on

end