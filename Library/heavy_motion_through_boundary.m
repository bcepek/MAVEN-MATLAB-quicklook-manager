function f1 = heavy_motion_through_boundary(products_H, products_O, products_O2, epoch_d1_disp)

verify = products_H.epoch>epoch_d1_disp(1) & products_H.epoch<epoch_d1_disp(end);
epoch = products_H.epoch(verify);
vO = products_O.v_mso(verify,:);
vO2 = products_O2.v_mso(verify,:);

%normal = repmat(normal, length(epoch), 1);
normal = products_O.pos_sc_mso(verify,:);
normal = normal./repmat(sqrt(sum(normal.^2, 2)),1,3);
angleO = acos( sum(vO.*normal, 2)./sqrt(sum(vO.^2, 2)) )*180/pi;
angleO2 = acos( sum(vO2.*normal, 2)./sqrt(sum(vO2.^2, 2)) )*180/pi;

f1 = figure(3);
plot(epoch, angleO, 'color', 'red')
hold on
plot(epoch, angleO2, 'color', 'blue')
hold off
datetick('x')
legend('O^+', 'O_2^+')
grid on
xlabel('UT, HH:MM')
ylabel({'angle between ion bulk velocity'; 'and magnetopause, degrees'})
set(f1, 'pos', [-960 316 847 390])
end