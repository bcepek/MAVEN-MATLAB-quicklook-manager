function plot_dayside_orbit_MSE(x1, products_H)

epoch_d1 = products_H.epoch;
x = products_H.pos_sc_mso(:, 1);
y = products_H.pos_sc_mso(:, 2);
z = products_H.pos_sc_mso(:, 3);

verify = products_H.epoch >= x1(1) & products_H.epoch <= x1(2);
ind = find(abs(products_H.epoch - x1(3)) == min(abs(products_H.epoch - x1(3))));
B_avg = mean(products_H.magf(verify, :));

x = x(ind);
y = y(ind);
z = z(ind);

E = -cross([-1, 0, 0], B_avg);
E = E/sqrt(sum(E.^2));
rotmax = [[1;0;0], cross([1;0;0],E'), E'];
XYZ_mse = inv(rotmax)*[x;y;z];
y = XYZ_mse(2);
z = XYZ_mse(3);

r = sqrt(x.^2+y.^2+z.^2);
SZA = acos(x./r)*180/pi;
y_sza = y*SZA/sqrt(y^2+z^2);
z_sza = z*SZA/sqrt(y^2+z^2);

circle = zeros(2000,2);
circle(1:1000, 1) = linspace(-1,1,1000);
circle(1001:2000, 1) = linspace(1,-1,1000);
circle(1:1000, 2) = sqrt(1 - circle(1:1000, 1).^2);
circle(1001:2000, 2) = -sqrt(1 - circle(1001:2000, 1).^2);

szas = (10:10:90)/100;

figure()
plot(circle(:,1), circle(:,2), 'linewidth', 2, 'color', 'red')
hold on
for i = 1:length(szas)
    cir = circle*szas(i);
    plot(cir(:,1), cir(:,2), 'linewidth', 1, 'color', [0.8 0.8 0.8])
    text(szas(i)*cos(pi/4), -szas(i)*cos(pi/4), [num2str(10*i) '°'])
    text(szas(i)*cos(pi/4), szas(i)*cos(pi/4), [num2str(10*i) '°'])
    text(-szas(i)*cos(pi/4), -szas(i)*cos(pi/4), [num2str(10*i) '°'])
    text(-szas(i)*cos(pi/4), szas(i)*cos(pi/4), [num2str(10*i) '°'])
end
axis equal
plot(y_sza/100, z_sza/100, 'marker', '.', 'markersize', 20, 'color', 'black')
xlabel('Y_{MSE}, R_M')
ylabel('Z_{MSE}, R_M')
xlim([-1 1])
ylim([-1 1])
end