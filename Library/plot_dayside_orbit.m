function plot_dayside_orbit(epoch_d1, x, y, z)

verify = x>=0;
x = x(verify);
y = y(verify);
z = z(verify);

r = sqrt(x.^2+y.^2+z.^2);
%x_s = x./r;
y_s = y./r;
z_s = z./r;

SZA = acos(x./r)*180/pi;

height = 3390*(r-1);

circle = zeros(2000,2);
circle(1:1000, 1) = linspace(-1,1,1000);
circle(1001:2000, 1) = linspace(1,-1,1000);
circle(1:1000, 2) = sqrt(1 - circle(1:1000, 1).^2);
circle(1001:2000, 2) = -sqrt(1 - circle(1001:2000, 1).^2);

szas = (10:10:80)/90;

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
scatter(y_s.*SZA/90, z_s.*SZA/90, 25, height, 'filled')
xlabel('Y_{MSO}, R_M')
ylabel('Z_{MSO}, R_M')
cbar_handle = colorbar;
ylabel(cbar_handle, 'Height, km')
colormap hsv
%disp([y_s(1), z_s(1)])
%text(y_s(1), z_s(1), datestr(epoch_d1(1), 'HH:MM:SS'))
xlim([-1 1])
ylim([-1 1])

interesting_altit = get(cbar_handle, 'ytick');

[Xtest, Ytest] = meshgrid(height, interesting_altit);
[~, I] = min(abs(Ytest - Xtest), [],2);

plot(y_s(I).*SZA(I)/90, z_s(I).*SZA(I)/90, 'color', 'black', 'marker', '*', 'linestyle', 'none', 'markersize', 10)
end