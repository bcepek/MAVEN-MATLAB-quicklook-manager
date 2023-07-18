function Box_draw(x1,products_H, products_O)

ind_start = find(abs(products_H.epoch - x1(3)) == min(abs(products_H.epoch - x1(3))));
ind_end = find(abs(products_H.epoch - x1(4)) == min(abs(products_H.epoch - x1(4))));
verify = ind_start:ind_end;
B_avg = mean(products_H.magf(verify, :));

x = products_H.pos_sc_mso(verify,1);
y = products_H.pos_sc_mso(verify,2);
z = products_H.pos_sc_mso(verify,3);

X = x/3390;
Y = y/3390;
Z = z/3390;

X_mso = X;
Y_mso = Y;
Z_mso = Z;

height = sqrt(x.^2 + y.^2 + z.^2);
x = x./height;
y = y./height;
z = z./height;

E = -cross([-1, 0, 0], B_avg);
E = E/sqrt(sum(E.^2));
rotmax = [[1;0;0], cross([1;0;0],E'), E'];

XYZ_mse = inv(rotmax)*[X';Y';Z'];
X = XYZ_mse(1,:);
Y = XYZ_mse(2,:);
Z = XYZ_mse(3,:);

dt = second(products_H.epoch(end) - products_H.epoch(end-1));
sc_v(1) = X(end) - X(end-1);
sc_v(2) = Y(end) - Y(end-1);
sc_v(3) = Z(end) - Z(end-1);
sc_v = 3390*sc_v/dt;
sc_v_norm = sc_v/(2*norm(sc_v));

E = inv(rotmax)*E';

magf_mse = [inv(rotmax)*products_H.magf(verify,:)']';
magf_norm = sqrt(sum(magf_mse.^2,2));
magf_mse = magf_mse/max(magf_norm);
v_p_mse = [inv(rotmax)*products_H.v_mso(verify,:)']';
v_p_norm = sqrt(sum(v_p_mse.^2,2));
v_p_mse = v_p_mse/max(v_p_norm);
v_O_mse = [inv(rotmax)*products_O.v_mso(verify,:)']';
v_O_norm = sqrt(sum(v_O_mse.^2,2));
v_O_mse = v_O_mse/max(v_O_norm);

f1_h = figure();
% Mars normal
Mars_norm = [X(end),Y(end),Z(end)]/norm([X(end),Y(end),Z(end)]);
h1 = quiver3(X(end),Y(end),Z(end), Mars_norm(1), Mars_norm(2), Mars_norm(3),...
    'color', 'black', 'linewidth', 2);
hold on
% electric field
h2 = quiver3(X(end),Y(end),Z(end), E(1),E(2),E(3), 'color', 'green', 'linewidth', 2);
% sun
h3 = quiver3(X(end),Y(end),Z(end), 1,0,0, 'color', 'red', 'linewidth', 2);
% s/c velocity
h4 = quiver3(X(end),Y(end),Z(end), sc_v_norm(1),sc_v_norm(2),sc_v_norm(3),...
    'color', 'blue', 'linewidth', 2);
% MF
h5 = quiver3(X',Y',Z',magf_mse(:,1),magf_mse(:,2),magf_mse(:,3), 0,...
    'color', [0.5,0.5,0], 'linewidth',1);
% v_p
h6 = quiver3(X',Y',Z',v_p_mse(:,1),v_p_mse(:,2),v_p_mse(:,3), 0,...
    'color', [0,0.8,0.8], 'linewidth',1);
h7 = quiver3(X',Y',Z',v_O_mse(:,1),v_O_mse(:,2),v_O_mse(:,3), 0,...
    'color', [0.8,0,0.8], 'linewidth',1);

xlabel('X_{MSE}')
ylabel('Y_{MSE}')
zlabel('Z_{MSE}')

[X_sp, Y_sp, Z_sp] = sphere(50);
surf(X_sp, Y_sp, Z_sp,'FaceColor',[1 0.41 0.16],'EdgeColor','none')
alpha(0.8)
lightangle(90,0)

subsat_p = [X(end), Y(end), Z(end)]/norm([X(end), Y(end), Z(end)]);
plot3([X(end), subsat_p(1)], [Y(end), subsat_p(2)], [Z(end), subsat_p(3)], 'color', 'black', 'linewidth', 1.5);
for SZA = 10:10:170
    radius = sin(SZA*pi/180);
    y_mars = linspace(-radius, radius, 100);
    z_mars = sqrt(radius^2-y_mars.^2);
    x_mars = cos(SZA*pi/180)*ones(size(y_mars));
    if(mod(SZA, 30) == 0)
        plot3(x_mars, y_mars, z_mars, 'color', [1 0 0], 'linewidth', 2)
        plot3(x_mars, y_mars, -z_mars, 'color', [1 0 0], 'linewidth', 2)
    else
        plot3(x_mars, y_mars, z_mars, 'color', [1 0 0])
        plot3(x_mars, y_mars, -z_mars, 'color', [1 0 0])
    end
end
for lat = -80:10:80
    radius = cos(lat*pi/180);
    y_mars = linspace(-radius, radius, 100);
    x_mars = sqrt(radius^2-y_mars.^2);
    z_mars = sin(lat*pi/180)*ones(size(y_mars));
    if(mod(lat, 30)==0)
        plot3(x_mars, y_mars, z_mars, 'color', 'black', 'linewidth', 2)
        plot3(-x_mars, y_mars, z_mars, 'color', 'black', 'linewidth', 2)
    else
        plot3(x_mars, y_mars, z_mars, 'color', 'black')
        plot3(-x_mars, y_mars, z_mars, 'color', 'black')
    end
end
for th_e = 0:15:345
    x_mars = linspace(-1,1,100);
    y_mars = sqrt(1-x_mars.^2);
    z_mars = zeros(size(x_mars));
    if(mod(th_e, 45)==0)
        h = plot3(x_mars, y_mars, z_mars, 'color', 'green', 'linewidth', 2);
    else
        h = plot3(x_mars, y_mars, z_mars, 'color', 'green');
    end
    rotate(h, [1 0 0], th_e, [0 0 0])
end

hold off

axis equal

legend([h1,h2,h3,h4,h5,h6,h7],...
    {'Mars normal', 'Electric field', 'Sun', 'S/C velocity', 'B', 'v_p', 'V_{O^+}'})

annotation('textbox', [0.6911,0.0333,0.2536,0.7001]...
    , 'String', ["MSO, R_M:" num2str([X_mso(end) Y_mso(end) Z_mso(end)], '%10.2f') newline...
    "MSE, R_M:" num2str(XYZ_mse(:,end)', '%10.2f') newline...
    "height = " num2str(3390*sqrt(X(end)^2+Y(end)^2+Z(end)^2)-3390, "%10.0f") "km" newline...
    "Spacecraft velocity, km/s:" num2str(sc_v, '%10.2f')])

pos = get(f1_h, 'position');
pos(1) = pos(1) + round(pos(3)/2);
pos(2) = pos(2) - round(pos(4)/2);

f2_h = figure();

h1 = quiver3(0,0,0, Mars_norm(1), Mars_norm(2), Mars_norm(3),...
    'color', 'black', 'linewidth', 2);
hold on
% electric field
h2 = quiver3(0,0,0, E(1),E(2),E(3), 'color', 'green', 'linewidth', 2);
% sun
h3 = quiver3(0,0,0, 1,0,0, 'color', 'red', 'linewidth', 2);
% s/c velocity
h4 = quiver3(0,0,0, 2*sc_v_norm(1),2*sc_v_norm(2),2*sc_v_norm(3),...
    'color', 'blue', 'linewidth', 2);
% MF
magf_mse_norm = magf_mse./sqrt(sum(magf_mse.^2,2));
h5 = quiver3(0*X',0*Y',0*Z',magf_mse_norm(:,1),magf_mse_norm(:,2),magf_mse_norm(:,3), 0,...
    'color', [0.5,0.5,0], 'linewidth',1);
% v_p
v_p_mse_norm = v_p_mse./sqrt(sum(v_p_mse.^2,2));
h6 = quiver3(0*X',0*Y',0*Z',v_p_mse_norm(:,1),v_p_mse_norm(:,2),v_p_mse_norm(:,3), 0,...
    'color', [0,0.8,0.8], 'linewidth',1);
v_O_mse_norm = v_O_mse./sqrt(sum(v_O_mse.^2,2));
h7 = quiver3(0*X',0*Y',0*Z',v_O_mse_norm(:,1),v_O_mse_norm(:,2),v_O_mse_norm(:,3), 0,...
    'color', [0.8,0,0.8], 'linewidth',1);

xlabel('X_{MSE}')
ylabel('Y_{MSE}')
zlabel('Z_{MSE}')

for lat = -80:10:80
    radius = cos(lat*pi/180);
    y_mars = linspace(-radius, radius, 100);
    x_mars = sqrt(radius^2-y_mars.^2);
    z_mars = sin(lat*pi/180)*ones(size(y_mars));
    if(mod(lat, 30) == 0)
        plot3(x_mars, y_mars, z_mars, 'color', 'black', 'linewidth', 2)
        plot3(-x_mars, y_mars, z_mars, 'color', 'black', 'linewidth', 2)
    else
        plot3(x_mars, y_mars, z_mars, 'color', 'black')
        plot3(-x_mars, y_mars, z_mars, 'color', 'black')
    end
end
for lon = 0:15:345
    z_mars = linspace(-1,1,100);
    x_mars = sqrt(1-z_mars.^2);
    y_mars = zeros(size(x_mars));
    if(mod(lon, 45) == 0)
        h = plot3(x_mars, y_mars, z_mars, 'color', [0 0.5 0], 'linewidth', 2);
    else
        h = plot3(x_mars, y_mars, z_mars, 'color', [0 0.5 0]);
    end
    rotate(h, [0 0 1], lon, [0 0 0])
end

hold off

axis equal

legend([h1,h2,h3,h4,h5,h6,h7],...
    {'Mars normal', 'Electric field', 'Sun', 'S/C velocity', 'B', 'v_p','v_{O^+}'})
set(f2_h, 'position', pos)

pos(1) = pos(1) + round(pos(3)/2);
pos(2) = pos(2) - round(pos(4)/2);
f3_h = figure();
set(f3_h, 'position', pos)

[x, y] = vec2merc(Mars_norm(1), Mars_norm(2), Mars_norm(3));
h1 = plot(x,y,'marker', '.', 'markersize', 20, 'color', 'black', 'linestyle', 'none');
hold on
% electric field
[x, y] = vec2merc(E(1),E(2),E(3));
h2 = plot(x,y,'marker', '.', 'markersize', 20, 'color', 'green', 'linestyle', 'none');
% sun
[x, y] = vec2merc(1,0,0);
h3 = plot(x,y,'marker', '.', 'markersize', 20, 'color', 'red', 'linestyle', 'none');
% s/c velocity
[x, y] = vec2merc(2*sc_v_norm(1),2*sc_v_norm(2),2*sc_v_norm(3));
h4 = plot(x,y,'marker', '.', 'markersize', 20, 'color', 'blue', 'linestyle', 'none');
% MF
[x, y] = vec2merc(magf_mse_norm(:,1),magf_mse_norm(:,2),magf_mse_norm(:,3));
h5 = plot(x,y,'marker', '.', 'markersize', 12, 'color', [0.5,0.5,0], 'linestyle', 'none');
% v_p
[x, y] = vec2merc(v_p_mse_norm(:,1),v_p_mse_norm(:,2),v_p_mse_norm(:,3));
h6 = plot(x,y,'marker', '.', 'markersize', 12, 'color', [0,0.8,0.8], 'linestyle', 'none');
% v_O
[x, y] = vec2merc(v_O_mse_norm(:,1),v_O_mse_norm(:,2),v_O_mse_norm(:,3));
h7 = plot(x,y,'marker', '.', 'markersize', 12, 'color', [0.8,0,0.8], 'linestyle', 'none');

legend([h1,h2,h3,h4,h5,h6,h7],...
    {'Mars normal', 'Electric field', 'Sun', 'S/C velocity', 'B', 'v_p','v_{O^+}'})

axis equal
hold off
xlim([-180 180])
ylim([-90 90])
grid on
set(gca, 'xtick', [-180:30:180], 'ytick', [-90:30:90])
xlabel('Azimuth angle in MSE, deg')
ylabel('Polar angle in MSE, deg')

    function [phi, th] = vec2merc(x,y,z)
        th = asin(z)*180/pi;
        phi = zeros(size(x));
        v = x > 0 & y <= 0;
        phi(v) = atan(y(v)./x(v));
        v = x <= 0 & y < 0;
        phi(v) = -pi + atan(y(v)./x(v));
        v = x < 0 & y >= 0;
        phi(v) = pi + atan(y(v)./x(v));
        v = x >= 0 & y > 0;
        phi(v) = pi + atan(y(v)./x(v));
        phi = phi*180/pi;
        
        phi(phi>180) = phi(phi>180) - 360;
        phi(phi<-180) = phi(phi<-180) + 360;
    end

end