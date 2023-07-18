function plot_Orbits_and_vectors_MSE(products_O, products_H, products_O2, products_mag,...
    time_input, altit_low, altit_high)

TimeStart = time_input(3);
TimeEnd = time_input(4);

elapsed = TimeEnd-TimeStart;
elapsed_sec = hour(elapsed)*3600+minute(elapsed)*60+second(elapsed);
if(elapsed_sec<4)
    TimeEnd = TimeEnd + datenum('00-01-0000 00:00:06', 'dd-mm-yyyy HH:MM:SS');
end

[x, y, z, ...
    n_p,  ~, ~, v_x_p,  v_y_p,  v_z_p,  ~,...
    n_O,  ~, ~, v_x_O,  v_y_O,  v_z_O,  ~,...
    n_O2, ~, ~, v_x_O2, v_y_O2, v_z_O2, ~,...
    Blx, Bly, Blz, ~,...
    ~, epoch,...
    ~, Bx, By, Bz, ~, x_mag,y_mag,z_mag] =...
    set_time_range_for_variables(products_O, products_H, products_O2,...
    products_mag, TimeStart, TimeEnd, altit_low, altit_high );

%==============Calculate rotation matrix to MSE=========================
[~, ~, ~, ...
    n_p_u,  ~, ~, v_x_p_u,  v_y_p_u,  v_z_p_u,  ~,...
    n_O_u,  ~, ~, v_x_O_u,  v_y_O_u,  v_z_O_u,  ~,...
    n_O2_u, ~, ~, v_x_O2_u, v_y_O2_u, v_z_O2_u, ~,...
    Bx_upstr, By_upstr, Bz_upstr, ~,...
    ~, ~,...
    ~, ~, ~, ~, ~, ~,~,~] =...
    set_time_range_for_variables(products_O, products_H, products_O2,...
    products_mag, time_input(1), time_input(2), altit_low, altit_high );

%**********************  velocities upstream ***********
v_x_bulk_u = (v_x_p_u.*n_p_u + v_x_O_u.*n_O_u*16 + v_x_O2_u.*n_O_u*32 )./(n_p_u + n_O_u*16 + n_O2_u*32 );
v_y_bulk_u = (v_y_p_u.*n_p_u + v_y_O_u.*n_O_u*16 + v_y_O2_u.*n_O_u*32 )./(n_p_u + n_O_u*16 + n_O2_u*32 );
v_z_bulk_u = (v_z_p_u.*n_p_u + v_z_O_u.*n_O_u*16 + v_z_O2_u.*n_O_u*32 )./(n_p_u + n_O_u*16 + n_O2_u*32 );

E_upstr = -cross([v_x_bulk_u, v_y_bulk_u, v_z_bulk_u],...
    [Bx_upstr, By_upstr, Bz_upstr]);
E_upstr = mean(E_upstr(~isnan(E_upstr(:,1)),:), 1);
E_upstr = E_upstr./norm(E_upstr);

% Calculation of rotation angle and rotation matrix from MSO to MSE
rt_YZ = atan(E_upstr(2)/E_upstr(3));
if(E_upstr(3)<0)
    rt_YZ = rt_YZ+pi;
end
rt_M = [1 0 0;...
    0 cos(rt_YZ) -sin(rt_YZ);...
    0 sin(rt_YZ) cos(rt_YZ)];
%========================================================================

%valid_p] = set_time_range_for_variables(products_O, products_H, products_O2, products_mag, TimeStart, TimeEnd, altit_low, altit_high );

%*****calc_sc_velocity*********
[~, ~, ~, ~, ~, ~, ~, ~, Rm] = constants;
cur_sec = hour(epoch)*3600+minute(epoch)*60+second(epoch);
dx = Rm*[diff(x), diff(y), diff(z)];
dt = diff(cur_sec);
sc_vel = dx./dt;
sc_vel(end+1,:) = sc_vel(end,:);

v_x_p = v_x_p + sc_vel(:,1);
v_y_p = v_y_p + sc_vel(:,2);
v_z_p = v_z_p + sc_vel(:,3);
v_x_O = v_x_O + sc_vel(:,1);
v_y_O = v_y_O + sc_vel(:,2);
v_z_O = v_z_O + sc_vel(:,3);
v_x_O2 = v_x_O2 + sc_vel(:,1);
v_y_O2 = v_y_O2 + sc_vel(:,2);
v_z_O2 = v_z_O2 + sc_vel(:,3);
%*****************************

step = 2;
R_cyl_O = sqrt(y.^2 + z.^2);

%**********************  velocities ***********
v_x_bulk = (v_x_p.*n_p + v_x_O.*n_O*16 + v_x_O2.*n_O*32 )./(n_p + n_O*16 + n_O2*32 );
v_y_bulk = (v_y_p.*n_p + v_y_O.*n_O*16 + v_y_O2.*n_O*32 )./(n_p + n_O*16 + n_O2*32 );
v_z_bulk = (v_z_p.*n_p + v_z_O.*n_O*16 + v_z_O2.*n_O*32 )./(n_p + n_O*16 + n_O2*32 );
%v_bulk = sqrt(v_x_bulk.^2 + v_y_bulk.^2 + v_z_bulk.^2);

%*********************  Electric fiel calculation ***********
Ex = v_z_bulk.*Bly - v_y_bulk.*Blz;
Ey = v_x_bulk.*Blz - v_z_bulk.*Blx;
Ez = v_y_bulk.*Blx - v_x_bulk.*Bly;

%=========================Rotate values to MSE======================

pos_sc_mse = rt_M*[x';y';z'];
x = pos_sc_mse(1, :)';
y = pos_sc_mse(2, :)';
z = pos_sc_mse(3, :)';
v_p_mse = rt_M*[v_x_p'; v_y_p'; v_z_p'];
v_x_p = v_p_mse(1, :)';
v_y_p = v_p_mse(2, :)';
v_z_p = v_p_mse(3, :)';
v_O_mse = rt_M*[v_x_O'; v_y_O'; v_z_O'];
v_x_O = v_O_mse(1, :)';
v_y_O = v_O_mse(2, :)';
v_z_O = v_O_mse(3, :)';
v_O2_mse = rt_M*[v_x_O2'; v_y_O2'; v_z_O2'];
v_x_O2 = v_O2_mse(1, :)';
v_y_O2 = v_O2_mse(2, :)';
v_z_O2 = v_O2_mse(3, :)';
E_mse = rt_M*[Ex'; Ey'; Ez'];
Ex = E_mse(1, :)';
Ey = E_mse(2, :)';
Ez = E_mse(3, :)';
B_mse = rt_M*[Bx'; By'; Bz'];
Bx = B_mse(1, :)';
By = B_mse(2, :)';
Bz = B_mse(3, :)';

%===================================================================

f = figure;

scrsz=get(0,'ScreenSize');
set(f,'OuterPosition',[0 scrsz(4)*0.01 scrsz(3)/2 scrsz(4)])


%************************* Figure(1.1) Cylindrical system ************
f1 = subplot(3,3,1);
%***************  Bow shock ************
yb=0:0.01:4.5;
xb=(-exp(yb).^0.2)*0.2+0.04*yb-0.2.*yb.^2.15+1.6;
plot1 = plot(xb,yb,'k'); %******** plot bow shock ******
hold on

%************   Magnetosphere boundary ******
ymb=0:0.01:2.5;
xmb=(-exp(ymb).^2)*0.05+0.1*ymb-0.2.*ymb.^2+1.25;
%rmb=sqrt(xmb.^2 + ymb.^2);
plot2 = plot(xmb,ymb,'k'); %***** plot magnetosphere boundary ***
hold on
axis equal tight
t=(0:0.001:2*pi);

%*************** plot Mars ****
plot_3 = plot(cos(t),sin(t),'k');  %*************** plot Mars

%**************  Plot properties *************
k1 = xlabel('X_{MSE}, R_M');
ylabel('R, R_M')
ylim([0 2])
xlim([-0.5 2])
%set(k1,'position', [0.802175019101182,-0.185223688875722,-1])


%************* plot orbit *********
plot4 = plot(x,R_cyl_O);

%************  Calculate
v_phi_O = atan_angle(v_z_O, v_y_O);
v_phi_p = atan_angle(v_z_p, v_y_p);
phi = atan_angle(z,y);
E_phi = atan_angle(Ez,Ey); %Electric field

% Scale
scale_O = max(sqrt(v_x_O(1:step:end).^2 + (sqrt( v_y_O(1:step:end).^2 + v_z_O(1:step:end).^2 ).*cos( v_phi_O (1:step:end) - phi(1:step:end) )).^2));
scale_p = max(sqrt(v_x_p(1:step:end).^2 + (sqrt( v_y_p(1:step:end).^2 + v_z_p(1:step:end).^2 ).*cos( v_phi_p (1:step:end) - phi(1:step:end) )).^2));
scale_E = max(sqrt(Ex(1:step:end).^2 + (sqrt( Ey(1:step:end).^2 + Ez(1:step:end).^2 ).*cos( E_phi (1:step:end) - phi(1:step:end))).^2));

%*********************  Plot velocity and electric field vectors in
%cylindrical system **************************
vectors1 = quiver(x(1:step:end),sqrt(y(1:step:end).^2+z(1:step:end).^2),...
    v_x_O(1:step:end)/scale_O, sqrt(v_y_O(1:step:end).^2 + v_z_O(1:step:end).^2 ).*cos( v_phi_O (1:step:end) - phi(1:step:end))/scale_O,...
    'Autoscale','off','Color','r');
vectors2 = quiver(x(1:step:end),sqrt(y(1:step:end).^2+z(1:step:end).^2),...
    v_x_p(1:step:end)/scale_p, sqrt(v_y_p(1:step:end).^2 + v_z_p(1:step:end).^2 ).*cos( v_phi_p (1:step:end) - phi(1:step:end))/scale_p,...
    'Autoscale','off','Color','b');
vectors3 = quiver(x(1:step:end),sqrt(y(1:step:end).^2+z(1:step:end).^2),...
    Ex(1:step:end)/scale_E, sqrt(Ey(1:step:end).^2 + Ez(1:step:end).^2 ).*cos( E_phi (1:step:end) - phi(1:step:end))/scale_E,...
    'Autoscale','off','Color','g');
grid on

O_val = sqrt(v_x_O(1:step:end).^2+sqrt( v_y_O(1:step:end).^2 + v_z_O(1:step:end).^2 ).*cos( v_phi_O (1:step:end) - phi(1:step:end) ).^2);
p_val = sqrt(v_x_p(1:step:end).^2+sqrt( v_y_p(1:step:end).^2 + v_z_p(1:step:end).^2 ).*cos( v_phi_p (1:step:end) - phi(1:step:end) ).^2);
E_val = sqrt(Ex(1:step:end).^2+sqrt( Ey(1:step:end).^2 + Ez(1:step:end).^2 ).*cos( E_phi (1:step:end) - phi(1:step:end) ).^2);

% Изображение текста
annotation('textbox',[0.64 0.15 0.1 0.1],'String',['v_O is from ' num2str(min(O_val)) ' to ' num2str(max(O_val)) ' km/s' newline ...
    ['v_p is from ' num2str(min(p_val)) ' to ' num2str(max(p_val))  ' km/s'] newline ...
    ['E is from ' num2str(min(E_val)) ' to ' num2str(max(E_val))  ' mV/m']], 'fontsize', 14);

%**********  Legend  **************
vectors1.DisplayName = 'O^+';
vectors2.DisplayName = 'p';
vectors3.DisplayName = '- [V_{bulk}xB]';
plot1.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot2.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot_3.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot4.Annotation.LegendInformation.IconDisplayStyle = 'off';
legend show
legend('location','northeastoutside')
set(gca, 'FontSize',14)

prop_subplot331 = get(f1, 'position');
disp(prop_subplot331)
%set(f1,'position',[0.6*prop_subplot331(1) prop_subplot331(2) 1*prop_subplot331(3) 1.2*prop_subplot331(4)]);
set(f1,'position',[0.6*prop_subplot331(1) 0.95*prop_subplot331(2) 1.4*prop_subplot331(3) 1.4*prop_subplot331(4)]);

plot_title = strcat(datestr(epoch(1),'yyyy-mm-dd'), {' '}, datestr(epoch(1),'HH:MM:SS'),{'-'}, datestr(epoch(end),'HH:MM:SS'));
ht = title(plot_title);
title_pos = get(ht, 'position');
set(ht, 'position', [title_pos(1)+2 1.1*title_pos(2)])

%******************* figure(1.2) Orbit and velocities vectors in YZ plane  **************
f4 = subplot(3,3,4);
t=(0:0.001:2*pi);
plot7 = plot3(zeros(size(t)),cos(t),sin(t),'k');
grid on
xlabel('X_{MSE}, R_M')
ylabel('Y_{MSE}, R_M')
zlabel('Z_{MSE}, R_M')
axis equal tight
hold on
%plot8 = plot3(zeros(size(y)),y,z);

% Scale_YZ
scale_O_YZ = max(sqrt(v_y_O(1:step:end).^2 + v_z_O(1:step:end).^2));
scale_p_YZ = max(sqrt(v_y_p(1:step:end).^2 + v_z_p(1:step:end).^2));
scale_E_YZ = max(sqrt(Ey(1:step:end).^2 + Ez(1:step:end).^2));

vectors4 = quiver3(x(1:step:end),y(1:step:end),z(1:step:end),v_x_O(1:step:end)/scale_O_YZ,v_y_O(1:step:end)/scale_O_YZ, v_z_O(1:step:end)/scale_O_YZ,'Autoscale','off','DisplayName','O','Color','r');
vectors5 = quiver3(x(1:step:end),y(1:step:end),z(1:step:end),v_x_p(1:step:end)/scale_p_YZ,v_y_p(1:step:end)/scale_p_YZ,v_z_p(1:step:end)/scale_p_YZ,'Autoscale','off','Color','b');
vectors6 = quiver3(x(1:step:end),y(1:step:end),z(1:step:end),Ex(1:step:end)/scale_E_YZ,Ey(1:step:end)/scale_E_YZ,Ez(1:step:end)/scale_E_YZ,'Autoscale','off','Color','g');
view(90,0)



%**********  Plot properties *************
xlim([-2 2])
ylim([-2 2])
zlim([-2 2])
set(gca, 'FontSize',14)
prop_subplot334 = get(f4, 'position');
set(f4,'position',[0.5*prop_subplot334(1) prop_subplot334(2) 1.1*prop_subplot334(3) 1.1*prop_subplot334(4)]);

%******************* figure(1.2) Orbit and velocities vectors in XY plane  **************
f5 = subplot(3,3,5);
t=(0:0.001:2*pi);
plot3(cos(t),sin(t), zeros(size(t)), 'k')
grid on
xlabel('X_{MSE}, R_M')
ylabel('Y_{MSE}, R_M')
zlabel('Z_{MSE}, R_M')
axis equal tight
hold on
%plot(x,y)

% Scale_XY
scale_O_XY = max(sqrt(v_x_O(1:step:end).^2 + v_y_O(1:step:end).^2));
scale_p_XY = max(sqrt(v_x_p(1:step:end).^2 + v_y_p(1:step:end).^2));
scale_E_XY = max(sqrt(Ex(1:step:end).^2 + Ey(1:step:end).^2));

quiver3(x(1:step:end),y(1:step:end),z(1:step:end),v_x_O(1:step:end)/scale_O_XY,v_y_O(1:step:end)/scale_O_XY,v_z_O(1:step:end)/scale_O_XY,'Autoscale','off','DisplayName','O','Color','r');
quiver3(x(1:step:end),y(1:step:end),z(1:step:end),v_x_p(1:step:end)/scale_p_XY,v_y_p(1:step:end)/scale_p_XY,v_z_p(1:step:end)/scale_p_XY,'Autoscale','off','Color','b');
quiver3(x(1:step:end),y(1:step:end),z(1:step:end),Ex(1:step:end)/scale_E_XY,Ey(1:step:end)/scale_E_XY,Ez(1:step:end)/scale_E_XY,'Autoscale','off','Color','g');
view(0,90)

%**********  Plot properties *************
xlim([-2 2])
ylim([-2 2])
zlim([-2 2])
set(gca, 'FontSize',14)
prop_subplot335 = get(f5, 'position');
set(f5,'position',[0.9*prop_subplot335(1) prop_subplot335(2) 1.1*prop_subplot335(3) 1.1*prop_subplot335(4)]);


%******************* figure(1.2) Orbit and velocities vectors in XZ plane  **************
f6 = subplot(3,3,6);
t=(0:0.001:2*pi);
plot3(cos(t),zeros(size(t)),sin(t),'k')
grid on
xlabel('X_{MSE}, R_M')
ylabel('Y_{MSE}, R_M')
zlabel('Z_{MSE}, R_M')
axis equal tight
hold on

% Scale_XZ
scale_O_XZ = max(sqrt(v_x_O(1:step:end).^2 + v_z_O(1:step:end).^2));
scale_p_XZ = max(sqrt(v_x_p(1:step:end).^2 + v_z_p(1:step:end).^2));
scale_E_XZ = max(sqrt(Ex(1:step:end).^2 + Ez(1:step:end).^2));

quiver3(x(1:step:end),y(1:step:end),z(1:step:end),v_x_O(1:step:end)/scale_O_XZ,v_y_O(1:step:end)/scale_O_XZ,v_z_O(1:step:end)/scale_O_XZ,'Autoscale','off','DisplayName','O','Color','r');
quiver3(x(1:step:end),y(1:step:end),z(1:step:end),v_x_p(1:step:end)/scale_p_XZ,v_y_p(1:step:end)/scale_p_XZ,v_z_p(1:step:end)/scale_p_XZ,'Autoscale','off','Color','b');
quiver3(x(1:step:end),y(1:step:end),z(1:step:end),Ex(1:step:end)/scale_E_XZ,Ey(1:step:end)/scale_E_XZ,Ez(1:step:end)/scale_E_XZ,'Autoscale','off','Color','g');
view(0,0)

%**********  Plot properties *************
xlim([-2 2])
ylim([-2 2])
zlim([-2 2])
set(gca, 'FontSize',14)
prop_subplot336 = get(f6, 'position');
set(f6,'position',[prop_subplot336(1) prop_subplot336(2) 1.1*prop_subplot336(3) 1.1*prop_subplot336(4)]);


%*******************   figure(1.2) **************
f3 = subplot(3,3,3);
t=(0:0.001:2*pi);
plot_3 = plot3(zeros(size(t)),cos(t),sin(t),'k');
hold on
xlabel('X_{MSE}, R_M')
ylabel('Y_{MSE}, R_M')
zlabel('Z_{MSE}, R_M')
axis equal tight
hold on
grid on
step_B = 100;
scale_B = sqrt(By(1:step_B:end).^2 + Bz(1:step_B:end).^2);
vectors1 = quiver3(x_mag(1:step_B:end),y_mag(1:step_B:end),z_mag(1:step_B:end),Bx(1:step_B:end)./scale_B,By(1:step_B:end)./scale_B, Bz(1:step_B:end)./scale_B,'Autoscale','off','Color',[1 0 1]);
view(90,0)
plot_3.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot4.Annotation.LegendInformation.IconDisplayStyle = 'off';
vectors1.DisplayName = 'B';
legend show
legend('location','northeast')

xlim([-2 2])
ylim([-2 2])
zlim([-2 2])
vectors1.DisplayName = 'B';
set(gca, 'FontSize',16)
set(f3,'position',[0.5703    0.7200    0.3347    0.2113]);

%*******************   figure(1.4) Magnetic field angular distribution in YZ plan **************
subplot(3,3,8)
atan_Bz_to_By = [];
for i = 1:length(Bz)
    if (Bz(i)>=0 && By(i) >0) || (Bz(i)>=0 && By(i)<0)
        atan_Bz_to_By(i) = atan2(Bz(i),By(i))*180/pi;
    elseif(Bz(i)<0 && By(i)<=0) || (Bz(i)<0 && By(i)>=0)
        atan_Bz_to_By (i) = 360 + atan2(Bz(i),By(i))*180/pi;
    end
end
polarhistogram(atan_Bz_to_By*pi/180,12, 'facecolor', [1 0 1])

set(gca, 'FontSize',13)

%*******************   figure(1.5) Electric field angular distribution in YZ plane **************
subplot(3,3,7)
polarhistogram(atan_angle(Ez,Ey),12, 'facecolor', 'green')
set(gca, 'FontSize',16)
end