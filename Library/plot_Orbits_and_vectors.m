function plot_Orbits_and_vectors(products_O, products_H, products_O2, products_mag,...
    time_input, altit_low, altit_high)

TimeStart = time_input(1);
TimeEnd = time_input(2);

elapsed = TimeEnd-TimeStart;
elapsed_sec = hour(elapsed)*3600+minute(elapsed)*60+second(elapsed);
if(elapsed_sec<4)
    TimeEnd = TimeEnd + datenum('00-01-0000 00:00:06', 'dd-mm-yyyy HH:MM:SS');
end

[x, y, z, ...
    n_p,  T_p, T_p_energy, v_x_p,  v_y_p,  v_z_p,  v_p,...
    n_O,  T_O, T_O_energy, v_x_O,  v_y_O,  v_z_O,  v_O,...
    n_O2, T_O2, T_O2_energy, v_x_O2, v_y_O2, v_z_O2, v_O2,...
    Blx, Bly, Blz, Bl,...
    altit, epoch,...
    mf_epoch, Bx, By, Bz, B, x_mag,y_mag,z_mag] = set_time_range_for_variables(products_O, products_H, products_O2, products_mag, TimeStart, TimeEnd, altit_low, altit_high );
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
SZA = 180 / pi * atan ( R_cyl_O./x);
SZA (x<0) = SZA (x<0)  + 180;

%**********************  velocities ***********
v_x_bulk = (v_x_p.*n_p + v_x_O.*n_O*16 + v_x_O2.*n_O*32 )./(n_p + n_O*16 + n_O2*32 );
v_y_bulk = (v_y_p.*n_p + v_y_O.*n_O*16 + v_y_O2.*n_O*32 )./(n_p + n_O*16 + n_O2*32 );
v_z_bulk = (v_z_p.*n_p + v_z_O.*n_O*16 + v_z_O2.*n_O*32 )./(n_p + n_O*16 + n_O2*32 );
v_bulk = sqrt(v_x_bulk.^2 + v_y_bulk.^2 + v_z_bulk.^2);

%*********************  Electric fiel calculation ***********
Ex = v_z_bulk.*Bly - v_y_bulk.*Blz;
Ey = v_x_bulk.*Blz - v_z_bulk.*Blx;
Ez = v_y_bulk.*Blx - v_x_bulk.*Bly;

%epoch = epoch( epoch >= TimeStart & epoch <= TimeEnd);
%}

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
rmb=sqrt(xmb.^2 + ymb.^2);
plot2 = plot(xmb,ymb,'k'); %***** plot magnetosphere boundary ***
hold on
axis equal tight
t=(0:0.001:2*pi);

%*************** plot Mars ****
plot_3 = plot(cos(t),sin(t),'k');  %*************** plot Mars

%**************  Plot properties *************
k1 = xlabel('X_{MSO}, R_M');
ylabel('R_{MSO}, R_M')
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
xlabel('X_{MSO}, R_M')
ylabel('Y_{MSO}, R_M')
zlabel('Z_{MSO}, R_M')
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

%*************  Legend *************

% vectors4.DisplayName = 'O^+';
% vectors5.DisplayName = 'p';
% vectors6.DisplayName = '- [V_{bulk}xB]';
% plot7.Annotation.LegendInformation.IconDisplayStyle = 'off';
% plot8.Annotation.LegendInformation.IconDisplayStyle = 'off';
% legend show
% legend('location','northeastoutside')

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
xlabel('X_{MSO}, R_M')
ylabel('Y_{MSO}, R_M')
zlabel('Z_{MSO}, R_M')
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
xlabel('X_{MSO}, R_M')
ylabel('Y_{MSO}, R_M')
zlabel('Z_{MSO}, R_M')
axis equal tight
hold on
%plot(x,z)

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
%f3 = subplot(3,3,2);
f3 = subplot(3,3,3);
t=(0:0.001:2*pi);
plot_3 = plot3(zeros(size(t)),cos(t),sin(t),'k');
hold on
%plot4 = plot(y,z);
xlabel('X_{MSO}, R_M')
ylabel('Y_{MSO}, R_M')
zlabel('Z_{MSO}, R_M')
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
%prop_subplot332 = get(f3, 'position');
set(f3,'position',[0.5703    0.7200    0.3347    0.2113]);

%*******************   figure(1.4) Magnetic field angular distribution in YZ plan **************
%subplot(3,3,5)
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
% h.Color = [1 0 1];
% ylabel('-B_y')
% xlabel('-B_z')
% vec_pos = get(get(gca, 'XLabel'), 'Position');
% [N,edges]=histcounts(atan_Bz_to_By,12);
% edges_new = edges(1:end-1) + diff(edges)/2;
set(gca, 'FontSize',13)


%*******************   figure(1.5) Electric field angular distribution in YZ plane **************
%subplot(3,3,6)
subplot(3,3,7)
polarhistogram(atan_angle(Ez,Ey),12, 'facecolor', 'green')
%rose(atan_angle(Ez,Ey),12, 'color', 'green')
% h.Color='green';
% ylabel('-E_y')
% xlabel('-E_z')
% vec_pos = get(get(gca, 'XLabel'), 'Position');
% [N,edges]=histcounts(atan_Bz_to_By,12);
% edges_new = edges(1:end-1) + diff(edges)/2;
set(gca, 'FontSize',16)

% %--------------------------------%
% f4 = subplot(3,3,7);
%
% plot_Mars;
% hold on
% plot(x, z, 'k')
% scale_XZ = max(sqrt(v_x_p(1:step:end).^2 + v_z_p(1:step:end).^2));
% quiver(x(1:step:end), z(1:step:end), v_x_p(1:step:end)/scale_XZ, v_z_p(1:step:end)/scale_XZ,'Autoscale','off','Color','b');
% axis equal tight
% grid on
% ylim([-2 2])
% xlim([-2 2])
% xlabel('X_{MSO}, R_M')
% ylabel('Z_{MSO}, R_M')
% plot(x(1), z(1),'*g')
% plot(x(end), z(end),'*r')
% set(gca, 'FontSize',16)
% prop_subplot337 = get(f4, 'position');
% set(f4,'position',[prop_subplot337(1) prop_subplot337(2) prop_subplot337(3)*1.1 prop_subplot337(4)*1.1]);
% %------------------------------%
%
% %------------------------------%
% f5=subplot(3,3,8);
% plot_Mars;
% hold on
% plot(x, y, 'k')
% scale_XY = max(sqrt(v_x_p(1:step:end).^2 + v_y_p(1:step:end).^2));
% quiver(x(1:step:end), y(1:step:end), v_x_p(1:step:end)/scale_XY, v_y_p(1:step:end)/scale_XY,'Autoscale','off','Color','b');
% axis equal tight
% grid on
% ylim([-2 2])
% xlim([-2 2])
% xlabel('X_{MSO}, R_M')
% ylabel('Y_{MSO}, R_M')
% plot(x(1), y(1),'*g')
% plot(x(end), y(end),'*r')
% % plot_title = strcat(datestr(epoch(1),'yyyy-mm-dd'), {' '}, datestr(epoch(1),'HH:MM:SS'),{'-'}, datestr(epoch(end),'HH:MM:SS'));
% % title(plot_title)
% set(gca, 'FontSize',16)
% prop_subplot338 = get(f5, 'position');
% set(f5,'position',[prop_subplot338(1) prop_subplot338(2) prop_subplot338(3)*1.1 prop_subplot338(4)*1.1]);
% %------------------------------%
%
% %------------------------------%
% f6=subplot(3,3,9);
% plot_Mars;
% hold on
% plot(y, z, 'k')
% scale_YZ = max(sqrt(v_y_p(1:step:end).^2 + v_z_p(1:step:end).^2));
% quiver(y(1:step:end), z(1:step:end), v_y_p(1:step:end)/scale_YZ, v_z_p(1:step:end)/scale_YZ,'Autoscale','off','Color','b');
% axis equal tight
% grid on
% ylim([-2 2])
% xlim([-2 2])
% xlabel('Y_{MSO}, R_M')
% ylabel('Z_{MSO}, R_M')
% plot(y(1), z(1),'*g')
% plot(y(end), z(end),'*r')
% set(gca, 'FontSize',16)
% prop_subplot339 = get(f6, 'position');
% set(f6,'position',[prop_subplot339(1) prop_subplot339(2) prop_subplot339(3)*1.1 prop_subplot339(4)*1.1]);
% %------------------------------%

% %****** The second figure ******
%
% f = figure;
%
% set(f,'OuterPosition',[0 scrsz(4)*0.01 scrsz(3) scrsz(4)])
% subplot(2, 3, 1)
% plot_Mars;
% hold on
% plot(x, z, 'k')
% vectors = quiver(x(1:step:end), z(1:step:end), v_x_p(1:step:end), v_z_p(1:step:end),'Autoscale','off','Color','b');
% % hU = get(vectors,'UData') ;
% % hV = get(vectors,'VData') ;
% % set(vectors,'UData',scale_velocity*hU,'VData',scale_velocity*hV)
% axis equal tight
% grid on
% ylim([-2 2])
% xlim([-2 2])
% xlabel('X_{MSO}, R_M')
% ylabel('Z_{MSO}, R_M')
% plot(x(1), z(1),'*g')
% plot(x(end), z(end),'*r')
% set(gca, 'FontSize',16)
%
%
% subplot(2, 3, 2)
% plot_Mars;
% hold on
% plot(x, y, 'k')
% vectors = quiver(x(1:step:end), y(1:step:end), v_x_p(1:step:end), v_y_p(1:step:end),'Autoscale','off','Color','b');
% % hU = get(vectors,'UData') ;
% % hV = get(vectors,'VData') ;
% % set(vectors,'UData',scale_velocity*hU,'VData',scale_velocity*hV)
% axis equal tight
% grid on
% ylim([-2 2])
% xlim([-2 2])
% xlabel('X_{MSO}, R_M')
% ylabel('Y_{MSO}, R_M')
% plot(x(1), y(1),'*g')
% plot(x(end), y(end),'*r')
% plot_title = strcat(datestr(epoch(1),'yyyy-mm-dd'), {' '}, datestr(epoch(1),'HH:MM:SS'),{'-'}, datestr(epoch(end),'HH:MM:SS'));
% title(plot_title)
% set(gca, 'FontSize',16)
%
%
% subplot(2, 3, 3)
% plot_Mars;
% hold on
% plot(y, z, 'k')
% vectors = quiver(y(1:step:end), z(1:step:end), v_y_p(1:step:end), v_z_p(1:step:end),'Autoscale','off','Color','b');
% % hU = get(vectors,'UData') ;
% % hV = get(vectors,'VData') ;
% % set(vectors,'UData',scale_velocity*hU,'VData',scale_velocity*hV)
% axis equal tight
% grid on
% ylim([-2 2])
% xlim([-2 2])
% xlabel('Y_{MSO}, R_M')
% ylabel('Z_{MSO}, R_M')
% plot(y(1), z(1),'*g')
% plot(y(end), z(end),'*r')
% set(gca, 'FontSize',16)
%
%
%
% subplot(2, 3, 4)
% plot_Mars;
% hold on
% plot(x, z, 'k')
% vectors = quiver(x(1:step:end), z(1:step:end), v_x_O(1:step:end), v_z_O(1:step:end),'Autoscale','off','Color','r');
% % hU = get(vectors,'UData') ;
% % hV = get(vectors,'VData') ;
% % set(vectors,'UData',scale_velocity*hU,'VData',scale_velocity*hV)
% axis equal tight
% grid on
% ylim([-2 2])
% xlim([-2 2])
% xlabel('X_{MSO}, R_M')
% ylabel('Z_{MSO}, R_M')
% plot(x(1), z(1),'*g')
% plot(x(end), z(end),'*r')
% set(gca, 'FontSize',16)
%
% subplot(2, 3, 5)
% plot_Mars;
% hold on
% plot(x, y, 'k')
% vectors = quiver(x(1:step:end), y(1:step:end), v_x_O(1:step:end), v_y_O(1:step:end),'Autoscale','off','Color','r');
% % hU = get(vectors,'UData') ;
% % hV = get(vectors,'VData') ;
% % set(vectors,'UData',scale_velocity*hU,'VData',scale_velocity*hV)
% axis equal tight
% grid on
% ylim([-2 2])
% xlim([-2 2])
% xlabel('X_{MSO}, R_M')
% ylabel('Y_{MSO}, R_M')
% plot(x(1), y(1),'*g')
% plot(x(end), y(end),'*r')
% title('O^+')
% set(gca, 'FontSize',16)
%
%
% subplot(2, 3, 6)
% plot_Mars;
% hold on
% plot(y, z, 'k')
% vectors = quiver(y(1:step:end), z(1:step:end), v_y_O(1:step:end), v_z_O(1:step:end),'Autoscale','off','Color','r');
% % hU = get(vectors,'UData') ;
% % hV = get(vectors,'VData') ;
% % set(vectors,'UData',scale_velocity*hU,'VData',scale_velocity*hV)
% axis equal tight
% grid on
% ylim([-2 2])
% xlim([-2 2])
% xlabel('Y_{MSO}, R_M')
% ylabel('Z_{MSO}, R_M')
% plot(y(1), z(1),'*g')
% plot(y(end), z(end),'*r')
% set(gca, 'FontSize',16)

end