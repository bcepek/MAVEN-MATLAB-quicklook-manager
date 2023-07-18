function xticks = bin_dir(timetick,mass_num,epoch_d1,eflux_d1_cleaned,swp_ind_d1,quat_mso_d1,...
    energy,theta,phi,x,y,z, mag_filename)

% Calculation of upstream B and E vectors
mag = load(mag_filename);
[~, mag_timenid(1)] = min(abs(mag.mf_epoch - timetick(1)));
[~, mag_timenid(2)] = min(abs(mag.mf_epoch - timetick(2)));
B_upstr = mean([mag.Bx(mag_timenid(1):mag_timenid(2))...
    mag.By(mag_timenid(1):mag_timenid(2))...
    mag.Bz(mag_timenid(1):mag_timenid(2))]);
B_upstr = B_upstr./sqrt(sum(B_upstr.^2, 2));
E_upstr = -cross([-1 0 0], B_upstr);
E_upstr = E_upstr./norm(E_upstr);

% Calculation of rotation angle and rotation matrix from MSO to MSE
rt_YZ = atan(E_upstr(2)/E_upstr(3));
if(E_upstr(3)<0)
    rt_YZ = rt_YZ+pi;
end
rt_M = [1 0 0;...
    0 cos(rt_YZ) -sin(rt_YZ);...
    0 sin(rt_YZ) cos(rt_YZ)];

TimeStart = timetick(3);
TimeEnd = timetick(4);
RGBmap = load('newColormap.mat');
RGBmap.cmap(1,:) = [0 0 0.5625];
aem = 1.66*10e-27;
q = 1.60218e-19;
choose_ind = find(epoch_d1 >= TimeStart & epoch_d1 <= TimeEnd);
%*****calc_sc_velocity*********
[~, ~, ~, ~, ~, ~, ~, ~, Rm] = constants;
dr = Rm*[diff(x) diff(y) diff(z)];
epoch1 = datevec(epoch_d1);
epoch1(:,1:3) = [];
epoch1(:,1) = epoch1(:,1)*60*60;
epoch1(:,2) = epoch1(:,2)*60;
dt = diff(epoch1(:,1)+epoch1(:,2)+epoch1(:,3));
sc_vel = dr./dt;
sc_vel(end+1,:) = sc_vel(end,:);

switch mass_num
    case 1
        mass = aem;
    case 2
        mass = 4*aem;
    case 5
        mass = 16*aem;
    case 6
        mass = 32*aem;
end

figure()
subplot(1,3,1)
[X_sp, Y_sp, Z_sp] = sphere(50);
surf(X_sp, Y_sp, Z_sp,'FaceColor',[1 0.41 0.16],'EdgeColor','none')
lightangle(90,0)
alpha(0.8)
hold on
subplot(1,3,2)
surf(X_sp, Y_sp, Z_sp,'FaceColor',[1 0.41 0.16],'EdgeColor','none')
lightangle(90,0)
alpha(0.8)
hold on
subplot(1,3,3)
surf(X_sp, Y_sp, Z_sp,'FaceColor',[1 0.41 0.16],'EdgeColor','none')
lightangle(90,0)
alpha(0.8)
hold on

coords = [x,y,z];
for vecnum = 1:size(coords, 1)
    coords(vecnum,:) = (rt_M*coords(vecnum,:)')';
end
x = coords(:,1);
y = coords(:,2);
z = coords(:,3);
for i = 1:2:length(choose_ind)
    
    % Нужный поток в момент времени
    eflux = eflux_d1_cleaned(:,:,mass_num,choose_ind(i));
    swp_ind = swp_ind_d1(choose_ind(i))+1;
    [bins, energies] = find(eflux>1e7);
    table = [bins, energies];
    table = sortrows(table,2, 'descend');
    bins = table(:,1);
    energies = table(:,2);
    %         bins = 1:64;
    %         energies = 1:32;
    
    for k = 1:length(bins)
        % Углы прилета в приборной системе координат
        theta_m(k) = theta(energies(k),swp_ind,bins(k),mass_num);
        phi_m(k) = phi(energies(k),swp_ind,bins(k),mass_num);
        
        % Вектор прилета в приборной системе координат
        FOV_x(k) = sqrt(2*energy(energies(k),swp_ind,bins(k),mass_num)*q/mass)*double(cos(theta_m(k))*cos(phi_m(k)));
        FOV_y(k) = sqrt(2*energy(energies(k),swp_ind,bins(k),mass_num)*q/mass)*double(cos(theta_m(k))*sin(phi_m(k)));
        FOV_z(k) = sqrt(2*energy(energies(k),swp_ind,bins(k),mass_num)*q/mass)*double(sin(theta_m(k)));
        colormap(k) = eflux(bins(k),energies(k));
        
    end
    
    % Переход вектора поля зрения в MSO
    FOV_mso = quatrotate(quatinv(quat_mso_d1(choose_ind(i),:)), [FOV_x;FOV_y;FOV_z]' );
    % Добавляем скорость спутника
    FOV_mso = FOV_mso + sc_vel(choose_ind(i),:);
    % Переход вектора поля зрения в MSE
    for vecnum = 1:size(FOV_mso, 1)
        FOV_mso(vecnum,:) = (rt_M*FOV_mso(vecnum,:)')';
    end
    FOV = max(sqrt(FOV_mso(:,1).^2+FOV_mso(:,2).^2+FOV_mso(:,3).^2));
    
    %% --- FIELD OF VIEW PROJECTION ON XY MSO --- %%
    subplot(1,3,1)
    [~,~,~,~] = plot_MPB_bowshock;
    xlabel('X_{MSE}, Rm')
    ylabel('Y_{MSE}, Rm')
    zlabel('Z_{MSE}, Rm')
    axis equal
    hold on
    grid on
    box on
    xlim([-2 2])
    ylim([-2 2])
    zlim([-2 2])
    for j = 1:size(FOV_mso,1)
        mapind = ceil(64*colormap(j)/max(colormap));
        quiver3(x(choose_ind(i)),y(choose_ind(i)),z(choose_ind(i)),...
            FOV_mso(j,1)/FOV,FOV_mso(j,2)/FOV,FOV_mso(j,3)/FOV,'Color',...
            [RGBmap.cmap(mapind, 1) RGBmap.cmap(mapind, 2) RGBmap.cmap(mapind, 3)]);
        %[1-colormap(j)/max(colormap) 1-colormap(j)/max(colormap) 1-colormap(j)/max(colormap)]);
    end
    view(0,90)
    
    %% --- FIELD OF VIEW PROJECTION ON XZ MSO --- %%
    subplot(1,3,2)
    
    [plot1,plot2,~,~] = plot_MPB_bowshock;
    rotate(plot1,[1 0 0],-90)
    rotate(plot2,[1 0 0],-90)
    %     rotate(plot_3,[1 0 0],-90)
    
    t=(pi/2:0.001:3*pi/2);
    t1=(3*pi/2:0.001:5*pi/2);
    plot3(cos(t1),zeros(length(t1),1),sin(t1),'k')
    fill3(cos(t),zeros(length(t),1),sin(t),'k')
    xlabel('X_{MSE}, Rm')
    ylabel('Y_{MSE}, Rm')
    zlabel('Z_{MSE}, Rm')
    grid on
    box on
    axis equal
    hold on
    title(['FOV projection on XY, XZ and YZ planes in MSE '...
        datestr(TimeStart, 'yyyy-mm-dd HH:MM:SS') ' - ' datestr(TimeEnd, 'HH:MM:SS')])
    xlim([-2 2])
    ylim([-2 2])
    zlim([-2 2])
    for j = 1:size(FOV_mso,1)
        mapind = ceil(64*colormap(j)/max(colormap));
        quiver3(x(choose_ind(i)),y(choose_ind(i)),z(choose_ind(i)),...
            FOV_mso(j,1)/FOV,FOV_mso(j,2)/FOV,FOV_mso(j,3)/FOV,'Color',...
            [RGBmap.cmap(mapind, 1) RGBmap.cmap(mapind, 2) RGBmap.cmap(mapind, 3)]);
        %[1-colormap(j)/max(colormap) 1-colormap(j)/max(colormap) 1-colormap(j)/max(colormap)]);
    end
    view(0,0)
    
    %% --- FIELD OF VIEW PROJECTION ON YZ MSO --- %%
    subplot(1,3,3)
    t1=(0:0.001:2*pi);
    plot_3 = plot(cos(t1),sin(t1),'k');
    rotate(plot_3,[0 1 0], -90)
    hold on
    
    xlabel('X_{MSE}, Rm')
    ylabel('Y_{MSE}, Rm')
    zlabel('Z_{MSE}, Rm')
    axis equal
    box on
    grid on
    hold on
    xlim([-2 2])
    ylim([-2 2])
    zlim([-2 2])
    for j = 1:size(FOV_mso,1)
        mapind = ceil(64*colormap(j)/max(colormap));
        quiver3(x(choose_ind(i)),y(choose_ind(i)),z(choose_ind(i)),...
            FOV_mso(j,1)/FOV,FOV_mso(j,2)/FOV,FOV_mso(j,3)/FOV,'Color',...
            [RGBmap.cmap(mapind, 1) RGBmap.cmap(mapind, 2) RGBmap.cmap(mapind, 3)]);
        %[1-colormap(j)/max(colormap) 1-colormap(j)/max(colormap) 1-colormap(j)/max(colormap)]);
    end
    view(90,0)
    
    %     % Stop signal
    %     pause(1)
    %     delete(h1)
    %     delete(h2)
    %     delete(h3)
end

end