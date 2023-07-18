function bin_dir(TimeStart,TimeEnd,mass_num,epoch_d1,eflux_d1_cleaned,swp_ind_d1,quat_mso_d1,...
                          energy,theta,phi,x,y,z)
    RGBmap = load('newColormap.mat');
    RGBmap.cmap(1,:) = [0 0 0.5625];
    aem = 1.66*10e-27;                  
    q = 1.60218e-19;
    choose_ind = find(epoch_d1 >= TimeStart & epoch_d1 <= TimeEnd);
    %*****calc_sc_velocity*********
    % [~, ~, ~, ~, ~, ~, ~, ~, Rm] = constants;
    % dr = Rm*[diff(x) diff(y) diff(z)];
    % epoch1 = datevec(epoch_d1);
    % epoch1(:,1:3) = [];
    % epoch1(:,1) = epoch1(:,1)*60*60;
    % epoch1(:,2) = epoch1(:,2)*60;   
    % dt = diff(epoch1(:,1)+epoch1(:,2)+epoch1(:,3));
    % sc_vel = zeros(size(dr,1)+1,3);
    % sc_vel(1:end-1,:) = dr./dt;
    % sc_vel(end,:) = sc_vel(end-1,:);
    %***** END calc_sc_velocity*********
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
    alpha(0.8)
    lightangle(90,0)
    hold on
    subplot(1,3,2)
    surf(X_sp, Y_sp, Z_sp,'FaceColor',[1 0.41 0.16],'EdgeColor','none')
    alpha(0.8)
    lightangle(90,0)
    hold on
    subplot(1,3,3)
    surf(X_sp, Y_sp, Z_sp,'FaceColor',[1 0.41 0.16],'EdgeColor','none')
    alpha(0.8)
    lightangle(90,0)
    hold on
    
    for i = 1:2:length(choose_ind)
        
        % Нужный поток в момент времени
        eflux = eflux_d1_cleaned(:,:,mass_num,choose_ind(i));
        swp_ind = swp_ind_d1(choose_ind(i))+1;
        [bins, energies] = find(eflux>1e6);
        table = [bins, energies];
        table = sortrows(table,2, 'descend');
        bins = table(:,1);
        energies = table(:,2);

        theta_m = zeros(1, length(bins));
        phi_m = zeros(1, length(bins));
        FOV_x = zeros(1, length(bins));
        FOV_y = zeros(1, length(bins));
        FOV_z = zeros(1, length(bins));
        colormap = zeros(1, length(bins));
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
            FOV = max(sqrt(FOV_mso(:,1).^2+FOV_mso(:,2).^2+FOV_mso(:,3).^2));
           
    %% --- FIELD OF VIEW PROJECTION ON XY MSO --- %%
    subplot(1,3,1)
    [~,~,~,~] = plot_MPB_bowshock;
    xlabel('X, Rm')
    ylabel('Y, Rm')
    zlabel('Z, Rm')
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
                   [RGBmap.cmap(mapind, 1) RGBmap.cmap(mapind, 2) RGBmap.cmap(mapind, 3)],...
                   'linewidth', 1.5);
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
    xlabel('X, Rm')
    ylabel('Y, Rm')
    zlabel('Z, Rm')
    grid on 
    box on
    axis equal 
    hold on

    title(['FOV projection on XY, XZ and YZ planes in MSO '...
        datestr(TimeStart, 'yyyy-mm-dd HH:MM:SS')...
        ' - '...
        datestr(TimeEnd, 'yyyy-MM-dd HH:mm:ss')])

    % title(['FOV projection on XY, XZ and YZ planes in MSO '...
    %     string(datetime(TimeStart, 'convertfrom', 'datenum', 'Format', 'yyyy-MM-dd HH:mm:ss'))...
    %     ' - '...
    %     string(datetime(TimeEnd, 'convertfrom', 'datenum', 'Format', 'yyyy-MM-dd HH:mm:ss'))])
    xlim([-2 2])
    ylim([-2 2])
    zlim([-2 2])
        for j = 1:size(FOV_mso,1)
                   mapind = ceil(64*colormap(j)/max(colormap));
                   quiver3(x(choose_ind(i)),y(choose_ind(i)),z(choose_ind(i)),...
                   FOV_mso(j,1)/FOV,FOV_mso(j,2)/FOV,FOV_mso(j,3)/FOV,'Color',...
                   [RGBmap.cmap(mapind, 1) RGBmap.cmap(mapind, 2) RGBmap.cmap(mapind, 3)],...
                   'linewidth', 1.5);
                   %[1-colormap(j)/max(colormap) 1-colormap(j)/max(colormap) 1-colormap(j)/max(colormap)]);
        end
    view(0,0)
    
%% --- FIELD OF VIEW PROJECTION ON YZ MSO --- %%
    subplot(1,3,3)
    t1=(0:0.001:2*pi);
    plot_3 = plot(cos(t1),sin(t1),'k');  
    rotate(plot_3,[0 1 0], -90)
    hold on 
    
    xlabel('X, Rm')
    ylabel('Y, Rm')
    zlabel('Z, Rm')
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
                   [RGBmap.cmap(mapind, 1) RGBmap.cmap(mapind, 2) RGBmap.cmap(mapind, 3)],...
                   'LineWidth', 1.5);
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