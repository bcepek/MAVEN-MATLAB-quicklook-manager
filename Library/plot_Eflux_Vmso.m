filename = '\\193.232.6.100\Data\Maven\Static\d1_32e4d16a8m\2016\01\mvn_sta_l2_d1-32e4d16a8m_20160101_v01_r05.cdf';

epoch = spdfcdfread(filename, 'variables', 'epoch');
eflux = spdfcdfread(filename, 'variables', 'eflux');
energy = spdfcdfread(filename, 'variables', 'energy');
theta = spdfcdfread(filename, 'variables', 'theta');
phi = spdfcdfread(filename, 'variables', 'phi');
quat_mso = spdfcdfread(filename, 'variables', 'quat_mso');
swp_ind = spdfcdfread(filename, 'variables', 'swp_ind');
mass_arr = spdfcdfread(filename, 'variables', 'mass_arr');

start_time = '19:20:00';
stop_time = '19:21:00';
%time 
difference = abs(datenum(datestr(epoch, 'HH:MM:SS')) - datenum(start_time));
start_num = find (difference == min(difference), 1);
difference = abs(datenum(datestr(epoch, 'HH:MM:SS')) - datenum(stop_time));
stop_num = find (difference == min(difference), 1);

theta = theta*pi/180;
phi = phi*pi/180;

tft = [start_num stop_num];

N = 100;
    
k_ind = 0;

for timenum=tft(1):tft(1)%tft(2) 
    eflux_V = zeros(N, N, N, 8);
    for m_num=1%[1,5,6] 
        
        k_ind = k_ind + 1;
        h = figure(k_ind);
        set(h, 'pos', [7 0 900 900])
        
        % find V_min and V_max
        %{
        V_x_min = realmax;
        V_x_max = -realmax;
        V_y_min = realmax;
        V_y_max = -realmax;
        V_z_min = realmax;
        V_z_max = -realmax;
        for enum=1:32
            for nphi=1:16
                for ntheta=1:4
                    tht = theta(enum, swp_ind(timenum)+1, ntheta, m_num);
                    ph = phi(enum, swp_ind(timenum)+1, 4*nphi - 3, m_num);
                    m = mass_arr(enum, swp_ind(timenum)+1, 4*nphi - 3 + ntheta - 1, m_num);
                    m = 1.66e-27*m;
                    E = energy(enum, swp_ind(timenum)+1, 4*nphi - 3 + ntheta - 1, m_num);
                    E = 1.6e-19*E;
                    V_x = sqrt(2*E/m)*cos(tht)*cos(ph);
                    V_y = sqrt(2*E/m)*cos(tht)*sin(ph);
                    V_z = sqrt(2*E/m)*sin(tht); 
                    
                    V_mso = quatrotate(quatinv(quat_mso(timenum,:)), [V_x, V_y, V_z]);
                    V_x = V_mso(1);
                    V_y = V_mso(2);
                    V_z = V_mso(3);
                    
                    if V_x<V_x_min
                        V_x_min = V_x;
                    end
                    if V_x>V_x_max
                        V_x_max = V_x;
                    end
                    
                    if V_y<V_y_min
                        V_y_min = V_y;
                    end
                                  
                    if V_y>V_y_max
                        V_y_max = V_y;
                    end
                    
                    if V_z<V_z_min
                        V_z_min = V_z;
                    end
                    if V_z>V_z_max
                        V_z_max = V_z;
                    end
                    
                end
            end
        end
        %}
        
        energy_range = squeeze(energy(:, swp_ind(timenum)+1, 1, m_num));
        theta_range = squeeze(theta(16, swp_ind(timenum)+1, 1:4, m_num));
        phi_range = squeeze(phi(16, swp_ind(timenum)+1, 1:4:64, m_num));
        
        E_max = max(energy_range);
        E_max = 1.6e-19*E_max;
        m = 1.66e-27;
        V_x_min = sqrt(2*E_max/m)*cos(min(min(abs(theta(:,swp_ind(timenum)+1, :, 1)))));
        
        dV_x = (V_x_max - V_x_min)/N;
        dV_y = (V_y_max - V_y_min)/N;
        dV_z = (V_z_max - V_z_min)/N;
             
        % fill in eflux_V(:,:,:,m_nim);
        for i=1:N
            for j=1:N
                for k=1:N
%                     tht = theta(enum, swp_ind(timenum)+1, ntheta, m_num);
%                     ph = phi(enum, swp_ind(timenum)+1, 4*nphi - 3, m_num);
%                     E = energy(enum, swp_ind(timenum)+1, 4*nphi - 3 + ntheta - 1, m_num);
%                     E = 1.6e-19*E;
%                     V_x = sqrt(2*E/m)*cos(tht)*cos(ph);
%                     V_y = sqrt(2*E/m)*cos(tht)*sin(ph);
%                     V_z = sqrt(2*E/m)*sin(tht);
%                     
%                     V_mso = quatrotate(quatinv(quat_mso(timenum,:)), [V_x, V_y, V_z]);
%                     V_x = V_mso(1);
%                     V_y = V_mso(2);
%                     V_z = V_mso(3);
%                     
%                     i = ceil((V_x-V_x_min)/dV_x);
%                     j = ceil((V_y-V_y_min)/dV_y);
%                     k = ceil((V_z-V_z_min)/dV_z);
%                     if i==0, i=1; end; if i==N+1, i=N; end 
%                     if j==0, j=1; end; if j==N+1, j=N; end
%                     if k==0, k=1; end; if k==N+1, k=N; end
                    %eflux_V(i,j,k,m_num) = eflux_V(i,j,k,m_num) + eflux(4*nphi - 3 + ntheta - 1, enum, m_num, timenum);
                    
                    V_x = V_x_min + 0.5*dV_x*(2*i-1);
                    V_y = V_y_min + 0.5*dV_y*(2*j-1);
                    V_z = V_z_min + 0.5*dV_z*(2*k-1);
                    
                    m = mass_arr(16, swp_ind(timenum)+1, 1, m_num);
                    m = m*1.66e-27;
                    V = sqrt(V_x^2 + V_y^2 +V_z^2);
                    E = 0.5*m*V^2/(1.6e-19);
                    tht = asin(V_z/V);  
                    %if (tht<-90)||(tht>90), disp('tht'); disp(V_z); disp(V); pause;  end
					
					% угол phi определяется как arctg(V_y/V_x), но он покрывает лишь половину промежутка от -180 до 180, поэтому нужны условия, описанные ниже
                    if (V_x==0)
                        if (V_y>0)
                            ph = pi/2;
                        elseif (V_y<0)
                            ph = -pi/2;
                        else
                            ph = 0;
                        end
                    else
                        ph = atan (V_y/V_x);
                        if ((V_y>=0)&&(V_x<0))
                            ph = ph + pi;
                        end
                        if ((V_y<0)&&(V_x<0))
                            ph = ph - pi;
                        end
                    end
                    %ph = ph*180/pi;
                    %if (ph<-180)||(ph>180), disp('ph'); disp(V_x); disp(V_y); pause;  end
                    
					% получившиеся углы theta и phi не равны точным значениям theta и phi в eflux. нужно найти ближайшие подходящие значения theta и phi.
					% для этого я нахожу минимаьную разницу между theta получившимися по Vx, Vy, Vz и теми, что есть в eflux. 
                    min_theta_difference = min(abs(theta_range - tht));
					% а потом нахожу индекс theta, при котором это достигается. 
                    tht_ind = find(abs(theta_range-tht)==min_theta_difference);
					% то же самое с phi 
                    min_phi_difference = min(abs(phi_range - ph));
                    phi_ind = find(abs(phi_range-ph)==min_phi_difference);
					% так же и для энергии
                    min_energy_difference = min(abs(energy_range - E));
                    en_ind = find(abs(energy_range-E)==min_energy_difference);  
                    
                    eflux_V(i,j,k,m_num) = eflux_V(i,j,k,m_num) + eflux(4*phi_ind - 3 + tht_ind - 1, en_ind, m_num, timenum);
                end
            end
        end
     

        V_x_label = (V_x_min+0.5*dV_x):dV_x:(V_x_max-0.5*dV_x); V_x_label=1e-3*V_x_label;
        V_y_label = (V_y_min+0.5*dV_y):dV_y:(V_y_max-0.5*dV_y); V_y_label=1e-3*V_y_label;
        V_z_label = (V_z_min+0.5*dV_z):dV_z:(V_z_max-0.5*dV_z); V_z_label=1e-3*V_z_label;

%         eflux_V_disp = squeeze(sum(eflux_V(:,:,1:N,m_num),3));
        %eflux_V_disp = squeeze(sum(eflux_V(:,:,round(N/2):round(N/2)+1,m_num),3));
        eflux_V_disp = squeeze(sum(eflux_V(:,:,round(N/2),m_num),3));
        
        verify = eflux_V_disp ~=0;
        colordata = [min(min(log10(eflux_V_disp(verify)))) max(max(log10(eflux_V_disp(verify))))];
        pcolor(V_x_label, V_y_label, log10(eflux_V_disp));
        if ~isempty(colordata)
            if colordata(1)~=colordata(2)
                caxis(colordata)
            end
        end    
        shading flat
        
        ttl = ['Mass = ', num2str(round(squeeze(mass_arr(16, swp_ind(timenum)+1, 1, m_num)))), 'a.u., Time = ', datestr(epoch(timenum))];
        title(ttl)
        xlabel('V_x, km/s (MSO)')
        ylabel('V_y, km/s (MSO)')

        bar_handle = colorbar;
        labels = get(bar_handle, 'yticklabel');
        barlabels = cell(size(labels, 1), 1);
        for i=1:size(labels, 1)
            barlabels{i} = ['10^{', labels{i}, '}'];
        end
        set(bar_handle, 'yticklabel', char(barlabels), 'FontWeight', 'bold', 'fontsize', 9)
        ylabel(bar_handle, 'Differential energy flux')
        
        %name = ['./Eflux_Vmso/Eflux_Vmso_', 'M=',num2str(round(squeeze(mass_arr(16, swp_ind(timenum)+1, 1, m_num)))), '_t=', num2str(datestr(epoch(timenum), 'yyyy.mm.dd')),'-',num2str(datestr(epoch(timenum), 'HH.MM.SS')), '.png'];
        %saveas (figure(k_ind), name);
        %delete(figure(k_ind))
    end
end