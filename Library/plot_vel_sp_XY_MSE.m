function xticks = plot_vel_sp_XY_MSE(TimeStart,mass_num,x,y,z,epoch,eflux,mass_arr,energy,swp_ind,quat_mso,phi,theta,mag_filename)
%(x1,mass_num,x,y,z,epoch_d1,eflux_d1_cleaned,mass_arr_d1,energy_d1,swp_ind_d1,magf_d1, quat_mso_d1,phi_d1,theta_d1,filename_mag)

search_time = TimeStart(3);
Rm = 3390;
aem = 1.66e-27;
q = 1.602e-19;

% Spacecraft velocity
dx = Rm*[diff(x), diff(y), diff(z)];
epoch1 = datevec(epoch);
epoch1(:,1:3) = [];
epoch1(:,1) = epoch1(:,1)*60*60;
epoch1(:,2) = epoch1(:,2)*60;
dt = diff(epoch1(:,1)+epoch1(:,2)+epoch1(:,3));
sc_vel = dx./dt;
sc_vel(end+1,:) = sc_vel(end,:);

% Calculation of upstream B and E vectors
mag = load(mag_filename);
[~, mag_timenid(1)] = min(abs(mag.mf_epoch - TimeStart(1)));
[~, mag_timenid(2)] = min(abs(mag.mf_epoch - TimeStart(2)));
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

h = figure();
for time_ind_add = 0:2
    
    mass_arr = 16*ones(size(mass_arr));
    phsdensity = zeros(size(eflux));
    for i = 1:size(eflux, 4)
        m_sq = aem*permute(squeeze(mass_arr(:, swp_ind(i)+1, :, :)), [2 1 3]);
        e_sq = permute(squeeze(energy(:, swp_ind(i)+1, :, :)), [2 1 3]);
        phsdensity(:, :, :, i) = 1e4*0.5*eflux(:, :, :, i).*(m_sq./(e_sq*q)).^2;
    end
    
    [~, time_ind] = min(abs(epoch - search_time));
    
    time_ind = time_ind+time_ind_add;
    
    epoch2 = epoch(time_ind);
    eflux1 = eflux(:,:,:, time_ind);
    phsdensity1 = phsdensity(:,:,:, time_ind);
    swp_ind1 = swp_ind(time_ind);
    quat_mso1 = quat_mso(time_ind, :);
    sc_vel1 = sc_vel(time_ind,:);
    sc_vel1 = (rt_M*(sc_vel1'))';   %rotate sc_vel to MSE
    
    % ВЫБИРАТЬ МАССУ ТУТ;
    timenum = 1;
    if(sum(sum(eflux1(:,:,mass_num,timenum))) ~=0 )
        vec_seq = zeros(64*32, 5);
        for enum = 1:size(eflux1,2)
            for dirnum = 1:size(eflux1, 1)
                if(eflux1(dirnum, enum, mass_num, timenum)~=0)
                    vec_num = (enum-1)*size(eflux1, 1) + dirnum;
                    if(dirnum>=17 && dirnum <=36)
                        vec_seq(vec_num, 5)=1;  %marking ionosphere
                    end
                    v = 1e-3*sqrt((2*q*energy(enum,swp_ind1(timenum)+1,dirnum,mass_num))/(aem*mass_arr(enum,swp_ind1(timenum)+1,dirnum,mass_num)));
                    vec_seq(vec_num, 1:4) = v*[cos(phi(enum,swp_ind1(timenum)+1,dirnum,mass_num))*cos(theta(enum,swp_ind1(timenum)+1,dirnum,mass_num)),...
                        sin(phi(enum,swp_ind1(timenum)+1,dirnum,mass_num))*cos(theta(enum,swp_ind1(timenum)+1,dirnum,mass_num)),...
                        sin(theta(enum,swp_ind1(timenum)+1,dirnum,mass_num)),...
                        (1/v)*(phsdensity1(dirnum,enum,mass_num,timenum))];
                end
            end
        end
        rows2delete = all(vec_seq(:,1:4)'==0)';
        vec_seq(rows2delete, :) = [];
    end
    vec_mso = quatrotate(quatinv(quat_mso1(end,:)), vec_seq(:,[1 2 3]));
    vec_mse = zeros(size(vec_mso));
    for i = 1:size(vec_mse, 1)
        vec_mse(i,:) = rt_M*(vec_mso(i,:)');
    end
    
    l = time_ind_add*3 + 1;
    scrsz=get(0,'ScreenSize');
    set(h,'OuterPosition',[0 scrsz(4)*0.01 scrsz(3)/2 scrsz(4)])
    
    subplot(3,3,l)
    
    %vec_seq(:,4) = rand(size(vec_seq,1),1);
    norm_coef = 476.444/max(vec_seq(:,4));
    ion_mark = vec_seq(:, 5)==1;
    scatter3(vec_mse(ion_mark,1)+double(repmat(sc_vel1(1),length(vec_mse(ion_mark,1)),1)),...
        vec_mse(ion_mark,2)+double(repmat(sc_vel1(2),length(vec_mse(ion_mark,2)),1)),...
        vec_mse(ion_mark,3)+double(repmat(sc_vel1(3),length(vec_mse(ion_mark,3)),1)),...
        vec_seq(ion_mark,4)*norm_coef,...
        'markeredgecolor', [0 0 0.5])%,'MarkerFaceColor',[0 0 0.5])
    alpha(0.5)
    hold on
    scatter3(vec_mse(~ion_mark,1)+double(repmat(sc_vel1(1),length(vec_mse(~ion_mark,1)),1)),...
        vec_mse(~ion_mark,2)+double(repmat(sc_vel1(2),length(vec_mse(~ion_mark,2)),1)),...
        vec_mse(~ion_mark,3)+double(repmat(sc_vel1(3),length(vec_mse(~ion_mark,3)),1)),...
        vec_seq(~ion_mark,4)*norm_coef,...
        'markeredgecolor', [0 0 0.5])%,'MarkerFaceColor',[0 0 0.5])
    alpha(0.5)
    
    %plot3(vec_seq(:,1), vec_seq(:,2), vec_seq(:,3), '.', 'markersize', 14)
    xlabel('Vx MSE, km/s')
    ylabel('Vy MSE, km/s')
    zlabel('Vz MSE, km/s')
    box on
    grid on
    axis equal
    xlim([-25 25])
    ylim([-25 25])
    zlim([-25 25])
    
    
    %mag_file = find_mag_file(search_time);
    B_ind = find(abs(mag.mf_epoch-epoch(time_ind)) == min(abs(mag.mf_epoch-epoch(time_ind))));
    B = [mag.Bx(B_ind), mag.By(B_ind), mag.Bz(B_ind)];
    B = B/sqrt(sum(B.^2,2));
    max_len = max(sqrt(sum(vec_mse.^2,2)));
    B = B*max_len/2;
    B = (rt_M*(B'))';   % rotate to MSE
    
    % Plot magnetic field vector
    hold on
    plot3([0 B(1)], [0 B(2)], [0 B(3)], 'linestyle', '-', 'color', 'red', 'linewidth', 2)
    hold off
    switch mass_num
        case 1
            mass_name = 'H^+';
        case 2
            mass_name = 'He^+';
        case 5
            mass_name = 'O^+';
        case 6
            mass_name = 'O_2^+';
    end
    view(0,90)
    
    
    subplot(3,3,l+1)
    scatter3(vec_mse(ion_mark,1)+double(repmat(sc_vel1(1),length(vec_mse(ion_mark,1)),1)),...
        vec_mse(ion_mark,2)+double(repmat(sc_vel1(2),length(vec_mse(ion_mark,2)),1)),...
        vec_mse(ion_mark,3)+double(repmat(sc_vel1(3),length(vec_mse(ion_mark,3)),1)),...
        vec_seq(ion_mark,4)*norm_coef,...
        'markeredgecolor', [0 0 0.5])%,'MarkerFaceColor',[0 0 0.5])
    alpha(0.5)
    hold on
    scatter3(vec_mse(~ion_mark,1)+double(repmat(sc_vel1(1),length(vec_mse(~ion_mark,1)),1)),...
        vec_mse(~ion_mark,2)+double(repmat(sc_vel1(2),length(vec_mse(~ion_mark,2)),1)),...
        vec_mse(~ion_mark,3)+double(repmat(sc_vel1(3),length(vec_mse(~ion_mark,3)),1)),...
        vec_seq(~ion_mark,4)*norm_coef,...
        'markeredgecolor', [0 0 0.5])%,'MarkerFaceColor',[0 0 0.5])
    alpha(0.5)
    xlabel('Vx MSE, km/s')
    ylabel('Vy MSE, km/s')
    zlabel('Vz MSE, km/s')
    box on
    grid on
    axis equal
    
    hold on
    plot3([0 B(1)], [0 B(2)], [0 B(3)], 'linestyle', '-', 'color', 'red', 'linewidth', 2)
    hold off
    view(0,0)
    
    xlim([-25 25])
    ylim([-25 25])
    zlim([-25 25])
    title([mass_name, ' ', datestr(epoch2)])
    
    subplot(3,3,l+2)
    scatter3(vec_mse(ion_mark,1)+double(repmat(sc_vel1(1),length(vec_mse(ion_mark,1)),1)),...
        vec_mse(ion_mark,2)+double(repmat(sc_vel1(2),length(vec_mse(ion_mark,2)),1)),...
        vec_mse(ion_mark,3)+double(repmat(sc_vel1(3),length(vec_mse(ion_mark,3)),1)),...
        vec_seq(ion_mark,4)*norm_coef,...
        'markeredgecolor', [0 0 0.5])%,'MarkerFaceColor',[0 0 0.5])
    alpha(0.5)
    hold on
    scatter3(vec_mse(~ion_mark,1)+double(repmat(sc_vel1(1),length(vec_mse(~ion_mark,1)),1)),...
        vec_mse(~ion_mark,2)+double(repmat(sc_vel1(2),length(vec_mse(~ion_mark,2)),1)),...
        vec_mse(~ion_mark,3)+double(repmat(sc_vel1(3),length(vec_mse(~ion_mark,3)),1)),...
        vec_seq(~ion_mark,4)*norm_coef,...
        'markeredgecolor', [0 0 0.5])%,'MarkerFaceColor',[0 0 0.5])
    alpha(0.5)
    xlabel('Vx MSE, km/s')
    ylabel('Vy MSE, km/s')
    zlabel('Vz MSE, km/s')
    box on
    grid on
    axis equal
    
    hold on
    plot3([0 B(1)], [0 B(2)], [0 B(3)], 'linestyle', '-', 'color', 'red', 'linewidth', 2)
    hold off
    view(90,0)
    xlim([-25 25])
    ylim([-25 25])
    zlim([-25 25])
end

end