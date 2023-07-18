function JxB_acceleration(x1, products_O2, epoch_d1_disp,...
    eflux_d1_cleaned, quat_mso_d1, theta_d1, phi_d1, swp_ind_d1)

fsize = 14;

if(x1(2) < x1(1))
    x1([1 2]) = x1([2 1]);
end
if(x1(4) < x1(3))
    x1([3 4]) = x1([4 3]);
end

mass_ind = 6;   %O2+

verify_pre = epoch_d1_disp > x1(1) & epoch_d1_disp < x1(2);
verify_post = epoch_d1_disp > x1(3) & epoch_d1_disp < x1(4);
[~, verify_pt] = min(abs(epoch_d1_disp - x1(5)));

[~, verify_pos] = min(abs(products_O2.epoch - x1(5)));
pos = products_O2.pos_sc_mso(verify_pos, :)/3390;

verify_magf = products_O2.epoch>=epoch_d1_disp(1) & products_O2.epoch<=epoch_d1_disp(end);
magf_mso = quatrotate(quatinv(quat_mso_d1), products_O2.magf(verify_magf,:));

B_pre = mean(magf_mso(verify_pre,:));
B_post = mean(magf_mso(verify_post,:));

[bin_ind, e_ind] = find(eflux_d1_cleaned(:,:,mass_ind,verify_pt) ==...
                max(max(eflux_d1_cleaned(:,:,mass_ind,verify_pt))));

theta = theta_d1(e_ind, swp_ind_d1(verify_pt)+1, bin_ind, mass_ind);
phi = phi_d1(e_ind, swp_ind_d1(verify_pt)+1, bin_ind, mass_ind);
x_sta = cos(theta)*cos(phi);
y_sta = cos(theta)*sin(phi);
z_sta = sin(theta);

pts_dir = quatrotate(quatinv(quat_mso_d1(verify_pt,:)), [x_sta, y_sta, z_sta]);

B_pre = B_pre/norm(B_pre);
B_post = B_post/norm(B_post);
pts_dir = pts_dir/norm(pts_dir);

n = cross(B_pre, B_post);
n = n/norm(n);

v_angle_with_norm = acos(dot(n, pts_dir));
if(v_angle_with_norm > pi/2)
    v_angle_with_norm = pi - v_angle_with_norm;
end
v_angle_with_mf_plane = pi/2 - v_angle_with_norm;

disp(v_angle_with_mf_plane*180/pi)

figure()

[X_sp, Y_sp, Z_sp] = sphere(100);
axes = surf(X_sp, Y_sp, Z_sp,'FaceColor',[1 0.41 0.16],'EdgeColor','none');
set(gca, 'fontsize', fsize)
alpha(1)
lightangle(90,0)
hold on

l1 = quiver3(pos(1)-B_pre(1),pos(2)-B_pre(2),pos(3)-B_pre(3),...
    B_pre(1),B_pre(2),B_pre(3),'off','color', 'red', 'LineWidth', 2);
l2 = quiver3(pos(1),pos(2),pos(3), B_post(1),B_post(2),B_post(3),'off','color', 'blue', 'LineWidth', 2);
l3 = quiver3(pos(1),pos(2),pos(3), pts_dir(1),pts_dir(2),pts_dir(3),'off','color', 'green', 'LineWidth', 2);
hold off
axis equal
box on
grid on
xlabel('X_{MSO}')
ylabel('Y_{MSO}')
zlabel('Z_{MSO}')

legend([l1,l2,l3], {'MF_{pre}','MF_{post}','{V_{O2^+}}'}, 'FontSize', fsize)

end