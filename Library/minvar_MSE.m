% The script performs the minimun variance analysis based on MAVEN/MAG
% measurements.
% 
% ====Key variables to be specified====
% timefrom - beginning of the analyzed interval
% timeto - end of the analyzed interval
% readpath - path to a folder with MAG data stored in *.mat files
%
% ====OUTPUT====
% Figure 1: Bx, By, Bz plot for the selected interval
% Figure 3: B godogram in (i, k) and (j, k) projections
% Figure 4: 3D godogram of B in minvar coordinates
% Figure 5: 3D godogram of B in MSO coordinates
% Console: Eigenvalues and eigenvectors of minvar matrix
%
%           Author: Sergey Shuvalov (shuvalovsergei@gmail.com,
%                                    +79153481805)
%           Last change date: 18.06.2018 16:06

function minvar_MSE(t,x,y,z)
% timefrom = datenum('2015-07-17 11:08:59');
% timeto = datenum('2015-07-17 11:09:23');

% readpath = '\\193.232.6.100\Data\Maven\data for MATLAB\MAG\';
% filename = dir([readpath 'mvn_mag_l2_*' datestr(t(3), 'yyyymmdd') '*.sts.mat']);
% if length(filename)<1
%     disp('MAG file not found')
%     return
% end
% load([readpath filename.name])

filename = find_mag_file(t(1));
load(filename);

B = [Bx, By, Bz];
epoch = mf_epoch;

choose_ind = epoch>t(1) & epoch<t(2);
Bsw = B(choose_ind, :);

choose_ind = find(epoch>t(3) & epoch<t(4));
epoch = epoch(choose_ind);
B = B(choose_ind, :);

figure
plot(epoch, B(:, 1), epoch, B(:, 2), epoch, B(:, 3))
hold on
plot(epoch, sqrt(sum(B.^2, 2)), 'color', 'black')
datetick('x','HH:MM:SS');
xlim([epoch(1) epoch(end)])
grid on
legend('Bx', 'By', 'Bz', '|B|')
ylabel('B, nT')

M = zeros(3);
for i = 1:3
    for j = i:3
        M(i, j) = mean(B(:, i).*B(:, j)) - mean(B(:, i))*mean(B(:, j));
        M(j, i) = M(i, j);
    end
end

[vec, val] = eig(M);
B_ = zeros(size(B));
for i = 1:size(B, 1)
    B_(i, :) = vec\B(i, :)';
end

% ===== Convert to MSE===============================================

% Calculation of upstream B and E vectors
B_upstr = mean(Bsw,1);
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

vec_MSE = rt_M*vec;

% ==== END MSE calc ===================================================

figure
subplot(1, 2, 1)
%===
p1 = [B_(1:end-1, 1), B_(1:end-1, 3)];
p2 = [B_(2:end, 1), B_(2:end, 3)];
dp = p2-p1;
quiver(p1(:, 1), p1(:, 2), dp(:, 1), dp(:, 2), 0)
xlabel('B_3, nT')
ylabel('B_1, nT')
hold on
%===
%plot(B_(:, 1), B_(:, 3))
hold on
plot(B_(1, 1), B_(1, 3), 'o', 'markerfacecolor', 'green', 'markersize', 12)
plot(B_(end, 1), B_(end, 3), 'o', 'markerfacecolor', 'red', 'markersize', 12)
hold off
axis equal
grid on
ylim_left = get(gca, 'ylim');
%xlim([min([B_(:, 1); B_(:, 2)]) max([B_(:, 1); B_(:, 2)])])

subplot(1, 2, 2)
%===
p1 = [B_(1:end-1, 2), B_(1:end-1, 3)];
p2 = [B_(2:end, 2), B_(2:end, 3)];xlabel('B, nT')
dp = p2-p1;
quiver(p1(:, 1), p1(:, 2), dp(:, 1), dp(:, 2), 0)
xlabel('B_2, nT')
ylabel('B_1, nT')
name = [datestr(t(3)) ' - ' datestr(t(4))];
name = name([1:23 36:43]);
title(name)
hold on
%===
%plot(B_(:, 2), B_(:, 3))
hold on
plot(B_(1, 2), B_(1, 3), 'o', 'markerfacecolor', 'green', 'markersize', 12)
plot(B_(end, 2), B_(end, 3), 'o', 'markerfacecolor', 'red', 'markersize', 12)
hold off
axis equal
grid on
ylim_right = get(gca, 'ylim');
ylim_fig3 = [min([ylim_left(1) ylim_right(1)]) max([ylim_left(2) ylim_right(2)])];
ylim(ylim_fig3)
subplot(1, 2, 1)
ylim(ylim_fig3)
%xlim([min([B_(:, 1); B_(:, 2)]) max([B_(:, 1); B_(:, 2)])])

figure
plot3(B_(:, 1), B_(:, 2), B_(:, 3))
hold on
plot3(B_(1, 1), B_(1, 2), B_(1, 3), 'o', 'markerfacecolor', 'green', 'markersize', 12)
plot3(B_(end, 1), B_(end, 2), B_(end, 3), 'o', 'markerfacecolor', 'red', 'markersize', 12)
xlabel('i')
ylabel('j')
zlabel('k')
hold off
axis equal
grid on

figure
plot3(B(:, 1), B(:, 2), B(:, 3))
hold on
plot3(B(1, 1), B(1, 2), B(1, 3), 'o', 'markerfacecolor', 'green', 'markersize', 12)
plot3(B(end, 1), B(end, 2), B(end, 3), 'o', 'markerfacecolor', 'red', 'markersize', 12)
xlabel('X_{MSO}')
ylabel('Y_{MSO}')
zlabel('Z_{MSO}')
hold off
axis equal
grid on
hold on
pos_sc_mso = mean([x y z]);
pos_sc_mso = max(max(B))*pos_sc_mso/norm(pos_sc_mso);
meanB = mean(B);
quiver3(meanB(1),meanB(2),meanB(3),pos_sc_mso(1),pos_sc_mso(2),pos_sc_mso(3))
hold off

disp('Eigenvalues are:')
disp(val)
disp('Eigenvectors MSO are:')
disp(vec)
%annotation('textbox', [0.007142857142857,0.55952380952381,0.253571428571429,0.428571428571433]...
annotation('textbox', [0.0071,0.2238,0.2536,0.7643]...
    , 'String', ["Eigenvectors MSO are:" num2str(vec(1,:), '%10.2f')...
    num2str(vec(2,:), '%10.2f') num2str(vec(3,:), '%10.2f') newline...
    "Eigenvectors MSE are:" num2str(vec_MSE(1,:), '%10.2f')...
    num2str(vec_MSE(2,:), '%10.2f') num2str(vec_MSE(3,:), '%10.2f') newline...
    "Eigenvalues are:" num2str(val(1,:), '%10.2f')...
    num2str(val(2,:), '%10.2f') num2str(val(3,:), '%10.2f')])

angle_ = acos(sum(pos_sc_mso'.*vec(:,1))/norm(pos_sc_mso))*180/pi;
disp(num2str(angle_))
end