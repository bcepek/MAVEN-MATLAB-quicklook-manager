function minvar_mf (timefrom, timeto, mf_epoch, Bx, By, Bz)


% timefrom = datenum(timefrom);
% timeto = datenum(timeto);

% load('mf_2015-01-04.mat')
% 
% 
% mf_epoch = mf_epoch;
% 

B = [Bx, By, Bz];
choose_ind = find(mf_epoch>timefrom & mf_epoch<timeto);
mf_epoch = mf_epoch(choose_ind);
B = B(choose_ind, :);



figure
plot(mf_epoch, B(:, 1), mf_epoch, B(:, 2), mf_epoch, B(:, 3))
datetick('x','HH:MM:SS');
xlim([mf_epoch(1) mf_epoch(end)])
grid on

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
%xlim([min([B_(:, 1); B_(:, 2)]) max([B_(:, 1); B_(:, 2)])])

subplot(1, 2, 2)
%===
p1 = [B_(1:end-1, 2), B_(1:end-1, 3)];
p2 = [B_(2:end, 2), B_(2:end, 3)];xlabel('B, nT')
dp = p2-p1;
quiver(p1(:, 1), p1(:, 2), dp(:, 1), dp(:, 2), 0)
xlabel('B_2, nT')
ylabel('B_1, nT')
name = [datestr(timefrom) ' - ' datestr(timeto)];
name = name([13:23 36:43]);
title(name)
hold on
%===
%plot(B_(:, 2), B_(:, 3))
hold on
plot(B_(1, 2), B_(1, 3), 'o', 'markerfacecolor', 'white', 'markersize', 12)
plot(B_(end, 2), B_(end, 3), 'o', 'markerfacecolor', 'black', 'markersize', 12)
hold off
axis equal
grid on
%xlim([min([B_(:, 1); B_(:, 2)]) max([B_(:, 1); B_(:, 2)])])

figure
plot3(B_(:, 1), B_(:, 2), B_(:, 3))
hold on
plot3(B_(1, 1), B_(1, 2), B_(1, 3), 'o', 'markerfacecolor', 'green', 'markersize', 12)
plot3(B_(end, 1), B_(end, 2), B_(end, 3), 'o', 'markerfacecolor', 'red', 'markersize', 12)
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


end