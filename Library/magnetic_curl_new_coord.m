function magnetic_curl_new_coord(x1, time_from, time_end, windowpos)
smooth_par = 128;
mu0 = 1.256637e-6;

mag_fname = find_mag_file(time_from);
mag = load(mag_fname);
verify = mag.mf_epoch > time_from & mag.mf_epoch < time_end;

mf_epoch = mag.mf_epoch(verify);
pos_sc_mso = mag.pos_sc_mso(verify, :);
Bx = mag.Bx(verify);
By = mag.By(verify);
Bz = mag.Bz(verify);

lpw_fname = find_lpw_we12(time_from);
lpw_data = spdfcdfread(lpw_fname, 'variables', 'data');
lpw_epoch = spdfcdfread(lpw_fname, 'variables', 'epoch');
verify = lpw_epoch < time_end & lpw_epoch > time_from;
lpw_epoch = lpw_epoch(verify);
lpw_data = lpw_data(verify);

E_dir_mso = zeros(size(mf_epoch, 1), 3);
lpw_vec = [0 1 0];
et = zeros(size(mf_epoch));
%pos_sc_mso = zeros(size(mag.pos_sc_mso(verify, :)));
for i = 1:length(mf_epoch)
    et(i) = cspice_str2et(datestr(mf_epoch(i), 'yyyy-mm-dd HH:MM:SS'));
    matrix = cspice_sxform('MAVEN_SPACECRAFT', 'MAVEN_MSO', et(i) );
    E_dir_mso(i, :) = (matrix(4:6, 4:6)*lpw_vec')';
    %[state, time] = cspice_spkezr('MAVEN', et(i), 'MAVEN_MSO', 'NONE', 'MARS');
    %pos_sc_mso(i, :) = state(1:3)';
end

% a = diff(pos_sc_mso(:, 1));
% verify = [1; find(a~=0)];
% pos_sc_mso_interp(:,1) = interp1(mf_epoch(verify), pos_sc_mso(verify, 1), ...
%                                  mf_epoch);
% pos_sc_mso_interp(:,2) = interp1(mf_epoch(verify), pos_sc_mso(verify, 2), ...
%                                  mf_epoch);
% pos_sc_mso_interp(:,3) = interp1(mf_epoch(verify), pos_sc_mso(verify, 3), ...
%                                  mf_epoch);
% dx = diff(pos_sc_mso_interp(:, 1));
% dy = diff(pos_sc_mso_interp(:, 2));
% dz = diff(pos_sc_mso_interp(:, 3));


dBx = diff(Bx);
dBy = diff(By);
dBz = diff(Bz);
dx = diff(pos_sc_mso(:,1));
dy = diff(pos_sc_mso(:,2));
dz = diff(pos_sc_mso(:,3));

J_vec = [dBz./dy - dBy./dz,...
         dBx./dz - dBz./dx,...
         dBy./dx - dBx./dy] * mu0;  % in pico-amperes/m^2
J = sqrt(sum(J_vec.^2, 2));

J_mvn = zeros(size(J));
for i = 1:length(J_mvn)
    J_mvn(i) = dot(J_vec(i,:), E_dir_mso(i, :));
end

% searching some Edir vector for MVA interval:
mva_startind = find(abs(mf_epoch-x1(1))==min(abs(mf_epoch-x1(1))));
mva_endind = find(abs(mf_epoch-x1(2))==min(abs(mf_epoch-x1(2))));
Edir_mva = E_dir_mso(mva_startind,:);
% searching mva vector:
B = [Bx(mva_startind:mva_endind), By(mva_startind:mva_endind), Bz(mva_startind:mva_endind)];
M = zeros(3);
for i = 1:3
    for j = i:3
        M(i, j) = mean(B(:, i).*B(:, j)) - mean(B(:, i))*mean(B(:, j));
        M(j, i) = M(i, j);
    end
end
[vec, ~] = eig(M);
mva_dir = vec(:,3);
% calculating Z vector:
Z_new = cross(Edir_mva, mva_dir);
Z_new = Z_new./norm(Z_new);
% calculating Y vector:
Y_new = cross(Z_new, Edir_mva);
% transform matrix:
transform_matrix = [Edir_mva', Y_new', Z_new'];
J_vec_new = zeros(size(J_vec));
for i = 1:size(J_vec_new, 1)
    J_vec_new(i, :) = (transform_matrix\J_vec(i, :)')';
end
B_new = zeros(size(Bx,1), 3);
for i = 1:length(Bx)
    B_new(i, :) = (transform_matrix\[Bx(i); By(i); Bz(i)])';
end

h = figure();
windowpos(4) = floor(windowpos(4)/2);
windowpos(4) = windowpos(4)*3;
set(h, 'position', windowpos);
subplot(3,1,1)
yyaxis left
%plot(mf_epoch(1:end-1), smooth(J_mvn, 320))
plot(mf_epoch(1:end-1), smooth(J_mvn, smooth_par))
ylim([-1e-5 1e-5])
ylabel('Current, pA/m^2')
datetick('x')
yyaxis right
plot(lpw_epoch, lpw_data)
ylabel('LPW U12 [V]')
lim = get(gca, 'ylim');
lim = max(abs(lim));
set(gca, 'ylim', [-lim lim])
yline(0)
xline(x1(1))
xline(x1(2))

subplot(3,1,2)
plot(mf_epoch, B_new)
ylabel('B in new c/s')
legend('x', 'y', 'z', 'AutoUpdate','off')
datetick('x')
grid on
xline(x1(1))
xline(x1(2))

subplot(3,1,3)
J_vec_new(:,1) = smooth(J_vec_new(:,1), smooth_par);
J_vec_new(:,2) = smooth(J_vec_new(:,2), smooth_par);
J_vec_new(:,3) = smooth(J_vec_new(:,3), smooth_par);
plot(mf_epoch(1:end-1), J_vec_new)
ylabel('J in new c/s')
legend('x', 'y', 'z', 'AutoUpdate','off')
datetick('x')
grid on
xline(x1(1))
xline(x1(2))

end