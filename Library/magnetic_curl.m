function magnetic_curl(time_from, time_end, windowpos)

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

% ====Calculating Magnetic curl START=======================
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
% ====Calculating Magnetic curl END=========================

% ====Projecting magnetic curl on MVN_Y axis=======
J_mvn = zeros(size(J));
for i = 1:length(J_mvn)
    J_mvn(i) = dot(J_vec(i,:), E_dir_mso(i, :));
end
% =================================================

h = figure();
windowpos(4) = floor(windowpos(4)/2);
set(h, 'position', windowpos);
yyaxis left
%plot(mf_epoch(1:end-1), smooth(J_mvn, 320))
plot(mf_epoch(1:end-1), smooth(J_mvn, 128))
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
end