function plume_calc2(x1, epoch_d1,eflux_d1_cleaned_sum,swp_ind_d1,energy_d1,x,y,z,filename_mag, products_H)

% plume interval
TimeStart = x1(3);
TimeEnd = x1(4);

[amu, anu, c, q, m_p, m_O, m_O2, m_CO2, Rm] = constants;

% calculating average upstream Vsw
verify = products_H.epoch >= x1(1) & products_H.epoch <= x1(2);
V_sw = mean(products_H.v_mso(verify,:),1);
V_sw = norm(V_sw)*1e3;  % in m/s

% calculating upstream magnetic/electric field from STS
mag = load(filename_mag);
[~, mag_timenid(1)] = min(abs(mag.mf_epoch - x1(1)));
[~, mag_timenid(2)] = min(abs(mag.mf_epoch - x1(2)));
B_upstr = mean([mag.Bx(mag_timenid(1):mag_timenid(2))...
    mag.By(mag_timenid(1):mag_timenid(2))...
    mag.Bz(mag_timenid(1):mag_timenid(2))]);
%B_upstr = B_upstr./sqrt(sum(B_upstr.^2, 2));
E_upstr = -cross([-V_sw 0 0], B_upstr)*1e-9;    % in V/m
%E_upstr = E_upstr./norm(E_upstr);

choose_ind = epoch_d1 >= TimeStart & epoch_d1 <= TimeEnd;


eflux = eflux_d1_cleaned_sum(:,:,choose_ind);
epoch = epoch_d1(choose_ind);
swp_ind = swp_ind_d1(choose_ind)+1;
    % coordinates in km
x = Rm*x(choose_ind);
y = Rm*y(choose_ind);
z = Rm*z(choose_ind);

    % calculating energy of maximum plume intensity
[~,energy_ind] = max(squeeze(eflux(:,5,:)),[],1);
plume_energy = energy_d1(energy_ind,swp_ind,1,5);
plume_energy = plume_energy(:,1);


for R = 400:20:500
voltage = zeros(size(x));
for i = 1:size(x,1)
    a = norm(E_upstr)^2;
    b = 2*(E_upstr*[x(i); y(i); z(i)]);
    c = x(i)^2 + y(i)^2 + z(i)^2 - (R+Rm)^2;
    D = b^2 - 4*a*c;
    if(D<0)
        voltage(i) = nan;
    else
        t(1) = (-b+sqrt(D))/(2*a);
        t(2) = (-b-sqrt(D))/(2*a);
        intersec = [E_upstr(1)*t(1) + x(i), E_upstr(2)*t(1) + y(i), E_upstr(3)*t(1) + z(i);...
                    E_upstr(1)*t(2) + x(i), E_upstr(2)*t(2) + y(i), E_upstr(3)*t(2) + z(i)];
        dist = [norm([x(i) y(i) z(i)]-intersec(1,:));...
                norm([x(i) y(i) z(i)]-intersec(2,:))];
        voltage(i) = 1000*min(dist)*norm(E_upstr);
    end
end

% p = polyfit(epoch, plume_energy, 2);
% y1 = polyval(p, epoch);

%for alpha = 0.003:0.001:0.01
figure()
semilogy(epoch, plume_energy)
datetick('x')
xlabel('UT, HH:MM')
ylabel('Energy, eV')
hold on
%semilogy(epoch, y1, 'color', 'red')
%semilogy(epoch(1:end-1), voltage(1:end-1)+A-min(A))
semilogy(epoch, voltage)
hold off
title(['R = ' num2str(R)])
legend('Plume measured energy', 'Plume predicted energy')

end
end