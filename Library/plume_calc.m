function plume_calc(x1, epoch_d1,eflux_d1_cleaned_sum,swp_ind_d1,energy_d1,x,y,z,filename_mag, products_H)
TimeStart = x1(3);
TimeEnd = x1(4);

[amu, anu, c, q, m_p, m_O, m_O2, m_CO2, Rm] = constants;

verify = products_H.epoch >= x1(1) & products_H.epoch <= x1(2);
V_sw = mean(products_H.v_mso(verify,:),1);
V_sw = norm(V_sw)*1e3;  % in m/s

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
x = Rm*x(choose_ind);
y = Rm*y(choose_ind);
z = Rm*z(choose_ind);

[~,energy_ind] = max(squeeze(eflux(:,5,:)),[],1);
plume_energy = energy_d1(energy_ind,swp_ind,1,5);
plume_energy = plume_energy(:,1);

dl = 1e3*[diff(x), diff(y), diff(z)];   %in meters
a = dot(repmat(E_upstr,size(dl,1),1),dl, 2);    %in Volts


figure()
semilogy(epoch, plume_energy)
datetick('x')
xlabel('UT, HH:MM')
ylabel('Energy, eV')
hold on
semilogy(epoch(1:end-1), 5+cumsum(a)-min(cumsum(a)))
hold off


end