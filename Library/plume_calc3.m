function plume_calc3(x1, epoch_d1,eflux_d1_cleaned_sum,swp_ind_d1,energy_d1,x,y,z,filename_mag, products_H)

[amu, anu, c, q, m_p, m_O, m_O2, m_CO2, Rm] = constants;

verify = products_H.epoch >= x1(1) & products_H.epoch <= x1(2);
V_sw = products_H.v_mso(verify,:)*1e3;  % in m/s;
sta_epoch = products_H.epoch(verify);

choose_ind = epoch_d1 >= x1(1) & epoch_d1 <= x1(2);

eflux = eflux_d1_cleaned_sum(:,:,choose_ind);
swp_ind = swp_ind_d1(choose_ind)+1;
    % coordinates in km
x = Rm*x(choose_ind)*1e3;
y = Rm*y(choose_ind)*1e3;
z = Rm*z(choose_ind)*1e3;

    % calculating energy of maximum plume intensity
[~,energy_ind] = max(squeeze(eflux(:,5,:)),[],1);
plume_energy = energy_d1(energy_ind,swp_ind,1,5);
plume_energy = plume_energy(:,1);

mag = load(filename_mag);
verify = mag.mf_epoch >= x1(1) & mag.mf_epoch <= x1(2);
mf_epoch = mag.mf_epoch(verify);
B_hc = [mag.Bx(verify), mag.By(verify), mag.Bz(verify)];
B_lc = zeros(length(sta_epoch)-2, 3);
for i = 2:length(sta_epoch)-1
    ltime = (sta_epoch(i-1)+sta_epoch(i))/2;
    rtime = (sta_epoch(i)+sta_epoch(i+1))/2;
    
    [~, lind] = min(abs(mf_epoch - ltime));
    [~, rind] = min(abs(mf_epoch - rtime));
    
    B_lc(i-1, :) = mean(B_hc(lind:rind,:),1);
end

V_sw([1 end], :) = [];
sta_epoch([1 end]) = [];
x([1 end]) = [];
y([1 end]) = [];
z([1 end]) = [];
plume_energy([1 end]) = [];
E = -cross(V_sw, B_lc)*1e-9;    % in V/m

voltage = zeros(size(sta_epoch));
voltage(1) = plume_energy(1);
for i = 1:size(E,1)-1
    dx = [x(i+1)-x(i), y(i+1)-y(i), z(i+1)-z(i)];
    dE = dot(E(i,:), dx);
    voltage(i+1) = voltage(i)+dE;
end

figure()
semilogy(sta_epoch, plume_energy)
datetick('x')
xlabel('UT, HH:MM')
ylabel('Energy, eV')
hold on
semilogy(sta_epoch, voltage)
hold off
legend('Plume measured energy', 'Plume predicted energy')

end