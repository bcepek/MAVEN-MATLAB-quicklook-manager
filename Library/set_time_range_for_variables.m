function [x, y, z, ... 
    n_p,  T_p, T_p_energy, v_x_p,  v_y_p,  v_z_p,  v_p,...
    n_O,  T_O, T_O_energy, v_x_O,  v_y_O,  v_z_O,  v_O,...
    n_O2, T_O2, T_O2_energy, v_x_O2, v_y_O2, v_z_O2, v_O2,...
    Blx, Bly, Blz, Bl,...
    altit, epoch,...
    mf_epoch, Bx, By, Bz, B, x_mag,y_mag,z_mag] = set_time_range_for_variables(products_O, products_H, products_O2, products_mag, TimeStart, TimeEnd, altit_low, altit_high )
    %valid_p] = set_time_range_for_variables(products_O, products_H, products_O2, products_mag, TimeStart, TimeEnd, altit_low, altit_high )

[~, ~, ~, ~, ~, ~, ~, ~, Rm] = constants;

epoch = products_O.epoch;
x = products_O.pos_sc_mso(:, 1);
y = products_O.pos_sc_mso(:, 2);
z = products_O.pos_sc_mso(:, 3);
altit = sqrt(x.^2 + y.^2 + z.^2) - Rm;

x = products_O.pos_sc_mso(( epoch >= TimeStart & epoch <= TimeEnd & altit >= altit_low & altit < altit_high),1)/Rm;
y = products_O.pos_sc_mso(( epoch >= TimeStart & epoch <= TimeEnd & altit >= altit_low & altit < altit_high),2)/Rm;
z = products_O.pos_sc_mso(( epoch >= TimeStart & epoch <= TimeEnd & altit >= altit_low & altit < altit_high),3)/Rm;


%********************* p *************
n_p = products_H.concentration(( epoch >= TimeStart & epoch <= TimeEnd & altit >= altit_low & altit < altit_high));
T_p = products_H.temp(( epoch >= TimeStart & epoch <= TimeEnd & altit >= altit_low & altit < altit_high));
T_p_energy = 0; % products_H.temp_from_energy(( epoch >= TimeStart & epoch <= TimeEnd & altit >= altit_low & altit < altit_high));
v_x_p = products_H.v_mso(( epoch >= TimeStart & epoch <= TimeEnd & altit >= altit_low & altit < altit_high),1);
v_y_p = products_H.v_mso(( epoch >= TimeStart & epoch <= TimeEnd & altit >= altit_low & altit < altit_high),2);
v_z_p = products_H.v_mso(( epoch >= TimeStart & epoch <= TimeEnd & altit >= altit_low & altit < altit_high),3);
v_p = sqrt (v_x_p.^2 + v_y_p.^2 + v_z_p.^2);
%valid_p = products_H.valid(( epoch >= TimeStart & epoch <= TimeEnd & altit >= altit_low & altit < altit_high));

%********************* O plus ***********
n_O = products_O.concentration(( epoch >= TimeStart & epoch <= TimeEnd & altit >= altit_low & altit < altit_high));
T_O = products_O.temp(( epoch >= TimeStart & epoch <= TimeEnd & altit >= altit_low & altit < altit_high));
T_O_energy = 0; % products_O.temp_from_energy(( epoch >= TimeStart & epoch <= TimeEnd & altit >= altit_low & altit < altit_high));
v_x_O = products_O.v_mso(( epoch >= TimeStart & epoch <= TimeEnd & altit >= altit_low & altit < altit_high),1);
v_y_O = products_O.v_mso(( epoch >= TimeStart & epoch <= TimeEnd & altit >= altit_low & altit < altit_high),2);
v_z_O = products_O.v_mso(( epoch >= TimeStart & epoch <= TimeEnd & altit >= altit_low & altit < altit_high),3);
v_O = sqrt (v_x_O.^2 + v_y_O.^2 + v_z_O.^2);


%********************* O2 *************
n_O2 = products_O2.concentration(( epoch >= TimeStart & epoch <= TimeEnd & altit >= altit_low & altit < altit_high));
T_O2 = products_O2.temp(( epoch >= TimeStart & epoch <= TimeEnd & altit >= altit_low & altit < altit_high));
T_O2_energy = 0; % products_O2.temp_from_energy(( epoch >= TimeStart & epoch <= TimeEnd & altit >= altit_low & altit < altit_high));
v_x_O2 = products_O2.v_mso(( epoch >= TimeStart & epoch <= TimeEnd & altit >= altit_low & altit < altit_high),1);
v_y_O2 = products_O2.v_mso(( epoch >= TimeStart & epoch <= TimeEnd & altit >= altit_low & altit < altit_high),2);
v_z_O2 = products_O2.v_mso(( epoch >= TimeStart & epoch <= TimeEnd & altit >= altit_low & altit < altit_high),3);
v_O2 = sqrt (v_x_O2.^2 + v_y_O2.^2 + v_z_O2.^2);

%******************** B field *****************
%***************  low temporal resolution from STATIC data *******
Blx = products_O.magf(( epoch >= TimeStart & epoch <= TimeEnd & altit >= altit_low & altit < altit_high),1);
Bly = products_O.magf(( epoch >= TimeStart & epoch <= TimeEnd & altit >= altit_low & altit < altit_high),2);
Blz = products_O.magf(( epoch >= TimeStart & epoch <= TimeEnd & altit >= altit_low & altit < altit_high),3);
Bl = sqrt(Blx.^2 + Bly.^2 + Blz.^2);

%***************  high temporal resolution from MAG data *******
mf_epoch = products_mag.mf_epoch;
Bx = products_mag.Bx(( mf_epoch >= TimeStart & mf_epoch <= TimeEnd));
By = products_mag.By(( mf_epoch >= TimeStart & mf_epoch <= TimeEnd));
Bz = products_mag.Bz(( mf_epoch >= TimeStart & mf_epoch <= TimeEnd));
x_mag = products_mag.pos_sc_mso(( mf_epoch >= TimeStart & mf_epoch <= TimeEnd),1)/Rm;
y_mag = products_mag.pos_sc_mso(( mf_epoch >= TimeStart & mf_epoch <= TimeEnd),2)/Rm;
z_mag = products_mag.pos_sc_mso(( mf_epoch >= TimeStart & mf_epoch <= TimeEnd),3)/Rm;
B = sqrt(Bx.^2 + By.^2 + Bz.^2);
mf_epoch = mf_epoch( mf_epoch >= TimeStart & mf_epoch <= TimeEnd);


altit_tmp = altit( epoch >= TimeStart & epoch <= TimeEnd & altit >= altit_low & altit < altit_high);
epoch = epoch( epoch >= TimeStart & epoch <= TimeEnd & altit >= altit_low & altit < altit_high);
altit = altit_tmp;

end
