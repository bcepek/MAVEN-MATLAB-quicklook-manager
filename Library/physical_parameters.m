function [epoch, x, y, z, altit, R_cyl_O, SZA, Bx, By, Bz, B,... 
          n_p,  v_p,  v_p_x,  v_p_y,  v_p_z, T_p,...
          n_O,  v_O,  v_O_x,  v_O_y,  v_O_z, T_O,...
          n_O2, v_O2, v_O2_x, v_O2_y, v_O2_z, T_O2,...
          flux_p, flux_O, flux_O2, magnetic_to_kinetic_en_ratio,...
          att_ind, swp_ind] = physical_parameters(products_H, products_O, products_O2)

[~, ~, ~, ~, m_p, m_O, m_O2, ~, Rm] = constants;

epoch = products_O.epoch;
att_ind = products_O.att_ind;

x = products_O.pos_sc_mso(:,1);
y = products_O.pos_sc_mso(:,2);
z = products_O.pos_sc_mso(:,3);

Bx = products_O.magf(:, 1);
By = products_O.magf(:, 2);
Bz = products_O.magf(:, 3);
B = sqrt(Bx.^2 + By.^2 + Bz.^2);

n_p = products_H.concentration;
v_MSO_p = products_H.v_mso;
v_p = sqrt(v_MSO_p(:, 1).^2 + v_MSO_p(:, 2).^2 + v_MSO_p(:, 3).^2);
v_p_x = v_MSO_p(:, 1);
v_p_y = v_MSO_p(:, 2);
v_p_z = v_MSO_p(:, 3);
T_p = products_H.temp;

n_O = products_O.concentration;
v_MSO_O = products_O.v_mso;
v_O = sqrt(v_MSO_O(:, 1).^2 + v_MSO_O(:, 2).^2 + v_MSO_O(:, 3).^2);
v_O_x = v_MSO_O(:, 1);
v_O_y = v_MSO_O(:, 2);
v_O_z = v_MSO_O(:, 3);
T_O = products_O.temp;

n_O2 = products_O2.concentration;
v_MSO_O2 = products_O2.v_mso;
v_O2 = sqrt(v_MSO_O2(:, 1).^2 + v_MSO_O2(:, 2).^2 + v_MSO_O2(:, 3).^2);
v_O2_x = v_MSO_O2(:, 1);
v_O2_y = v_MSO_O2(:, 2);
v_O2_z = v_MSO_O2(:, 3);
T_O2 = products_O2.temp;

swp_ind = products_H.swp_ind;

flux_p = n_p*10^5.*v_p;
flux_O = n_O*10^5.*v_O;
flux_O2 = n_O2*10^5.*v_O2;
magnetic_to_kinetic_en_ratio = (B.^2/(8*pi)*10^-10) ./ ((n_p.*v_p.^2*m_p + n_O.*v_O.^2*m_O + n_O2.*v_O2.^2*m_O2)/2*10^10) ;


altit = sqrt(x.^2 + y.^2 + z.^2) - Rm;
R_cyl_O = sqrt(y.^2 + z.^2);
SZA = 180 / pi * atan ( R_cyl_O./x);
SZA(x < 0) = SZA(x < 0)  + 180;


end