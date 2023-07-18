function AdditionalPlots_Callback(source, eventdata, epoch_d1_disp, energy_d1, swp_ind_d1, eflux_d1_cleaned_sum, quat_mso_d1, magf_d1, eflux_d1_cleaned,...
    products_O, products_H, products_O2, products_mag,...
    altit_low, altit_high, x, y, z, pos_sc_mso_d1_disp, quat_mso_d1_disp,...
    mass_arr_d1, phi_d1, theta_d1, filename_mag, windowpos)
str = source.String;
val = source.Value;
switch str{val}
    case 'Zoom on time axis (2)'
        n = 2;
        [x1, ~] = ginput(n);
        %MVN_STA_MAG_parameters_interactive_func(x1(1), x1(2), windowpos)
        MVN_main(x1(1), x1(2));
        
    case 'Time stamp (1)'
        n = 1;
        [x1, ~] = ginput(n);
        time_stamp(x1)

    case 'Ion velocity MSO (0)'
        plot_ion_velocity_MSO(epoch_d1_disp, products_O, products_H, products_O2, windowpos);

    case 'n-T diagram (2)'
        n = 2;
        [x1, ~] = ginput(n);
        plot_n_T (x1(1), x1(2), products_O, products_H, products_O2);

    case 'Box (4)'
        n = 4;
        [x1, ~] = ginput(n);
        Box_draw(x1,products_H,products_O);

    case 'Energy Spectrum AVG 5 WINDOW (1)'
        n = 1;
        [x1, ~] = ginput(n);
        plot_energy_spectrum_avg_5_wind (x1, epoch_d1_disp, energy_d1, swp_ind_d1, eflux_d1_cleaned_sum)

    case 'Energy Spectra Set (2)'
        n = 2;
        [x1, ~] = ginput(n);
        plot_energy_spectra_set (x1, epoch_d1_disp, energy_d1, swp_ind_d1, eflux_d1_cleaned_sum)%,epoch_svyspec,diff_en_fluxes,energy)

    case 'Minimum Variance MSE (4)'
        n = 4;
        [x1, ~] = ginput(n);
        verify = epoch_d1_disp>x1(1) & epoch_d1_disp<x1(2);
        minvar_MSE(x1,x(verify),y(verify),z(verify))

    case 'Minimum Variance (2)'
        n = 2;
        [x1, ~] = ginput(n);
        verify = epoch_d1_disp>x1(1) & epoch_d1_disp<x1(2);
        minvar_MSO(x1,x(verify),y(verify),z(verify))

    case 'Orbit and Vectors (2)'
        n = 2;
        [x1, ~] = ginput(n);
        plot_Orbits_and_vectors(products_O, products_H, products_O2, products_mag,...
            x1, altit_low, altit_high);
    case 'Orbit and Vectors MSE (4)'
        n = 4;
        [x1, ~] = ginput(n);
        plot_Orbits_and_vectors_MSE(products_O, products_H, products_O2, products_mag,...
            x1, altit_low, altit_high);

    case '3D Distribution Function H (low E) (2)'
        n = 2;
        [x1, ~] = ginput(n);
        mass_num = 1;
        energies_to_plot = 17:32;
        d1_hammer(x1(1), x1(2), mass_num, epoch_d1_disp, eflux_d1_cleaned, energy_d1, swp_ind_d1, magf_d1, quat_mso_d1, energies_to_plot);

    case '3D Distribution Function O (low E) (2)'
        n = 2;
        [x1, ~] = ginput(n);
        mass_num = 5;
        energies_to_plot = 17:32;
        d1_hammer(x1(1), x1(2), mass_num, epoch_d1_disp, eflux_d1_cleaned, energy_d1, swp_ind_d1, magf_d1, quat_mso_d1, energies_to_plot);

    case '3D Distribution Function O2 (low E) (2)'
        n = 2;
        [x1, ~] = ginput(n);
        mass_num = 6;
        energies_to_plot = 17:32;
        d1_hammer(x1(1), x1(2), mass_num, epoch_d1_disp, eflux_d1_cleaned, energy_d1, swp_ind_d1, magf_d1, quat_mso_d1, energies_to_plot);

    case '3D Distribution Function H (high E) (2)'
        n = 2;
        [x1, ~] = ginput(n);
        mass_num = 1;
        energies_to_plot = 1:16;
        d1_hammer(x1(1), x1(2), mass_num, epoch_d1_disp, eflux_d1_cleaned, energy_d1, swp_ind_d1, magf_d1, quat_mso_d1, energies_to_plot);

    case '3D Distribution Function O (high E) (2)'
        n = 2;
        [x1, ~] = ginput(n);
        mass_num = 5;
        energies_to_plot = 1:16;
        d1_hammer(x1(1), x1(2), mass_num, epoch_d1_disp, eflux_d1_cleaned, energy_d1, swp_ind_d1, magf_d1, quat_mso_d1, energies_to_plot);

    case '3D Distribution Function O2 (high E) (2)'
        n = 2;
        [x1, ~] = ginput(n);
        mass_num = 6;
        energies_to_plot = 1:16;
        d1_hammer(x1(1), x1(2), mass_num, epoch_d1_disp, eflux_d1_cleaned, energy_d1, swp_ind_d1, magf_d1, quat_mso_d1, energies_to_plot);

    case '3D Distribution Function H (average) (2)'
        n = 2;
        [x1, ~] = ginput(n);
        mass_num = 1;
        d1_mercator_av(x1(1), x1(2), mass_num, epoch_d1_disp, eflux_d1_cleaned, energy_d1, swp_ind_d1);
    case '3D Distribution Function O (average) (2)'
        n = 2;
        [x1, ~] = ginput(n);
        mass_num = 5;
        d1_mercator_av(x1(1), x1(2), mass_num, epoch_d1_disp, eflux_d1_cleaned, energy_d1, swp_ind_d1);
    case '3D Distribution Function O2 (average) (2)'
        n = 2;
        [x1, ~] = ginput(n);
        mass_num = 6;
        d1_mercator_av(x1(1), x1(2), mass_num, epoch_d1_disp, eflux_d1_cleaned, energy_d1, swp_ind_d1);

    case 'Velocity space XYZ H+ (1)'
        n = 1;
        [x1, ~] = ginput(n);
        mass_num = 1;
        plot_vel_sp_XY(x1,mass_num,x,y,z,epoch_d1_disp,eflux_d1_cleaned,mass_arr_d1,energy_d1,swp_ind_d1,...
            quat_mso_d1,phi_d1,theta_d1,filename_mag);

    case 'Velocity space XYZ H+ MSE (3)'
        n = 3;
        [x1, ~] = ginput(n);
        mass_num = 1;
        plot_vel_sp_XY_MSE(x1,mass_num,x,y,z,epoch_d1_disp,eflux_d1_cleaned,mass_arr_d1,energy_d1,swp_ind_d1,...
            quat_mso_d1,phi_d1,theta_d1,filename_mag);

    case 'Velocity space XYZ O+ (1)'
        n = 1;
        [x1, ~] = ginput(n);
        mass_num = 5;
        plot_vel_sp_XY(x1,mass_num,x,y,z,epoch_d1_disp,eflux_d1_cleaned,mass_arr_d1,energy_d1,swp_ind_d1,...
            quat_mso_d1,phi_d1,theta_d1,filename_mag);

    case 'Velocity space XYZ O+ MSE (3)'
        n = 3;
        [x1, ~] = ginput(n);
        mass_num = 5;
        plot_vel_sp_XY_MSE(x1,mass_num,x,y,z,epoch_d1_disp,eflux_d1_cleaned,mass_arr_d1,energy_d1,swp_ind_d1,...
            quat_mso_d1,phi_d1,theta_d1,filename_mag);

    case 'Velocity space XYZ O2+ (1)'
        n = 1;
        [x1, ~] = ginput(n);
        mass_num = 6;
        plot_vel_sp_XY(x1,mass_num,x,y,z,epoch_d1_disp,eflux_d1_cleaned,mass_arr_d1,energy_d1,swp_ind_d1,...
            quat_mso_d1,phi_d1,theta_d1,filename_mag);

    case 'Velocity space XYZ O2+ MSE (3)'
        n = 3;
        [x1, ~] = ginput(n);
        mass_num = 6;
        plot_vel_sp_XY_MSE(x1,mass_num,x,y,z,epoch_d1_disp,eflux_d1_cleaned,mass_arr_d1,energy_d1,swp_ind_d1,...
            quat_mso_d1,phi_d1,theta_d1,filename_mag);

    case 'View dayside orbit (0)'
        plot_dayside_orbit(epoch_d1_disp, x, y, z);

    case 'View dayside orbit MSE (3)'
        n = 3;
        [x1, ~] = ginput(n);
        plot_dayside_orbit_MSE(x1, products_H);

    case 'Plot sub-satellite point (0)'
        [~,~,~] = plot_subsat_point(epoch_d1_disp(1), epoch_d1_disp(end));

    case 'Orientation (0)'
        plot_STA_orientation(epoch_d1_disp, pos_sc_mso_d1_disp, quat_mso_d1_disp)

    case 'FOV bin direction H+ (2)'
        n = 2;
        [x1, ~] = ginput(n);
        mass_num = 1;
        bin_dir(x1(1), x1(2), mass_num,epoch_d1_disp,eflux_d1_cleaned,swp_ind_d1,quat_mso_d1,energy_d1,theta_d1,phi_d1,x,y,z)

    case 'FOV bin direction O+ (2)'
        n = 2;
        [x1, ~] = ginput(n);
        mass_num = 5;
        bin_dir(x1(1), x1(2), mass_num,epoch_d1_disp,eflux_d1_cleaned,swp_ind_d1,quat_mso_d1,energy_d1,theta_d1,phi_d1,x,y,z)

    case 'FOV bin direction O2+ (2)'
        n = 2;
        [x1, ~] = ginput(n);
        mass_num = 6;
        bin_dir(x1(1), x1(2), mass_num,epoch_d1_disp,eflux_d1_cleaned,swp_ind_d1,quat_mso_d1,energy_d1,theta_d1,phi_d1,x,y,z)

    case 'FOV bin direction H+ MSE (4)'
        n = 4;
        [x1, ~] = ginput(n);
        mass_num = 1;
        bin_dir_MSE(x1, mass_num,epoch_d1_disp,eflux_d1_cleaned,swp_ind_d1,quat_mso_d1,energy_d1,theta_d1,phi_d1,x,y,z,filename_mag)

    case 'FOV bin direction O+ MSE (4)'
        n = 4;
        [x1, ~] = ginput(n);
        mass_num = 5;
        bin_dir_MSE(x1, mass_num,epoch_d1_disp,eflux_d1_cleaned,swp_ind_d1,quat_mso_d1,energy_d1,theta_d1,phi_d1,x,y,z,filename_mag)

    case 'FOV bin direction O2+ MSE (4)'
        n = 4;
        [x1, ~] = ginput(n);
        mass_num = 6;
        bin_dir_MSE(x1, mass_num,epoch_d1_disp,eflux_d1_cleaned,swp_ind_d1,quat_mso_d1,energy_d1,theta_d1,phi_d1,x,y,z,filename_mag)

    case 'O, O2 angle with surface (0)'
        heavy_motion_through_boundary(products_H, products_O, products_O2, epoch_d1_disp);

    case 'Wavelet (0)'
        plot_wavelet_B(epoch_d1_disp(1), epoch_d1_disp(end))

    case 'JxB acceleration (5)'
        n = 5;
        [x1, ~] = ginput(n);
        JxB_acceleration(x1, products_O2, epoch_d1_disp, eflux_d1_cleaned, quat_mso_d1, theta_d1, phi_d1, swp_ind_d1);

    case 'Plume calc'
        n = 4;
        [x1, ~] = ginput(n);
        plume_calc(x1, epoch_d1_disp,eflux_d1_cleaned_sum,swp_ind_d1,energy_d1,x,y,z,filename_mag, products_H);

    case 'Plume calc 2'
        n = 4;
        [x1, ~] = ginput(n);
        plume_calc2(x1, epoch_d1_disp,eflux_d1_cleaned_sum,swp_ind_d1,energy_d1,x,y,z,filename_mag, products_H);

    case 'Plume calc 3'
        n = 2;
        [x1, ~] = ginput(n);
        plume_calc3(x1, epoch_d1_disp,eflux_d1_cleaned_sum,swp_ind_d1,energy_d1,x,y,z,filename_mag, products_H);

end

end