function AdditionalPlots_Callback(source, eventdata, epoch_d1, energy_d1, swp_ind_d1, eflux_d1_cleaned_sum, quat_mso_d1, magf_d1, eflux_d1_cleaned,...
    products_O, products_H, products_O2, products_mag,...
    altit_low, altit_high, x, y, z, epoch_svyspec, diff_en_fluxes,energy,mf_epoch,B,...
    eflux_d1,domega_d1,domega_d1_sum,choose_ind_d1,theta_d1,phi_d1,mass_arr_d1,filename_mag)
str = source.String;
val = source.Value;
switch str{val}
    case 'Zoom on time axis'
        n = 2;
        [x1, ~] = ginput(n);
        MVN_STA_MAG_parameters_interactive_func(x1(1), x1(2))
    case 'Energy Spectrum AVG 3 WINDOW'
        n = 1;
        [x1, ~] = ginput(n);
        plot_energy_spectrum_avg_3_wind (x1, epoch_d1, energy_d1, swp_ind_d1, eflux_d1_cleaned_sum)
    case 'Energy Spectrum AVG 5 WINDOW'
        n = 1;
        [x1, ~] = ginput(n);
        plot_energy_spectrum_avg_5_wind (x1, epoch_d1, energy_d1, swp_ind_d1, eflux_d1_cleaned_sum)
    case 'Minimum Variance'
        n = 2;
        [x1, ~] = ginput(n);
        minvar_3_14(x1(1), x1(2))
    case 'Orbit and Vectors'
        n = 2;
        [x1, ~] = ginput(n);
        TimeStart = x1(1);
        TimeEnd   = x1(2);
        plot_Orbits_and_vectors(products_O, products_H, products_O2, products_mag,...
            TimeStart, TimeEnd, altit_low, altit_high);
    case 'Energy Spectra Set'
        n = 2;
        [x1, ~] = ginput(n);
        plot_energy_spectra_set (x1, epoch_d1, energy_d1, swp_ind_d1, eflux_d1_cleaned_sum,epoch_svyspec,diff_en_fluxes,energy)
        
    case '3D Distribution Function H (low E)'
        n = 2;
        [x1, ~] = ginput(n);
        mass_num = 1;
        energies_to_plot = 17:32;
        h = d1_mercator(x1(1), x1(2), mass_num, epoch_d1, eflux_d1_cleaned, energy_d1, swp_ind_d1, magf_d1, quat_mso_d1, energies_to_plot);
        
    case '3D Distribution Function O (low E)'
        n = 2;
        [x1, ~] = ginput(n);
        mass_num = 5;    
        energies_to_plot = 17:32;
        h = d1_mercator(x1(1), x1(2), mass_num, epoch_d1, eflux_d1_cleaned, energy_d1, swp_ind_d1, magf_d1, quat_mso_d1, energies_to_plot);
        
    case '3D Distribution Function O2 (low E)'
        n = 2;
        [x1, ~] = ginput(n);
        mass_num = 6;
        energies_to_plot = 17:32;
        h = d1_mercator(x1(1), x1(2), mass_num, epoch_d1, eflux_d1_cleaned, energy_d1, swp_ind_d1, magf_d1, quat_mso_d1, energies_to_plot);
        
    case '3D Distribution Function H (high E)'
        n = 2;
        [x1, ~] = ginput(n);
        mass_num = 1;
        energies_to_plot = 1:16;
        h = d1_mercator(x1(1), x1(2), mass_num, epoch_d1, eflux_d1_cleaned, energy_d1, swp_ind_d1, magf_d1, quat_mso_d1, energies_to_plot);
        
    case '3D Distribution Function O (high E)'
        n = 2;
        [x1, ~] = ginput(n);
        mass_num = 5;    
        energies_to_plot = 1:16;
        h = d1_mercator(x1(1), x1(2), mass_num, epoch_d1, eflux_d1_cleaned, energy_d1, swp_ind_d1, magf_d1, quat_mso_d1, energies_to_plot);
        
    case '3D Distribution Function O2 (high E)'
        n = 2;
        [x1, ~] = ginput(n);
        mass_num = 6;
        energies_to_plot = 1:16;
        h = d1_mercator(x1(1), x1(2), mass_num, epoch_d1, eflux_d1_cleaned, energy_d1, swp_ind_d1, magf_d1, quat_mso_d1, energies_to_plot); 
    case '3D Distribution Function H (average)'
        n = 2;
        [x1, ~] = ginput(n);
        mass_num = 1;
        h = d1_mercator_av(x1(1), x1(2), mass_num, epoch_d1, eflux_d1_cleaned, energy_d1, swp_ind_d1);
    case '3D Distribution Function O (average)'
        n = 2;
        [x1, ~] = ginput(n);
        mass_num = 5;
        h = d1_mercator_av(x1(1), x1(2), mass_num, epoch_d1, eflux_d1_cleaned, energy_d1, swp_ind_d1);
    case '3D Distribution Function O2 (average)'
        n = 2;
        [x1, ~] = ginput(n);
        mass_num = 6;
        h = d1_mercator_av(x1(1), x1(2), mass_num, epoch_d1, eflux_d1_cleaned, energy_d1, swp_ind_d1);
        
    case 'View dayside orbit'
        plot_dayside_orbit(epoch_d1, x, y, z);
        
    case 'Plot sub-satellite point'
        [~,~,~] = plot_subsat_point(epoch_d1(1), epoch_d1(end));
        
    case 'Velocity-density ratio'
        n = 4;
        [x1, ~] = ginput(n);
        plot_densities_velosity_ratio(x1(1),x1(2),x1(3),x1(4), products_O, products_H, products_O2,altit_low, altit_high);
        
    case 'Wavelet analysis and velocity'
        n = 2;
        [x1, ~] = ginput(n);
        plot_wavelet_B(xticks,mf_epoch,B,x1(1),x1(2), products_O, products_H, products_O2,altit_low, altit_high);
        
    case 'Velocity space XYZ H+'
        n = 1;
        [x1, ~] = ginput(n);
        mass_num = 1;
        plot_vel_sp_XY(x1,mass_num,x,y,z,epoch_d1,eflux_d1_cleaned,mass_arr_d1,energy_d1,swp_ind_d1,magf_d1,...
                       quat_mso_d1,phi_d1,theta_d1,filename_mag);
       
    case 'Velocity space XYZ O+'
        n = 1;
        [x1, ~] = ginput(n);
        mass_num = 5;
        plot_vel_sp_XY(x1,mass_num,x,y,z,epoch_d1,eflux_d1_cleaned,mass_arr_d1,energy_d1,swp_ind_d1,magf_d1,...
                       quat_mso_d1,phi_d1,theta_d1,filename_mag);
        
    case 'Velocity space XYZ O2+'
        n = 1;
        [x1, ~] = ginput(n);
        mass_num = 6;
        plot_vel_sp_XY(x1,mass_num,x,y,z,epoch_d1,eflux_d1_cleaned,mass_arr_d1,energy_d1,swp_ind_d1,magf_d1,...
                       quat_mso_d1,phi_d1,theta_d1,filename_mag);
    
    case 'Pressure analysis'
        n = 2;
        [x1, ~] = ginput(n);
        plot_energies(x1(1),x1(2), xticks, products_O, products_H, products_O2, products_mag, altit_low, altit_high)
        
    case 'FOV bin direction H+'
        n = 2;
        [x1, ~] = ginput(n);
        mass_num = 1;
        bin_dir(x1(1), x1(2), mass_num,epoch_d1,eflux_d1_cleaned,swp_ind_d1,quat_mso_d1,energy,theta_d1,phi_d1,x,y,z)
        
    case 'FOV bin direction O+'
        n = 2;
        [x1, ~] = ginput(n);
        mass_num = 5;
        bin_dir(x1(1), x1(2), mass_num,epoch_d1,eflux_d1_cleaned,swp_ind_d1,quat_mso_d1,energy,theta_d1,phi_d1,x,y,z)
        
    case 'FOV bin direction O2+'
        n = 2;
        [x1, ~] = ginput(n);
        mass_num = 6;
        bin_dir(x1(1), x1(2), mass_num,epoch_d1,eflux_d1_cleaned,swp_ind_d1,quat_mso_d1,energy,theta_d1,phi_d1,x,y,z)
        
end

end