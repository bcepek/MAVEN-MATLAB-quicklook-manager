function [epoch_swia_spec, num_accum_swia,decom_flag_swia,spectra_counts_swia,...
          spectra_diff_en_fluxes_swia,geom_factor_swia,de_over_e_spectra_swia,...
          accum_time_spectra_swia,energy_spectra_swia,num_spec_swia] = read_variables_from_cdf_swia_spec(filename_swia_spec)


data = spdfcdfread(filename_swia_spec, 'variables', {
                                              'num_accum',...
                                              'decom_flag',...
                                              'spectra_counts',...
                                              'spectra_diff_en_fluxes',...
                                              'geom_factor',...
                                              'de_over_e_spectra',...
                                              'accum_time_spectra',...
                                              'energy_spectra',...
                                              'num_spec',...
                                              'epoch'
                                             });
epoch_swia_spec = data{10};
num_accum_swia = data{1};
decom_flag_swia = data{2};
spectra_counts_swia = data{3};
spectra_diff_en_fluxes_swia = data{4};
geom_factor_swia = data{5};
de_over_e_spectra_swia = data{6};
accum_time_spectra_swia = data{7};
energy_spectra_swia = data{8};
num_spec_swia = data{9};



