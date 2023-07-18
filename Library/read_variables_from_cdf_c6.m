function [epoch_c6, eflux_c6, pos_sc_mso_c6, energy_c6, nenergy_c6, denergy_c6, theta_c6, phi_c6, domega_c6, nbins_c6,...
          swp_ind_c6, mass_arr_c6, nanode_c6, ndef_c6, quat_mso_c6, magf_c6,...
          sc_pot_c6, att_ind_c6, quality_flag_c6, valid_c6 ] = read_variables_from_cdf_c6(filename_c6)

%{
epoch_c6 = spdfcdfread(filename_c6, 'variables', 'epoch');
eflux_c6 = spdfcdfread(filename_c6, 'variables', 'eflux');
pos_sc_mso_c6 = spdfcdfread(filename_c6, 'variables', 'pos_sc_mso');
energy_c6 = spdfcdfread(filename_c6, 'variables', 'energy');
nenergy_c6 = spdfcdfread(filename_c6, 'variables', 'nenergy');
denergy_c6 = spdfcdfread(filename_c6, 'variables', 'denergy');
theta_c6 = spdfcdfread(filename_c6, 'variables', 'theta');
phi_c6 = spdfcdfread(filename_c6, 'variables', 'phi');
domega_c6 = spdfcdfread(filename_c6, 'variables', 'domega');
nbins_c6 = spdfcdfread(filename_c6, 'variables', 'nbins');
swp_ind_c6 = spdfcdfread(filename_c6, 'variables', 'swp_ind');
mass_arr_c6 = spdfcdfread(filename_c6, 'variables', 'mass_arr');
nanode_c6 = spdfcdfread(filename_c6, 'variables', 'nanode');
ndef_c6 = spdfcdfread(filename_c6, 'variables', 'ndef');
quat_mso_c6 = spdfcdfread(filename_c6, 'variables', 'quat_mso');
magf_c6 = spdfcdfread(filename_c6, 'variables', 'magf');
sc_pot_c6 = spdfcdfread(filename_c6, 'variables', 'sc_pot');
att_ind_c6 = spdfcdfread(filename_c6, 'variables', 'att_ind');
%}

data = spdfcdfread(filename_c6, 'variables', {'epoch',...
                                           'eflux',...
                                           'pos_sc_mso',...
                                           'energy',...
                                           'nenergy',...
                                           'denergy',...
                                           'theta',...
                                           'phi',...
                                           'domega',...
                                           'nbins',...
                                           'swp_ind',...
                                           'mass_arr',...
                                            'nanode',...
                                            'ndef',...
                                            'quat_mso',...
                                            'magf',...
                                            'sc_pot',...
                                            'att_ind',...
                                            'quality_flag',...
                                            'valid',...
                                            'tof_arr',...
                                            'twt_arr'});
epoch_c6        = data{1};
eflux_c6        = data{2};
pos_sc_mso_c6   = data{3};
energy_c6       = data{4};
nenergy_c6      = data{5};
denergy_c6      = data{6};
theta_c6        = data{7};
phi_c6          = data{8};
domega_c6       = data{9};
nbins_c6        = data{10};
swp_ind_c6      = data{11};
mass_arr_c6     = data{12};
nanode_c6       = data{13};
ndef_c6         = data{14};
quat_mso_c6     = data{15};
magf_c6         = data{16};
sc_pot_c6       = data{17};
att_ind_c6      = data{18};


%{
1 = electrostatic attenuation, 
2 = mechanical attenuation, 
3 = mechanical and electrostatic attenuation). 
%}

quality_flag_c6 = data{19};
%{
Quality flag,
Bit 0 – test pulser on
Bit 1 – diagnostic mode
Bit 2 - dead time correction >2 flag
Bit 3 – detector droop correction >2 flag
Bit 4 – dead time correction not at event time
Bit 5 – electrostatic attenuator problem
Bit 6 – attenuator change during accumulation
Bit 7 – mode change during accumulation
Bit 8 – LPW interference with data
Bit 9 – high background
Bit 10 – no background subtraction array
Bit 11 – missing spacecraft potential
Bit 12 – inflight calibration incomplete
Bit 13 – geometric factor problem
Bit 14 – ion suppression problem
Bit 15 – 0
%}

valid_c6 = data{20};
%{
Validity flag, 1 for valid data, 0 for invalid data
%}

tof_arr_c6 = data{21}; %Gives average TOF value for mass bins.
twt_arr_c6 = data{22}; %Gives number of TOF bins in a given mass bin. Used for normalizing a mass spectra.


theta_c6 = theta_c6*pi/180;
phi_c6 = phi_c6*pi/180;
magf_c6 = quatrotate(quatinv(quat_mso_c6), magf_c6);






