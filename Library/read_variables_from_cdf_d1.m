function [epoch_d1, eflux_d1, pos_sc_mso_d1, energy_d1, nenergy_d1, denergy_d1, theta_d1, phi_d1, domega_d1, nbins_d1,...
          swp_ind_d1, mass_arr_d1, nanode_d1, ndef_d1, quat_mso_d1, magf_d1,...
          sc_pot_d1, att_ind_d1, quality_flag_d1, valid_d1 ] = read_variables_from_cdf_d1(filename_d1)

%{
epoch_d1 = spdfcdfread(filename_d1, 'variables', 'epoch');
eflux_d1 = spdfcdfread(filename_d1, 'variables', 'eflux');
pos_sc_mso_d1 = spdfcdfread(filename_d1, 'variables', 'pos_sc_mso');
energy_d1 = spdfcdfread(filename_d1, 'variables', 'energy');
nenergy_d1 = spdfcdfread(filename_d1, 'variables', 'nenergy');
denergy_d1 = spdfcdfread(filename_d1, 'variables', 'denergy');
theta_d1 = spdfcdfread(filename_d1, 'variables', 'theta');
phi_d1 = spdfcdfread(filename_d1, 'variables', 'phi');
domega_d1 = spdfcdfread(filename_d1, 'variables', 'domega');
nbins_d1 = spdfcdfread(filename_d1, 'variables', 'nbins');
swp_ind_d1 = spdfcdfread(filename_d1, 'variables', 'swp_ind');
mass_arr_d1 = spdfcdfread(filename_d1, 'variables', 'mass_arr');
nanode_d1 = spdfcdfread(filename_d1, 'variables', 'nanode');
ndef_d1 = spdfcdfread(filename_d1, 'variables', 'ndef');
quat_mso_d1 = spdfcdfread(filename_d1, 'variables', 'quat_mso');
magf_d1 = spdfcdfread(filename_d1, 'variables', 'magf');
sc_pot_d1 = spdfcdfread(filename_d1, 'variables', 'sc_pot');
att_ind_d1 = spdfcdfread(filename_d1, 'variables', 'att_ind');
%}

data = spdfcdfread(filename_d1, 'variables', {'epoch',...
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
epoch_d1        = data{1};
eflux_d1        = data{2};
pos_sc_mso_d1   = data{3};
energy_d1       = data{4};
nenergy_d1      = data{5};
denergy_d1      = data{6};
theta_d1        = data{7};
phi_d1          = data{8};
domega_d1       = data{9};
nbins_d1        = data{10};
swp_ind_d1      = data{11};
mass_arr_d1     = data{12};
nanode_d1       = data{13};
ndef_d1         = data{14};
quat_mso_d1     = data{15};
magf_d1         = data{16};
sc_pot_d1       = data{17};
att_ind_d1      = data{18};


%{
1 = electrostatic attenuation, 
2 = mechanical attenuation, 
3 = mechanical and electrostatic attenuation). 
%}

quality_flag_d1 = data{19};
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

valid_d1 = data{20};
%{
Validity flag, 1 for valid data, 0 for invalid data
%}

tof_arr_d1 = data{21}; %Gives average TOF value for mass bins.
twt_arr_d1 = data{22}; %Gives number of TOF bins in a given mass bin. Used for normalizing a mass spectra.


theta_d1 = theta_d1*pi/180;
phi_d1 = phi_d1*pi/180;
magf_d1 = quatrotate(quatinv(quat_mso_d1), magf_d1);






