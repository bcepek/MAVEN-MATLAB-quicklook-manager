function [epoch_swia,time_met_swia,time_unix_swia,atten_state_swia,telem_mode_swia,quality_flag_swia,...
         decom_flag_swia,density_swia,pressure_swia,velcoity_swia,velocity_mso_swia,...
         temperature_swia,temperature_mso_swia,pindex_swia,vindex_swia,tindex_swia,p_label_swia,...
         v_label_swia,t_label_swia,num_mom_swia] = read_variables_from_cdf_swia_mom(filename_swia)


data = spdfcdfread(filename_swia, 'variables', {'epoch',...
                                              'time_met',...
                                              'time_unix',...
                                              'atten_state',...
                                              'telem_mode',...
                                              'quality_flag',...
                                              'decom_flag',...
                                              'density',...
                                              'pressure',...
                                              'velocity',...
                                              'velocity_mso',...
                                              'temperature',...
                                              'temperature_mso',...
                                              'pindex',...
                                              'vindex',...
                                              'tindex',...
                                              'p_label',...
                                              'v_label',...
                                              't_label',...
                                              'num_mom',...
                                             });
epoch_swia = data{1};
time_met_swia = data{2};
time_unix_swia = data{3};
atten_state_swia = data{4};
telem_mode_swia = data{5};
quality_flag_swia = data{6};
decom_flag_swia = data{7};
density_swia = data{8};
pressure_swia = data{9};
velcoity_swia = data{10};
velocity_mso_swia = data{11};
temperature_swia = data{12};
temperature_mso_swia = data{13};
pindex_swia = data{14};
vindex_swia = data{15};
tindex_swia = data{16};
p_label_swia = data{17};
v_label_swia = data{18};
t_label_swia = data{19};
num_mom_swia = data{20};


