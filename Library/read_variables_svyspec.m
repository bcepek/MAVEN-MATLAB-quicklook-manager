function [epoch_svyspec,num_accum,counts,diff_en_fluxes,weight_factor,...
          geom_factor, g_engy, de_over_e, accum_time, energy,num_spec] = read_variables_svyspec(filename, TimeStart, TimeEnd)
% read physical paramteres, measured by SWEA MAVEN 
% Written by Konstantkin Kim

epoch_svyspec = spdfcdfread(filename, 'variables','epoch');
num_accum = spdfcdfread(filename, 'variables','num_accum');
counts = spdfcdfread(filename, 'variables','counts');
diff_en_fluxes = spdfcdfread(filename, 'variables','diff_en_fluxes');
weight_factor = spdfcdfread(filename, 'variables','weight_factor');
geom_factor = spdfcdfread(filename, 'variables','geom_factor');
g_engy = spdfcdfread(filename, 'variables','g_engy');
de_over_e = spdfcdfread(filename, 'variables','de_over_e');
accum_time = spdfcdfread(filename, 'variables','accum_time');
energy = spdfcdfread(filename, 'variables','energy');
num_spec = spdfcdfread(filename, 'variables','num_spec');
% Отбор времен
choose_ind_svyspec = find(epoch_svyspec >= TimeStart & epoch_svyspec <= TimeEnd);

epoch_svyspec = epoch_svyspec(choose_ind_svyspec);
num_accum = num_accum(choose_ind_svyspec);
counts = counts(choose_ind_svyspec,:);
diff_en_fluxes = diff_en_fluxes(choose_ind_svyspec,:);

end 