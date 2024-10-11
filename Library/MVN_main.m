function [f, sub_handles] = MVN_main(TimeStart, TimeEnd)

%----!!! Please adjust this set of parameters as you wish !!!---
caxis_lims = []; % set two values as in example "caxis_lims = [2 4]" or delete values to use automatic limits
LineWidth = 1;
FontSize = 12;
plotsGap = 0.010;
GraphicsSmoothing = 'on';
subPlotAmount = 8;
%---------------------------------------------------

sub_handles = zeros(subPlotAmount, 1);

altit_low = 0; %350; % [km]
altit_high = Inf; %2600;  % [km]
% step = 5;
% step_text = 20;
% 
% [amu, anu, c, q, m_p, m_O, m_O2, m_CO2, Rm] = constants;


filename_d1 = find_CDF_file_d1_Version2(TimeStart);
if filename_d1 == -1
    disp("mvn_sta_d1 file is not found on local machine. Trying to download it from LASP SDC...")
    filename_d1 = download_sta_d1(TimeStart);
    disp("mvn_sta_d1 downloaded")
end

global epoch_d1 eflux_d1 pos_sc_mso_d1 energy_d1 nenergy_d1 denergy_d1 theta_d1 phi_d1 domega_d1 nbins_d1...
    swp_ind_d1 mass_arr_d1 nanode_d1 ndef_d1 quat_mso_d1 magf_d1...
    sc_pot_d1 att_ind_d1 quality_flag_d1
%sc_pot_d1 att_ind_d1 quality_flag_d1 valid_d1
[epoch_d1, eflux_d1, pos_sc_mso_d1, energy_d1, nenergy_d1, denergy_d1, theta_d1, phi_d1, domega_d1, nbins_d1,...
    swp_ind_d1, mass_arr_d1, nanode_d1, ndef_d1, quat_mso_d1, magf_d1,...
    sc_pot_d1, att_ind_d1, quality_flag_d1] = read_variables_from_cdf_d1(filename_d1);
%sc_pot_d1, att_ind_d1, quality_flag_d1, valid_d1 ] = read_variables_from_cdf_d1(filename_d1);

% mass_num = 1; H+
% mass_num = 5; O+
% mass_num = 6; O2+

choose_ind_d1 = find(epoch_d1 >= TimeStart & epoch_d1 <= TimeEnd);

if isempty(epoch_d1(choose_ind_d1))
    disp('ATTENTION!!!')
    disp(['There is no d1 data for time range ', datestr(TimeStart, 'yyyy-mm-dd HH:MM:SS'), ' - ', datestr(TimeEnd, 'yyyy-mm-dd HH:MM:SS')])
    disp(['Time range for d1 data is from ', datestr(epoch_d1(1)), ' to ', datestr(epoch_d1(end)) ])
    epoch_d1_disp = epoch_d1(choose_ind_d1);
else
    epoch_d1_disp = epoch_d1(choose_ind_d1);
end
%height = height(choose_ind);

eflux_d1_disp = eflux_d1(:, :, :, choose_ind_d1);
swp_ind_d1_disp = swp_ind_d1(choose_ind_d1);
quat_mso_d1_disp = quat_mso_d1(choose_ind_d1,:);
%valid_d1_disp = valid_d1(choose_ind_d1);
%quality_flag_d1_disp = quality_flag_d1(choose_ind_d1);
pos_sc_mso_d1_disp = pos_sc_mso_d1(choose_ind_d1, :);
magf_d1_disp = magf_d1(choose_ind_d1, :);
%sc_pot_d1_disp = sc_pot_d1(choose_ind_d1);
%att_ind_d1_disp = att_ind_d1(choose_ind_d1);


%---BEGIN--- Calculate STATIC full field of view from d1(domega_d1_sum) ------------

domega_d1_sum = zeros(nenergy_d1, size(eflux_d1_disp, 4), size(eflux_d1_disp, 3));
for en = 1:nenergy_d1
    domega_d1_sum(en, :, :) = sum(domega_d1(en, swp_ind_d1_disp +1, :, :),3); %swp_ind_d1+1
    
end
domega_d1_sum = permute(domega_d1_sum, [1 3 2]);

%---END--- Calculate STATIC full field of view from d1 (domega_d1_sum) ------------

%---BEGIN--- 3D Cleaning O+ from d1 product---------
eflux_d1_cleaned = zeros(size(eflux_d1_disp));
eflux_d1_cleaned(:, :, 1, :) = eflux_d1_disp(:, :, 1, :);
temp = eflux_d1_disp(:,:, 5,:) - 0.08*eflux_d1_disp(:,:, 1,:) >= 0;
eflux_O = eflux_d1_disp(:,:, 5,:);
eflux_p = eflux_d1_disp(:,:, 1,:);
eflux_O(temp) = eflux_O(temp) - 0.08*eflux_p(temp);
clear temp;
temp = eflux_d1_disp(:,:, 5,:) - 0.08*eflux_d1_disp(:,:, 1,:) < 0;
eflux_O(temp) = 0;
eflux_d1_cleaned(:,:, 5,:) = eflux_O;
%---END--- 3D Cleaning O+ from d1 product---------



%---BEGIN--- 3D Cleaning O2+ from d1 product ---------
clear temp
temp = eflux_d1_disp(:,:, 6,:) - 0.08*eflux_d1_disp(:,:, 1,:) >= 0;
eflux_O2 = eflux_d1_disp(:,:, 6,:);
eflux_O2(temp) = eflux_O2(temp) - 0.08*eflux_p(temp);
clear temp;
temp = eflux_d1_disp(:,:, 6,:) - 0.08*eflux_d1_disp(:,:, 1,:) < 0;
eflux_O2(temp) = 0;
eflux_d1_cleaned(:,:, 6,:) = eflux_O2;
%---END--- 3D Cleaning O2+ from d1 product ---------



%---BEGIN--- Integrating d1 product by angles ---------
eflux_d1_cleaned_sum = zeros(size(eflux_d1_disp, 2), size(eflux_d1_disp, 3), size(eflux_d1_disp, 4));
for time = 1:size(eflux_d1_cleaned, 4)
    eflux_d1_cleaned_sum( :, :, time) = squeeze(sum(squeeze(    eflux_d1_cleaned(:, :, :, time) ).*...
        permute(    squeeze(        domega_d1(:, swp_ind_d1_disp(time)+1, :, :)       ), [2, 1, 3]),...
        1));
end
eflux_d1_cleaned_sum = eflux_d1_cleaned_sum./domega_d1_sum;
%---END--- Integrating d1 product by angles ---------


caxis_flux = [eflux_d1_cleaned_sum(:, 1, :)  eflux_d1_cleaned_sum(:, 5, :)  eflux_d1_cleaned_sum(:, 6, :)];
%{
caxis_lims = [min(min(log10(caxis_flux(caxis_flux > 0)))) max(max(max(max(log10(caxis_flux)))))];  % this way doesn't take into outliers account, so the new one below
%}
% Delete outliers
[bins, bin_centers] = hist(log10(caxis_flux(caxis_flux > 0)), 1000);
bins( bins <= 10 ) = 0;
bin_centers = bin_centers(logical(bins));
if isempty(caxis_lims)
    caxis_lims = [min(bin_centers) max(bin_centers)];
end
%}

[products_H, products_O, products_O2] = find_STA_mat_file_cleaned(TimeStart);
if(any([isempty(products_H), isempty(products_O), isempty(products_O2)]))
    disp("at least one of pre-calculated moments of mvn_sta_d1 is missing. Pre-calculating now...")
    calc_sta_d1_moments_cleaned(filename_d1)
    disp("mvn_sta_d1 moments pre-calculated")
    [products_H, products_O, products_O2] = find_STA_mat_file_cleaned(TimeStart);
end

filename_mag = find_mag_file(TimeStart);
if(filename_mag == -1)
    disp("mag mat file not found. trying to download sts file from LASP SDC...")
    filename_mag_sts = download_mag_ss(TimeStart);
    if filename_mag_sts ~= -1
        disp("mag sts file downloaded. Converting to mat...")
        mag_sts2mat(filename_mag_sts)
        filename_mag = find_mag_file(TimeStart);
    end
end
if(filename_mag~=-1)
    products_mag = load (filename_mag);
else
    products_mag = [];
end


% if isempty(products_H) || isempty(products_O) || isempty(products_O2) || isempty(products_mag)
%     disp('There is not some product file')
%     return
% end


[x, y, z, ...
    n_p,  T_p, T_p_energy, v_x_p,  v_y_p,  v_z_p,  v_p,...
    n_O,  T_O, T_O_energy, v_x_O,  v_y_O,  v_z_O,  v_O,...
    n_O2, T_O2, T_O2_energy, v_x_O2, v_y_O2, v_z_O2, v_O2,...
    Blx, Bly, Blz, Bl,...
    altit, epoch,...
    mf_epoch, Bx, By, Bz, B] = set_time_range_for_variables(products_O, products_H, products_O2, products_mag, TimeStart, TimeEnd, altit_low, altit_high );
%valid_p] = set_time_range_for_variables(products_O, products_H, products_O2, products_mag, TimeStart, TimeEnd, altit_low, altit_high );


R_O = sqrt(x.^2 + y.^2+z.^2);
R_cyl_O = sqrt(y.^2 + z.^2);
SZA = 180 / pi * atan ( R_cyl_O./x);
SZA (x<0) = SZA (x<0)  + 180;


v_x_bulk = v_x_p;
v_y_bulk = v_y_p;
v_z_bulk = v_z_p;
v_bulk = sqrt(v_x_bulk.^2 + v_y_bulk.^2 + v_z_bulk.^2);


%*********************************
Ex = v_z_bulk.*Bly - v_y_bulk.*Blz;
Ey = v_x_bulk.*Blz - v_z_bulk.*Blx;
Ez = v_y_bulk.*Blx - v_x_bulk.*Bly;
%****** components of normal to external boundary of barrier ****


%% *************** ELECTRON DISTRIBUTION FUNCTIONS ******************%
filename_svyspec = find_svyspec(TimeStart);
if filename_svyspec==-1
    disp("mvn_swe_svyspec not found. Trying to download from LASP SDC...")
    filename_svyspec = download_swe_svyspec(TimeStart);
    %filename_svyspec = find_svyspec(TimeStart);
end
if(filename_svyspec ~= -1)
    [epoch_svyspec,num_accum,counts,diff_en_fluxes,weight_factor,...
        geom_factor, g_engy, de_over_e, accum_time, energy,num_spec] = read_variables_svyspec(filename_svyspec, TimeStart, TimeEnd);
end
%% --- SWIA ENERGY-TIME SPECTROGRAM AND PRESSURE --- %%
%Root for SWIA moments data
filename_swia_mom = find_CDF_file_swia_moments(TimeStart);
if filename_swia_mom==-1
    disp("mvn_swi_svymom not found. Trying to download from LASP SDC...")
    download_swi_svymom(TimeStart);
    filename_swia_mom = find_CDF_file_swia_moments(TimeStart);
end
% Read moments variables
if(length(filename_swia_mom) > 1)
    [epoch_swia,time_met_swia,time_unix_swia,atten_state_swia,telem_mode_swia,quality_flag_swia,...
        decom_flag_swia,density_swia,pressure_swia,velcoity_swia,velocity_mso_swia,...
        temperature_swia,temperature_mso_swia,pindex_swia,vindex_swia,tindex_swia,p_label_swia,...
        v_label_swia,t_label_swia,num_mom_swia] = read_variables_from_cdf_swia_mom(filename_swia_mom);
    
    choose_ind_swia = find(epoch_swia >= TimeStart & epoch_swia <= TimeEnd);
end

% Root for SWIA spectra data
filename_swia_spec = find_CDF_file_swia_spec(TimeStart);
if filename_swia_spec==-1
    disp("mvn_swi_svyspec not found. Trying to download from LASP SDC...")
    filename_swia_spec = download_swi_svyspec(TimeStart);
    if(filename_swia_spec ~= -1)
        filename_swia_spec = find_CDF_file_swia_spec(TimeStart);
    end
end

% Read moments variables
if(filename_swia_spec ~= -1)
    [epoch_swia_spec, num_accum_swia, decom_flag_swia, spectra_counts_swia,...
        spectra_diff_en_fluxes_swia, geom_factor_swia, de_over_e_spectra_swia,...
        accum_time_spectra_swia, energy_spectra_swia, num_spec_swia] = read_variables_from_cdf_swia_spec(filename_swia_spec);
end
%%

%***************** source parameters figure(1) *******************
f = figure();
f.GraphicsSmoothing = GraphicsSmoothing;
%set(f,'defaultLegendAutoUpdate','off', 'position', [9 49 944 948])
%set(f,'defaultLegendAutoUpdate','off', 'position', [-1279 57 1280 948])
set(f, 'position', [-1279 57 1280 948])



%******************* main settings for plot positions ****************

scrsz=get(0,'ScreenSize');
set(1,'OuterPosition',[0 scrsz(4)*0.01 scrsz(3)*3/4 scrsz(4)]) % first argument corresponds to figure number
%clf('reset')
topGap=0.04;
bottomGap=0.1;
%plotsGap=plotsGap;
%plotLeftGap = 0.09;
plotLeftGap = 0.12;
%plotLength = 0.79;
plotLength = 0.7;
plotHight=(1-topGap-bottomGap-(subPlotAmount-1)*plotsGap)/subPlotAmount;
%******************* main settings for plot positions ****************

%{
hControl    = uicontrol('Style','popupmenu',...
           'String',{'On', 'Off'},...
    'Units','normalized', 'Position',[0.75,0.96,0.05,0.03],...
    'Callback',{@Control_Callback});
%}


%---BEGIN--- plot 0 ---------------
xticks = plot_south_north_hemisph_label ( bottomGap,  plotsGap, plotLeftGap, plotLength, FontSize, epoch, z, altit);
xlim([TimeStart TimeEnd])
%---END--- plot 0 ---------------



plotNumber = 1;


%---BEGIN---plot 1 ---------------
%plotNumber = plotNumber + 1;
%[~, sub_handles(1)] = plot_Bx (plotNumber, bottomGap, plotHight, plotsGap, plotLeftGap, plotLength, FontSize, xticks, caxis_lims,...
if(~isempty(mf_epoch))
    [~, sub_handles(1)] = plot_B_components_and_magnitude (plotNumber, bottomGap, plotHight, plotsGap, plotLeftGap, plotLength, FontSize, xticks, caxis_lims,...
        epoch_d1_disp, nenergy_d1, energy_d1, swp_ind_d1_disp, eflux_d1_cleaned_sum,...
        mf_epoch, Bx, By, Bz, LineWidth, Blx, Bly, Blz,x,y,z);
    xlim([TimeStart TimeEnd])
end
%----END----plot 1 ---------------

%---BEGIN---plot 2 ---------------
plotNumber = plotNumber + 1;
if(~isempty(epoch))
    [~, sub_handles(2)] = plot_velocities (plotNumber, bottomGap, plotHight, plotsGap, plotLeftGap, plotLength, FontSize, xticks,...
        epoch,...
        v_p,  v_O,  v_O2);
    xlim([TimeStart TimeEnd])
end
%----END----plot 2 ---------------

%---BEGIN---plot 3 ---------------
plotNumber = plotNumber + 1;
if(~isempty(epoch))
    [~, sub_handles(3)] = plot_densities_and_ratios (plotNumber, bottomGap, plotHight, plotsGap, plotLeftGap, plotLength, FontSize, xticks,...
        epoch,...
        n_p,  n_O,  n_O2, LineWidth);
    xlim([TimeStart TimeEnd])
end
%----END----plot 3 ---------------

%---BEGIN---plot 4 ---------------
plotNumber = plotNumber + 1;
if(filename_svyspec ~= -1)
    [~, sub_handles(4)] = plot_electron_spectrogram(plotNumber, bottomGap, plotHight, plotsGap, plotLeftGap, plotLength, FontSize, xticks, caxis_lims,...
        epoch_svyspec,energy, diff_en_fluxes);
    xlim([TimeStart TimeEnd])
end

%----END----plot 4 ---------------

%---BEGIN---plot 5 ---------------
plotNumber = plotNumber + 1;
if(~isempty(epoch))
    [~, sub_handles(5)] = plot_O2_plus_spectrogram (plotNumber, bottomGap, plotHight, plotsGap, plotLeftGap, plotLength, FontSize, xticks, caxis_lims,...
        epoch_d1_disp,nenergy_d1, energy_d1, swp_ind_d1_disp, eflux_d1_cleaned_sum,...
        mf_epoch, Bx, By, Bz, n_O, n_O2);
    xlim([TimeStart TimeEnd])
end
%----END----plot 5 ---------------


%---BEGIN---plot 6 ---------------
plotNumber = plotNumber + 1;
if(~isempty(epoch))
    [~, sub_handles(6)] = plot_O_plus_spectrogram (plotNumber, bottomGap, plotHight, plotsGap, plotLeftGap, plotLength, FontSize, xticks, caxis_lims,...
        epoch_d1_disp,nenergy_d1, energy_d1, swp_ind_d1_disp, eflux_d1_cleaned_sum);
    xlim([TimeStart TimeEnd])
end
%----END----plot 6 ---------------


%
%---BEGIN---plot 7 ---------------
plotNumber = plotNumber + 1;
if(~isempty(epoch))
    % [~] = plot_proton_spectrogram (plotNumber, bottomGap, plotHight, plotsGap, plotLeftGap, plotLength, FontSize, xticks, caxis_lims,...
    %     epoch_d1_disp, nenergy_d1, energy_d1, swp_ind_d1_disp, eflux_d1_cleaned_sum, altit, SZA,...
    %     mf_epoch, B);
    [~, sub_handles(7)] = plot_proton_spectrogram_withPressure (plotNumber, bottomGap, plotHight, plotsGap, plotLeftGap, plotLength, FontSize, xticks, caxis_lims,...
        epoch_d1_disp, nenergy_d1, energy_d1, swp_ind_d1_disp, eflux_d1_cleaned_sum, altit, SZA,...
        n_p, v_x_p);
    xlim([TimeStart TimeEnd])
end
%----END----plot 7 ---------------

%---BEGIN---plot 8 ---------------
plotNumber = plotNumber + 1;
if(length(filename_swia_mom) > 1 && filename_swia_spec ~= -1)
    [~, sub_handles(8)] = plot_proton_spectrogram_swia(plotNumber, bottomGap, plotHight, plotsGap, plotLeftGap, plotLength, FontSize, xticks, caxis_lims,...
        epoch_swia_spec,energy_spectra_swia, spectra_diff_en_fluxes_swia,density_swia,velocity_mso_swia,choose_ind_swia);
    xlim([TimeStart TimeEnd])
end
%----END----plot 8 ---------------


plot_title = strcat('MAVEN', {' '}, datestr(TimeStart,'yyyy-mm-dd'), {' '}, datestr(TimeStart,'HH:MM:SS'),{' - '},datestr(TimeEnd,'HH:MM:SS'));
title(plot_title)

% %----!!! Please adjust this set of parameters as you wish !!!---
% caxis_lims = []; % set two values as in example "caxis_lims = [2 4]" or delete values to use automatic limits
% LineWidth = 1;
% FontSize = 12;
% plotsGap = 0.010;
% GraphicsSmoothing = 'on';
% subPlotAmount = 4;
% f1 = figure(2);
% f1.GraphicsSmoothing = GraphicsSmoothing;
% set(f1,'defaultLegendAutoUpdate','off', 'position', [9 49 944 948])
%
% %******************* main settings for plot positions ****************
%
% scrsz=get(0,'ScreenSize');
% set(1,'OuterPosition',[0 scrsz(4)*0.01 scrsz(3)*3/4 scrsz(4)]) % first argument corresponds to figure number
% topGap=0.04;
% bottomGap=0.1;
% plotsGap=plotsGap;
% plotLeftGap = 0.12;
% plotLength = 0.7;
% plotHight=(1-topGap-bottomGap-(subPlotAmount-1)*plotsGap)/subPlotAmount;
%
% %---BEGIN--- plot 0 ---------------
% plot_south_north_hemisph_label ( bottomGap,  plotsGap, plotLeftGap, plotLength, FontSize, epoch, z, altit)
% %---END--- plot 0 ---------------
%
%
% plotNumber = 4;
% [~] = plot_wavelet_B(plotNumber, bottomGap, plotHight, plotsGap, plotLeftGap, plotLength, FontSize, xticks,...
%                                      mf_epoch,B);
%
% plotNumber = 3;
% [~] = plot_velocities_components_p(plotNumber, bottomGap, plotHight, plotsGap, plotLeftGap, plotLength, FontSize, xticks,...
%                                      epoch,v_x_p,  v_y_p,  v_z_p,  v_p);
%
% plotNumber = 2;
% [~] = plot_velocities_components_O(plotNumber, bottomGap, plotHight, plotsGap, plotLeftGap, plotLength, FontSize, xticks,...
%                                      epoch,v_x_O,  v_y_O,  v_z_O,  v_O);
% plotNumber = 1;
% [~] = plot_velocities_components_O2(plotNumber, bottomGap, plotHight, plotsGap, plotLeftGap, plotLength, FontSize, xticks,...
%                                      epoch,v_x_O2,  v_y_O2,  v_z_O2,  v_O2);

windowpos = get(f, 'position');

cButton    = uicontrol('Style','popupmenu', 'String',...
    {'Zoom on time axis (2)', 'Time stamp (1)',...
    'Magnetic Curl (0)' 'Magnetic Curl new coord (2)',...
    'Ion velocity MSO (0)', 'n-T diagram (2)', 'Box (4)','Energy Spectrum AVG 5 WINDOW (1)',...
    'Energy Spectra Set (2)', 'Minimum Variance MSE (4)', 'Minimum Variance (2)', 'Orbit and Vectors (2)', 'Orbit and Vectors MSE (4)',...
    '3D Distribution Function H (low E) (2)',...
    '3D Distribution Function O (low E) (2)',...
    '3D Distribution Function O2 (low E) (2)',...
    '3D Distribution Function H (high E) (2)',...
    '3D Distribution Function O (high E) (2)',...
    '3D Distribution Function O2 (high E) (2)',...
    '3D Distribution Function H (average) (2)',...
    '3D Distribution Function O (average) (2)',...
    '3D Distribution Function O2 (average) (2)',...
    'Velocity space XYZ H+ (1)', 'Velocity space XYZ H+ MSE (3)',...
    'Velocity space XYZ O+ (1)', 'Velocity space XYZ O+ MSE (3)',...
    'Velocity space XYZ O2+ (1)', 'Velocity space XYZ O2+ MSE (3)',...
    'View dayside orbit (0)', 'View dayside orbit MSE (3)', 'Plot sub-satellite point (0)', 'Orientation (0)',...
    'FOV bin direction H+ (2)', 'FOV bin direction O+ (2)', 'FOV bin direction O2+ (2)',...
    'FOV bin direction H+ MSE (4)', 'FOV bin direction O+ MSE (4)', 'FOV bin direction O2+ MSE (4)',...
    'O, O2 angle with surface (0)',...
    'Wavelet (0)', 'JxB acceleration (5)'},...
    'Units','normalized', 'Position', [0.8, 0.98, 0.08, 0.02],...
    'Callback',{@AdditionalPlots_Callback, epoch_d1_disp, energy_d1, swp_ind_d1_disp, eflux_d1_cleaned_sum, quat_mso_d1_disp, magf_d1_disp, eflux_d1_cleaned,...
    products_O, products_H, products_O2, products_mag,...
    altit_low, altit_high, x, y, z,...
    pos_sc_mso_d1_disp, quat_mso_d1_disp,...
    mass_arr_d1, phi_d1, theta_d1, filename_mag, windowpos});

end