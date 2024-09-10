function [xticks, sub_handle] = plot_proton_spectrogram_swia(plotNumber, bottomGap, plotHight, plotsGap, plotLeftGap, plotLength, FontSize, xticks, caxis_lims,...
                                     epoch_swia,energy_spectra_swia, spectra_diff_en_fluxes_swia,density_swia,velocity_mso_swia,choose_ind_swia)
[amu, anu, c, q, m_p, m_O, m_O2, m_CO2, Rm] = constants;
                                 
plotBottomPosition=bottomGap+(plotHight+plotsGap)*(plotNumber-1);
positionVector2=[plotLeftGap, plotBottomPosition, plotLength, plotHight];
sub_handle = subplot('Position',positionVector2);


pcolor(repmat(epoch_swia', 48, 1), repmat(energy_spectra_swia,1,length(epoch_swia)), log10(squeeze(spectra_diff_en_fluxes_swia')))
shading flat
%datetick('x',13)
datetick('x')
xlim([min(epoch_swia(choose_ind_swia)) max(epoch_swia(choose_ind_swia))])
set(gca, 'YScale', 'log', 'ycolor', 'black', 'ytick', [1 10 100 1000 10000]);
ylabel({'(1)'; 'SWIA E'; '[eV]'}, 'color', 'black');
grid on
set(gca, 'layer', 'top')
cbar_handle = colorbar;

ylabels = get(cbar_handle, 'ticklabels');
for i=1:length(ylabels)
    ylabels{i} = ['10^{' num2str(ylabels{i}) '}'];
end
set(cbar_handle, 'ticklabels', ylabels, 'fontsize', FontSize)

p6 = get(sub_handle, 'position');
p6(3) = plotLength; %p6(3)*1.06;
set(sub_handle, 'pos', p6)

xlabel(cbar_handle, {'eflux', '[eV/(eV cm^2 s sr)]'}, 'FontSize', FontSize) 
set(gca,'xtick',xticks)
set(gca,'xticklabel',[])
set(gca, 'FontSize', FontSize)
hold on
%limits = ylim;
load('newColormap.mat');
colormap(gca, cmap);
set(gca, 'gridalpha', 0.5, 'minorgridalpha', 0.5)

yyaxis right
pressure = 1e9*(m_p*1e-3*density_swia(choose_ind_swia,:)*1e6.*sum((velocity_mso_swia(choose_ind_swia,:)*1e3).^2,2))/2;  % in nPa
% pressure = pressure/q;  % in eV/m3
% pressure = pressure*1e6;    % in eV/cm3
plot(epoch_swia(choose_ind_swia), pressure*1e-8,...
     'Color',[0.82 0.14 0.32],'LineWidth',1.5)
ylabel({'Pressure,'; 'dyn cm^{-2}'}, 'Color',[0.82 0.14 0.32])
set(gca, 'yscale', 'log', 'ytick', [1e-16 1e-15 1e-14 1e-13 1e-12 1e-11 1e-10 1e-9],...
    'ycolor', [0.82 0.14 0.32])
ylim([1e-13 1e-8])
grid on
%grid minor 
box on 
xlim([min(epoch_swia(choose_ind_swia)) max(epoch_swia(choose_ind_swia))])


end