function [xticks, sub_handle] = plot_electron_spectrogram(plotNumber, bottomGap, plotHight, plotsGap, plotLeftGap, plotLength, FontSize, xticks, caxis_lims,...
                                     epoch_svyspec,energy, diff_en_fluxes)

plotBottomPosition=bottomGap+(plotHight+plotsGap)*(plotNumber-1);
positionVector2=[plotLeftGap, plotBottomPosition, plotLength, plotHight];
sub_handle = subplot('Position',positionVector2);


pcolor(repmat(epoch_svyspec', 64, 1), repmat(energy,1,length(epoch_svyspec)), log10(squeeze(diff_en_fluxes')))
%caxis(caxis_lims)
shading flat
datetick('x')
xlim([epoch_svyspec(1) epoch_svyspec(end)])
set(gca, 'YScale', 'log', 'ycolor', 'black', 'ytick', [1 10 100 1000 10000]);
ylabel({'(5)'; 'SWEA E'; '[eV]'}, 'color', 'black');
grid on
set(gca, 'layer', 'top')
cbar_handle = colorbar;
p6 = get(sub_handle, 'position');
p6(3) = plotLength; %p6(3)*1.06;
set(sub_handle, 'pos', p6)

xlabel(cbar_handle, {'log10(eflux)', '[eV/(eV cm^2 s sr)]'}, 'FontSize', FontSize) %'FontWeight', 'bold')
set(gca,'xticklabel',[])
set(gca, 'FontSize', FontSize)
set(gca, 'gridalpha', 0.5, 'minorgridalpha', 0.5)

hold on

load('newColormap.mat')
colormap(gca, cmap);

v = ~isnan(diff_en_fluxes(1,:))';
ylim([min(energy(v)), max(energy(v))])

%colormap jet

end