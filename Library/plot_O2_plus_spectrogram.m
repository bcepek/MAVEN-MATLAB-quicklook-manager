function [xticks, sub_handle] = plot_O2_plus_spectrogram (plotNumber, bottomGap, plotHight, plotsGap, plotLeftGap, plotLength, FontSize, xticks, caxis_lims,...
                                     epoch_d1,nenergy_d1, energy_d1, swp_ind_d1, eflux_d1_cleaned_sum,...
                                     mf_epoch, Bx, By, Bz, n_O, n_O2)


plotBottomPosition = bottomGap + (plotHight+plotsGap) * (plotNumber-1);
positionVector2 = [plotLeftGap, plotBottomPosition, plotLength, plotHight];
sub_handle = subplot('Position',positionVector2);

% set(sub_handle6, 'pos', p6)
% pcolor(repmat(epoch_c6', nenergy_c6, 1), energy_c6(:,swp_ind_c6+1,1),  eflux_O2_c6);
pcolor(repmat(epoch_d1', nenergy_d1, 1), (energy_d1(:,swp_ind_d1 + 1,2,5)),  log10(squeeze(eflux_d1_cleaned_sum(:, 6,:))) );
%caxis(caxis_lims)
shading flat
% datetick('x')
xlim([epoch_d1(1) epoch_d1(end)])
set(gca, 'YScale', 'log', 'ycolor', 'black', 'ytick', [1 10 100 1000 10000]);
ylabel({'(4)'; 'STATIC E'; 'O_2^+ [eV]'});
grid on
set(gca, 'layer', 'top')
set(gca,'xticklabel',[])
set(gca, 'FontSize', FontSize)
set(gca,'xtick',xticks)
% colormap(gca, jet)
load('newColormap.mat')
colormap(gca, cmap);
% plot([epoch(high_boundary_index) epoch(high_boundary_index)], [limits(1) limits(2)], 'color', 'white' ,'LineStyle','--')
% plot([epoch(low_boundary_index) epoch(low_boundary_index)], [limits(1) limits(2)], 'color', 'white', 'LineStyle','--' )
% caxis([ -1 6])

set(gca, 'gridalpha', 0.5, 'minorgridalpha', 0.5)

yyaxis right
p1 = plot(epoch_d1, n_O2./n_O, 'color', 'black');
p1.Color(4) = 0.5;
ylabel('n_{O2}/n_O', 'color', 'black')
set(gca, 'ycolor', 'black')
ylim([0 10])

end