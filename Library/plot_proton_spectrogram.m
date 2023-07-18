function xticks = plot_proton_spectrogram (plotNumber, bottomGap, plotHight, plotsGap, plotLeftGap, plotLength, FontSize, xticks, caxis_lims,...
                                     epoch_d1,nenergy_d1, energy_d1, swp_ind_d1, eflux_d1_sum, altit, SZA,...
                                     mf_epoch, B)



plotBottomPosition=bottomGap+(plotHight+plotsGap)*(plotNumber-1);
positionVector2=[plotLeftGap, plotBottomPosition, plotLength, plotHight];
sub_handle7 = subplot('Position',positionVector2);
% pcolor(repmat(epoch_c6', nenergy_c6, 1), energy_c6(:,swp_ind_c6+1,1),  eflux_p_c6);
pcolor(repmat(epoch_d1', nenergy_d1, 1), (energy_d1(:,swp_ind_d1 + 1,2,5)),  log10(squeeze(eflux_d1_sum(:, 1, :))) );
%caxis(caxis_lims)
shading flat
datetick('x')
xlim([epoch_d1(1) epoch_d1(end)])
set(gca, 'YScale', 'log', 'ycolor', 'black', 'ytick', [1 10 100 1000 10000]);
ylabel({'p energy'; '[eV]'}, 'color', 'black');
grid on
set(gca, 'layer', 'top')

p7 = get(sub_handle7, 'position');
p7(3) = plotLength; %p6(3)*1.06;
set(sub_handle7, 'pos', p7)
plot_title = strcat('MAVEN STATIC', {' '}, datestr(epoch_d1(1),'yyyy-mm-dd'), {' '}, datestr(epoch_d1(1),'HH:MM:SS'),{' - '},datestr(epoch_d1(end),'HH:MM:SS'));
%{
plot_title = strcat(datestr(epoch_d1(1),'yyyy-mm-dd'), {' '}, datestr(epoch_d1(1),'HH:MM:SS'),{' - '},datestr(epoch_d1(end),'HH:MM:SS'),...
    {' '}, 'Altit = ', num2str(round(altit(1))), {'-'},num2str(round(altit(end))),{'; '},...
    'SZA = ', num2str(round(SZA(find(min(abs(600 - altit)) == abs(600 - altit))))),...
    '(Altit= ', num2str(round(altit(find(min(abs(600 - altit)) == abs(600 - altit))))), ')');
%}
title(plot_title)
hold on
limits = ylim;
set(gca,'xticklabel',[])
set(gca, 'FontSize', FontSize)
% plot([epoch(high_boundary_index) epoch(high_boundary_index)], [limits(1) limits(2)], 'color', 'white' ,'LineStyle','--')
% plot([epoch(low_boundary_index) epoch(low_boundary_index)], [limits(1) limits(2)], 'color', 'white', 'LineStyle','--' )
% colormap(gca, jet)
load('newColormap.mat')
colormap(gca, cmap);
% caxis([ -1 6])
set(gca, 'gridalpha', 0.5, 'minorgridalpha', 0.5)

end