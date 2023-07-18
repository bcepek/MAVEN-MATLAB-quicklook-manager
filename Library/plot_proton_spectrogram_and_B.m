function xticks = plot_proton_spectrogram_and_B (plotNumber, bottomGap, plotHight, plotsGap, plotLeftGap, plotLength, FontSize, xticks, caxis_lims,...
                                     epoch_d1,nenergy_d1, energy_d1, swp_ind_d1, eflux_d1_sum, altit, SZA,...
                                     mf_epoch, B)



plotBottomPosition=bottomGap+(plotHight+plotsGap)*(plotNumber-1);
positionVector2=[plotLeftGap, plotBottomPosition, plotLength, plotHight];
sub_handle7 = subplot('Position',positionVector2);
% pcolor(repmat(epoch_c6', nenergy_c6, 1), energy_c6(:,swp_ind_c6+1,1),  eflux_p_c6);
pcolor(repmat(epoch_d1', nenergy_d1, 1), (energy_d1(:, swp_ind_d1 + 1, 2, 5)),  log10(squeeze(eflux_d1_sum(:, 1,:))) );
caxis(caxis_lims)
shading flat
datetick('x')
xlim([epoch_d1(1) epoch_d1(end)])
set(gca, 'YScale', 'log', 'ycolor', 'black', 'ytick', [1 10 100 1000 10000]);
ylabel({'p energy'; '[eV]'}, 'color', 'black');
grid on
set(gca, 'layer', 'top')
plot_title = strcat(datestr(epoch_d1(1),'yyyy-mm-dd'), {' '}, datestr(epoch_d1(1),'HH:MM:SS'),{' - '},datestr(epoch_d1(end),'HH:MM:SS'),...
    {' '}, 'Altit = ', num2str(round(altit(1))), {'-'},num2str(round(altit(end))),{'; '},...
    'SZA = ', num2str(round(SZA(find(min(abs(600 - altit)) == abs(600 - altit))))),...
    '(Altit= ', num2str(round(altit(find(min(abs(600 - altit)) == abs(600 - altit))))), ')');
title(plot_title)
hold on
limits = ylim;
% plot([epoch(high_boundary_index) epoch(high_boundary_index)], [limits(1) limits(2)], 'color', 'white' ,'LineStyle','--')
% plot([epoch(low_boundary_index) epoch(low_boundary_index)], [limits(1) limits(2)], 'color', 'white', 'LineStyle','--' )

yyaxis right
%plot(epoch_cdf, sqrt(sum(magf.^2, 2)), 'color', 'red')
plot(mf_epoch,B, 'color', 'red')
ylim([0 inf])
ylabel('B [nT]')
set(gca,'xticklabel',[])
p7 = get(sub_handle7, 'position');
set(gca, 'FontSize', FontSize)
set(gca,'YColor','red')
if isempty(xticks )
    xticks = get(gca,'xtick');
end
set(gca,'xtick',xticks)
colormap(gca, jet)

end