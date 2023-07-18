function xticks = plot_O_plus_spectrogram_and_B (plotNumber, bottomGap, plotHight, plotsGap, plotLeftGap, plotLength, FontSize, xticks, caxis_lims,...
                                     epoch_d1,nenergy_d1, energy_d1, swp_ind_d1, eflux_d1_cleaned_sum,...
                                        mf_epoch, B)
                                 

plotBottomPosition=bottomGap+(plotHight+plotsGap)*(plotNumber-1);
positionVector2=[plotLeftGap, plotBottomPosition, plotLength, plotHight];
sub_handle6 = subplot('Position',positionVector2);

% pcolor(repmat(epoch_cdf', nenergy, 1), energy(:,swp_ind+1,1),  flux_O);
pcolor(repmat(epoch_d1', nenergy_d1, 1), (energy_d1(:,swp_ind_d1 + 1, 2, 5)),  log10(squeeze(eflux_d1_cleaned_sum(:, 5,:))) );
caxis(caxis_lims)
shading flat
datetick('x')
xlim([epoch_d1(1) epoch_d1(end)])
set(gca, 'YScale', 'log', 'ycolor', 'black', 'ytick', [1 10 100 1000 10000]);
ylabel({'O^+ energy'; '[eV]'}, 'color', 'black');
grid on
set(gca, 'layer', 'top')
cbar_handle = colorbar;
p6 = get(sub_handle6, 'position');
p6(3) = plotLength; %p6(3)*1.06;
set(sub_handle6, 'pos', p6)
% cbar_tick = get(cbar_handle, 'ytick');
% barlabels = cell(size(cbar_tick, 1), 1);
% for j=1:size(cbar_tick, 2)
%     barlabels{j} = ['10^{', num2str(cbar_tick(j)), '}'];
% end
% set(cbar_handle, 'xticklabel', barlabels)
xlabel(cbar_handle, {'log10(eflux)', '[1/(eV cm^2 s sr)]'}, 'FontSize', FontSize) %'FontWeight', 'bold')
set(gca,'xticklabel',[])
set(gca, 'FontSize', FontSize)
set(gca,'xtick',xticks)
hold on
limits = ylim;
colormap(gca, parula)
% plot([epoch(high_boundary_index) epoch(high_boundary_index)], [limits(1) limits(2)], 'color', 'white' ,'LineStyle','--')
% plot([epoch(low_boundary_index) epoch(low_boundary_index)], [limits(1) limits(2)], 'color', 'white', 'LineStyle','--' )
%---END--- plot 6 ---------------
% caxis([ -1 6])
% cbar_handle.Ticks = [-1, 0 ,1, 2, 3, 4, 5, 6];

yyaxis right
%plot(epoch_cdf, sqrt(sum(magf.^2, 2)), 'color', 'red')
plot(mf_epoch,B, 'color', 'red')
ylim([0 inf])
ylabel('B [nT]')
set(gca,'xticklabel',[])
p7 = get(sub_handle6, 'position');
set(gca, 'FontSize', FontSize)
set(gca,'YColor','red')
if isempty(xticks )
    xticks = get(gca,'xtick');
end
set(gca,'xtick',xticks)
colormap(gca, parula)


end