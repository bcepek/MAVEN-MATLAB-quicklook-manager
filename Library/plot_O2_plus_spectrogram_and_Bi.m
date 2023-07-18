function xticks = plot_O2_plus_spectrogram_and_Bi (plotNumber, bottomGap, plotHight, plotsGap, plotLeftGap, plotLength, FontSize, xticks, caxis_lims,...
                                     epoch_d1,nenergy_d1, energy_d1, swp_ind_d1, eflux_d1_cleaned_sum,...
                                     mf_epoch, Bx, By, Bz)


plotBottomPosition = bottomGap + (plotHight+plotsGap) * (plotNumber-1);
positionVector2 = [plotLeftGap, plotBottomPosition, plotLength, plotHight];
sub_handle = subplot('Position',positionVector2);

% set(sub_handle6, 'pos', p6)
% pcolor(repmat(epoch_c6', nenergy_c6, 1), energy_c6(:,swp_ind_c6+1,1),  eflux_O2_c6);
pcolor(repmat(epoch_d1', nenergy_d1, 1), (energy_d1(:,swp_ind_d1 + 1, 2, 5)),  log10(squeeze(eflux_d1_cleaned_sum(:, 6,:))) );
caxis(caxis_lims)
shading flat
% datetick('x')
xlim([epoch_d1(1) epoch_d1(end)])
set(gca, 'YScale', 'log', 'ycolor', 'black', 'ytick', [1 10 100 1000 10000]);
ylabel({'O_2^+ energy'; '[eV]'}, 'color', 'yellow');
grid on
set(gca, 'layer', 'top')
set(gca,'xticklabel',[])
set(gca, 'FontSize', FontSize)
set(gca,'xtick',xticks)

yyaxis right
a = plot(mf_epoch,Bx,  'r');
legend('Bx')
hold on 
b = plot(mf_epoch,By,  'g');
c = plot(mf_epoch,Bz,  'color', 'white');
ylabel('B [nT]')
set(gca,'xticklabel',[])
lgd = legend([a,b,c] , {'B_x', 'B_y', 'B_z'});

lgd.Color = 'black';
lgd.TextColor = 'white';
lgd.Location = 'west';
ax = gca;
yyaxis(ax, 'right') 
set(gca, 'FontSize', FontSize)
set(gca,'YColor','black')
set(gca,'xtick',xticks)
% set(sub_handle6, 'pos', p6)
% p6 = get(sub_handle6, 'position');
% p6(3) = p6(3)*1.06;
B_lim_arr = 10:10:100;
B_max =  round(max([max(abs(Bx)), max(abs(By)), max(abs(Bz))]));
B_lim_arr = B_lim_arr(B_lim_arr - B_max >=0);
B_lim_ind = find( (abs(B_lim_arr - B_max) == min(abs(B_lim_arr - B_max)))  );
B_lim = B_lim_arr(B_lim_ind);
ylim([-B_lim B_lim])
grid on
set(gca,'GridColor','white')
limits = ylim;
colormap(gca, jet)
% plot([epoch(high_boundary_index) epoch(high_boundary_index)], [limits(1) limits(2)], 'color', 'white' ,'LineStyle','--')
% plot([epoch(low_boundary_index) epoch(low_boundary_index)], [limits(1) limits(2)], 'color', 'white', 'LineStyle','--' )


end