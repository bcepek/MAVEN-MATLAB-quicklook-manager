function xticks = plot_Bx (plotNumber, bottomGap, plotHight, plotsGap, plotLeftGap, plotLength, FontSize, xticks, caxis_lims,...
                                     epoch_d1,nenergy_d1, energy_d1, swp_ind_d1, eflux_d1_cleaned_sum,...
                                     mf_epoch, Bx, By, Bz, LineWidth, Blx, Bly, Blz,x,y,z)

plotBottomPosition = bottomGap + (plotHight+plotsGap) * (plotNumber-1);
positionVector2 = [plotLeftGap, plotBottomPosition, plotLength, plotHight];
sub_handle = subplot('Position',positionVector2);

positive = Bx>0;
Bx_pos = Bx;
Bx_pos(~positive) = nan;
Bx_neg = Bx;
Bx_neg(positive) = nan;

plot(mf_epoch,Bx_pos, 'color', [1 0 0], 'LineWidth', LineWidth);
hold on
plot(mf_epoch,Bx_neg, 'color', [0 0 1], 'LineWidth', LineWidth);
hold off

grid on
ylabel({'(8)'; 'B_X [nT]'})
set(gca,'xticklabel',[])

xlim([epoch_d1(1) epoch_d1(end)])
set(gca,'xtick',xticks)

set(gca, 'FontSize', FontSize)
set(gca,'YColor','black')
set(gca,'xtick',xticks)

% B_lim_arr = 10:10:2000;
% B_max =  round(max([max(abs(Bx)), max(abs(By)), max(abs(Bz))]));
% B_lim_arr = B_lim_arr(B_lim_arr - B_max >=0);
% B_lim_ind = find( (abs(B_lim_arr - B_max) == min(abs(B_lim_arr - B_max)))  );
% B_lim = B_lim_arr(B_lim_ind);
% ylim([-B_lim B_lim])
grid on

colormap(gca, jet)

set(gca, 'gridalpha', 0.5, 'minorgridalpha', 0.5)

yyaxis right
sza = acos(x./sqrt(x.^2+y.^2+z.^2))*180/pi;
plot(epoch_d1,sza, 'linestyle', '--', 'color', 'black')
ylabel('SZA, deg', 'color', 'black')
set(gca, 'ycolor', 'black', 'ytick', [0 45 90 135 180])
ylim([0 180])

end