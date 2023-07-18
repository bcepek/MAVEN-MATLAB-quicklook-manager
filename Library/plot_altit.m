function xticks = plot_altit (plotNumber, bottomGap, plotHight, plotsGap, plotLeftGap, plotLength, FontSize, xticks,...
                                     epoch_d1, altit, valid_d1, quality_flag_d1, swp_ind_d1 )


plotBottomPosition=bottomGap+(plotHight+plotsGap)*(plotNumber-1);
positionVector2=[plotLeftGap, plotBottomPosition, plotLength, plotHight];
subplot('Position',positionVector2);
%*********************************************


% plot(epoch_d1, valid_d1, 'b','LineWidth',1.5)
% grid on
% hold on
% plot(epoch_d1, quality_flag_d1)

plot(epoch_d1, swp_ind_d1 )
ylabel('swp ind')
% grid on
% hold on
% plot(epoch_d1, altit)
% legend('att', 'quality')


datetick('x','HH:MM:SS');
% ylabel ({'attenuator', 'state'})
% xlabel ('time')

set(gca,'xticklabel',[])
if plotNumber == 1
    datetick
    xlabel('UT')
end
set(gca, 'FontSize',FontSize)
xlim([epoch_d1(1) epoch_d1(end)])
% set(gca,'YTick',[1 10 10^2 10^3 10^4]);
if isempty(xticks )
    xticks = get(gca,'xtick');
end
set(gca,'xtick',xticks)
% ylim([0 3])
% plot([epoch(high_boundary_index) epoch(high_boundary_index)], [limits(1) limits(2)], 'color', 'black' ,'LineStyle','--')
% plot([epoch(low_boundary_index) epoch(low_boundary_index)], [limits(1) limits(2)], 'color', 'black', 'LineStyle','--' )




yyaxis right
plot(epoch_d1, altit)
ylabel('altit, km')
% plot(epoch_d1, sc_pot_d1, 'r','LineWidth',1.5)
% l = legend('0 - closed; 1 - elect; 2 - mech; 3 - (elect + mech)', 'valid',  'sc pot');
% l.Location = 'west';


end