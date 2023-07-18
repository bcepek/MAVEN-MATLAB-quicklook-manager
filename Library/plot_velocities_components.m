function xticks = plot_velocities_components (plotNumber, bottomGap, plotHight, plotsGap, plotLeftGap, plotLength, FontSize, xticks,...
                                     epoch,...
                                     v_x_p,  v_y_p,  v_z_p,  v_p,...
                                     v_x_O,  v_y_O,  v_z_O,  v_O,...
                                     v_x_O2, v_y_O2, v_z_O2, v_O2,...
                                     valid_p);

                                                                                                
plotBottomPosition=bottomGap+(plotHight+plotsGap)*(plotNumber-1);
positionVector=[plotLeftGap, plotBottomPosition, plotLength, plotHight];
subplot('Position',positionVector);
%*********************************************
plot(epoch, v_x_p,'k','LineWidth',1);
grid on
hold on
plot(epoch, v_x_O,'LineWidth',1);
%plot(epoch, v_z_p,'r','LineWidth',1);
% semilogy(epoch, v_O,'LineWidth',1.5) %,'LineStyle','--');
% semilogy(epoch, v_O2,'r','LineWidth',1.5) %,'LineStyle',':');

%datetick('x','HH:MM:SS');
set(gca,'xticklabel',[])
ylabel ({'velocity', 'components', '[km/s]'})
% xlabel ('time')
l = legend('p_{Vx}','O_{Vx}');
l.Location = 'west';
set(gca,'xticklabel',[])

set(gca, 'FontSize',FontSize)
xlim([epoch(1) epoch(end)])
% ylim([4 500])

if isempty(xticks )
    xticks = get(gca,'xtick');
end
set(gca,'xtick',xticks)
set(gca,'xticklabel',[])
ylim([-100 20])
% plot([epoch(high_boundary_index) epoch(high_boundary_index)], [limits(1) limits(2)], 'color', 'black' ,'LineStyle','--')
% plot([epoch(low_boundary_index) epoch(low_boundary_index)], [limits(1) limits(2)], 'color', 'black', 'LineStyle','--' )

%{
yyaxis right
plot(epoch, valid_p, 'LineWidth', 1.5 )
ax = gca;
%ax.YScale = 'log';
ylim([0.1 10])
% ylim([10^-2 10^3])
%set(gca,'YTick',[ 10^-1 10^0 10^1]);
%set(gca,'YColor',[0 0.498 0])
ylabel('valid')
xlim([epoch(1) epoch(end)])
% plot_title = strcat(datestr(epoch(1),'yyyy-mm-dd'), {' '}, datestr(epoch(1),'HH:MM:SS'),{' - '}, datestr(epoch(end),'HH:MM:SS'),{' Altitude range: '},num2str(altit_low),{' - '},num2str(altit_high), ' km');
set(gca,'xtick',xticks)
limits = ylim;
%}


end