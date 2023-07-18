function xticks = plot_velocities (plotNumber, bottomGap, plotHight, plotsGap, plotLeftGap, plotLength, FontSize, xticks,...
                                     epoch,...
                                     v_p,  v_O,  v_O2)


                                                                                                
plotBottomPosition=bottomGap+(plotHight+plotsGap)*(plotNumber-1);
positionVector=[plotLeftGap, plotBottomPosition, plotLength, plotHight];
subplot('Position',positionVector);
%*********************************************
semilogy(epoch, v_p,'k','LineWidth',1.5);
grid on
hold on
semilogy(epoch, v_O,'LineWidth',1.5) %,'LineStyle','--');
semilogy(epoch, v_O2,'r','LineWidth',1.5) %,'LineStyle',':');

datetick('x','HH:MM:SS');
ylabel ({'(7)'; 'STATIC v'; '[km/s]'})
% xlabel ('time')
l = legend('p','O^+','O_2^+');
l.Location = 'west';
set(gca,'xticklabel',[])
if plotNumber == 1
    datetick
    xlabel('UT')
end
set(gca, 'FontSize',FontSize)
xlim([epoch(1) epoch(end)])
ylim([4 500])
set(gca,'YTick',[1 10 10^2 10^3 10^4]);
if isempty(xticks )
    xticks = get(gca,'xtick');
end
set(gca,'xtick',xticks)
% plot([epoch(high_boundary_index) epoch(high_boundary_index)], [limits(1) limits(2)], 'color', 'black' ,'LineStyle','--')
% plot([epoch(low_boundary_index) epoch(low_boundary_index)], [limits(1) limits(2)], 'color', 'black', 'LineStyle','--' )

set(gca, 'gridalpha', 0.5, 'minorgridalpha', 0.5)


end