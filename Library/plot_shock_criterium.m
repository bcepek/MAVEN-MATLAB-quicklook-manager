function xticks = plot_shock_criterium (plotNumber, bottomGap, plotHight, plotsGap, plotLeftGap, plotLength, FontSize, xticks,...
                                     epoch_d1, n_p, n_O, n_O2)


                          
                                 
plotBottomPosition=bottomGap+(plotHight+plotsGap)*(plotNumber-1);
positionVector=[plotLeftGap, plotBottomPosition, plotLength, plotHight];
subplot('Position',positionVector);
%*********************************************
plot(epoch_d1, (n_p + n_O*16 + n_O2*32)./(n_p + n_O + n_O2),'k','LineWidth',1.5);
hold on
y1 = repmat(4/3,   size(epoch_d1));
plot(epoch_d1, y1,'r','LineWidth',1)
y2 = repmat(25/16, size(epoch_d1));
plot(epoch_d1, y2,'g','LineWidth',1)
y3 = repmat(49/24, size(epoch_d1));
plot(epoch_d1, y3 ,'b','LineWidth',1)
grid on
% hold on
% plot(epoch_d1, valid_d1, 'b','LineWidth',1.5)
% plot(epoch_d1, quality_flag_d1)



ylabel('M/M_{\infty}')

% xlabel ('time')
datetick('x','HH:MM:SS');
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

legend('M/M_{\infty}', 'y = 2','y = 5/3', 'y = 7/5' )

%{
yyaxis right
plot(epoch_d1, sc_pot_d1, 'r','LineWidth',1.5)
l = legend('0 - closed; 1 - elect; 2 - mech; 3 - (elect + mech)', 'valid',  'sc pot');
l.Location = 'west';
%}
ylim([0 4])

end