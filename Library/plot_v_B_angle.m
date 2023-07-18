function xticks = plot_v_B_angle (plotNumber, bottomGap, plotHight, plotsGap, plotLeftGap, plotLength, FontSize, xticks,...
                                     epoch,...
                                     v_x_p,  v_y_p,  v_z_p, v_p,...
                                     v_x_O,  v_y_O,  v_z_O, v_O,...
                                     v_x_O2, v_y_O2, v_z_O2, v_O2,...
                                     Blx, Bly, Blz, Bl, LineWidth)

plotBottomPosition=bottomGap+(plotHight+plotsGap)*(plotNumber-1);
positionVector2=[plotLeftGap, plotBottomPosition, plotLength, plotHight];
subplot('Position',positionVector2);
%*************************************************
plot(epoch, acos( (v_x_p.*Blx + v_y_p.*Bly + v_z_p.*Blz)./(v_p.*Bl) )*180/pi ,'k', 'LineWidth', LineWidth)
grid on
hold on
plot(epoch, acos( (v_x_O.*Blx + v_y_O.*Bly + v_z_O.*Blz)./(v_O.*Bl) )*180/pi , 'LineWidth', LineWidth) %,'LineStyle','--')
plot(epoch, acos( (v_x_O2.*Blx + v_y_O2.*Bly + v_z_O2.*Blz)./(v_O2.*Bl) )*180/pi ,'r', 'LineWidth', LineWidth) %,'LineStyle',':')
% plot(epoch, acos( (v_x_CO2.*Blx + v_y_CO2.*Bly + v_z_CO2.*Blz)./(v_CO2.*Bl) )*180/pi ,'g', 'LineWidth', 2)
l = legend('p','O^+','O_2^+');
l.Location = 'west';
ylabel({'V - B angle [^o]'})
% datetick %('x','HH:MM');
set(gca,'XTickLabel',[])
xlim([epoch(1) epoch(end)])
% xticks = get(gca, 'xtick');
set(gca, 'FontSize', FontSize)
% set(gca,'XTickLabel',[])
set(gca,'YTick',[0 45 90 135 180]);
set(gca,'xtick',xticks)
ylim([0 180])
limits = ylim;
% plot([epoch(high_boundary_index) epoch(high_boundary_index)], [limits(1) limits(2)], 'color', 'black' ,'LineStyle','--')
% plot([epoch(low_boundary_index) epoch(low_boundary_index)], [limits(1) limits(2)], 'color', 'black', 'LineStyle','--' )

end