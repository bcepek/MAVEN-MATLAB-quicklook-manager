function xticks = plot_temp_from_energy (plotNumber, bottomGap, plotHight, plotsGap, plotLeftGap, plotLength, FontSize, xticks,...
                                     epoch,...
                                     T_p,  T_O,  T_O2, T_p_energy,  T_O_energy,  T_O2_energy, v_O,...
                                     mf_epoch, B, LineWidth)

[~, ~, ~, ~, m_p, m_O, m_O2, ~, ~] = constants;                                 

plotBottomPosition=bottomGap+(plotHight+plotsGap)*(plotNumber-1);
positionVector2=[plotLeftGap, plotBottomPosition, plotLength, plotHight];
subplot('Position',positionVector2);
%*****************************************

semilogy(datetime(datevec(epoch)), T_p_energy,...
         datetime(datevec(epoch)), T_O_energy,...
         datetime(datevec(epoch)), T_O2_energy)%,...
 %        datetime(datevec(epoch)), movstd(T_O_energy, 15)) %,'LineStyle',':'   )
l = legend ('T_p energy','T_O energy','T_{O2} energy');
l.Location = 'west';
grid on
% ylabel (strcat(char(425),' eV cm^{-3}')) %n_im_iV_i^2/2,
ylabel('[eV cm^{-3}]')
datetick('x','HH:MM:SS');
% set(gca,'XTickLabel',[])
ylim([0.1 10^4])
set(gca,'YTick',[10^0 10^1 10^2 10^3 10^4 ]);
%set(gca,'XTickLabel',[])
set(gca, 'FontSize',FontSize)
xlim([datetime(datevec(epoch(1))) datetime(datevec(epoch(end)))])
%set(gca,'xtick',xticks)
limits = ylim;
% plot([epoch(high_boundary_index) epoch(high_boundary_index)], [limits(1) limits(2)], 'color', 'black' ,'LineStyle','--')
% plot([epoch(low_boundary_index) epoch(low_boundary_index)], [limits(1) limits(2)], 'color', 'black', 'LineStyle','--' )

% yyaxis right
% semilogy(datetime(datevec(epoch)), 1e-3*1e6/1.6e-19* m_O * v_O.^2  / 2 )

end