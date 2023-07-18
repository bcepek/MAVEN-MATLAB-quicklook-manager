function xticks = plot_energies (plotNumber, bottomGap, plotHight, plotsGap, plotLeftGap, plotLength, FontSize, xticks,...
                                     epoch,...
                                     n_p, v_p, T_p,  n_O, v_O, T_O, n_O2, v_O2, T_O2,...
                                     mf_epoch, B, LineWidth)

[~, ~, ~, ~, m_p, m_O, m_O2, ~, ~] = constants;                                 

plotBottomPosition=bottomGap+(plotHight+plotsGap)*(plotNumber-1);
positionVector2=[plotLeftGap, plotBottomPosition, plotLength, plotHight];
subplot('Position',positionVector2);
%*****************************************
semilogy(epoch, (n_p.*v_p.^2*m_p + n_O.*v_O.^2*m_O + n_O2.*v_O2.^2*m_O2)/2*10^10/(1.6*10^-12),'k', 'LineWidth', LineWidth) %,'LineStyle','--')
hold on
semilogy(mf_epoch, B.^2 /(8 * pi)*10^-10/(1.6*10^-12), 'LineWidth', LineWidth)
semilogy(epoch, ( n_p.*T_p + n_O.*T_O + n_O2.*T_O2 ),'r', 'LineWidth', LineWidth) %,'LineStyle',':'   )
l = legend ('nmV^2/2','B^2/(8pi)','nkT');
l.Location = 'west';
grid on
% ylabel (strcat(char(425),' eV cm^{-3}')) %n_im_iV_i^2/2,
ylabel('[eV cm^{-3}]')
datetick('x','HH:MM:SS');
% set(gca,'XTickLabel',[])
ylim([10^1 10^4])
set(gca,'YTick',[ 10^2 10^3 10^4]);
set(gca,'XTickLabel',[])
set(gca, 'FontSize',FontSize)
xlim([epoch(1) epoch(end)])
set(gca,'xtick',xticks)
limits = ylim;
% plot([epoch(high_boundary_index) epoch(high_boundary_index)], [limits(1) limits(2)], 'color', 'black' ,'LineStyle','--')
% plot([epoch(low_boundary_index) epoch(low_boundary_index)], [limits(1) limits(2)], 'color', 'black', 'LineStyle','--' )

set(gca, 'gridalpha', 0.5, 'minorgridalpha', 0.5)

% yyaxis right
% n0 = 3000; 
% v0 = 5;
% T0 = 0.5;
% plot(epoch, (1 - exp(-(v0*T0/n0)*n_O./(v_O.*T_O))),'Color',[0 0.498 0], 'LineWidth', 1.5 )
% ax = gca;
% %ax.YScale = 'log';
% %ylim([0 1])
% % ylim([10^-2 10^3])
% %set(gca,'YTick',[ 10^-1 10^0 10^1]);
% set(gca,'YColor',[0 0.498 0])
% ylabel('1-exp(-n_O*v_0*T_0/(n_0*v_O*T_O))')
xlim([epoch(1) epoch(end)])
% plot_title = strcat(datestr(epoch(1),'yyyy-mm-dd'), {' '}, datestr(epoch(1),'HH:MM:SS'),{' - '}, datestr(epoch(end),'HH:MM:SS'),{' Altitude range: '},num2str(altit_low),{' - '},num2str(altit_high), ' km');
set(gca,'xtick',xticks)

end