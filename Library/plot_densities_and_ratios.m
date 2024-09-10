function [xticks, sub_handle] = plot_densities_and_ratios (plotNumber, bottomGap, plotHight, plotsGap, plotLeftGap, plotLength, FontSize, xticks,...
                                     epoch,...
                                     n_p,  n_O,  n_O2, LineWidth)

plotBottomPosition=bottomGap+(plotHight+plotsGap)*(plotNumber-1);
positionVector2=[plotLeftGap, plotBottomPosition, plotLength, plotHight];
sub_handle = subplot('Position',positionVector2);
l_p = semilogy(epoch, n_p,'k','LineWidth',LineWidth);
hold on
l_o = semilogy(epoch, n_O,'LineWidth',LineWidth); %,'LineStyle','--');
l_o2 = semilogy(epoch, n_O2,'r','LineWidth',LineWidth); %,'LineStyle',':');
%semilogy(epoch, (n_O + n_O2)./n_p,'Color', [0.9 0.5 0], 'LineWidth', LineWidth )
% semilogy(epoch, n_CO2,'g','LineWidth',2);
datetick('x','HH:MM:SS');
ylabel ({'(6)', 'STATIC n', '[cm^{-3}]'})
ylim([10^-2 10^3])
set(gca,'YTick',[ 10^-2 10^-1 10^0 10^1 10^2 10^3 ]);
%l = legend ([l_p, l_o, l_o2],'p','O^+','O_2^+'); %, 'n_h/n_p');
% l.Location = 'west';
set(gca,'xticklabel',[])
set(gca, 'FontSize', FontSize)
grid on

yyaxis right
%plot(epoch, n_p./(n_p + n_O + n_O2),'Color',[0 0.498 0], 'LineWidth', LineWidth)
plot(epoch, n_p./(n_p + n_O + n_O2),'Color',[0 0 0.8], 'LineWidth', LineWidth)
ax = gca;
ax.YScale = 'lin';
ylim([0 1])
% ylim([10^-2 10^3])
% set(gca,'YTick',[ 10^-2 10^-1 10^0 10^1 10^2 10^3 ]);
%set(gca,'YColor',[0 0.498 0])
set(gca,'YColor',[0 0 0.8])
ylabel('n_p/(n_p + n_h)')
xlim([epoch(1) epoch(end)])
% plot_title = strcat(datestr(epoch(1),'yyyy-mm-dd'), {' '}, datestr(epoch(1),'HH:MM:SS'),{' - '}, datestr(epoch(end),'HH:MM:SS'),{' Altitude range: '},num2str(altit_low),{' - '},num2str(altit_high), ' km');
set(gca,'xtick',xticks)
%limits = ylim;
% plot([epoch(high_boundary_index) epoch(high_boundary_index)], [limits(1) limits(2)], 'color', 'black' ,'LineStyle','--')
% plot([epoch(low_boundary_index) epoch(low_boundary_index)], [limits(1) limits(2)], 'color', 'black', 'LineStyle','--' )

set(gca, 'gridalpha', 0.5, 'minorgridalpha', 0.5)
l = legend ([l_p, l_o, l_o2],'p','O^+','O_2^+', 'AutoUpdate','off');
l.Location = 'west';
end