function xticks = plot_velocities_and_ionosph_criterion (plotNumber, bottomGap, plotHight, plotsGap, plotLeftGap, plotLength, FontSize, xticks,...
                                     epoch,...
                                     v_p, n_O, v_O, T_O,  v_O2, LineWidth)


plotBottomPosition=bottomGap+(plotHight+plotsGap)*(plotNumber-1);
positionVector2=[plotLeftGap, plotBottomPosition, plotLength, plotHight];
subplot('Position',positionVector2);
%*********************************************
h1 = semilogy(epoch, v_p,'k','LineWidth', LineWidth);
grid on
hold on
h2 = semilogy(epoch, v_O,'LineWidth', LineWidth); %,'LineStyle','--');
h3 = semilogy(epoch, v_O2,'r','LineWidth', LineWidth); %,'LineStyle',':');
datetick('x','HH:MM:SS');
ylabel ({'velocity', '[km/s]'})
% xlabel ('time')
l = legend([h1 h2 h3], {'p','O^+','O_2^+'});
l.Location = 'west';
set(gca,'xticklabel',[])
set(gca, 'FontSize',FontSize)
xlim([epoch(1) epoch(end)])
ylim([4 500])
set(gca,'YTick',[1 10 10^2 10^3 10^4]);
if isempty(xticks )
    xticks = get(gca,'xtick');
end
set(gca,'xtick',xticks)
limits = ylim;
% plot([epoch(high_boundary_index) epoch(high_boundary_index)], [limits(1) limits(2)], 'color', 'black' ,'LineStyle','--')
% plot([epoch(low_boundary_index) epoch(low_boundary_index)], [limits(1) limits(2)], 'color', 'black', 'LineStyle','--' )


yyaxis right
plot(epoch, (1 - exp(-n_O./(v_O.*T_O))),'Color',[0 0.498 0], 'LineWidth', 1.5 )
ax = gca;
%ax.YScale = 'log';
ylim([0 1])
% ylim([10^-2 10^3])
%set(gca,'YTick',[ 10^-1 10^0 10^1]);
set(gca,'YColor',[0 0.498 0])
ylabel('1-exp(-n_O/(v_O*T_O))')
xlim([epoch(1) epoch(end)])
% plot_title = strcat(datestr(epoch(1),'yyyy-mm-dd'), {' '}, datestr(epoch(1),'HH:MM:SS'),{' - '}, datestr(epoch(end),'HH:MM:SS'),{' Altitude range: '},num2str(altit_low),{' - '},num2str(altit_high), ' km');
set(gca,'xtick',xticks)
limits = ylim;

end