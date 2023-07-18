function xticks = plot_densities (plotNumber, bottomGap, plotHight, plotsGap, plotLeftGap, plotLength, FontSize, xticks,...
                                     epoch,...
                                     n_p,  n_O,  n_O2)



plotBottomPosition=bottomGap+(plotHight+plotsGap)*(plotNumber-1);
positionVector=[plotLeftGap, plotBottomPosition, plotLength, plotHight];
subplot('Position',positionVector);
%************************************
semilogy(epoch, n_p,'k','LineWidth',1.5);
grid on
hold on
semilogy(epoch, n_O,'LineWidth',1.5) %,'LineStyle','--');
semilogy(epoch, n_O2,'r','LineWidth',1.5) %,'LineStyle',':');
% semilogy(epoch, n_CO2,'g','LineWidth',2);
datetick('x','HH:MM:SS');
ylabel ({'number', 'density [cm^{-3}]'})
ylim([10^-2 10^3])
set(gca,'YTick',[ 10^-2 10^-1 10^0 10^1 10^2 10^3 ]);
l = legend ('p','O^+','O_2^+');
l.Location = 'west';
set(gca,'xticklabel',[])
set(gca, 'FontSize', FontSize)
xlim([epoch(1) epoch(end)])
set(gca,'xtick',xticks)


end