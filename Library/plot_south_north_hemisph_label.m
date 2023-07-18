function xticks = plot_south_north_hemisph_label ( bottomGap,  plotsGap, plotLeftGap, plotLength, FontSize, epoch, z, altit)


plotHightSpecial = bottomGap / 5;
plotBottomPosition=bottomGap - plotHightSpecial - plotsGap;
positionVector2=[plotLeftGap, plotBottomPosition, plotLength, plotHightSpecial];
sub_1 = subplot('Position',positionVector2);

if(~isempty(epoch))
    Z_color = ones(size(epoch));
    neg_z = z<0;
    Z_color(neg_z) = - Z_color(neg_z);
    Z_color(end)=1;
    pcolor([epoch epoch]', [zeros(size(epoch)) ones(size(epoch))]', [Z_color Z_color]')
end

datetick
%xtickangle(30)
shading flat
c = colorbar;
colorbar('Ticks',[c.Limits(1) c.Limits(2)], 'TickLabels', {'Z_{MSO} < 0','Z_{MSO} > 0'});
set(gca, 'YTick',[])
set(gca, 'FontSize', FontSize)
xlabel('UT')
if(~isempty(epoch))
    xlim([epoch(1) epoch(end)])
end
% ylabel('sign(z)', 'Rotation', 0)

p6 = get(sub_1, 'position');
p6(3) = plotLength; %p6(3)*1.06;
set(sub_1, 'pos', p6)
xticks = get(gca, 'xtick');
set(gca, 'layer', 'top')



yyaxis right
interesting_altit = [100 200 300 400 500 600 700 800 900 1000 1500 2000 3000 4000 4500 5000 5500 6000];
[Xtest, Ytest] = meshgrid(altit, interesting_altit);
[M, I] = min(abs(Ytest - Xtest), [],2);
I = I(M < 50);
interesting_altit = interesting_altit(M < 50);
% I = find(min(abs(Ytest - Xtest), [],2) == abs(Ytest - Xtest), [], 2);
plot(epoch(I), ones(size(I))*1, 'w.', 'markersize', 12)
hold on
text(epoch(I), ones(size(I))*1, num2cell(interesting_altit),'Color','w', 'HorizontalAlignment', 'left', ...
    'fontsize', 12, 'fontweight', 'bold')

interesting_altit_minor = 100:10:1000;
interesting_altit_minor(rem(interesting_altit_minor, 100)==0)=[];
[Xtest, Ytest] = meshgrid(altit, interesting_altit_minor);
[M, I] = min(abs(Ytest - Xtest), [],2);
I = I(M < 3);
interesting_altit_minor = interesting_altit_minor(M < 50);
hold on
text(epoch(I), ones(size(I))*1, '|','Color','w', 'HorizontalAlignment', 'center', ...
    'fontsize', 12, 'fontweight', 'bold')

ax = gca;
yyaxis right
ax.YTick = [];
load('two_colors_colormap.mat');
colormap(two_colors_colormap);

end