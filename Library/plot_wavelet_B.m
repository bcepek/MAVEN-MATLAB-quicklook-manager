function plot_wavelet_B(TimeStart, TimeEnd)

load(find_mag_file(TimeStart), 'B', 'mf_epoch', 'Bx', 'By', 'Bz');

% Wavlet преобразование магнитного поля
figure()
subplot(1,2,1)
choose_ind = find(mf_epoch >= TimeStart & mf_epoch <= TimeEnd);
[wt,f] = cwt(B(choose_ind),32);
% [wtx,~] = cwt(Bx(choose_ind),32);
% [wty,~] = cwt(By(choose_ind),32);
% [wtz,~] = cwt(Bz(choose_ind),32);
% 
% wt = abs(wt)-abs(wtx)-abs(wty)-abs(wtz);

pcolor(repmat(mf_epoch(choose_ind),1,length(f))',repmat(f,1,length(mf_epoch(choose_ind))),abs(wt).^2);
mfepoch = mf_epoch(choose_ind);
shading flat
colormap jet
c = colorbar;
c.Label.String = 'Amplitude^2';
c.FontSize = 14; c.FontWeight = 'bold';
ylabel({'Frequency','Hz'},'fontsize',14,'fontweight','bold')
datetick('x','HH:MM');
set(gca, 'YScale', 'log', 'ycolor', 'black');
xlim([mfepoch(1) mfepoch(end)])
ylim([min(f) max(f)])

% Циклотронные частоты 
hold on
B_mean = smooth(B(choose_ind)*1e-9,500);
aem = 1.67e-27;
q = 1.6e-19;
mass = [1, 16, 32];
freq_c = zeros(size(B_mean, 1), length(mass));
for i = 1:length(mass)
    freq_c(:, i) = q*B_mean./(aem*mass(i));
end
col = ['r' 'g' 'b'];
for i=1:length(mass)
    plot(mf_epoch(choose_ind), freq_c(:, i), 'color', col(i))
end
hold off

subplot(1,2,2)
wt2 = abs(wt).^2;
dt = 1/32;
for i = 1:size(wt,1)
    for j = 1:(size(wt,2)-1)
        GWS(i,j) = (wt2(i,j+1)+wt2(i,j+1))/2*dt; 
    end
end
PS = [GWS(:,1) cumsum(GWS,2)];

plot(PS(:,end)/max(PS(:,end)),f,'k','LineWidth',2)
ylim([min(f) max(f)])
ylabel({'Frequency','Hz'},'fontsize',14,'fontweight','bold')
xlabel('Normalized GWS','fontsize',14,'fontweight','bold')
grid on 

%====накладываем гирочастоты
B_mean = mean(B(choose_ind))*1e-9;
% B_mean = norm(B_IMF)*1e-9;
aem = 1.67e-27;
q = 1.6e-19;
mass = [1, 16, 32];
freq_c = zeros(size(B_mean, 1), length(mass));
for i = 1:length(mass)
    freq_c(1, i) = q*B_mean./(aem*mass(i));
    freq_c(2, i) = q*std(B(choose_ind)*1e-9)./(aem*mass(i));
end

yline(freq_c(1,1),'r','H^{+}')
yline(freq_c(1,1)+freq_c(2,1),'--r',['\omega_{c,H^{+}} = ' num2str(freq_c(1,1)+freq_c(2,1)) ' Hz'],'LineWidth',2)
yline(freq_c(1,1)-freq_c(2,1),'--r',['\omega_{c,H^{+}} = ' num2str(freq_c(1,1)-freq_c(2,1)) ' Hz'],'LineWidth',2)
hold on
fill([0 0 1 1], [freq_c(1,1)-freq_c(2,1) freq_c(1,1)+freq_c(2,1) freq_c(1,1)+freq_c(2,1) freq_c(1,1)-freq_c(2,1)], [1 0 0], 'facealpha', 0.1, 'edgealpha', 0)

y1 = yline(freq_c(1,2),'-','O^{+}');
y1.Color = [0.4660 0.6740 0.1880];
y1 = yline(freq_c(1,2)+freq_c(2,2),'--',['\omega_{c,O^{+}} = ' num2str(freq_c(1,2)+freq_c(2,2)) ' Hz'],'LineWidth',2);
y1.Color = [0.4660 0.6740 0.1880];
yline(freq_c(1,2)-freq_c(2,2),'--',['\omega_{c,O^{+}} = ' num2str(freq_c(1,2)-freq_c(2,2)) ' Hz'],'LineWidth',2);
y1.Color = [0.4660 0.6740 0.1880];
fill([0 0 1 1], [freq_c(1,2)-freq_c(2,2) freq_c(1,2)+freq_c(2,2) freq_c(1,2)+freq_c(2,2) freq_c(1,2)-freq_c(2,2)], [0.4660 0.6740 0.1880], 'facealpha', 0.1, 'edgealpha', 0)


yline(freq_c(1,3),'b','O_2^{+}')
yline(freq_c(1,3)+freq_c(2,3),'--b',['\omega_{c,O_2^{+}} = ' num2str(freq_c(1,3)+freq_c(2,3)) ' Hz'],'LineWidth',2)
yline(freq_c(1,3)-freq_c(2,3),'--b',['\omega_{c,O_2^{+}} = ' num2str(freq_c(1,3)-freq_c(2,3)) ' Hz'],'LineWidth',2)
fill([0 0 1 1], [freq_c(1,3)-freq_c(2,3) freq_c(1,3)+freq_c(2,3) freq_c(1,3)+freq_c(2,3) freq_c(1,3)-freq_c(2,3)], [0 0 1], 'facealpha', 0.1, 'edgealpha', 0)

set(gca,'yscale','log')
end