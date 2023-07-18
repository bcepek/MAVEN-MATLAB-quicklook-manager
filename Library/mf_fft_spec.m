timefrom = datenum('2019-10-24 02:00:00');
timeto = datenum('2019-10-24 03:30:00');
window_size = 320;

mf_filename = '\\193.232.6.100\Data\Maven\data for MATLAB\MAG\mvn_mag_l2_2019297ss_20191024_v01_r01.sts.mat';

%mf_data = load(mf_filename);
load(mf_filename);
%mf_epoch = mf_data(:, 7) + datenum(['00-Jan-', num2str(mf_data(1, 1)), ' 00:00:00']);

choose_ind = find(mf_epoch>timefrom & mf_epoch<timeto);
mf_epoch = mf_epoch(choose_ind);
B = [Bx By Bz B];
B = B(choose_ind, :);


%B = [mf_data(:, [8 9 10]), sqrt(sum(mf_data(:, [8 9 10]).^2, 2))];
spectrogram = zeros(window_size/2+1, size(B, 1)-window_size);
for i=1:size(B, 1)-window_size
    Y = fft(B(i:i+window_size-1, 4));
    spectrogram(:, i) = Y(1:window_size/2+1);
end

x = mf_epoch(window_size/2:size(mf_epoch, 1)-window_size/2-1);
y = 32*(0:(window_size/2))/window_size;
pcolor(x, y, log(abs(spectrogram/window_size)))
datetick('x','HH:MM:SS');
xlim([x(1) x(end)])
shading flat
colorbar
colormap(jet)
xlabel('UT, HH:MM:SS')
ylabel('Frequency, Hz')
set(gca, 'yscale', 'log')
% 
% %====накладываем гирочастоты
% B_mean = smooth(B(:, 4), window_size)*1e-9;
% aem = 1.67e-27;
% q = 1.6e-19;
% mass = [1, 16, 32];
% freq_c = zeros(size(B_mean, 1), length(mass));
% for i = 1:length(mass)
%     freq_c(:, i) = q*B_mean./(aem*mass(i));
% end
% 
% hold on
% for i=1:length(mass)
%     plot(mf_epoch, freq_c(:, i), 'color', 'blue')
% end
% %=====================
% B_mean = smooth(B(:, 4), window_size)*1e-9;
% aem = 1.67e-27;
% q = 1.6e-19;
% mass = [1, 16, 32];
% freq_c = zeros(size(B_mean, 1), length(mass));
% for i = 1:length(mass)
%     freq_c(:, i) = q*B_mean./(2*pi*aem*mass(i));
% end
% 
% hold on
% for i=1:length(mass)
%     plot(mf_epoch, freq_c(:, i), 'color', 'green')
% end
% hold off