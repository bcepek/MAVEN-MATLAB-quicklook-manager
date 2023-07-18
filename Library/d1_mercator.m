% %filename = '\\193.232.6.100\Data\Maven\Static\d1_32e4d16a8m\2015\01\mvn_sta_l2_d1-32e4d16a8m_20150104_v01_r12.cdf';
% clear
% date = datenum('2017-11-10 21:10:01','yyyy-mm-dd HH:MM:SS');
% filename = find_CDF_file_d1(date);
% 
% epoch = spdfcdfread(filename, 'variables', 'epoch');
% eflux = spdfcdfread(filename, 'variables', 'eflux');
% energy = spdfcdfread(filename, 'variables', 'energy');
% theta = spdfcdfread(filename, 'variables', 'theta');
% dtheta = spdfcdfread(filename, 'variables', 'dtheta');
% phi = spdfcdfread(filename, 'variables', 'phi');
% dphi = spdfcdfread(filename, 'variables', 'dphi');
% swp_ind = spdfcdfread(filename, 'variables', 'swp_ind');
% mass_arr = spdfcdfread(filename, 'variables', 'mass_arr');
% magf = spdfcdfread(filename, 'variables', 'magf');
% quat_mso = spdfcdfread(filename, 'variables', 'quat_mso');
% k=0;

%eflux = rand(size(eflux));
%eflux([26 27 30 31 34 35],:,:,:) = 0;
function h_plot = d1_mercator(timestart, timeend, mass_num, epoch, eflux, energy, swp_ind, magf, quat_mso, energies_to_plot)
timelen = 10;
start_time = find(abs(epoch - timestart) == min(abs(epoch - timestart)));
stop_time = find(abs(epoch - timeend) == min(abs(epoch - timeend)));
if(stop_time-start_time < timelen-1)
    stop_time = start_time + timelen - 1;
end

load('newColormap.mat', 'cmap');
%===============================
Phi = zeros(5, 7);
Phi(:, 1) = [-8.1158724;-4.0579362;0;4.0579362;8.1158724]*pi/180;
Phi(:, 2) = [-11.25812535;-5.62906265;0;5.62906275;11.25812545]*pi/180;
Phi(:, 3) = [-15.6169737;-7.8084863;0;7.8084885;15.6169759]*pi/180;
Phi(:, 4) = [-21.6634545;-10.8317275;0;10.8317265;21.6634535]*pi/180;
Phi(:, 5) = [-30.050973;-15.025485;0;15.025491;30.050979]*pi/180;
Phi(:, 6) = [-41.685922;-20.842962;0;20.842958;41.685918]*pi/180;
Phi(:, 7) = [-45.866664;-22.933332;0;22.933332;45.866664]*pi/180;

Lambda = [-168.75000;-146.25000;-123.75000;-101.25000;-78.750000;-56.250000;-33.750000;-11.250000;11.250000;33.750000;56.250000;78.750000;101.25000;123.75000;146.25000;168.75000;191.25000]*pi/180;
X = zeros(size(Phi, 1), length(Lambda), 7);
Y = zeros(size(Phi, 1), length(Lambda), 7);
for en = 1:7
    for i = 1:length(Lambda)
        for j = 1:size(Phi, 1)
            X(j, i, en) = 2*sqrt(2)*cos(Phi(j, en))*sin(Lambda(i)/2)/sqrt(1+cos(Phi(j, en))*cos(Lambda(i)/2));
            Y(j, i, en) = sqrt(2)*sin(Phi(j, en))/sqrt(1+cos(Phi(j, en))*cos(Lambda(i)/2));
        end
    end
end
dotcont = zeros(100, 2);
Phi = linspace(-pi/2, pi/2, 100);
Lambda = [-168.75, 191.25]*pi/180;
c=1;
for i=1:length(Lambda)
    for j=1:length(Phi)
        dotcont(c, 1) = 2*sqrt(2)*cos(Phi(j))*sin(Lambda(i)/2)/sqrt(1+cos(Phi(j))*cos(Lambda(i)/2));
        dotcont(c, 2) = sqrt(2)*sin(Phi(j))/sqrt(1+cos(Phi(j))*cos(Lambda(i)/2));
        c=c+1;
    end
end
%clear Phi Lambda i j
%===============================

eflux_colordata = eflux(:, :, :, start_time:stop_time);
verify = eflux_colordata ~=0;
colordata = [min(min(min(min(min(log10(eflux_colordata(verify))))))), max(max(max(max(max(log10(eflux_colordata(verify)))))))];

%while(start_time+timelen-1 <= stop_time)
while(start_time <= stop_time)
    tft = [start_time, start_time+timelen-1];
    start_time = start_time + timelen;

    h_plot = figure();
    set(h_plot, 'pos', [0 0 1920 1280])
    
    eflux_disp = eflux(:, :, :, tft(1):1:tft(2));
    c = 1;
    
    for enum=energies_to_plot
        for timenum =1:timelen
            h = subplot(16, timelen, c);
            fdist = zeros(5, 17);
            fdist(1, 1:16) = eflux_disp(1:4:64, enum, mass_num, timenum);
            fdist(2, 1:16) = eflux_disp(2:4:64, enum, mass_num, timenum);
            fdist(3, 1:16) = eflux_disp(3:4:64, enum, mass_num, timenum);
            fdist(4, 1:16) = eflux_disp(4:4:64, enum, mass_num, timenum);
%             if(enum>7)
%                 en=7;
%             else
%                 en = enum;
%             end
            pcolor(log10(fdist));
            if(any(c==1:timelen:enum*timelen))
                hl = ylabel([num2str(round(energy(enum, swp_ind(timenum+tft(1)-1)+1, 1, mass_num))), 'eV'], 'fontweight', 'bold', 'rotation', 0);
                p = get(hl, 'pos');
                p(1) = p(1) - 2;
                set(hl, 'pos', p)
            end
            caxis(colordata)
            colormap(cmap)
            hold on
            hold off
            axis equal
            
            p = get(h, 'pos');
            p(2) = p(2) + 0.04;
            p(4) = p(4)*1.3;
            p(2) = p(2) - (enum-1)*0.0015;
            p(2) = p(2) + 0.006;
            p(3) = p(3)*1.3;
            set(h, 'pos', p)
            
            shading flat
            set(gca, 'xtick', [], 'ytick', [])
            if(enum == energies_to_plot(1))
                title(datestr(epoch(timenum+tft(1)-1), 'HH:MM:SS'), 'FontWeight', 'bold')
            end
            
            %                 if (enum == 32)
            %                     set(gca, 'xtick', 8.5, 'xticklabel', num2str(swp_ind(timenum+tft(1)-1)+1), 'fontsize', 6)
            %                 end
            
            
            
            hold off
            
            
            c = c+1;
        end
    end
    
    
    bar_handle = colorbar;
    set(bar_handle, 'pos', [0.918    0.11    0.027  0.867])
    labels = get(bar_handle, 'yticklabel');
    barlabels = cell(size(labels, 1), 1);
    for i=1:size(labels, 1)
        barlabels{i} = ['10^{', labels{i}, '}'];
    end
    set(bar_handle, 'yticklabel', char(barlabels), 'FontWeight', 'bold', 'fontsize', 10)
    switch mass_num
        case 1
            ylabel(bar_handle, 'H^+ Differential energy flux')
        case 5
            ylabel(bar_handle, 'O^+ Differential energy flux')
        case 6
            ylabel(bar_handle, 'O_2^+ Differential energy flux')
    end
end