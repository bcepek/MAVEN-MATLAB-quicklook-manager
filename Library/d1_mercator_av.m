
function h_plot = d1_mercator_av(timestart, timeend, mass_num, epoch, eflux, energy, swp_ind)

% timestart = datenum('2019-08-03 07:30:45', 'yyyy-mm-dd HH:MM:SS');
% timeend = datenum('2019-08-03 07:36:45', 'yyyy-mm-dd HH:MM:SS');

% timestart = datenum('2019-07-29 03:19:37', 'yyyy-mm-dd HH:MM:SS');
% timeend = datenum('2019-07-29 03:26:57', 'yyyy-mm-dd HH:MM:SS');

skip_par = 1;
energies_to_plot = 32; 
timelen = 7;
start_time = find(abs(epoch - timestart) == min(abs(epoch - timestart)));
stop_time = find(abs(epoch - timeend) == min(abs(epoch - timeend)));
if(stop_time-start_time < timelen*skip_par-1)
    stop_time = start_time + timelen*skip_par - 1;
    if(stop_time>length(epoch))
        stop_time = length(epoch);
    end
end

load('newColormap.mat', 'cmap');
% ЦВЕТОВАЯ ШКАЛА
eflux_colordata = eflux(:, 1:2:end, :, start_time:stop_time)+eflux(:, 2:2:end, :, start_time:stop_time);
verify = eflux_colordata ~=0;
colordata = [min(min(min(min(min(log10(eflux_colordata(verify))))))), max(max(max(max(max(log10(eflux_colordata(verify)))))))];

% ПОСТРОЕНИЕ ГРАФИКОВ 
%while(start_time+timelen*skip_par-1 <= stop_time)
while(start_time <= stop_time)
    tft = [start_time, start_time+timelen*skip_par-1];
    start_time = start_time + timelen*skip_par;

    h_plot = figure();
    %set(h_plot, 'pos', [0 0 1920 1280])
    set(h_plot, 'pos', [193 253 1176 626])
    
    eflux_disp = eflux(:, :, :, tft(1):1:tft(2));
    c = 1;
    
    %for enum = 1:2:energies_to_plot
    for enum = 3:2:28
        for timenum =1:skip_par:timelen*skip_par
            %h = subplot(16, timelen, c);
            h = subplot(13, timelen, c);
            fdist = zeros(5, 17);
            fdist(1, 1:16) = eflux_disp(1:4:64, enum, mass_num, timenum) + eflux_disp(1:4:64, enum+1, mass_num, timenum);
            fdist(2, 1:16) = eflux_disp(2:4:64, enum, mass_num, timenum) + eflux_disp(2:4:64, enum+1, mass_num, timenum);
            fdist(3, 1:16) = eflux_disp(3:4:64, enum, mass_num, timenum) + eflux_disp(3:4:64, enum+1, mass_num, timenum);
            fdist(4, 1:16) = eflux_disp(4:4:64, enum, mass_num, timenum) + eflux_disp(4:4:64, enum+1, mass_num, timenum);

            pcolor(log10(fdist));
            if(any(c==1:timelen:enum*timelen))
                energy_av = round((energy(enum, swp_ind(timenum+tft(1)-1)+1, 1, mass_num) + energy(enum+1, swp_ind(timenum+tft(1)-1)+1, 1, mass_num))/2);
                hl = ylabel([num2str(energy_av) , ' eV'], 'fontweight', 'bold', 'rotation', 0);
                p = get(hl, 'pos');
                p(1) = p(1) - 2;
                set(hl, 'pos', p)
            end
            %caxis([colordata(1) 10^(7.7)])
            %caxis(colordata)
            caxis([colordata(1) 7.5])
            colormap(cmap)
            hold on
            hold off
            axis equal
            box on
            
            ylim([1 5])
            
            p = get(h, 'pos');
%             p(2) = p(2) + 0.04;
%             p(4) = p(4)*1.3;
%             p(2) = p(2) - (enum-1)*0.0015;
%             p(2) = p(2) + 0.006;
%             p(3) = p(3)*1.3;
%             p(4) = p(4)*2;

            p(2) = p(2)*1;
            p(3) = p(3)*1.2;
            p(4) = p(4)*2;
            
            set(h, 'pos', p)
            
            
            
            shading flat
            set(gca, 'xtick', [], 'ytick', [])
            if (enum == 3)
                title(datestr(epoch(timenum+tft(1)-1), 'HH:MM:SS'), 'FontWeight', 'bold')
            end
            if(enum==energies_to_plot-1 && timenum==-skip_par+1+round((timelen*skip_par)/2))
                xlabel(datestr(epoch(timenum+tft(1)-1), 'yyyy-mmm-dd'), 'FontWeight', 'bold', 'fontsize', 14)
            end
                
                     
            hold off
                       
            c = c+1;
        end
    end
    
    % ПОДПИСЬ ЦВЕТОВОЙ ШКАЛЫ СПРАВА
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

