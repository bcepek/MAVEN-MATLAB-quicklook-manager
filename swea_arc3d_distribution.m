timestart = datenum('2020-12-07 10:15:00','yyyy-mm-dd HH:MM:SS');
timeend = datenum('2020-12-07 10:30:00','yyyy-mm-dd HH:MM:SS');

filename = find_swea_arc3d(timestart);
epoch = spdfcdfread(filename, 'variables', 'epoch');
eflux = spdfcdfread(filename, 'variables', 'diff_en_fluxes');
energy = spdfcdfread(filename, 'variables', 'energy');
azim = spdfcdfread(filename, 'variables', 'azim');
elev = spdfcdfread(filename, 'variables', 'elev');

verify_time = epoch>timestart & epoch<timeend;
epoch = epoch(verify_time);
eflux = eflux(:, :, :, verify_time);    % Azim, Energy, Elev, Time

energies_to_plot = 19:44;
times_to_plot = 51:71;
timelen = length(times_to_plot);
fsize = 12;

eflux_disp = eflux(:,energies_to_plot,:,times_to_plot);
caxis_lims = [min(min(min(min(log10(eflux_disp(eflux_disp~=0)))))),...
    max(max(max(max(log10(eflux_disp(eflux_disp~=0))))))];

f = figure(1);
set(f, 'position', [2628, 54, 2120, 1251])
plot_num = 0;
for enum = energies_to_plot
    for timenum = times_to_plot
        plot_num = plot_num + 1;
        ax = subplot(length(energies_to_plot), timelen, plot_num);
        plotdist = squeeze(log10(eflux(:,enum,:,timenum)))';
        plotdist(end+1, end+1) = 0;
        p = pcolor(plotdist);
        clim(caxis_lims)
        colormap jet
        ax.XTick = [];
        ax.YTick = [];
        ax.Position(3) = ax.Position(3)*1.1;
        ax.Position(4) = ax.Position(4)*1.1;
        shading flat
        if(enum == energies_to_plot(1))
            title(datestr(epoch(timenum), 'HH:MM:SS'), 'FontWeight', 'bold', 'fontsize', fsize)
        end
        if(timenum == times_to_plot(1))
            hl = ylabel([num2str(round(energy(enum))), 'eV'], 'fontweight', 'bold', 'rotation', 0,...
                    'FontSize', fsize, 'color', 'black');
        end
    end
end