timestart = datenum('2020-12-07 10:15:00','yyyy-mm-dd HH:MM:SS');
timeend = datenum('2020-12-07 10:30:00','yyyy-mm-dd HH:MM:SS');

filename = find_swea_arcpad(timestart);
epoch = spdfcdfread(filename, 'variables', 'epoch');
eflux = spdfcdfread(filename, 'variables', 'diff_en_fluxes');
energy = spdfcdfread(filename, 'variables', 'energy');
pa = spdfcdfread(filename, 'variables', 'pa');
d_pa = spdfcdfread(filename, 'variables', 'd_pa');
g_pa = spdfcdfread(filename, 'variables', 'g_pa');

verify_time = epoch>timestart & epoch<timeend;
epoch = epoch(verify_time);
eflux = eflux(:, :, verify_time);    % Pitch, Energy, Time
pa = pa(:,:,verify_time);

pa_lims = [0, 60; 60, 120; 120, 180];
%pa_lims = [0, 30; 30, 60; 60, 90; 90, 120; 120, 150; 150, 180];
f = figure();
set(f, 'Position', [2648 1563 520 619])
spectra = zeros(length(pa_lims), length(energy), length(epoch));
for pa_int_num = 1:size(pa_lims, 1)
    for timenum = 1:length(epoch)
        for enum = 1:length(energy)
            choose_pa = pa(:, enum, timenum) > pa_lims(pa_int_num, 1) &...
                pa(:, enum, timenum) < pa_lims(pa_int_num, 2);
            spectra(pa_int_num, enum, timenum) = spectra(pa_int_num, enum, timenum) +...
                sum(eflux(choose_pa, enum, timenum))/sum(choose_pa);
        end
    end
end

caxis_lim = log10([min(min(min(spectra(spectra~=0)))),...
    max(max(max(spectra)))]);

for pa_int_num = 1:size(pa_lims, 1)
    h = subplot(size(pa_lims, 1),1,pa_int_num);
    h.Position(4) = 1.1*h.Position(4);
    h.Position(1) = 0.5*h.Position(1);
    pcolor(epoch, energy, log10(squeeze(spectra(pa_int_num, :, :))))
    set(gca, 'gridalpha', 0.5, 'minorgridalpha', 0.5)
    ylabel('Energy [eV]')
    shading flat
    load('newColormap.mat')
    colormap(gca, cmap);
    set(gca, 'YScale', 'log')
    datetick('x')
    title(['Pitch-angles ' num2str(pa_lims(pa_int_num, 1)) '° - ' num2str(pa_lims(pa_int_num, 2)) '°'])
    ylim([27 4700])
    caxis(caxis_lim)
    grid on
    set(gca, 'layer', 'top')
    if(pa_int_num ~= size(pa_lims,1))
        set(gca, 'XTickLabel', [])
    else
        cbar_handle = colorbar;
        set(cbar_handle, 'Position', [0.8519, 0.1099, 0.0346, 0.8239])
        xlabel(cbar_handle, 'Differential energy flux [eV/(eV cm^2 s sr)]')
        ylabels = get(cbar_handle, 'ticklabels');
        for i=1:length(ylabels)
            ylabels{i} = ['10^{' num2str(ylabels{i}) '}'];
        end
        set(cbar_handle, 'ticklabels', ylabels)
    end
end