function plot_energy_spectra_set (x1, epoch_d1, energy_d1, swp_ind_d1, eflux_d1_cleaned_sum)

h = figure;
%set(h, 'position', [-1003 147 905 803])
epoch_ind_start = find (abs(x1(1) - epoch_d1) == min(abs(x1(1) - epoch_d1)));
epoch_ind_end   = find (abs(x1(2) - epoch_d1) == min(abs(x1(2) - epoch_d1)));


subplot(1, 3, 1)
for i = epoch_ind_start : epoch_ind_end
    
    semilogx( (energy_d1(:,swp_ind_d1(i) + 1, 2, 1)),  log10(squeeze(eflux_d1_cleaned_sum(:, 1, i))) + ...
        (i - epoch_ind_start )* 1, 'k','LineWidth', 1);
    hold on
end
xlabel('Energy [eV]')
ylabel('Log differential energy flux [arb. units]')
xlim([1 1e4])
grid on
legend('p')

subplot(1, 3, 2)
for i = epoch_ind_start : epoch_ind_end
    semilogx( (energy_d1(:,swp_ind_d1(i) + 1, 2, 5)),  log10(squeeze(eflux_d1_cleaned_sum(:, 5, i))) + ...
        (i - epoch_ind_start )* 1, 'b', 'LineWidth', 1  );
    hold on
end
xlim([10^(-0.5) 10^(3.5)])
grid on
legend('O^+')
xlabel('Energy [eV]')
title('')

subplot(1, 3, 3)
for i = epoch_ind_start : epoch_ind_end
    semilogx( log10(energy_d1(:,swp_ind_d1(i) + 1, 2, 6)),  log10(squeeze(eflux_d1_cleaned_sum(:, 6, i))) + ...
        (i - epoch_ind_start )* 1, 'r','LineWidth',1 );
    hold on
end
xlim([10^(-0.5) 10^(3.5)])
grid on
legend('{O_2}^+')
xlabel('Energy [eV]')

end