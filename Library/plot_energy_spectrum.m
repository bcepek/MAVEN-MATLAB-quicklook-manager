function plot_energy_spectrum (x1, epoch_d1, energy_d1, swp_ind_d1, eflux_d1_cleaned_sum)
            
            figure
%             subplot(1, 2, 1)
            epoch_ind = find (abs(x1 - epoch_d1) == min(abs(x1 - epoch_d1)));
            %loglog((energy_d1(:,swp_ind_d1(epoch_ind) + 1, 2, 1)),  (squeeze(eflux_d1_cleaned_sum(:, 1, epoch_ind))),'k','LineWidth',1.5);
            loglog(energy_d1(:,swp_ind_d1(epoch_ind) + 1, 2, 1),...
                squeeze(sum(eflux_d1_cleaned_sum(:, 1, epoch_ind - 1 : epoch_ind + 1), 3)),...
                'k','LineWidth',1.5);
            hold on
            loglog(energy_d1(:,swp_ind_d1(epoch_ind) + 1, 2, 5),...
                squeeze(sum(eflux_d1_cleaned_sum(:, 5, epoch_ind -1 : epoch_ind + 1), 3)),...
                'LineWidth',1.5);
            loglog(energy_d1(:,swp_ind_d1(epoch_ind) + 1, 2, 6),...
                squeeze(sum(eflux_d1_cleaned_sum(:, 6, epoch_ind - 1 : epoch_ind + 1), 3)),...
                'r','LineWidth',1.5 );
            legend ('p','O^+','O_2^+')
            xlabel('energy, eV')
            ylabel('flux')
            grid on
            title(datestr(x1))
            ylim([10^3 10^8])
            xlim([1 10^5])

%             subplot(1, 2, 2)
%             plot()
            
end
