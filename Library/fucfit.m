function fucfit (x1, epoch_d1, energy_d1, swp_ind_d1, eflux_d1_cleaned_sum)
            
            figure
%             subplot(1, 2, 1)
            epoch_ind = find (abs(x1 - epoch_d1) == min(abs(x1 - epoch_d1)));
            %loglog((energy_d1(:,swp_ind_d1(epoch_ind) + 1, 2, 1)),  (squeeze(eflux_d1_cleaned_sum(:, 1, epoch_ind))),'k','LineWidth',1.5);
            
%             loglog(energy_d1(:,swp_ind_d1(epoch_ind) + 1, 2, 1),...
%                 squeeze(sum(eflux_d1_cleaned_sum(:, 1, epoch_ind - 1 : epoch_ind + 1), 3)),...
%                 'k','LineWidth',1.5);
            x = energy_d1(:,swp_ind_d1(epoch_ind) + 1, 2, 1);
            y = squeeze(sum(eflux_d1_cleaned_sum(:, 1, epoch_ind - 1 : epoch_ind + 1), 3));
            f = fit(x,y,'gauss1');
            loglog(x,y,'k','marker', '.', 'markersize', 14, 'linestyle', 'none')
            hold on
            locX = linspace(x(1), x(end), 10000);
            p1 = loglog(locX,f(locX), 'k', 'linewidth', 1.5, 'marker', 'none');
            
%             loglog(energy_d1(:,swp_ind_d1(epoch_ind) + 1, 2, 5),...
%                 squeeze(sum(eflux_d1_cleaned_sum(:, 5, epoch_ind - 2 : epoch_ind + 2), 3)),...
%                 'LineWidth',1.5);
            x = energy_d1(:,swp_ind_d1(epoch_ind) + 1, 2, 5);
            y = squeeze(sum(eflux_d1_cleaned_sum(:, 5, epoch_ind - 2 : epoch_ind + 2), 3));
            f = fit(x,y,'gauss1');
            loglog(x,y,'marker', '.', 'markersize', 14, 'linestyle', 'none', 'color', [0.85 0.33 0.10])
            locX = linspace(x(1), x(end), 10000);
            p2 = loglog(locX,f(locX), 'linewidth', 1.5, 'marker', 'none', 'color', [0.85 0.33 0.10]);
            
%             loglog(energy_d1(:,swp_ind_d1(epoch_ind) + 1, 2, 6),...
%                 squeeze(sum(eflux_d1_cleaned_sum(:, 6, epoch_ind - 2 : epoch_ind + 2), 3)),...
%                 'r','LineWidth',1.5 );
            x = energy_d1(:,swp_ind_d1(epoch_ind) + 1, 2, 6);
            y = squeeze(sum(eflux_d1_cleaned_sum(:, 6, epoch_ind - 2 : epoch_ind + 2), 3));
            f = fit(x,y,'gauss1');
            loglog(x,y,'marker', '.', 'markersize', 14, 'linestyle', 'none', 'color', [1 0 0])
            locX = linspace(x(1), x(end), 10000);
            p3 = loglog(locX,f(locX), 'linewidth', 1.5, 'marker', 'none', 'color', [1 0 0]);
            
            legend ([p1 p2 p3], {'p','O^+','O_2^+'})
            xlabel('energy, eV')
            ylabel('eflux, eV/(eV cm^2 s sr)')
            grid on
            title(datestr(x1))
            ylim([10^3 10^8])
            xlim([1 10^5])

%             subplot(1, 2, 2)
%             plot()
            
end
