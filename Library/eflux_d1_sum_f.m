function [eflux_d1_sum, domega_d1_sum] = eflux_d1_sum_f (eflux_d1, nenergy_d1, nbins_d1, domega_d1, swp_ind_d1 )

eflux_d1_sum = zeros(size(eflux_d1, 2), size(eflux_d1, 3), size(eflux_d1, 4));
domega_d1_sum = zeros(nenergy_d1, size(eflux_d1, 4), size(eflux_d1, 3));
for en = 1:nenergy_d1
    domega_d1_sum(en, :, :) = sum(domega_d1(en, swp_ind_d1+1, :, :),3);
    
    for time = 1:size(eflux_d1, 4)
        for bin = 1:nbins_d1
            eflux_d1_sum(en,:,time) = eflux_d1_sum(en,:,time)...
                + (  squeeze(eflux_d1(bin,en,:,time)).*squeeze(domega_d1(en,swp_ind_d1(time)+1,bin, :))  )';
        end
    end
end
domega_d1_sum = permute(domega_d1_sum, [1 3 2]);
eflux_d1_sum = eflux_d1_sum./domega_d1_sum;

end