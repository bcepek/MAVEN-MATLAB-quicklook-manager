function calc_date_moments_cleaned_2(filename, writepath)
% written by Sergey and Vladimir
% clc
% clear
mass = [1 16 32];
titles = {'H^+', 'O^+', 'O_2^+', 'CO_2^+'};

q = 1.602177335e-19;
aem = 1.66054021010e-27;
k = 1.38064852e-23;

data = spdfcdfread(filename, 'variables', {'epoch',...
    'eflux',...
    'pos_sc_mso',...
    'energy',...
    'nenergy',...
    'denergy',...
    'theta',...
    'phi',...
    'domega',...
    'nbins',...
    'swp_ind',...
    'mass_arr',...
    'nanode',...
    'ndef',...
    'quat_mso',...
    'magf',...
    'sc_pot',...
    'att_ind',...
    'quality_flag',...
    'valid'});
epoch = data{1};
eflux  = data{2};
pos_sc_mso = data{3};
energy = data{4};
nenergy = data{5};
denergy = data{6};
theta = data{7}*pi/180;
phi = data{8}*pi/180;
domega = data{9};
swp_ind = data{11};
mass_arr = data{12};
nanode = data{13};
ndef = data{14};
quat_mso = data{15};
magf = data{16};
sc_pot = data{17};
att_ind = data{18};
quality_flag = data{19};
valid = data{20};

% eflux = load("eflux_mat.mat");
% eflux = eflux.eflux_new;

disp(['cdf ' filename ' read successfully'])

%---BEGIN---clean O+ mass channel----------
temp = eflux(:,:, 5,:) - 0.08*eflux(:,:, 1,:) >= 0;
eflux_O = eflux(:,:, 5,:);
eflux_p = eflux(:,:, 1,:);
eflux_O(temp) = eflux_O(temp) - 0.08*eflux_p(temp);
clear temp;
temp = eflux(:,:, 5,:) - 0.08*eflux(:,:, 1,:) < 0;
eflux_O(temp) = 0;
eflux(:,:, 5,:) = eflux_O;
clear temp
%---END---clean O+ mass channel----------

%---BEGIN---clean O2+ mass channel----------
temp = eflux(:,:, 6,:) - 0.08*eflux(:,:, 1,:) >= 0;
eflux_O2 = eflux(:,:, 6,:);
eflux_O2(temp) = eflux_O2(temp) - 0.08*eflux_p(temp);
clear temp;
temp = eflux(:,:, 6,:) - 0.08*eflux(:,:, 1,:) < 0;
eflux_O2(temp) = 0;
eflux(:,:, 6,:) = eflux_O2;
%---END---clean O2+ mass channel----------

for mass_id = 1:length(mass)
    switch mass_id
        case 1
            mass_num_range = 1;
            mass_num = 1;
        case 2
            mass_num_range = 5;
            mass_num = 5;
        case 3
            mass_num_range = 6;
            mass_num = 6;
            %         case 4
            %             mass_num_range = 7;
            %             mass_num = 7;
            %         case 5
            %             mass_num_range = 2;
            %             mass_num = 2;
    end
    %for mass_id = 1:1
    %     mass_num_range = 1:2;
    %     mass_num = mass_id;
    
    onemass_eflux = reshape( sum(eflux(:, :, mass_num_range, :), 3), [size(eflux, 1) size(eflux, 2) size(eflux, 4)] );
    v = 1*sqrt(2*q*energy./(aem*mass_arr));
    
    
    phsdensity = zeros(size(onemass_eflux));
    for i = 1:size(onemass_eflux, 3)
        m_sq = aem*permute(squeeze(mass_arr(:, swp_ind(i)+1, :, mass_num)), [2 1]);
        e_sq = permute(squeeze(energy(:, swp_ind(i)+1, :, mass_num)), [2 1]);
        phsdensity(:, :, i) = 1e4*0.5*onemass_eflux(:, :, i).*(m_sq./(e_sq*q)).^2;
    end
    
    %disp('phsdensity calculated')
    
    concentration = zeros(length(epoch), 1);
    v_st = zeros(length(epoch), 3);
    v_par = zeros(length(epoch), 1);
    v_perp = zeros(length(epoch), 1);
    temp = zeros(length(epoch), 1);
    Tpar = zeros(length(epoch), 1);
    Tperp = zeros(length(epoch), 1);
    %disp('start moments calculation')
    for timenum = 1:length(epoch)
        %          en_int = find(squeeze(energy(:, swp_ind(timenum)+1, 1, mass_num))>=30);
        %          for en = en_int(1):en_int(end)
        for en = 1:nenergy
            for nphi = 1:nanode
                for ntheta = 1:ndef
                    bin = ndef*(nphi-1)+ntheta; %CHECK
                    %                     i = [en, swp_ind(timenum)+1, bin, mass_num];
                    %                     volume = q*v(i(1),i(2),i(3),i(4))*domega(i(1),i(2),i(3),i(4))*denergy(i(1),i(2),i(3),i(4))/(aem*mass_arr(i(1),i(2),i(3),i(4)));
                    swp_i = swp_ind(timenum);
                    v_ccl = v(en, swp_i+1, bin, mass_num);
                    volume = q*v_ccl*domega(en,swp_i+1, bin, mass_num)*denergy(en, swp_i+1,bin, mass_num)/(aem*mass_arr(en, swp_i+1, bin, mass_num));
                    concentration(timenum) = concentration(timenum) + volume*phsdensity(bin, en, timenum);
                    v_st(timenum, 1) = v_st(timenum, 1) + volume*v_ccl*cos(  phi(en,swp_i+1,bin,mass_num))*cos(theta(en,swp_i+1,bin,mass_num))*phsdensity(bin, en, timenum);
                    v_st(timenum, 2) = v_st(timenum, 2) + volume*v_ccl*sin(  phi(en,swp_i+1,bin,mass_num))*cos(theta(en,swp_i+1,bin,mass_num))*phsdensity(bin, en, timenum);
                    v_st(timenum, 3) = v_st(timenum, 3) + volume*v_ccl*sin(theta(en,swp_i+1,bin,mass_num))*phsdensity(bin, en, timenum);
                    % {
                    % pitch calculation takes a lot of time, so it can be
                    % switched on when necessary
                    pitch = acos((...
                        cos(  phi(en, swp_i+1, bin, mass_num))*cos(theta(en, swp_i+1, bin, mass_num))*magf(timenum, 1) + ...
                        sin(  phi(en, swp_i+1, bin, mass_num))*cos(theta(en, swp_i+1, bin, mass_num))*magf(timenum, 2) + ...
                        sin(theta(en,swp_i+1,bin,mass_num))*magf(timenum, 3)...
                        )/sqrt(sum(magf(timenum, :).^2, 2)));
                        v_par(timenum) = v_par(timenum) + volume*v(en,swp_i+1,bin,mass_num)*phsdensity(bin, en, timenum)*cos(pitch);
                        v_perp(timenum) = v_perp(timenum) + volume*v(en,swp_i+1,bin,mass_num)*phsdensity(bin, en, timenum)*sin(pitch);
                    % }
                end
            end
        end
    end
    
    for i=1:3
        v_st(:, i) = v_st(:, i)./concentration;
    end
    v_mso = quatrotate(quatinv(quat_mso), v_st);
    v_par = v_par./concentration;
    v_perp = v_perp./concentration;
    
    for timenum = 1:length(epoch)
        %          en_int = find(squeeze(energy(:, swp_ind(timenum)+1, 1, mass_num))>=30);
        %          for en = en_int(1):en_int(end)
        for en = 1:nenergy
            for nphi = 1:nanode
                for ntheta = 1:ndef
                    bin = ndef*(nphi-1)+ntheta;
                    %                     i = [en, swp_ind(timenum)+1, bin, mass_num];
                    v_ccl = v(en,swp_i+1,bin,mass_num);
                    swp_i = swp_ind(timenum);
                    volume = q*v_ccl*domega(en,swp_i+1,bin,mass_num)*denergy(en,swp_i+1,bin,mass_num)/(aem*mass_arr(en,swp_i+1,bin,mass_num));
                    cur_v = [v_ccl*cos(phi(en,swp_i+1,bin,mass_num))*cos(theta(en,swp_i+1,bin,mass_num)),...
                        v_ccl*sin(phi(en,swp_i+1,bin,mass_num))*cos(theta(en,swp_i+1,bin,mass_num)),...
                        v_ccl*sin(theta(en,swp_i+1,bin,mass_num))];
                    % {
                    % pitch calculation takes a lot of time, so it can be
                    % switched on when necessary
                    pitch = acos((...
                        cos(phi(en,swp_i+1,bin,mass_num))*cos(theta(en,swp_i+1,bin,mass_num))*magf(timenum, 1) + ...
                        sin(phi(en,swp_i+1,bin,mass_num))*cos(theta(en,swp_i+1,bin,mass_num))*magf(timenum, 2) + ...
                        sin(theta(en,swp_i+1,bin,mass_num))*magf(timenum, 3)...
                        )./sqrt(sum(magf(timenum, :).^2, 2)));
                        % }
                        %temp(timenum) = temp(timenum) + aem*mass_arr(en,swp_i+1,bin,mass_num)*sum((cur_v-v_st(timenum)).^2)*phsdensity(bin, en, timenum)*volume;
                        temp(timenum) = temp(timenum) + aem*mass_arr(en,swp_i+1,bin,mass_num)*sum((cur_v-v_st(timenum,:)).^2,2)*phsdensity(bin, en, timenum)*volume;
                        Tpar(timenum) = Tpar(timenum) + aem*mass_arr(en,swp_i+1,bin,mass_num)*(v_ccl*cos(pitch)-v_par(timenum)).^2*phsdensity(bin, en, timenum)*volume;
                        Tperp(timenum) = Tperp(timenum) + aem*mass_arr(en,swp_i+1,bin,mass_num)*(v_ccl*sin(pitch)-v_perp(timenum)).^2*phsdensity(bin, en, timenum)*volume;
                end
            end
        end
    end
    temp = temp./(3*q*concentration);
    Tpar = Tpar./(3*q*concentration);
    Tperp = Tperp./(3*q*concentration);
    
    %{
            figure(1)
            subplot(3, 1, 1)
            semilogy(epoch, concentration/1e6, 'color', 'black', 'linewidth', 2)
            datetick('x','HH:MM:SS');
            set(gca, 'xticklabel', [])
            grid on
            ylabel('n, cm^{-3}')
            title(titles{mass_id})
    
            subplot(3, 1, 2)
            plot(epoch, v_mso(:, 1)/1e3, 'color', 'red', 'linewidth', 2)
            hold on
            plot(epoch, v_mso(:, 2)/1e3, 'color', 'green', 'linewidth', 2)
            plot(epoch, v_mso(:, 3)/1e3, 'color', 'blue', 'linewidth', 2)
            plot(epoch, sqrt(v_mso(:, 1).^2+v_mso(:, 2).^2+v_mso(:, 3).^2)/1e3, 'color', 'black', 'linewidth', 2)
            legend('V_x MSO', 'V_y MSO', 'V_z MSO', 'V')
            datetick('x','HH:MM:SS');
            hold off
            grid on
            ylabel('V, km/s')
    
            subplot(3, 1, 3)
            plot(epoch, temp, 'color', 'black', 'linewidth', 2)
            hold off
            datetick('x','HH:MM:SS');
            grid on
            ylabel('T, eV')
            xlabel('UT, HH:MM:SS')
    %}
    %data = [datevec(epoch), double(magf), double(pos_sc_mso), concentration/1e6, v_mso/1e3, temp, Tpar, Tperp];
    savename = [num2str(mass_id), 'mass'];
    %save('lowE', 'data')
    concentration = concentration/1e6;
    v_mso = v_mso/1e3;
    switch mass_id
        case 1
            save([writepath, '/', filename(53:end), '_H.mat'], 'epoch', 'magf', 'pos_sc_mso', 'concentration', 'v_mso', 'temp', 'Tpar', 'Tperp',...
                'sc_pot', 'att_ind', 'quality_flag', 'valid')
        case 2
            save([writepath, '/', filename(53:end), '_O.mat'], 'epoch', 'magf', 'pos_sc_mso', 'concentration', 'v_mso', 'temp', 'Tpar', 'Tperp',...
                'sc_pot', 'att_ind', 'quality_flag', 'valid')
        case 3
            save([writepath, '/', filename(53:end), '_O2.mat'], 'epoch', 'magf', 'pos_sc_mso', 'concentration', 'v_mso', 'temp', 'Tpar', 'Tperp',...
                'sc_pot', 'att_ind', 'quality_flag', 'valid')
            %         case 4
            %             save('060115_CO2', 'epoch', 'magf', 'pos_sc_mso', 'concentration', 'v_mso', 'temp', 'Tpar', 'Tperp')
            %         case 5
            %             save('050115_He',  'epoch', 'magf', 'pos_sc_mso', 'concentration', 'v_mso', 'temp', 'Tpar', 'Tperp')
    end
end


