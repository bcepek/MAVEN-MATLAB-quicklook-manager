function h_plot = d1_hammer(timestart, timeend, mass_num, epoch, eflux, energy, swp_ind, magf, quat_mso, energies_to_plot)

fsize = 10;
mf_markersize = 7;
sw_markersize = 15;
antisw_markersize = 5;

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

Phi = flipud(Phi);

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
dotcont(100:200,:) = flipud(dotcont(100:200,:));
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
    set(h_plot, 'pos', [0 0 1417 928])

    eflux_disp = eflux(:, :, :, tft(1):1:tft(2));
    quat_mso_disp = quat_mso(tft(1):tft(2), :);
    c = 1;

    for enum=energies_to_plot
        for timenum =1:timelen
            h = subplot(16, timelen, c);
            fdist = zeros(5, 17);
            fdist(1, 1:16) = eflux_disp(1:4:64, enum, mass_num, timenum);
            fdist(2, 1:16) = eflux_disp(2:4:64, enum, mass_num, timenum);
            fdist(3, 1:16) = eflux_disp(3:4:64, enum, mass_num, timenum);
            fdist(4, 1:16) = eflux_disp(4:4:64, enum, mass_num, timenum);
            if(enum>7)
                en=7;
            else
                en = enum;
            end
            pcolor(squeeze(X(:, :, en)), squeeze(Y(:, :, en)), log10(fdist));
            shading flat
            ax = gca;
            set(ax, 'xtick', [], 'ytick', [], 'color', 'none')
            ax.XColor = 'none';
            ax.YColor = 'none';
            if(any(c==1:timelen:enum*timelen))
                hl = ylabel([num2str(round(energy(enum, swp_ind(timenum+tft(1)-1)+1, 1, mass_num))), 'eV'], 'fontweight', 'bold', 'rotation', 0,...
                    'FontSize', fsize, 'color', 'black');
                p = get(hl, 'pos');
                p(1) = p(1) - 2;
                set(hl, 'pos', p)
            end
            caxis(colordata)
            colormap(cmap)
            hold on
            plot(dotcont(1:end-1, 1), dotcont(1:end-1, 2), 'color', 'black', 'LineWidth',1.5)
            hold off
            axis equal

            p = get(h, 'pos');
            % p(2) = p(2) + 0.04;
            % p(4) = p(4)*1.3;
            % p(2) = p(2) - (enum-1)*0.0015;
            % p(2) = p(2) + 0.006;
            % p(3) = p(3)*1.3;

            p(2) = p(2) + 0.04;
            p(2) = p(2) - (enum-1)*0.004;
            p(2) = p(2) + 0.006;
            p(3) = p(3)*1.3;
            p(4) = p(4)*1.3;
            
            
            set(h, 'pos', p)
            %ylim(1.0005*[min(dotcont(:,2)), max(dotcont(:,2))])



            if(enum == energies_to_plot(1))
                title(datestr(epoch(timenum+tft(1)-1), 'HH:MM:SS'), 'FontWeight', 'bold', 'fontsize', fsize)
            end

            %Magnetic field
            mf = magf (tft(1):tft(2), :);
            Bx = mf(timenum, 1); By = mf(timenum, 2); Bz = mf(timenum, 3);
            %B = quatrotate(angle2quat(pi, 0, 0, 'ZYX'), [Bx, By, Bz]);
            %Bx = -B(1); By = -B(2); Bz = -B(3);
            B = sqrt(Bx^2 + By^2 + Bz^2);

            thetam = asin(Bz/B);
            if (Bx==0)
                if (By>0)
                    lambdam = pi/2;
                elseif (By<0)
                    lambdam = -pi/2;
                else
                    lambdam = 0; %!
                end
            else
                lambdam = atan (By/Bx);
                if ((By>=0)&&(Bx<0))
                    lambdam = lambdam + pi;
                end
                if ((By<0)&&(Bx<0))
                    lambdam = lambdam - pi;
                end
            end

            hold on
            Xmf = 2*sqrt(2)*cos(thetam)*sin(lambdam/2)/sqrt(1+cos(lambdam)*cos(lambdam/2));
            Ymf = sqrt(2)*sin(thetam)/sqrt(1+cos(thetam)*cos(lambdam/2));
            plot(Xmf, Ymf, 'pentagram', 'color', 'red', 'markersize', mf_markersize, 'linewidth', 1.5)
            %End.Magnetic field

            soldir = quatrotate(quat_mso_disp(timenum, :), [1 0 0]);
            %soldir = quatrotate(angle2quat(pi, 0, 0, 'ZYX'), soldir);
            thetas = asin(soldir(3));
            if (soldir(1)==0)
                if (soldir(2)>0)
                    lambdas = pi/2;
                elseif (soldir(2)<0)
                    lambdas = -pi/2;
                else
                    lambdas = 0; %!
                end
            else
                lambdas = atan (soldir(2)/soldir(1));
                if ((soldir(2)>=0)&&(soldir(1)<0))
                    lambdas = lambdas + pi;
                end
                if ((soldir(2)<0)&&(soldir(1)<0))
                    lambdas = lambdas - pi;
                end
            end
            Xs = 2*sqrt(2)*cos(thetas)*sin(lambdas/2)/sqrt(1+cos(lambdas)*cos(lambdas/2));
            Ys = sqrt(2)*sin(thetas)/sqrt(1+cos(thetas)*cos(lambdas/2));
            plot(Xs, Ys, '.', 'color', 'red', 'markersize', sw_markersize, 'linewidth', 2)

            soldir = quatrotate(quat_mso_disp(timenum, :), [-1 0 0]);
            %soldir = quatrotate(angle2quat(pi, 0, 0, 'ZYX'), soldir);
            thetas = asin(soldir(3));
            if (soldir(1)==0)
                if (soldir(2)>0)
                    lambdas = pi/2;
                elseif (soldir(2)<0)
                    lambdas = -pi/2;
                else
                    lambdas = 0; %!
                end
            else
                lambdas = atan (soldir(2)/soldir(1));
                if ((soldir(2)>=0)&&(soldir(1)<0))
                    lambdas = lambdas + pi;
                end
                if ((soldir(2)<0)&&(soldir(1)<0))
                    lambdas = lambdas - pi;
                end
            end
            Xs = 2*sqrt(2)*cos(thetas)*sin(lambdas/2)/sqrt(1+cos(lambdas)*cos(lambdas/2));
            Ys = sqrt(2)*sin(thetas)/sqrt(1+cos(thetas)*cos(lambdas/2));
            plot(Xs, Ys, 'o', 'color', 'red', 'markersize', antisw_markersize, 'linewidth', 1.5)



            %                 if (enum == 32)
            %                     set(gca, 'xtick', 8.5, 'xticklabel', num2str(swp_ind(timenum+tft(1)-1)+1), 'fontsize', 6)
            %                 end



            hold off

            c = c+1;
        end
    end


    bar_handle = colorbar;
    %set(bar_handle, 'pos', [0.918    0.11    0.027  0.867])
    set(bar_handle, 'pos', [0.929291460832745,0.132543103448276,0.014956951305575,0.851293103448276])
    labels = get(bar_handle, 'yticklabel');
    barlabels = cell(size(labels, 1), 1);
    for i=1:size(labels, 1)
        barlabels{i} = ['10^{', labels{i}, '}'];
    end
    set(bar_handle, 'yticklabel', char(barlabels), 'FontWeight', 'bold', 'fontsize', fsize)
    switch mass_num
        case 1
            ylabel(bar_handle, 'H^+ Differential energy flux')
        case 5
            ylabel(bar_handle, 'O^+ Differential energy flux')
        case 6
            ylabel(bar_handle, 'O_2^+ Differential energy flux')
    end
end