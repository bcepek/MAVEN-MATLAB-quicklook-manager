%load('mvn_mag_l2_2017309pc1s_20171105_v01_r01.sts.mat');

% TimeStart = datenum('2017-11-05 14:00:00', 'yyyy-mm-dd HH:MM:SS');
% TimeEnd = datenum('2017-11-05 15:00:00', 'yyyy-mm-dd HH:MM:SS');

function [mf_epoch, lon, lat]=plot_subsat_point(TimeStart, TimeEnd)
root = '/Users/sesh2112/data/Matlab_MAVEN/pre-calculated/MAG_PC';
files = dir([root '/mvn_mag_l2_*' datestr(TimeStart, 'yyyymmdd')...
    '_v*_r*.sts.mat']);
if(isempty(files))
    disp(['MAG PC file for ' datestr(TimeStart, 'yyyy-mm-dd')...
        ' not found'])
    return
else
    load([root '/' files(1).name], 'pos_sc_mso', 'mf_epoch')
end

verify = TimeStart<mf_epoch & mf_epoch<TimeEnd;
mf_epoch = mf_epoch(verify);
pos_sc_mso = pos_sc_mso(verify, :);

%-----LONGITUDE----------
lon = -atan(pos_sc_mso(:,2)./pos_sc_mso(:,1));

verify = pos_sc_mso(:,1)<0;
lon(verify) = lon(verify) + pi;

verify = pos_sc_mso(:,1)>0 & pos_sc_mso(:,2)>0;
lon(verify) = lon(verify) + 2*pi;

lon = lon*180/pi;
%------------------------

%-----LATIITUDE----------
lat = atan(pos_sc_mso(:, 3)./sqrt(sum(pos_sc_mso(:, [1 2]).^2,2)))*180/pi;
%------------------------

    % searching for markers every 5 min
time2sec = 60*minute(mf_epoch)+second(mf_epoch);
time2sec(second(mf_epoch)>=59.5) = time2sec(second(mf_epoch)>=59.5)-60;
searchsec = (0:5:55)*60;
tmp = abs(time2sec-searchsec);
[row, ~] = find(tmp<1);
row(diff(row)==1)=[];

crust = imread('Crust.png');
figure()
imshow(crust)
lon_im = lon*size(crust, 2)/360;
lat_im = -(lat-90)*size(crust, 1)/180;
hold on
plot(lon_im, lat_im, '.', 'color', 'blue', 'markersize', 14)
plot(lon_im(1), lat_im(1), '.', 'color', 'green', 'markersize', 50)
plot(lon_im(end), lat_im(end), '.', 'color', 'red', 'markersize', 50)
plot(lon_im(row), lat_im(row), '.', 'color', 'white', 'markersize', 30)
hold off

end