% This script calculates the velocity of DeHoffmann-Teller frame with the
% usage of MAG data and proton velocity taken from cleaned STATIC data.
%
% WARNING. Boundaries of time interval must belong to the same day.
%
%
% =====SYNTAX=====
% INPUT:
% 1) time_start [1x1] double
% 2) time_end [1x1] double
%
% OUTPUT:
% 1) DeHoffmann-Teller velocity in km/s [1x3] double
%
%           Author: Sergey Shuvalov (shuvalovsergei@gmail.com,
%                                    +79153481805)
%           Last change date: 21.11.2018 14:22

function V = Vdht(time_start, time_end)
% clear
% time_start = datenum('20160104 12:00:00', 'yyyymmdd HH:MM:SS');
% time_end = datenum('20160104 13:00:00', 'yyyymmdd HH:MM:SS');
root = '\\193.232.6.100\общие\Владимир\MAVEN_data_in_MatLab_format_cleaned\';

timestamp = datestr(time_start, 'yyyymmdd');
file = dir([root 'mvn_sta_l2_d1-32e4d16a8m_' timestamp '_*.cdf_H.mat']);

if(length(file)<1)
    disp(['No MAG files found for date ', datestr(time_start, 'yyyy-mm-dd')])   % if no files found
    return
end

data = load([root '\' file.name]);
verify = data.epoch>=time_start & data.epoch<=time_end;
B = data.magf(verify, :);
u = data.v_mso(verify, :)';

K = zeros(3, 3, size(B, 1));
K(1, 1, :) = B(:, 2).^2 + B(:, 3).^2;
K(1, 2, :) = -B(:, 1).*B(:, 2);
K(1, 3, :) = -B(:, 1).*B(:, 3);
K(2, 1, :) = K(1, 2, :);
K(2, 2, :) = B(:, 1).^2 + B(:, 3).^2;
K(2, 3, :) = -B(:, 2).*B(:, 3);
K(3, 1, :) = K(1, 3, :);
K(3, 2, :) = K(2, 3, :);
K(3, 3, :) = B(:, 1).^2 + B(:, 2).^2;

meanvel = zeros(3, size(K, 3));
for i = 1:size(K, 3)
    meanvel(:, i) = K(:, :, i)*u(:, i);
end

V = mean(K, 3)\mean(meanvel, 2);
end