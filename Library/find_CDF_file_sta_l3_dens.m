function filename = find_CDF_file_sta_l3_dens(cur_date)
% finds STA CDF file for d1 product by date and returns full
% filename;
% written by Sergey Shuvalov, created as a separate function by Vladimir
% Ermakov;
% last modified June 13, 2018;

load('paths.mat', 'paths', 'slash');
root = paths.sta_l3;

%root = '/Users/sesh2112/data/Matlab_MAVEN/STATIC_d1';

cur_year = datestr(cur_date, 'yyyy');
cur_month = datestr(cur_date, 'mm');
cur_day = datestr(cur_date, 'dd');

monthpath = [root slash cur_year slash cur_month slash];
month_filelist = dir([monthpath, '*.cdf']);

if isempty(month_filelist)
    filename = -1;
end

for i = 1:length(month_filelist)
    if(all(month_filelist(i).name(22:23) == cur_day))
        filename = [monthpath month_filelist(i).name];
        isfound = true;
        break
    elseif i == length(month_filelist)
        isfound = false;
        disp([datestr(cur_date) ' not found'])
        filename = -1;
    end
end


% isfound = false;
% disp([datestr(cur_date) ' not found'])
% filename_d1 = -1;
end
