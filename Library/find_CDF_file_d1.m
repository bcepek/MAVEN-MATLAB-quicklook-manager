function filename_d1 = find_CDF_file_d1(date)
% finds STA CDF file for d1 product by date and returns full
% filename;
% written by Sergey Shuvalov, created as a separate function by Vladimir
% Ermakov;
% last modified June 13, 2018;

root = '\\193.232.6.100\Data\Maven\Static\d1_32e4d16a8m_Version2';

cur_date = datenum(date);

cur_year = num2str(year(cur_date));
cur_month = num2str(month(cur_date));
if(length(cur_month) < 2)
    cur_month = ['0' cur_month];
end
cur_day = num2str(day(cur_date));
if(length(cur_day) < 2)
    cur_day = ['0' cur_day];
end
monthpath = [root '\' cur_year '\' cur_month '\'];
month_filelist = dir([monthpath, '*.cdf']);

for i = 1:length(month_filelist)
    if(all(month_filelist(i).name(32:33) == cur_day))
        filename_d1 = [monthpath month_filelist(i).name];
        isfound = true;
        break
    elseif i == length(month_filelist)
        isfound = false;
        disp([datestr(cur_date) ' not found'])
        filename_d1 = -1;
    end
end
end
