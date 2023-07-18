function filename_c6 = find_CDF_file_c6(date)
% finds STA CDF file for c6 product by date and returns full filename
% written by Sergey Shuvalov, created as a separate function by Vladimir
% Ermakov
% last modified June 13, 2018

% адрес каталога с cdf-ками
cat = '\\193.232.6.100\Data Raw\Maven\Static\c6_32e64m';

% Ищем файл нужного дня
timestr = datestr(date, 'yyyy-mm-dd HH:MM:SS');

cur_year = timestr(1:4);
cur_month = timestr(6:7);
cur_day = timestr(9:10);

catpath = [cat '\' cur_year '\' cur_month '\'];

filelist = dir([catpath, '*.cdf']);
for filenum = 1:size(filelist, 1)
    if filelist(filenum).name(28:29) == cur_day
        name = filelist(filenum).name;
        break
    end
    if filenum == size(filelist, 1)
        disp('File not found')
        return
    end
end

% Формируем имя файла нужного дня
filename_c6 = [catpath, name];

end