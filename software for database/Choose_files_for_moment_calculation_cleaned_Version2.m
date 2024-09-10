clear
datefrom = '2022-12-19';
dateto   = '2022-12-19';
writepath = '/Users/sesh2112/data/Matlab_MAVEN/pre-calculated/STATIC_moments';
root = '/Users/sesh2112/data/Matlab_MAVEN/STATIC_d1';

date_from_num = datenum(datefrom);
date_to_num = datenum(dateto);


cur_date = date_from_num;
while cur_date <= date_to_num
    disp(['STA moments start converting ' datestr(cur_date, 'yyyy-mm-dd')])
    cur_year = num2str(year(cur_date));
    cur_month = num2str(month(cur_date));
    if(length(cur_month) < 2)
        cur_month = ['0' cur_month];
    end
    cur_day = num2str(day(cur_date));
    if(length(cur_day) < 2)
        cur_day = ['0' cur_day];
    end
    monthpath = [root '/' cur_year '/' cur_month '/'];
    month_filelist = dir([monthpath, '*.cdf']);
    
    for i = 1:length(month_filelist)
        if(all(month_filelist(i).name(32:33) == cur_day))
            filename = [monthpath month_filelist(i).name];
            isfound = true;
%             if month_filelist(i).bytes > 2*10^9
%                  copyfile(filename, 'E:\MatLab\temp\')
%                  filename = ['E:\MatLab\temp\' month_filelist(i).name];
%             end
            break
        elseif i == length(month_filelist)
            isfound = false;
            disp([datestr(cur_date) ' not found'])
        end
    end
    
    if(isfound)
        
%         info = cdfinfo(filename);
%         if info.FileSize < 2*10^9
            calc_date_moments_cleaned_2(filename, writepath);
%         else
%             copyfile(filename, 'E:\MatLab\temp\')
%         end
    end
    
    cur_date = cur_date + datenum('01-Jan-0000');
end