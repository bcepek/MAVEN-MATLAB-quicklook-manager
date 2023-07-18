clear
datefrom = '2020-12-01';
dateto = '2020-12-31';
writepath = '/Users/sesh2112/data/Matlab_MAVEN/pre-calculated/MAG_PC';
root = '/Users/sesh2112/data/Matlab_MAVEN/MAG';

from = datenum(datefrom);
to = datenum(dateto);
i = 0;
cur_date = from;
while cur_date<=to
    disp(['MAG Start converting ' datestr(cur_date, 'yyyy-mm-dd')])
    isfound = 0;
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
    
    file_to_convert = dir([monthpath, '*_',datestr(cur_date, 'yyyymmdd'), '_v*_r*.sts']);
    if(isempty(file_to_convert))
        isfound = false;
        disp([datestr(cur_date) ' not found'])
    else
        filename = [monthpath file_to_convert.name];
        isfound = true;
    end
    
    if(isfound)
        Read_magnetic_field_data(filename, writepath);
    end
    
    cur_date = cur_date + datenum('01-Jan-0000');

    i = i+1;
end



