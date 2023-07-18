datefrom = datenum('2020-12-01', 'yyyy-mm-dd');
dateto = datenum('2020-12-31', 'yyyy-mm-dd');

writepath = '/Users/sesh2112/data/Matlab MAVEN/pre-calculated/MAG_SS';
root = '/Users/sesh2112/data/Matlab_MAVEN/MAG';

cur_date = datefrom;
while cur_date <= dateto
    filename = find_mag_file(cur_date);
    if(filename ~= -1)
        matfile = load(filename);
        if(~isfield(matfile, 'pos_sc_mso'))
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
            month_filelist = dir([monthpath, '*.sts']);
            %disp(i)
            for i = 1:length(month_filelist)
                if(all(month_filelist(i).name(28:29) == cur_day))
                    filename = [monthpath month_filelist(i).name];
                    isfound = true;
                    break
                elseif i == length(month_filelist)
                    isfound = false;
                    disp([datestr(cur_date) ' not found'])
                end
            end
            %disp(i)
            if(isfound)
                Read_magnetic_field_data(filename, writepath);
            end
        end
    end
    disp([datestr(cur_date, 'yyyy-mm-dd') ' processed'])
    cur_date = cur_date + datenum('01-Jan-0000');
end