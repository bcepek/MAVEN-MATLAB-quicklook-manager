function filename_svyspec = find_lpw_bursthf(cur_date) 

load('paths.mat', 'paths', 'slash')
root = paths.lpw_bursthf;
%root = '/Users/sesh2112/data/Matlab_MAVEN/LPW/we12';
cur_year = datestr(cur_date, 'yyyy');
cur_month = datestr(cur_date, 'mm');
cur_day = datestr(cur_date, 'dd');

monthpath = [root slash cur_year slash cur_month slash];
month_filelist = dir([monthpath, '*.cdf']);

if(isempty(month_filelist))
    isfound = false;
    disp(['lpw we12 data for ' datestr(cur_date) ' not found'])
    filename_svyspec = -1;
end

for i = 1:length(month_filelist)
    if(all(month_filelist(i).name(30:31) == cur_day))
        filename_svyspec = [monthpath month_filelist(i).name];
        isfound = true;
        break
    elseif i == length(month_filelist)
        isfound = false;
        disp(['lpw we12 data for ' datestr(cur_date) ' not found'])
        filename_svyspec = -1;
    end
end

end