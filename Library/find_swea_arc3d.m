function filename_svyspec = find_swea_arc3d(cur_date) 

%root = '\\193.232.6.100\Data\Maven\Static\d1_32e4d16a8m_Version2';
%root = '\\193.232.6.100\общие\Константин\Data\MAVEN\SWEA';
load('paths.mat', 'paths')
root = paths.swe_arc3d;
%root = '/Users/sesh2112/data/Matlab_MAVEN/SWEA/arc3d';
cur_year = datestr(cur_date, 'yyyy');
cur_month = datestr(cur_date, 'mm');
cur_day = datestr(cur_date, 'dd');

monthpath = [root '/' cur_year '/' cur_month '/'];
month_filelist = dir([monthpath, '*.cdf']);

if(isempty(month_filelist))
    isfound = false;
    disp([datestr(cur_date) ' not found'])
    filename_svyspec = -1;
end

for i = 1:length(month_filelist)
    if(all(month_filelist(i).name(24:25) == cur_day))
        filename_svyspec = [monthpath month_filelist(i).name];
        isfound = true;
        break
    elseif i == length(month_filelist)
        isfound = false;
        disp([datestr(cur_date) ' not found'])
        filename_svyspec = -1;
    end
end

end