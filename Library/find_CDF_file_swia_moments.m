
function filename_d1 = find_CDF_file_swia_moments(cur_date)
% finds STA CDF file for SWIA distribution function moments
% filename;

load('paths.mat', 'paths')
root = paths.swi_onboardsvymom;
%root = '/Users/sesh2112/data/Matlab_MAVEN/SWIA/onboard_svy_mom';

% cur_date = datenum(Date);

%cur_year = num2str(year(cur_date));
cur_year = datestr(cur_date, 'yyyy');
% cur_month = num2str(month(cur_date));
% if(length(cur_month) < 2)
%     cur_month = ['0' cur_month];
% end
cur_month = datestr(cur_date, 'mm');
% cur_day = num2str(day(cur_date));
% if(length(cur_day) < 2)
%     cur_day = ['0' cur_day];
%end
cur_day = datestr(cur_date, 'dd');
monthpath = [root '/' cur_year '/' cur_month '/'];
month_filelist = dir([monthpath, '*.cdf']);

if(length(month_filelist) == 0)
    isfound = false;
   filename_d1 = -1;
end

for i = 1:length(month_filelist)
    if(all(month_filelist(i).name(32:33) == cur_day))
        filename_d1 = [monthpath month_filelist(i).name];
        isfound = true;
        break
    elseif i == length(month_filelist)
        isfound = false;
        filename_d1 = -1;
    end
end



end
