function filename = find_mag_pc_file(cur_date)

datefrom = datestr(cur_date, 'yyyymmdd');	% '2015-01-03';

load('paths.mat', 'paths', 'slash')
monthpath_MAG = paths.mag_pc_1s;

%monthpath_MAG = '\\193.232.6.100\Data\Maven\data for MATLAB\MAG-pc';

month_filelist_MAG =  dir([monthpath_MAG slash 'mvn_mag_l2_*' datefrom '*sts.mat']);

if(length(month_filelist_MAG)<1)
    disp(['data for MAG for ' datefrom ' are not found'])
    filename = -1;
else
    filename = [monthpath_MAG slash month_filelist_MAG(1).name];
end
