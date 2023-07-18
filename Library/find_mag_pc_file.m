function filename = find_mag_pc_file(cur_date)

datefrom = datestr(cur_date, 'yyyymmdd');	% '2015-01-03';

monthpath_MAG = '\\193.232.6.100\Data\Maven\data for MATLAB\MAG-pc';

month_filelist_MAG =  dir([monthpath_MAG '\' 'mvn_mag_l2_*' datefrom '*sts.mat']);

if(length(month_filelist_MAG)<1)
    disp(['data for MAG for ' datefrom ' are not found'])
    filename = -1;
else
    filename = [monthpath_MAG '\' month_filelist_MAG(1).name];
end
