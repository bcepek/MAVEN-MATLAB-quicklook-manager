function filename = find_mag_file(cur_date)

datefrom = datestr(cur_date, 'yyyymmdd');	% '2015-01-03';

monthpath_MAG = '/Users/sesh2112/data/Matlab_MAVEN/pre-calculated/MAG_SS';

month_filelist_MAG =  dir([monthpath_MAG '/' 'mvn_mag_l2_*' datefrom '*sts.mat']);

if(length(month_filelist_MAG)<1)
    disp(['data for MAG for ' datefrom ' are not found'])
    filename = -1;
else
    filename = [monthpath_MAG '/' month_filelist_MAG(1).name];
end
